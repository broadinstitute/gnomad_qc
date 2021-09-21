from gnomad.sample_qc.relatedness import get_duplicated_samples, infer_families
from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources.sample_qc import *
from gnomad_qc.v2.resources import CURRENT_FAM, fam_path, get_gnomad_data
import numpy as np
import pandas as pd
import random
from datetime import datetime
from pyspark.sql import SQLContext
import argparse
import sys
import hail as hl
from typing import Dict, List, Set, Tuple
from collections import Counter
import logging

logger = logging.getLogger("create_fam")


class GnomADRelatedData:
    """
    This class encapsulates the necessary data to create trios.
    All its properties are lazily evaluated.

    """

    def __init__(self, data_type: str):
        self.data_type = data_type
        self._kin_ht = None
        self._meta_pd = None
        self._sample_project = None
        self._dups = None
        self._sample_to_dups = None


    @property
    def kin_ht(self) -> hl.Table:
        """
        Get the gnomAD kinship table computed by `pc_relate`

        :return: Kinship table
        :rtype: Table
        """
        if self._kin_ht is None:
            kin_ht = hl.read_table(relatedness_ht_path)
            kin_ht = kin_ht.filter((kin_ht.i.data_type == self.data_type) & (kin_ht.j.data_type == self.data_type))
            kin_ht = kin_ht.annotate(i=kin_ht.i.s, j=kin_ht.j.s).key_by('i', 'j')
            self._kin_ht = kin_ht
        return self._kin_ht

    @property
    def meta_pd(self) -> pd.DataFrame:
        """
        Get the metadata required for gnomAD trio inference

        :return: Metadata dataframe
        :rtype: DataFrame
        """
        if self._meta_pd is None:
            qc_ht = hl.read_table(qc_ht_path(self.data_type, part='hard_filters'))
            self._meta_pd = qc_ht.select(qc_ht.s,
                                   qc_ht.project_id,
                                   qc_ht.is_female
                                   ).to_pandas()
        return self._meta_pd

    @property
    def sample_project(self) -> Dict[str, str]:
        """
        Gets a dict mapping sample IDs with their project ID.

        :return: Dict of sample ID -> Project ID
        :rtype: dict of str -> str
        """
        if self._sample_project is None:
            self._sample_project = {row.s: row.project_id for row in self.meta_pd.itertuples()}
        return self._sample_project

    @property
    def dups(self) -> List[Set[str]]:
        """
        Gets the list of all sets of duplicated samples

        :return: Duplicated samples
        :rtype: list of set of str
        """
        if self._dups is None:
            self._dups = get_duplicated_samples(self.kin_ht)
        return self._dups


    @property
    def sample_to_dups(self) -> Dict[str, Set[str]]:
        """
        Gets a dict mapping each sample ID with the set of all sample IDs that are duplicates of that sample (based on `pc_relate` kinship)

        :return: Dict of sample -> duplicates
        :rtype: dict of str -> set of str
        """
        if self._sample_to_dups is None:
            self._sample_to_dups = dict()
            for dup in self.dups:
                for d in dup:
                    self._sample_to_dups[d] = dup.difference({d})
        return self._sample_to_dups


def get_dup_trios(
        pedigree: hl.Pedigree,
        sample_to_dups: Dict[str, Set[str]]
) -> pd.DataFrame:
    """
    Given a pedigree and a set of duplicate samples for each sample, creates a dataframe with one row per combination of duplicate samples.
    Rows are trios, columns are 'fam_id', 's', 'is_female', 'pat_id', 'mat_id'

    :param Pedigree pedigree: Input pedigree
    :param dict of str -> set of str sample_to_dups: for each sample, the set of all its duplicates
    :return: Dataframe with trios including combination of duplicates
    :rtype: DataFrame
    """

    dup_trios = []
    trio_id = 0
    for trio in pedigree.trios:
        trio_id += 1
        s_dup = [trio.s]
        if trio.s in sample_to_dups:
            s_dup.extend(sample_to_dups[trio.s])
        for s in s_dup:
            pat_dup = [trio.pat_id]
            if trio.pat_id in sample_to_dups:
                pat_dup.extend(sample_to_dups[trio.pat_id])
            for pat_id in pat_dup:
                mat_dup = [trio.mat_id]
                if trio.mat_id in sample_to_dups:
                    mat_dup.extend(sample_to_dups[trio.mat_id])
                for mat_id in mat_dup:
                    dup_trios.append((trio_id,
                                      s,
                                      trio.fam_id,
                                      pat_id,
                                      mat_id,
                                      trio.is_female))

    return pd.DataFrame(dup_trios, columns=['trio_id','s','fam_id','pat_id','mat_id','is_female'])


def pick_duplicates(
        dup_trios: pd.DataFrame,
        rank_pd: pd.DataFrame,
        sample_project: Dict[str, str]
) -> pd.DataFrame:
    """

    Picks duplicates to create trios according the following criteria:
    1) All samples (or as many as possible) within a trio should be in the same project
    2) As many trios within the same family should be in the same project
    3) Sum of the ranks of samples in the trio should be as low as possible

    :param DataFrame dup_trios: Input dataframe of trios, including all duplicate combinations
    :param DataFrame rank_pd: Dataframe with the rank of each sample
    :param dict of str -> str sample_project: mapping of sample ID to project ID
    :return: Dataframe of unique trios picking the optimal duplicate combination
    :rtype: DataFrame
    """

    # Add rank for each sample
    rank = rank_pd.set_index('s')
    dup_trios = dup_trios.set_index('s', drop=False).join(rank).rename(columns={'rank': 'rank_s'})
    dup_trios = dup_trios.set_index('pat_id', drop=False).join(rank).rename(columns={'rank': 'rank_pat'})
    dup_trios = dup_trios.set_index('mat_id', drop=False).join(rank).rename(columns={'rank': 'rank_mat'})
    dup_trios.reset_index()

    #Add trio-rank
    dup_trios['rank_trio'] = dup_trios.rank_s + dup_trios.rank_pat + dup_trios.rank_mat

    # Add project for each sample
    dup_trios = dup_trios.merge(
            dup_trios.apply(lambda row: pd.Series({
                'proj_s': sample_project[row.s],
                'proj_pat': sample_project[row.pat_id],
                'proj_mat': sample_project[row.mat_id]
            }), axis=1),
            left_index=True, right_index=True
        )

    # Annotate project with max number of samples in the trio
    def annotate_trio_project_stats(row):
        max_project = Counter([row.proj_s, row.proj_pat, row.proj_mat]).most_common(1)[0]
        return pd.Series({
            'max_project': max_project[0],
            'n_samples_project': max_project[1]
        })

    dup_trios = dup_trios.merge(
            dup_trios.apply(annotate_trio_project_stats, axis=1),
            left_index=True, right_index=True
        )

    # Annotate the family-level count of complete trios for each trio
    def annotate_proj_per_fam(df):
        df = df[df.n_samples_project == 3]
        return pd.DataFrame(list(Counter(df.max_project).items()), columns=['fam_project', 'n_fam_project'])

    proj_fam = dup_trios.groupby(dup_trios.fam_id).apply(annotate_proj_per_fam).reset_index()
    dup_trios = pd.merge(dup_trios, proj_fam, how='left', left_on=['fam_id','max_project'], right_on=['fam_id','fam_project'])

    # For each trio, keep a single set of dups based on
    # 1. Max number of samples in the trio in the same project
    # 2. Max number of trios in the family in the same project
    # 3. Rank of samples in the trio
    return dup_trios.groupby(dup_trios.trio_id).apply(lambda df: df.sort_values(by=['n_samples_project', 'n_fam_project', 'rank_trio'], ascending=[False, False, True])[0:1])


def ped_to_pandas(
        ped: hl.Pedigree,
        complete_trios: bool = True
) -> pd.DataFrame:
    """

    Creates a pandas dataframe from a hail pedigree.
    Rows are trios, columns are 'fam_id', 's', 'is_female', 'pat_id', 'mat_id'

    :param Pedigree ped: Input pedigree
    :param bool complete_trios: Whether to only output rows for complete trios
    :return: Dataframe
    :rtype: DataFrame
    """
    trios = ped.complete_trios() if complete_trios else ped.trios()
    return pd.DataFrame(
        [(trio.fam_id, trio.s, trio.is_female, trio.pat_id, trio.mat_id) for trio in trios],
        columns=['fam_id', 's', 'is_female', 'pat_id', 'mat_id']
    )


def pandas_to_ped(
        ped_pd: pd.DataFrame
):
    """
    Creates a Hail Pedigree object from trios stored as rows in a DataFrame.
    Input columns should contain 'fam_id', 's', 'is_female', 'pat_id', 'mat_id'

    :param DataFrame ped_pd: Input DataFrame
    :return: Pedigree
    :rtype: Pedigree
    """
    return hl.Pedigree([hl.Trio(s=row.s, is_female=row.is_female, pat_id=row.pat_id, mat_id=row.mat_id, fam_id=str(row.fam_id)) for row in ped_pd.itertuples()])


def merge_pedigree_pandas(
    pedigrees: List[Tuple[str, pd.DataFrame]],
    sample_to_dups: Dict[str, Set[str]],
    group_by_fam: bool = False
) -> pd.DataFrame:
    """

    Unionizes multiple pedigrees (in pandas format) into a single dataframe.
    Takes into consideration duplicates and will match individuals to a single sample name while
    keeping original sample name. The 'main' sample name is determined based on the order in the
    list of pedigrees (first pedigree has priority over 2nd, etc.)

    Output columns:
    's', 'pat_id', 'mat_id', 'is_female', 'fam_id', 'ped_name', 's_origin', 'pat_id_origin', 'mat_id_origin'

    Where the first 3 columns contain a unique name per sample taking duplicates into account.

    If `group_by_fam` is set, the output contains a single row per trio and the `*_origin`, `fam_id`, `is_female` and `ped_name` columns contain
    ordered comma-delimited `str` with values from all the different pedigrees where the trio appears.
    In addition an extra column `sex_inconsistent` is created containing a `bool` which is `True` if `is_female` was the same in all pedigrees with that trio.

    :param list of (str, DataFrame) pedigrees: Input pedigrees
    :param dict of str -> set of str sample_to_dups: sample to duplicates mapping
    :param bool group_by_fam: Whether output should be grouped by family
    :return: Dataframe with all pedigrees
    :rtype: DataFrame
    """

    def consolidate_names(ped_name, ped_pd, sample_main_name):
        new_names_pd = ped_pd.apply(lambda row: pd.Series({
            's': sample_main_name[row.s],
            'pat_id': sample_main_name[row.pat_id],
            'mat_id': sample_main_name[row.mat_id],
        }), axis=1)
        ped_pd = ped_pd.rename(columns={c: f'{c}_origin' for c in ['s', 'pat_id', 'mat_id']})
        ped_pd['ped_name'] = ped_name
        ped_pd = ped_pd.merge(new_names_pd, left_index=True, right_index=True)
        return ped_pd.set_index(['s', 'pat_id', 'mat_id'])

    # Create a directed samples -> dup dict
    sample_main_name = {}
    for ped_name, pedigree in pedigrees:
        for row in pedigree.itertuples():
            for s in [row.s, row.mat_id, row.pat_id]:
                if s not in sample_main_name:
                    s_names = [s]
                    if s in sample_to_dups:
                        s_names.extend(sample_to_dups[s])
                    for s_name in s_names:
                        sample_main_name[s_name] = s

    # Consolidate names and concatenate dataframes
    ped_pd = pd.concat([consolidate_names(ped_name, current_ped_pd, sample_main_name) for ped_name, current_ped_pd in pedigrees])

    if group_by_fam:
        ped_pd = ped_pd.groupby(['s', 'pat_id', 'mat_id']).apply(lambda df: pd.Series({
            'is_female': list(df.is_female)[0],
            'fam_id': ",".join(df.fam_id),
            'ped_name': ",".join(df.ped_name),
            's_origin': ",".join(df.s_origin),
            'pat_id_origin': ",".join(df.pat_id_origin),
            'mat_id_origin': ",".join(df.pat_id_origin),
            'sex_inconsistent': len(set(df.is_female)) > 1
        }))

    ped_pd = ped_pd.reset_index()
    return ped_pd


def add_pedigree_meta(
        ped_pd: pd.DataFrame,
        kin_ht: hl.Table = None,
        meta_pd: pd.DataFrame = None,
        mendel_per_sample_ht: hl.Table = None,
) -> pd.DataFrame:
    """

    Adds the following metadata columns to a pedigree dataframe if specified:
    * `kin_ht`: Adds parent/child and parent/parent kin and k2
    * `meta_pd`: Adds all sample-level metadata found in dataframe to s, pat, and mat
    * `mendel_per_sample`: Adds mendel errors

    :param DataFrame ped_pd: Input pedigree dataframe
    :param Table kin_ht: Kinship table (e.g. from `pc_relate`). Expects cols `i` and `j` as strings and `kin`, `ibd2` columns.
    :param DataFrame meta_pd: Sample-level metadata
    :param Table mendel_per_sample_ht: Per-sample mendel errors from hail `mendel_errors`
    :return: Dataframe annotated with given metadata
    :rtype: DataFrame
    """
    if mendel_per_sample_ht is not None:
        mv_per_sample_pd = mendel_per_sample_ht.to_pandas()
        ped_pd = ped_pd.set_index(['s', 'fam_id']).join(mv_per_sample_pd.set_index(['s', 'fam_id']))
        ped_pd = ped_pd.reset_index()

    if meta_pd is not None:
        meta_pd = meta_pd.set_index('s')
        ped_pd = ped_pd.set_index('s', drop=False).join(
            meta_pd.rename(columns={col: f'{col}_s' for col in meta_pd.columns }), rsuffix="_meta")
        ped_pd = ped_pd.set_index('mat_id', drop=False).join(
            meta_pd.rename(columns={col: f'{col}_mat' for col in meta_pd.columns }), rsuffix="_meta")
        ped_pd = ped_pd.set_index('pat_id', drop=False).join(
            meta_pd.rename(columns={col: f'{col}_pat' for col in meta_pd.columns }), rsuffix="_meta")
        ped_pd = ped_pd.reset_index(drop=True)

    if kin_ht is not None:
        ped_samples = hl.literal({s for s in list(ped_pd.s) + list(ped_pd.mat_id) + list(ped_pd.pat_id)})
        kin = {tuple(sorted((row.i, row.j))): (row.kin, row.ibd2) for row in kin_ht.filter(ped_samples.contains(kin_ht.i) & ped_samples.contains(kin_ht.j)).collect()}

        def get_kin(row):
            s_pat_kin = kin.get(tuple(sorted((row.s, row.pat_id))))
            s_mat_kin = kin.get(tuple(sorted((row.s, row.mat_id))))
            pat_mat_kin = kin.get(tuple(sorted((row.mat_id, row.pat_id))))

            return pd.Series({
                's_pat_kin': s_pat_kin[0] if s_pat_kin is not None else None,
                's_pat_ibd2': s_pat_kin[1] if s_pat_kin is not None else None,
                's_mat_kin': s_mat_kin[0] if s_mat_kin is not None else None,
                's_mat_ibd2': s_mat_kin[1] if s_mat_kin is not None else None,
                'pat_mat_kin': pat_mat_kin[0] if pat_mat_kin is not None else None,
                'pat_mat_ibd2': pat_mat_kin[1] if pat_mat_kin is not None else None
            })

        ped_pd = ped_pd.merge(
            ped_pd.apply(get_kin, axis=1),
            left_index=True, right_index=True
        )

    return ped_pd


def create_fake_pedigree(
        n: int,
        sample_list: List[str],
        real_pedigree: hl.Pedigree = None
) -> hl.Pedigree:
    """

    Generates a "fake" pedigree made of trios created by sampling 3 random samples in the sample list.
    If `real_pedigree` is given, then children from the real pedigrees won't be used as probands.
    This functions insures that:
    - All probands are unique
    - All individuals in a trio are different

    :param int n: Number of fake trios desired in the pedigree
    :param list of str sample_list: List of samples
    :param Pedigree real_pedigree: Optional pedigree to exclude children from
    :return: Fake pedigree
    :rtype: Pedigree
    """

    probands = set()
    if real_pedigree is not None:
        probands = {trio.s for trio in real_pedigree.trios}.intersection(set(sample_list))
        if len(probands) == len(sample_list):
            raise ValueError("Full sample list for fake trios generation needs to include samples that aren't probands in the real trios.")

    fake_trios = []
    for i in range(n):
        mat_id, pat_id = random.sample(sample_list, 2)
        s = random.choice(sample_list)
        while s in probands.union({mat_id, pat_id}):
            s = random.choice(sample_list)

        probands.add(s)

        fake_trios.append(hl.Trio(
            s=s,
            pat_id=pat_id,
            mat_id=mat_id,
            fam_id=str(str(i+1)),
            is_female=True
        ))

    return hl.Pedigree(fake_trios)


def infer_ped(related_data: GnomADRelatedData) -> None:
    """

    Infers trios based on `pc_relate` kinship output.
    Writes a CSV containing one row per trio.
    If there are duplicate samples, each combination of duplicate samples will be present in the output.

    :param GnomADRelatedData related_data: Input data for inference
    :return: Nothing
    :rtype: None
    """
    
    logger.info(f"Inferring pedigree for {related_data.data_type}")
    sex = {row.s: row.is_female for row in related_data.meta_pd.itertuples()}
    dups_to_remove = {s for d in related_data.dups for s in list(d)[1:]}

    raw_ped = infer_families(related_data.kin_ht, sex, dups_to_remove)
    logger.info(f"Found {len(raw_ped.complete_trios())} complete trios in {related_data.data_type}.")

    # Create dataframe with all combinations of dups
    dup_trios_pd = get_dup_trios(raw_ped, related_data.sample_to_dups)
    logger.info(f"Found {len(dup_trios_pd)} trios combinations with dups in {related_data.data_type}.")
    with hl.hadoop_open(dup_pedigree_tsv_path(related_data.data_type), 'w') as out:
        dup_trios_pd.to_csv(out, sep="\t", index=False)


def select_dups(related_data: GnomADRelatedData) -> None:
    """
    Loads gnomAD CSV file containing pedigrees with duplicates and writes the raw gnomAD ped file by selecting a single row per combination of duplicate samples maximizing:
    1) All samples (or as many as possible) within a trio should be in the same project
    2) As many trios within the same family should be in the same project
    3) Sum of the quality ranks (based on sample qc) of samples in the trio should be as low as possible

    :param GnomADRelatedData related_data: Input data
    :return: Nothing
    :rtype: None
    """
    dup_trios_pd = hl.import_table(dup_pedigree_tsv_path(related_data.data_type), impute=True).to_pandas()

    rank_ht = hl.import_table(rank_annotations_path('joint'), impute=True)
    rank_ht = rank_ht.filter(rank_ht.data_type == related_data.data_type)
    rank_pd = rank_ht.select(rank_ht.s, rank_ht.rank).to_pandas()

    raw_trios = pick_duplicates(dup_trios_pd, rank_pd, related_data.sample_project)
    raw_ped = pandas_to_ped(raw_trios)
    logger.info(f"Generated {len(raw_ped.complete_trios())} complete trios after optimizing families in {related_data.data_type}.")
    raw_ped.write(raw_fam_path(related_data.data_type))


def create_meta(related_data: GnomADRelatedData, fake_fam_prop: float, old_version: str, overwrite: bool) -> None:
    """
    Creates and writes a dataframe with metadata to evaluate gnomAD trios from the raw ped file.
    In order to compare the raw ped, metadata is also generated for:
    1) A number of fake families are generated
    2) The previous iteration of the ped file (old_version)

    :param GnomADRelatedData related_data: Input data
    :param float fake_fam_prop: Number of fake trios to generate as a proportion of the number of real families in the data
    :param str old_version: Version of previous iteration to load
    :param bool overwrite: Whether to overwrite previous data
    :return: Nothing
    :rtype: None
    """

    raw_ped = hl.Pedigree.read(raw_fam_path(related_data.data_type), delimiter="\\t")

    n_fake_trios = int(fake_fam_prop * len(raw_ped.complete_trios()))
    logger.info(f"Generating fake pedigree with {n_fake_trios} trios for {related_data.data_type}")
    fake_fams = create_fake_pedigree(n_fake_trios, list(related_data.meta_pd.s), raw_ped)

    fake_fams.write(fake_fam_path(related_data.data_type))

    logger.info(f"Running mendel_errors on {related_data.data_type}")

    # Run mendel errors on families made of random samples to establish expectation in non-trios:
    pedigrees = [('new', raw_ped),
                 ('old', hl.Pedigree.read(fam_path(related_data.data_type, version=old_version), delimiter="\\t")),
                 ('fake', hl.Pedigree.read(fake_fam_path(related_data.data_type), delimiter="\\t"))]

    ped_pd = merge_pedigree_pandas([(name, ped_to_pandas(ped)) for name, ped in pedigrees],
                                   related_data.sample_to_dups,
                                   True)

    # Run mendel_errors
    all_ped = pandas_to_ped(ped_pd)
    gnomad = get_gnomad_data(related_data.data_type)
    fam_samples = hl.literal({s for trio in all_ped.trios for s in [trio.s, trio.mat_id, trio.pat_id]})
    gnomad = gnomad.filter_cols(fam_samples.contains(gnomad.s))
    all_errors, per_fam, per_sample, _ = hl.mendel_errors(gnomad['GT'], all_ped)

    all_errors.write(sample_qc_mendel_ht_path(related_data.data_type, "all_errors"), overwrite=overwrite)
    per_fam.write(sample_qc_mendel_ht_path(related_data.data_type, "per_fam"), overwrite=overwrite)
    per_sample.write(sample_qc_mendel_ht_path(related_data.data_type, "per_sample"), overwrite=overwrite)

    # Merge all metadata
    ped_pd = add_pedigree_meta(ped_pd=ped_pd,
                               meta_pd=related_data.meta_pd,
                               kin_ht=related_data.kin_ht,
                               mendel_per_sample_ht=per_sample)

    # Write merged pedigrees as HT
    sql_context = SQLContext(hl.spark_context())
    hl.Table.from_spark(sql_context.createDataFrame(ped_pd)).write(merged_pedigrees_ht_path(related_data.data_type), overwrite=overwrite)
    

def create_ped(related_data: GnomADRelatedData, new_version: str, max_mv_z: int = 3):
    """

    Loads the raw gnomAD ped, applies Mendelian Violations cutoff in order to produce a final ped file and writes final gnomAD ped file.

    :param GnomADRelatedData related_data: Input data
    :param str new_version: String containing the new version name to write the data to
    :param int max_mv_z: Max number of std devs above the mean number of MVs in inferred trios to keep trio.
    :return: Nothing
    :rtype: None
    """
    raw_ped = hl.Pedigree.read(raw_fam_path(related_data.data_type), delimiter="\\t")
    logger.info(f"Loaded raw {related_data.data_type} pedigree containing {len(raw_ped.trios)} trios")

    # Filter families
    ped_meta = hl.read_table(merged_pedigrees_ht_path(related_data.data_type)).to_pandas()

    ped_meta = ped_meta[ped_meta.ped_name.str.contains('new')]
    mean_errors = np.mean(ped_meta.errors)
    std_errors = np.std(ped_meta.errors)

    filtered_s = set(ped_meta[ped_meta.errors > mean_errors + max_mv_z * std_errors].s)

    # Write final fam file
    final_ped = hl.Pedigree([trio for trio in raw_ped.trios if trio.s not in filtered_s])
    final_ped.write(fam_path(related_data.data_type, version=new_version))

    logger.info(f"Wrote final {related_data.data_type} pedigree with {len(final_ped.trios)} trios.")


def main(args):

    data_types = []
    if args.exomes:
        data_types.append('exomes')
    if args.genomes:
        data_types.append('genomes')

    for data_type in data_types:

        related_data = GnomADRelatedData(data_type)

        # Infer families from the data
        if args.infer_ped:
            infer_ped(related_data)

        #Select duplicates and write raw ped file
        if args.select_dups:
            select_dups(related_data)

        #Create trio metadata for evaluation
        if args.create_meta:
            create_meta(related_data, fake_fam_prop=args.fake_fam_prop, old_version=args.old_version, overwrite=args.overwrite)

        #Apply mendel errors threshold and create final ped file
        if args.create_ped:
            create_ped(related_data, max_mv_z=args.max_mv_z, new_version=args.new_version)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    parser.add_argument('--exomes', help='Use exomes data.', action='store_true')
    parser.add_argument('--genomes', help='Use genomes data.', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')

    infer_ped_grp = parser.add_argument_group('Pedigree inference')
    infer_ped_grp.add_argument('--infer_ped', help='Infer pedigrees and create a table with all trios with all duplicates.', action='store_true')

    select_dup_grp = parser.add_argument_group('Select duplicates')
    select_dup_grp.add_argument('--select_dups', help='Select a single duplicate sample for each trio, optimizing for (1) number of trio/family members within the same project, and (2) rank in sampleqc list.', action='store_true')

    meta_grp = parser.add_argument_group('Pedigree quality control (metadata creation)')
    meta_grp.add_argument('--fake_fam_prop', help='Number of fake trios to generate as a proportion of the total number of trios found in the data. (default=0.1)', default=0.1, type=float)
    meta_grp.add_argument('--create_meta', help='Creates a HT with metadata (mendel error, kinship, sample meta) about trios (old, new and fake to allow for comparison).', action='store_true')
    meta_grp.add_argument('--old_version', help='Version for the old data. (default=current version specified in resources)', default=CURRENT_FAM)

    final_ped_grp = parser.add_argument_group('Pedigree filtering (final pedigree generation)')
    final_ped_grp.add_argument('--create_ped', help='Creates filtered pedigree file for downstream usage', action='store_true')
    final_ped_grp.add_argument('--new_version', help="Version for the new data. (default='YYYY-MM-DD')", default=datetime.today().strftime('%Y-%m-%d'))
    final_ped_grp.add_argument('--max_mv_z', help='Max number of std devs above the mean number of MVs in inferred trios to keep trio. (default=3)', default=3, type=int)

    args = parser.parse_args()
    if not args.exomes and not args.genomes:
        sys.exit('Error: At least one of --exomes or --genomes must be specified.')

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)

