# gnomad_qc

This repo contains the complete set of scripts used to perform sample and variant QC for the gnomAD v2 and v3 releases. We will continue to update and improve upon the code to handle new releases as they grow in size and complexity and as they require increasingly sophisticated QC treatment. The current code therefore represents the most recent iteration of our pipelines and is guaranteed to change over time.

NB: The scripts make reference to gnomAD-related metadata files (not public) and may perform procedures that are not strictly necessary for quality control of all germline datasets. For example, the gnomAD v2 dataset comprises both exomes and genomes, and a substantial portion of the code is written to handle technical differences between those call sets, as well as to perform relevant joint analyses (such as inferring cryptically related individuals across exomes and genomes). These steps may not be relevant for all call sets.

We therefore encourage users to browse through the code and identify modules and functions that will be useful in their own pipelines, and to edit and reconfigure the gnomAD pipeline to suit their particular analysis and QC needs.

A more extensive overview and explanation of the gnomAD QC process is available on the [gnomAD browser](https://gnomad.broadinstitute.org/news/) and may help inform users’ design decisions for other pipelines.
  * [v2.1](https://gnomad.broadinstitute.org/news/2018-10-gnomad-v2-1/)
  * [v3](https://gnomad.broadinstitute.org/news/2019-10-gnomad-v3-0/)
  * [v3.1](https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/)

Note also that many basic functions and file paths used in the code are imported from a separate repo, [gnomad_methods](https://github.com/broadinstitute/gnomad_methods).
