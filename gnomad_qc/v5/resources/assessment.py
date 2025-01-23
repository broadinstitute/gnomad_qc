"""Script containing assessment resources."""


def get_indexed_array_for_missingness_path() -> str:
    """
    Return path to indexed array annotations to use for missingness checks for federated data.

    The first column contains the keys (names of the array annotations), and the second column contains the values (names of the globals containing
    the mapping of group name to index for that key).

    :return: Output path of file containing indexed array annotations to use for missingness checks.
    """
    return "gs://gnomad/v5.0/assessment/indexed_arrays.txt"
