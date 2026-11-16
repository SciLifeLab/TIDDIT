"""Validation helpers for names that end up in work file paths."""


def names_unusable_in_file_names(names):
    """Return the names that cannot be used as a single file name.

    TIDDIT names its work files after the sample and the contigs in the BAM
    header. A name containing a path separator would point the file into
    another directory instead.
    """

    return sorted(set(name for name in names if "/" in name))
