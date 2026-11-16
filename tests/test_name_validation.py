import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from tiddit.name_validation import names_unusable_in_file_names


def test_ordinary_names_are_accepted():
    assert names_unusable_in_file_names(
        ["sample1", "chr1", "HLA-A*01:01:01:01", "chrUn_KI270302v1", ".."]
    ) == []


def test_names_with_a_path_separator_are_reported_once_and_sorted():
    assert names_unusable_in_file_names(
        ["sample1", "batch1/sample2", "../../chr1", "chr2", "batch1/sample2"]
    ) == ["../../chr1", "batch1/sample2"]
