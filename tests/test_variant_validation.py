import sys
from pathlib import Path
from types import SimpleNamespace

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from tiddit.variant_validation import invalid_candidate_interval_diagnostics


def _write_empty_bam(pysam, bam_path, contigs=("chr1",)):
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": contig, "LN": 1000} for contig in contigs],
    }
    with pysam.AlignmentFile(bam_path, "wb", header=header):
        pass


def test_valid_candidate_intervals_are_unchanged():
    assert invalid_candidate_interval_diagnostics(
        "chr1", "chr2", 17, 100, 100, 200, 250, 6, 4, 1
    ) == []


def test_zero_width_candidate_interval_has_exact_diagnostic():
    assert invalid_candidate_interval_diagnostics(
        "chr1", "chr2", 17, 100, 150, 201, 200, 6, 4, 1
    ) == [
        "\t".join(
            (
                "TIDDIT_REJECT_INVALID_INTERVAL",
                "chrA=chr1",
                "chrB=chr2",
                "cluster=17",
                "side=B",
                "contig=chr2",
                "start=201",
                "end=200",
                "width=0",
                "N_discordants=6",
                "N_splits=4",
                "N_contigs=1",
            )
        )
    ]


@pytest.mark.parametrize(
    "start_a,end_a,start_b,end_b,expected_sides,expected_widths",
    [
        (100, 98, 200, 250, ["side=A"], ["width=-1"]),
        (100, 99, 251, 250, ["side=A", "side=B"], ["width=0", "width=0"]),
    ],
)
def test_all_invalid_sides_are_reported_without_swapping_bounds(
    start_a, end_a, start_b, end_b, expected_sides, expected_widths
):
    diagnostics = invalid_candidate_interval_diagnostics(
        "chr1", "chr2", "cluster-2", start_a, end_a, start_b, end_b, 3, 2, 0
    )

    assert len(diagnostics) == len(expected_sides)
    for diagnostic, side, width in zip(diagnostics, expected_sides, expected_widths):
        assert side in diagnostic
        assert width in diagnostic
    assert f"start={start_a}" in diagnostics[0]
    assert f"end={end_a}" in diagnostics[0]


def test_compiled_get_region_rejects_zero_width_before_fetch(tmp_path):
    pysam = pytest.importorskip("pysam")
    variant = pytest.importorskip("tiddit.tiddit_variant")

    bam_path = tmp_path / "empty.bam"
    _write_empty_bam(pysam, bam_path)
    pysam.index(str(bam_path))

    with pysam.AlignmentFile(bam_path, "rb") as alignment:
        with pytest.raises(
            ValueError,
            match=r"^TIDDIT_INVALID_INTERVAL contig=chr1 start=101 end=100 width=0$",
        ):
            variant.get_region(alignment, "chr1", 101, 100, 101, 5, 500, {"chr1": 0})


def test_invalid_candidate_is_rejected_and_logged_deterministically(tmp_path, capsys):
    pysam = pytest.importorskip("pysam")
    variant = pytest.importorskip("tiddit.tiddit_variant")

    bam_path = tmp_path / "empty.bam"
    _write_empty_bam(pysam, bam_path, contigs=("chr1", "chr2"))
    prefix = tmp_path / "run"
    index_dir = tmp_path / "run_tiddit"
    index_dir.mkdir()
    csi_path = index_dir / "sample.csi"
    pysam.samtools.index("-c", "-o", str(csi_path), str(bam_path), catch_stdout=False)

    args = SimpleNamespace(
        ref=None,
        o=str(prefix),
        p=3,
        r=3,
        skip_assembly=True,
        threads=2,
    )
    cluster = {
        "chr1": {
            "chr2": {
                17: {
                    "N_discordants": 6,
                    "N_splits": 4,
                    "N_contigs": 1,
                    "startA": 100,
                    "endA": 150,
                    "startB": 201,
                    "endB": 200,
                }
            }
        }
    }

    calls = variant.main(
        str(bam_path),
        cluster,
        args,
        {},
        5,
        ["sample"],
        {},
        {"chr1": 0, "chr2": 1},
        500,
        {},
    )

    assert calls == {"chr1": [], "chr2": []}
    assert capsys.readouterr().err == (
        "TIDDIT_REJECT_INVALID_INTERVAL\tchrA=chr1\tchrB=chr2\tcluster=17"
        "\tside=B\tcontig=chr2\tstart=201\tend=200\twidth=0"
        "\tN_discordants=6\tN_splits=4\tN_contigs=1\n"
    )
