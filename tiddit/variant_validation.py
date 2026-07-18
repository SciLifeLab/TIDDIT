"""Validation helpers for structural-variant candidate intervals."""


def invalid_candidate_interval_diagnostics(
    chr_a,
    chr_b,
    cluster,
    start_a,
    end_a,
    start_b,
    end_b,
    n_discordants,
    n_splits,
    n_contigs,
):
    """Return deterministic diagnostics for invalid inclusive intervals.

    TIDDIT candidate intervals are inclusive, so their width is
    ``end - start + 1``.  Any interval with ``end < start`` is invalid and
    must be rejected rather than reordered or assigned synthetic coverage.
    """

    diagnostics = []
    for side, contig, start, end in (
        ("A", chr_a, start_a, end_a),
        ("B", chr_b, start_b, end_b),
    ):
        width = end - start + 1
        if width > 0:
            continue
        diagnostics.append(
            "\t".join(
                (
                    "TIDDIT_REJECT_INVALID_INTERVAL",
                    f"chrA={chr_a}",
                    f"chrB={chr_b}",
                    f"cluster={cluster}",
                    f"side={side}",
                    f"contig={contig}",
                    f"start={start}",
                    f"end={end}",
                    f"width={width}",
                    f"N_discordants={n_discordants}",
                    f"N_splits={n_splits}",
                    f"N_contigs={n_contigs}",
                )
            )
        )
    return diagnostics
