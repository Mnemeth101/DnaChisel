from Bio.Data import CodonTable

from dnachisel import (
    AvoidFrameshiftStartStopCodons,
    DnaOptimizationProblem,
    EnforceTranslation,
    translate,
)


def has_frameshift_start_stop(sequence, genetic_table="Standard"):
    table = CodonTable.unambiguous_dna_by_name[genetic_table]
    start_codons = set(table.start_codons)
    stop_codons = set(table.stop_codons)
    for frame in (1, 2):
        for index in range(frame, len(sequence) - 2, 3):
            codon = sequence[index : index + 3]
            if (codon in start_codons) or (codon in stop_codons):
                return True
    return False


def test_AvoidFrameshiftStartStopCodons():
    sequence = "AATGAACTGCAAGCTGAA"
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            EnforceTranslation(),
            AvoidFrameshiftStartStopCodons(),
        ],
        logger=None,
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    assert not has_frameshift_start_stop(problem.sequence)
    assert translate(problem.sequence) == translate(sequence)
