"""Demo for AvoidFrameshiftStartStopCodons."""

from Bio.Data import CodonTable

from dnachisel import (
    AvoidFrameshiftStartStopCodons,
    DnaOptimizationProblem,
    EnforceTranslation,
    translate,
)


def find_frameshift_start_stop(sequence, genetic_table="Standard"):
    table = CodonTable.unambiguous_dna_by_name[genetic_table]
    start_codons = set(table.start_codons)
    stop_codons = set(table.stop_codons)
    hits = []
    for frame in (1, 2):
        for index in range(frame, len(sequence) - 2, 3):
            codon = sequence[index : index + 3]
            if codon in start_codons or codon in stop_codons:
                hits.append((frame, index, codon))
    return hits


def main():
    sequence = "AATGAACTGCAAGCTGAA"
    print("Initial translation (+0):", translate(sequence))
    print("Frameshift hits:", find_frameshift_start_stop(sequence))

    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            EnforceTranslation(),
            AvoidFrameshiftStartStopCodons(),
        ],
        logger=None,
    )

    print("Constraints pass before:", problem.all_constraints_pass())
    problem.resolve_constraints()
    print("Constraints pass after:", problem.all_constraints_pass())

    print("Optimized sequence:", problem.sequence)
    print("Optimized translation (+0):", translate(problem.sequence))
    print("Frameshift hits after optimization:", find_frameshift_start_stop(problem.sequence))


if __name__ == "__main__":
    main()
