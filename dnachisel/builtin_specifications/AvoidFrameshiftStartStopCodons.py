from Bio.Data import CodonTable

from ..Location import Location
from ..Specification import SpecEvaluation
from ..Specification.Specification import Specification


class AvoidFrameshiftStartStopCodons(Specification):
    """Avoid start/stop codons in the +1/+2 frames relative to the coding frame.

    This specification checks the alternative reading frames (+1 and +2) of the
    specified region (or the whole sequence if no location is provided) and
    flags any start or stop codons found there.
    """

    best_possible_score = 0

    def __init__(self, genetic_table="Standard", location=None, frames=(1, 2), boost=1.0):
        self.genetic_table = genetic_table
        self.location = Location.from_data(location)
        self.frames = frames
        self.boost = boost

    def initialized_on_problem(self, problem, role):
        return self._copy_with_full_span_if_no_location(problem)

    def localized(self, location, problem=None, with_righthand=True):
        if self.location is None:
            return self
        if self.location.overlap_region(location) is None:
            return None
        extended_location = location.extended(2, right=with_righthand)
        new_location = self.location.overlap_region(extended_location)
        return self.copy_with_changes(location=new_location)

    def _codon_location(self, location, codon_start_index):
        indices = location.indices
        codon_indices = indices[codon_start_index : codon_start_index + 3]
        if len(codon_indices) < 3:
            return None
        start = min(codon_indices)
        end = max(codon_indices) + 1
        return Location(start, end, location.strand)

    def evaluate(self, problem):
        location = (
            self.location
            if self.location is not None
            else Location(0, len(problem.sequence))
        )
        subsequence = location.extract_sequence(problem.sequence)
        table = CodonTable.unambiguous_dna_by_name[self.genetic_table]
        start_codons = set(table.start_codons)
        stop_codons = set(table.stop_codons)

        errors_locations = []
        for frame in self.frames:
            for index in range(frame, len(subsequence) - 2, 3):
                codon = subsequence[index : index + 3]
                if (codon in start_codons) or (codon in stop_codons):
                    location_match = self._codon_location(location, index)
                    if location_match is not None:
                        errors_locations.append(location_match)

        return SpecEvaluation(
            self,
            problem,
            score=-len(errors_locations),
            locations=errors_locations,
            message=(
                "All OK."
                if len(errors_locations) == 0
                else "Start/stop codons found at %s" % errors_locations
            ),
        )

    def __str__(self):
        return "AvoidFrameshiftStartStopCodons(%s)" % self.location
