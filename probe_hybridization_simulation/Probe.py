from collections import namedtuple

import numpy as np

import Target


class Probe:
    def __init__(self,
                 target,
                 lenght: int):
        self.target = target
        self.length = lenght
        assert self.length <= len(target.seq), "Probe is longer than target"
        assert self.length >= 4, "Probe is too small"


    def __generate_output(self,
                          sequence: str,
                          probe_type: str,
                          structure: dict):
        output = namedtuple("Probe", "seq type str")
        return output(sequence, probe_type, structure)


    def genertare_probe(self):
        probe_start_allow = len(self.target.seq) - self.length
        probe_start = np.random.random_integers(1, high=probe_start_allow)
        probe_end = probe_start + self.length
        probe_seq = self.target.seq[probe_start:probe_end]
        probe_info = Target.Target.get_fragment_info(target=self.target,
                                                     start=probe_start,
                                                     end=probe_end)
        probe_info_output = namedtuple("Probe",
                                       "seq start end len inum enum str")

        return probe_info_output(seq=probe_seq,
                                 start=probe_info.start,
                                 end=probe_info.end,
                                 len=len(probe_seq),
                                 inum=probe_info.inum,
                                 enum=probe_info.enum,
                                 str=probe_info.str)


    def mutate_probe(self, hamming_distance: int):
        pass
