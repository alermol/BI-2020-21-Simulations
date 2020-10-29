from collections import namedtuple

import numpy as np

import Target


class Probe:
    def __init__(self,
                 target,
                 lenght: int,
                 vector_len: int = 3000,
                 vector_gc_content: float = 0.5):
        self.target = target
        self.length = lenght
        self.vector_len = vector_len
        self.vector_gc_content = vector_gc_content
        self.ALPHABET = ["A", "T", "G", "C"]
        assert self.length <= len(target.seq), "Probe is longer than target"
        assert self.length >= 4, "Probe is too small"
        assert self.vector_len > 200, "Vector is too short"
        assert 0.0 <= vector_gc_content <= 1.0, "Wrong GC-content in vector"


    def __generate_output(self,
                          sequence: str,
                          probe_type: str,
                          structure: dict):
        output = namedtuple("Probe", "seq type str")
        return output(sequence, probe_type, structure)


    def generate_vector(self):
        gc_probability = self.vector_gc_content / 2
        at_probability = (1 - self.vector_gc_content) / 2
        nucl_probability = ([at_probability] * 2) + ([gc_probability] * 2)
        vector_seq = np.random.choice(self.ALPHABET,
                                      size=self.vector_len,
                                      p=nucl_probability)
        return "".join(vector_seq)


    def generate_probe(self):
        probe_start_allow = len(self.target.seq) - self.length
        probe_start = np.random.randint(0, high=probe_start_allow + 1)
        probe_end = probe_start + self.length
        probe_info = Target.Target.get_fragment_info(target=self.target,
                                                     start=probe_start,
                                                     end=(probe_end - 1))
        probe_seq = probe_info.seq + self.generate_vector()
        probe_start = probe_info.start # probe start on target sequence
        probe_end = probe_info.end # probe end on target sequence
        probe_structure = probe_info.adj_str
        probe_structure["vector"] = (self.length,
                                     self.vector_len + self.length - 1)
        probe_info_output = namedtuple("Probe",
                                       "seq start end ins_len str")
        return probe_info_output(seq=probe_seq,
                                 start=probe_start,
                                 end=probe_end,
                                 ins_len=self.length,
                                 str=probe_structure)
