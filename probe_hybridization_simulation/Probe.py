from collection import namedtuple

import nampy as np

class Probe:
    def __init__(self,
                 target,
                 lenght: int):
        self.target = target
        self.length = lenght
        assert self.length <= len(target.seq), "Probe is longer than target"


    def __generate_output(self,
                          sequence: str,
                          probe_type: str,
                          structure: dict):
        output = namedtuple("Probe", "seq type str")
        return output(sequence, probe_type, structure)


    def genertare_probe(self):
        pass


    def mutate_probe(self, hamming_distance: int):
        pass
