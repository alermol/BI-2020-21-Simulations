
from collections import namedtuple
import tqdm

import numpy as np
from Bio.Seq import Seq

import Target


class NickTranslation:
    """
    Generate nick-translation products from probe with uniformly
    distributed from a to b size
    """
    def __init__(self,
                 probe,
                 min_len: int = 590,
                 max_len: int = 600):
        self.probe = probe
        self.min_len = min_len
        self.max_len = max_len


    def __generate_output(self):
        pass


    def generate_possible_fragments(self):
        fragments_pool = []
        for k in tqdm.tqdm(range(self.min_len, self.max_len + 1),
                           desc="Fragments generation"):
            fragments = self.__sliding_window(sequence=self.probe.seq,
                                              win_size=k)
            fragments_pool += [i for i in fragments]
            fragments_pool_ss = [i[::-1] for i in fragments_pool]
        return tuple(fragments_pool + fragments_pool_ss)


    def __sliding_window(self,
                         sequence: str,
                         win_size: int):
        num_of_chunks = len(sequence) - win_size + 1
        for i in range(0, num_of_chunks):
            yield sequence[i:(i + win_size)]
