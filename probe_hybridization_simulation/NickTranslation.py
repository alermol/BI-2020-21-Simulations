
from collections import namedtuple
from itertools import islice

import numpy as np

import Target


class NickTranslation:
    """
    Generate nick-translation products from probe with uniformly
    distributed from a to b size
    """
    def __init__(self,
                 probe,
                 min_len: int = 200,
                 max_len: int = 600):
        self.probe = probe
        self.min_len = min_len
        self.max_len = max_len


    def __generate_output(self):
        pass


    def generate_possible_fragments(self):
        fragments_pool = []
        for k in range(self.min_len, self.max_len + 1):
            fragments = self.__sliding_window(sequence=self.probe.seq,
                                              win_size=k)
            fragments_pool += [i for i in fragments]
        return tuple(fragments_pool)


    def __sliding_window(self,
                         sequence: str,
                         win_size: int):
        num_of_chunks = len(sequence) - win_size + 1
        for i in range(0, num_of_chunks):
            yield sequence[i:(i + win_size)]
