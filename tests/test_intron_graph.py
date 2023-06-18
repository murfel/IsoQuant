import os
from collections import namedtuple

import gffutils
import pytest

from src.gene_info import GeneInfo
from src.intron_graph import IntronGraph, IntronCollector


class TestIntronGraph:
    source_dir = os.path.dirname(os.path.realpath(__file__))
    gffutils_db = gffutils.FeatureDB(os.path.join(source_dir, 'toy_data/synth.db'), keep_order=True)
    gene_db = gffutils_db['ENSMUSG00000020196.10']

    def __init__(self):
        self.gene_info = GeneInfo([self.gene_db], self.gffutils_db)

    def test_add_edge(self):
        ig = IntronGraph(params, self.gene_info, read_assignments)
