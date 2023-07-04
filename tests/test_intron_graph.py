import collections
import os
from collections import defaultdict

import gffutils

from src.gene_info import GeneInfo
from src.intron_graph import IntronGraph
from src.isoform_assignment import BasicReadAssignment
from src.serialization import read_int, TERMINATION_INT


class TestIntronGraph:
    source_dir = os.path.dirname(os.path.realpath(__file__))
    gffutils_db = gffutils.FeatureDB(os.path.join(source_dir, 'toy_data/synth.db'), keep_order=True)
    gene_db = gffutils_db['ENSMUSG00000020196.10']

    def get_read_assignments(self):
        dump_filename = '/home/natasha/PycharmProjects/IsoQuant/temp/OUT/aux/OUT.save'
        chr_id = 11

        multimapped_reads = defaultdict(list)

        with open(dump_filename + "_multimappers_chr" + str(chr_id), "rb") as multimap_loader:
            list_size = read_int(multimap_loader)
            while list_size != TERMINATION_INT:
                for i in range(list_size):
                    a = BasicReadAssignment.deserialize(multimap_loader)
                    if a.chr_id == chr_id:
                        multimapped_reads[a.read_id].append(a)
        return multimapped_reads.values()

    def test_add_edge(self):
        self.gene_info = GeneInfo([self.gene_db], self.gffutils_db)
        read_assignments = self.get_read_assignments()
        Params = collections.namedtuple('Params', 'delta, min_novel_intron_count, debug, min_novel_isolated_intron_abs')
        params = Params(1, 1, False, 1)
        ig = IntronGraph(params, self.gene_info, read_assignments)
