import collections
from collections import defaultdict

import gffutils
from pyfaidx import Fasta

from src.dataset_processor import ReadAssignmentLoader
from src.gene_info import GeneInfo
from src.intron_graph import IntronGraph
from src.isoform_assignment import BasicReadAssignment
from src.serialization import read_int, TERMINATION_INT


class TestIntronGraph:
    def get_read_assignments(self):
        #  ./isoquant.py --keep_tmp --reference tests/simple_data/chr9.4M.fa.gz  --genedb tests/simple_data/chr9.4M.gtf.gz --fastq tests/simple_data/chr9.4M.ont.sim.fq.gz  --data_type nanopore -p test_data -o temp
        dump_filename = '/home/natasha/PycharmProjects/IsoQuant/temp/test_data/aux/test_data.save'
        chr_id = 'chr9'
        gffutils_db = gffutils.FeatureDB('/home/natasha/PycharmProjects/IsoQuant/temp/chr9.4M.gtf.db', keep_order=True)
        chr_dump_file = dump_filename + "_" + chr_id  # dump_filename - путь к save файлам XXX_save
        current_chr_record = Fasta('/home/natasha/PycharmProjects/IsoQuant/tests/simple_data/chr9.4M.fa.gz')[chr_id]

        #(venv) $ ./isoquant.py --keep_tmp --reference tests/toy_data/MAPT.Mouse.reference.fasta  --genedb tests/toy_data/MAPT.Mouse.genedb.gtf --fastq tests/toy_data/MAPT.Mouse.ONT.simulated.fastq  --data_type nanopore -o test_data_toy
        # dump_filename = '/home/natasha/PycharmProjects/IsoQuant/temp/test_data/aux/test_data.save'

        multimapped_reads = defaultdict(list)
        multimappers_filename = dump_filename + "_multimappers_" + chr_id
        multimap_loader = open(multimappers_filename, "rb")
        list_size = read_int(multimap_loader)
        while list_size != TERMINATION_INT:
            for i in range(list_size):
                a = BasicReadAssignment.deserialize(multimap_loader)
                if a.chr_id == chr_id:
                    multimapped_reads[a.read_id].append(a)
            list_size = read_int(multimap_loader)

        # reading a 'genome chunk' / 'chunk' of a genome (possibly several genes) (~region = any consecutive sub-genome)
        loader = ReadAssignmentLoader(chr_dump_file, gffutils_db, current_chr_record, multimapped_reads)
        while loader.has_next():
            self.gene_info, self.assignment_storage = loader.get_next()
            break

    def test_add_edge(self):
        self.get_read_assignments()

        ModelConstructionStrategy = collections.namedtuple('ModelConstructionStrategy',
                                                           ('min_novel_intron_count',
                                                            'graph_clustering_ratio', 'graph_clustering_distance',
                                                            'min_novel_isolated_intron_abs',
                                                            'min_novel_isolated_intron_rel',
                                                            'terminal_position_abs', 'terminal_position_rel',
                                                            'terminal_internal_position_rel',
                                                            'min_known_count', 'min_nonfl_count',
                                                            'min_novel_count', 'min_novel_count_rel',
                                                            'min_mono_count_rel', 'singleton_adjacent_cov',
                                                            'fl_only', 'novel_monoexonic',
                                                            'delta', 'apa_delta', 'debug'))
        strategy = ModelConstructionStrategy(1, 0.5, 20, 3, 0.02, 1, 0.05, 0.05, 1, 3, 3, 0.02, 0.02, 10, False,
                                             False, 6, 50, False)  # 'default_ont' + delta + apa_delta + debug
        ig = IntronGraph(strategy, self.gene_info, self.assignment_storage)

        # take the first item in assignment storage, then dublicate 10 times into a new ass storage

        ig.print_graph()  # define logger
