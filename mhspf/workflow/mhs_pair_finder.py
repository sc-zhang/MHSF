import argparse
import os
from mhspf.utils.MHS import MHS, MHSTuple, MHSTable
from mhspf.utils.utils import MHSIO
from mhspf.utils.message import Message


def get_opts():
    group = argparse.ArgumentParser()
    group.add_argument('-g', '--genome', help="Input genome file", required=True)
    group.add_argument('--min_size', help="Minium size of MHS, default=5", type=int, default=5)
    group.add_argument('--max_size', help="Maximum size of MHS, default=7", type=int, default=7)
    group.add_argument('-c', '--count', help="Count of most frequency MHS for each length, "
                                             "default=20", type=int, default=20)
    group.add_argument('--min_pair_dist', help="Minimum distance between two MHSs tuples, default=500",
                       type=int, default=500)
    group.add_argument('--max_pair_dist', help="Maximum distance between two MHSs tuples, default=3000",
                       type=int, default=3000)
    group.add_argument('--max_element_dist', help="Maximum distance between two MHSs, default=50",
                       type=int, default=50)
    group.add_argument('-o', '--output', help="Output directory", required=True)
    group.add_argument('-t', '--threads', help="Threads count, default=10", type=int, default=10)
    return group.parse_args()


def main():
    opts = get_opts()
    genome = opts.genome
    min_element_size = opts.min_size
    max_element_size = opts.max_size
    most_freq_mhs_count = opts.count
    min_pair_dist = opts.min_pair_dist
    max_pair_dist = opts.max_pair_dist
    max_element_dist = opts.max_element_dist
    out_dir = opts.output
    threads = opts.threads

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    Message.info("Getting micro homologous sequences")
    mhs_dir = os.path.join(out_dir, "01.mhs")
    if not os.path.exists(mhs_dir):
        mhs_finder = MHS(genome, mhs_dir, min_element_size, max_element_size, threads)
        mhs_finder.get_mhs()
    else:
        Message.info("01.mhs found, skipping...")

    Message.info("Getting most frequent MHS sequence")
    freq_fn = os.path.join(out_dir, "02.most_freq_mhs.txt")
    if not os.path.exists(freq_fn):
        mhs_io = MHSIO(mhs_dir, most_freq_mhs_count, freq_fn)
        mhs_io.get_most_freq_mhs()
    else:
        Message.info("02.most_freq_mhs found, skipping...")

    Message.info("Getting MHS tuples")
    mhs_tuple_dir = os.path.join(out_dir, "03.mhs_tuples")
    if not os.path.exists(mhs_tuple_dir):
        mhs_tuple_finder = MHSTuple(mhs_dir, freq_fn, mhs_tuple_dir,
                                    min_pair_dist, max_pair_dist,
                                    max_element_dist, threads)
        mhs_tuple_finder.get_mhs_pairs()
    else:
        Message.info("03.mhs_tuples found, skipping...")

    Message.info("Constructing 3-tuple MHS pairs table")
    out_table = os.path.join(out_dir, "mhs_pairs.tsv")
    if not os.path.exists(out_table):
        mhs_table_constructor = MHSTable(mhs_tuple_dir, out_table)
        mhs_table_constructor.construct_table()
    else:
        Message.info("mhs_pairs.tsv found, skipping...")

    Message.info("Finished")
