import os
import pathos
from mhspf.utils.utils import Fasta
from mhspf.utils.message import Message
import itertools
import random
    

class MHS:
    def __init__(self, genome, out_dir, min_size, max_size, threads):
        self.__genome = genome
        self.__out_dir = out_dir
        self.__min_size = min_size
        self.__max_size = max_size
        self.__threads = threads
        if not os.path.exists(self.__out_dir):
            os.makedirs(self.__out_dir)

    def _sub_get_mhs(self, seq, fn):
        Message.info("\tPID: %d Starting" % os.getpid())
        MHS_db = {}
        for mhs_len in range(self.__min_size, self.__max_size + 1):
            Message.info("\tPID: %d, Getting with length: %d" % (os.getpid(), mhs_len))
            for _ in range(len(seq) - (mhs_len - 1)):
                mhs = seq[_: _ + mhs_len]
                if "N" in mhs or len(mhs) != mhs_len:
                    continue
                if mhs not in MHS_db:
                    MHS_db[mhs] = []
                MHS_db[mhs].append(_ + 1)

        Message.info("\tPID: %d, Writing MHS" % os.getpid())
        with open(fn, 'w') as fout:
            for mhs in sorted(MHS_db):
                if len(MHS_db[mhs]) < 2:
                    continue
                fout.write("%s\t%s\n" % (mhs, ','.join(map(str, MHS_db[mhs]))))
        Message.info("\tPID: %d Finished" % os.getpid())

    def get_mhs(self):
        Message.info("PID: %d, Loading genome" % os.getpid())
        fa_io = Fasta(self.__genome)
        fa_io.load()

        Message.info("PID: %d, Getting MHS" % os.getpid())
        pool = pathos.multiprocessing.Pool(processes=self.__threads)
        res = []
        for gid in fa_io.fa_db:
            fn = os.path.join(self.__out_dir, "%s.mhs" % gid)
            r = pool.apply_async(self._sub_get_mhs, (fa_io.fa_db[gid], fn,))
            res.append(r)
        pool.close()
        pool.join()

        for r in res:
            try:
                r.get()
            except Exception as e:
                Message.error("Error: {}".format(e))
        Message.info("PID: %d, Finished" % os.getpid())


class MHSTuple:
    def __init__(self, mhs_dir, freq_fn, out_dir, min_pair_dist, max_pair_dist, max_elem_dist, threads):
        self.__mhs_dir = mhs_dir
        self.__freq_fn = freq_fn
        self.__out_dir = out_dir
        self.__min_pair_dist = min_pair_dist
        self.__max_pair_dist = max_pair_dist
        self.__max_elem_dist = max_elem_dist
        self.__threads = threads
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
    
    def _get_mhs_pairs(self, mhs_file, out_file):
        Message.info("\tPID: %d, Loading most frequency MHS" % (os.getpid()))
        most_freq_mhs = []
        freq_db = {}
        with open(self.__freq_fn, 'r') as fin:
            for line in fin:
                if line[0] == '#':
                    continue
                data = line.strip().split()
                if len(set(list(data[1]))) == 1:
                    continue
                most_freq_mhs.append(data[1])
                freq_db[data[1]] = int(data[2])

        pos_db = {}
        with open(mhs_file, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                mhs = data[0]
                pos_db[mhs] = list(map(int, data[1].split(',')))

        all_mhs_pairs = [mhs_pairs for mhs_pairs in itertools.combinations(most_freq_mhs, 3)]
        random.shuffle(all_mhs_pairs)
        sel_mhs_pairs = 10 if len(all_mhs_pairs) // 10 < 10 else len(all_mhs_pairs) // 10
        sel_idx = 0
        with open(out_file, 'w') as fout:
            for mhs_pairs in all_mhs_pairs[:sel_mhs_pairs]:
                sel_idx += 1
                pmhs_pairs = list(mhs_pairs)
                random.shuffle(pmhs_pairs)
                pmhs_pairs = tuple(pmhs_pairs)
                pmhs_pairs_list = list(pmhs_pairs)
                Message.info("\tPID: %d, [%d/%d] Dealing %s with (%s)" % (
                    os.getpid(), sel_idx, sel_mhs_pairs, mhs_file, ','.join(pmhs_pairs_list)))
                tmp_mhs_pairs = []
                pos_size1 = len(pos_db[pmhs_pairs_list[0]])
                pos_start2 = 0
                pos_size2 = len(pos_db[pmhs_pairs_list[1]])
                pos_start3 = 0
                pos_size3 = len(pos_db[pmhs_pairs_list[2]])
                for idx1 in range(pos_size1):
                    pos1 = pos_db[pmhs_pairs_list[0]][idx1]
                    for idx2 in range(pos_start2, pos_size2):
                        pos2 = pos_db[pmhs_pairs_list[1]][idx2]
                        if pos2 - pos1 < len(pmhs_pairs_list[0]):
                            pos_start2 += 1
                            continue
                        if pos2 - pos1 > self.__max_elem_dist:
                            break
                        for idx3 in range(pos_start3, pos_size3):
                            pos3 = pos_db[pmhs_pairs_list[2]][idx3]
                            if pos3 - pos2 < len(pmhs_pairs_list[1]):
                                pos_start3 += 1
                                continue
                            if pos3 - pos2 > self.__max_elem_dist:
                                break
                            tmp_mhs_pairs.append([pos1, pos2, pos3])

                if len(tmp_mhs_pairs) == 0:
                    continue

                mhs_pairs = []
                mhs_pairs_regions_count = len(tmp_mhs_pairs)
                for i in range(mhs_pairs_regions_count - 1):
                    for j in range(i + 1, mhs_pairs_regions_count):
                        if tmp_mhs_pairs[j][0] - tmp_mhs_pairs[i][-1] < self.__min_pair_dist:
                            continue
                        if tmp_mhs_pairs[j][0] - tmp_mhs_pairs[i][-1] > self.__max_pair_dist:
                            break
                        mhs_pairs.append([tmp_mhs_pairs[i], tmp_mhs_pairs[j]])
                        break

                if len(mhs_pairs) == 0:
                    continue

                fout.write("# %s\n" % (','.join(pmhs_pairs_list)))
                fout.write("# %s\n" % (','.join(map(str, [freq_db[_] for _ in pmhs_pairs_list]))))
                for mhs_pairs1, mhs_pairs2 in mhs_pairs:
                    fout.write("%s\t%s\n" % (','.join(map(str, mhs_pairs1)), ','.join(map(str, mhs_pairs2))))

        Message.info("\tPID: %d, Finished" % os.getpid())
    
    def get_mhs_pairs(self):
        Message.info("PID: %d, Starting" % os.getpid())
        pool = pathos.multiprocessing.Pool(processes=self.__threads)
        res = []
        for fn in os.listdir(self.__mhs_dir):
            in_fn = os.path.join(self.__mhs_dir, fn)
            out_fn = os.path.join(self.__out_dir, "%ss" % fn)
            r = pool.apply_async(self._get_mhs_pairs, (in_fn, out_fn,))
            res.append(r)
        pool.close()
        pool.join()

        for r in res:
            try:
                r.get()
            except Exception as e:
                Message.error("Error: {}".format(e))
        Message.info("PID: %d, Finished" % os.getpid())


class MHSTable:
    def __init__(self, mhs_tuple_dir, out_table):
        self.__mhs_tuple_dir = mhs_tuple_dir
        self.__out_table = out_table

    def construct_table(self):
        freq_db = {}
        mhs_pairs_set = set()
        for fn in os.listdir(self.__mhs_tuple_dir):
            in_fn = os.path.join(self.__mhs_tuple_dir, fn)
            with open(in_fn, 'r') as fin:
                for line in fin:
                    if line[0] == '#':
                        data = line.strip().split()
                        vals = data[1].split(',')
                        if not vals[0].isdigit():
                            mhs_pairs = vals
                            mhs_pairs_set.add(tuple(mhs_pairs))
                        else:
                            freq = list(map(int, vals))
                            for idx in range(3):
                                mhs = mhs_pairs[idx]
                                if mhs not in freq_db:
                                    freq_db[mhs] = freq[idx]

        with open(self.__out_table, 'w') as fout:
            fout.write("ID\tMHS1\tMHS2\tMHS3\n")
            idx = 1
            for mhs_pairs in mhs_pairs_set:
                tmp = []
                for mhs in mhs_pairs:
                    tmp.append("%s,%s" % (mhs, freq_db[mhs]))
                fout.write("%d\t%s\n" % (idx, '\t'.join(tmp)))
                idx += 1
