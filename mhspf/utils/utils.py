import os

class Fasta:
    def __init__(self, fasta_file):
        self.__fasta_file = fasta_file
        self.fa_db = {}

    def load(self):
        with open(self.__fasta_file, 'r') as fin:
            for line in fin:
                if line[0] == '>':
                    gid = line.strip().split()[0][1:]
                    self.fa_db[gid] = []
                else:
                    self.fa_db[gid].append(line.strip().upper())

        for gid in self.fa_db:
            self.fa_db[gid] = ''.join(self.fa_db[gid])


class MHSIO:
    def __init__(self, mhs_dir, cnt, freq_fn):
        self.__mhs_dir = mhs_dir
        self.__cnt = cnt
        self.__freq_fn = freq_fn
        
    def get_most_freq_mhs(self):
        freq_db = {}
        for fn in os.listdir(self.__mhs_dir):
            full_fn = os.path.join(self.__mhs_dir, fn)
            with open(full_fn, 'r') as fin:
                for line in fin:
                    data = line.strip().split()
                    MHS = data[0]
                    mhs_len = len(data[0])
                    mhs_cnt = len(data[1].split(','))
                    if mhs_len not in freq_db:
                        freq_db[mhs_len] = {}
                    if MHS not in freq_db[mhs_len]:
                        freq_db[mhs_len][MHS] = 0
                    freq_db[mhs_len][MHS] += mhs_cnt

        with open(self.__freq_fn, 'w') as fout:
            fout.write("#MHS_len\tMHS_seq\tMHS_cnt\n")
            for mhs_len in sorted(freq_db):
                cur_cnt = 0
                for MHS in sorted(freq_db[mhs_len], key=lambda x: freq_db[mhs_len][x], reverse=True):
                    if len(set(list(MHS))) == 1:
                        continue
                    cur_cnt += 1
                    fout.write("%d\t%s\t%d\n" % (mhs_len, MHS, freq_db[mhs_len][MHS]))
                    if cur_cnt > self.__cnt:
                        break
