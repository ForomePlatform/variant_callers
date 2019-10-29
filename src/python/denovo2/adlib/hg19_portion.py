#  Copyright (c) 2019. Sergey Trifonov, Michael Bouzinier
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

#========================================
# Autonomous portion of hg19
# For debug purposed can be serialized
#========================================
class Hg19_Portion:
    PREFIX = "#Hg19_Block"
    def __init__(self, chrom, diap, letters):
        self.mChrom = chrom
        self.mDiap = diap
        self.mLetters = letters

    def getChrom(self):
        return self.mChrom

    def getDiap(self):
        return self.mDiap

    def getShift(self):
        return self.mDiap[0]

    def getLetters(self):
        return self.mLetters

    def getSize(self):
        return self.mDiap[1] - self.mDiap[0]

    def getLetter(self, pos):
        assert self.mDiap[0] <= pos < self.mDiap[1]
        return self.mLetters[pos - self.mDiap[0]]

    def save(self, fname):
        with open(fname, "w") as outp:
            print(self.PREFIX, file = outp)
            print("\t".join(map(str, [self.mChrom] + self.mDiap)), file = outp)
            shift, total = 0, self.mDiap[1] - self.mDiap[0]
            while total > 0:
                size = min(total, 50)
                print(self.mLetters[shift:shift + size], file = outp)
                shift += size
                total -= size

    @classmethod
    def load(cls, fname):
        with open(fname, "r") as inp:
            line = inp.readline()
            assert line.rstrip() == cls.PREFIX
            chrom, pos0, pos1 = map(int, inp.readline().rstrip().split())
            lines = []
            total = pos1 - pos0
            while total > 0:
                letters = inp.readline().rstrip()
                lines.append(letters)
                total -= len(letters)
            assert total == 0
            assert not inp.readline()
        return cls(chrom, [pos0, pos1], ''.join(lines))

#========================================
# __main__ forms searialization of some portions
# W/o locations work for predefined ones
# Format of location: <chrom>/<pos>, example: 8/46848652
#========================================
if __name__=="__main__":
    import sys
    from .hg19_rd import Hg19_Reader
    test_data = [
        [1, 228345458,  "CTATGA"],
        [1, 142569315,  "T"],
        [1, 152188319,  "C"],
        [1, 248605314,  "C"],
        [1, 248605316,  "C"],
        [1, 142569310,  "T"],
        [1, 248605354,  "C"],
        [1, 248605007,  "ACTT"],
        [1, 248604978,  "C"],
        [1, 144823897,  "T"],
        [1, 146696617,  "G"],
        [1, 146696610,  "TTTC"],
        [1, 146055394,  "A"],
        [1, 146696607,  "TC"],
        [1, 143232043,  "A"]]
    check_mode = True

    rd = Hg19_Reader(sys.argv[1], chrom_is_int = True, upper_case = True)

    if len(sys.argv) > 2:
        test_data = [list(map(int, arg.split('/')))
            for arg in sys.argv[2:]]
        check_mode = False

    min_chrom = min([dt[0] for dt in test_data])
    seq_pos = []
    for dt in test_data:
        if dt[0] == min_chrom:
            seq_pos.append(dt[1])
    min_pos = min(seq_pos)

    cnt = 0
    while rd.read():
        cnt += 1
        if cnt % 100 == 0:
            print (cnt, "portions...", file = sys.stderr)
        if rd.getCurChrom() != min_chrom:
            continue
        diap = rd.getCurDiap()
        if not (diap[0] <= min_pos < diap[1]):
            continue
        print ("Portion:", rd.getCurChrom(), diap, file = sys.stderr)
        hg19_block = Hg19_Portion(rd.getCurChrom(), diap, rd.getCurLetters())
        hg19_block.save("hg19_%d_%d.part" % (rd.getCurChrom(), diap[0]))
        for idx in range(len(test_data) -1, -1, -1):
            dt = test_data[idx]
            if dt[0] != min_chrom or not (diap[0] <= dt[1] < diap[1]):
                continue
            if check_mode:
                print ("\t".join(["P", str(min_chrom), str(dt[1]), dt[2],
                    hg19_block.getLetter(dt[1])]))
            del test_data[idx]
        if len(test_data) == 0:
            continue
        min_chrom = min([td[0] for td in test_data])
        seq_pos = []
        for dt in test_data:
            if dt[0] == min_chrom:
                seq_pos.append(dt[1])
        min_pos = min(seq_pos)

    rd.close()
    print("Done", file = sys.stderr)
