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

import sys, random, pysam, array
from .hg19_rd import Hg19_Reader
from .hg19_portion import Hg19_Portion
from .ad_miner import mineSamFilePos, AD_PortionMiner
from .ad_lib import AD_LibReader

#========================================
def selectSamFiles(lib_rd, sam_file):
    sam_list = []
    with open(sam_file, "r") as inp:
        for line in inp:
            sam_list.append(line.strip())
    samfiles = []
    for smp_name in lib_rd.iterSampleNames():
        sam_fnames = []
        for fname in sam_list:
            if fname.find(smp_name) >=0:
                sam_fnames.append(fname)
        print("For", smp_name, "we have", sam_fnames)
        assert len(sam_fnames) == 1
        samfiles.append(pysam.AlignmentFile(sam_fnames[0]))
    return samfiles


#========================================
def toArray(vv):
    return array.array('H',
        [int(vv[0][0]), int(vv[0][1]), int(vv[1][0]), int(vv[1][1])])


#========================================
# Check index library vs BAM-files
#> python3 -m adlib.check_lib random <input.ldx> /data/hg19.fasta <bam-file-list>
#
#> python3 -m adlib.check_lib portion <input.ldx> <hg19-portion> <bam-file-list> <smp-idx>
#
#========================================
if sys.argv[1] not in ("random", "portion"):
    print ("Mode should be random/portion:", sys.argv[1], file = sys.stderr)
    sys.exit()

if sys.argv[1] == "random":
    lib_rd = AD_LibReader(sys.argv[2])
    hg19_rd = Hg19_Reader(sys.argv[3], chrom_is_int = True, upper_case = True)
    samfiles = selectSamFiles(lib_rd, sys.argv[4])
    rH = random.Random(179)

    cnt = 0
    cnt_bad = 0
    while True:
        portion_info = lib_rd._nextPortions()
        if portion_info is None:
            break
        cnt += 1
        if cnt % 100 == 0:
            print("Checked %d portions" % cnt,  file = sys.stderr)
            sys.stderr.flush()
        lib_portions = lib_rd._getCurPortions()
        if any([portion.isComplex() for portion in lib_portions]):
            test_n = 50
            print("Complex portion", cnt, file = sys.stderr)
        else:
            test_n = 20
        pos_seq = sorted([lib_portions[0].getShift() +
                rH.randint(0, lib_portions[0].getSize() - 1)
                for idx in range(test_n)])
        pos_idx = 0
        while pos_idx < len(pos_seq):
            if hg19_rd.getCurChrom() != lib_portions[0].getChrom():
                hg19_rd.read()
                continue
            diap = hg19_rd.getCurDiap()
            if not (diap[0] <= pos_seq[pos_idx] < diap[1]):
                hg19_rd.read()
                continue
            ref_letter = hg19_rd.getCurLetters()[pos_seq[pos_idx] - diap[0]]
            sam_ad_seq = [mineSamFilePos(samfile, hg19_rd.getCurChrom(),
                pos_seq[pos_idx], ref_letter)
                for samfile in samfiles]
            lib_ad_seq = []
            for portion in lib_portions:
                lib_ad_seq.append(toArray(portion.getAD(pos_seq[pos_idx])))
            if sam_ad_seq != lib_ad_seq:
                print("Mismatch in chrom=%d pos=%d" %
                    (hg19_rd.getCurChrom(), pos_seq[pos_idx]))
                print("\tSAM:", sam_ad_seq)
                print("\tLDX:", lib_ad_seq)
                cnt_bad += 1
                sys.stdout.flush()
            pos_idx += 1

    print("Done:", cnt, "bad:", cnt_bad, file = sys.stderr)
    sys.exit()

assert sys.argv[1] == "portion"
lib_rd = AD_LibReader(sys.argv[2])
hg19_portion = Hg19_Portion.load(sys.argv[3])
samfiles = selectSamFiles(lib_rd, sys.argv[4])
idx_idx = int(sys.argv[5])
sam_file = samfiles[idx_idx]
miner = AD_PortionMiner(hg19_portion, sam_file)
print("Started", file = sys.stderr)
for pos in range(*hg19_portion.getDiap()):
    ldx_ad = toArray(lib_rd.getAD_seq(hg19_portion.getChrom(), pos)[idx_idx])
    mine_ad = array.array('H', miner.getAD(pos))
    sam_ad = mineSamFilePos(sam_file, hg19_portion.getChrom(), pos,
        hg19_portion.getLetter(pos))
    print("\t".join([str(pos),
        "ldx " + " ".join(map(str, ldx_ad)),
        "sam " + " ".join(map(str, sam_ad)),
        "mine " + " ".join(map(str, mine_ad))]))


