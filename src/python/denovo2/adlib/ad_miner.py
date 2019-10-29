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

import pysam
import array
from .ad_person import AD_Portion, AD_PersonData

#========================================
# Logic of pysam
#========================================
def _pileup(samfile, chrom, pos_from, pos_to):
    sam_chrom = str(chrom) if 0 < chrom <= 22 else {0: "M", 23: "X", 24: "Y"}[chrom]
    try:
        return samfile.pileup("chr" + sam_chrom, pos_from, pos_to)
    except ValueError:
        pass
    if chrom == 0:
        sam_chrom = "MT"
    return samfile.pileup(sam_chrom, pos_from, pos_to)

#========================================
MQ_thresh = -100.
BQ_thresh = -100.

def evalPileUpColumn(pileupcolumn, ref_letter, counts, pos_idx):
    global MQ_thresh, BQ_thresh
    for pileupread in pileupcolumn.pileups:
        if pileupread.is_del or pileupread.is_refskip:
            continue
        q_pos = pileupread.query_position
        MQ = pileupread.alignment.mapping_quality
        BQ = pileupread.alignment.query_qualities[q_pos]
        if MQ < MQ_thresh or BQ <BQ_thresh:
            continue
        is_ref =  (ref_letter ==
            pileupread.alignment.query_sequence[q_pos].upper())
        idx_rev = 1 if pileupread.alignment.is_reverse else 0
        counts[(pos_idx << 2) + (idx_rev << 1) + (0 if is_ref else 1)] += 1

def mineSamFilePos(samfile, chrom, pos, ref_letter):
    counts = array.array('H', [0] * 4)
    for pileupcolumn in _pileup(samfile, chrom, pos - 1, pos):
        if pos - 1 == pileupcolumn.pos:
            evalPileUpColumn(pileupcolumn, ref_letter, counts, 0)
    return counts

#========================================
# Miner: mines data using pysam 
#   and serializes result in array form
#========================================
class AD_PortionMiner:
    PILE_BLOCK = 100000

    def __init__(self, hg19_portion, samfile):
        pos_shift = hg19_portion.getShift()
        size = hg19_portion.getSize()
        
        self.mRefHg19 = hg19_portion
        self.mCounts = array.array('H', [0] * (size << 2))
        self.mReport = None


        for pos_min in range(pos_shift, pos_shift + size, self.PILE_BLOCK):
            pos_max = min(pos_min + self.PILE_BLOCK, pos_shift + size)
            for pileupcolumn in _pileup(samfile,
                    self.mRefHg19.getChrom(), pos_min - 1, pos_max - 1):
                if pos_min - 1 <= pileupcolumn.pos < pos_max - 1:
                    pos_idx = pileupcolumn.pos + 1 - pos_shift
                    evalPileUpColumn(pileupcolumn,
                        self.mRefHg19.getLetter(pos_shift + pos_idx),
                        self.mCounts, pos_idx)

    def getChrom(self):
        return self.mRefHg19.getChrom()

    def getShift(self):
        return self.mRefHg19.getShift()

    def getSize(self):
        return self.mRefHg19.getSize()

    def getAD(self, pos):
        idx0 = (pos - self.mRefHg19.getShift()) << 2
        return [self.mCounts[idx] for idx in range(idx0, idx0 + 4)]

    def report(self):
        return self.mReport

    def debugReport(self):
        for idx0 in range(0, len(self.mCounts), 4):
            print("\t".join(["Q"] +
                map(str, [self.mCounts[idx] for idx in range(idx0, idx0+4)])))

    def _flushData(self, prefix, head, pos_ref, ad_tab, bin_output):
        bin_output.write(prefix)
        head.tofile(bin_output)
        pos_ref.tofile(bin_output)
        ad_tab.tofile(bin_output)

    def toFile(self, bin_output):
        head = array.array('L', [self.mRefHg19.getChrom(), 
            self.mRefHg19.getShift(), self.mRefHg19.getSize(), 0])
        pos_ref = array.array('H')
        ad_tab = array.array('H')
        ad_dict = {(0, 0, 0, 0): 0}
        self.mReport = []
        for no in range(self.mRefHg19.getSize()):
            idx0 = no << 2
            pos_ad = tuple([self.mCounts[idx] 
                for idx in range(idx0, idx0 + 4)])
            if pos_ad not in ad_dict:
                if len(ad_dict) >= 64000:
                    head[3] = len(ad_dict)
                    self._flushData(AD_Portion.BLOCK_PRE, 
                        array.array('L',
                            [head[0], head[1], len(pos_ref), len(ad_dict)]),
                        pos_ref, ad_tab, bin_output)
                    head[1] += len(pos_ref)
                    head[2] -= len(pos_ref)
                    self.mReport += ["added", len(pos_ref), len(ad_dict)]
                    pos_ref = array.array('H')
                    ad_tab = array.array('H')
                    ad_dict = {(0, 0, 0, 0): 0}
                ad_dict[pos_ad] = len(ad_dict)
                ad_tab.extend(pos_ad)
            pos_ref.append(ad_dict[pos_ad])

        head[3] = len(ad_dict)
        assert head[3] < 64000
        self._flushData(AD_Portion.BLOCK_BASE, 
            head, pos_ref, ad_tab, bin_output)
        self.mReport += [self.getChrom(), self.getShift(), 
            self.getSize(), len(ad_dict)]

#========================================
class AD_PersonDataWriter:
    def __init__(self, fname):
        self.mOutput = open(fname, 'wb')
        self.mOutput.write(AD_PersonData.PREFIX)
        self.mRootPos = self.mOutput.tell()
        array.array('L', [0, 0]).tofile(self.mOutput)
        self.mTab = array.array('L')

    def addPortion(self, portion):
        self.mTab.extend([portion.getChrom(), portion.getShift(), 
            portion.getSize(), self.mOutput.tell()])
        portion.toFile(self.mOutput)

    def close(self):
        root_array = array.array('L', [self.mOutput.tell(), len(self.mTab) >> 2])
        self.mTab.tofile(self.mOutput)
        self.mOutput.seek(self.mRootPos)
        root_array.tofile(self.mOutput)
        self.mOutput.close()
        self.mOutput = None

#========================================
# __main__ is using for deep development
# Modes:
# create: create partial AD-data file for some positions and retrive AD:
# > python3 ad_miner.py create <input-sample.bam> <output.idx> <hg19-portions>
#
# read: retrive AD-data for some positions (if ad-data exists):
#  >python3 ad_miner.py read <input.idx>
#========================================
if __name__=="__main__":
    import sys, time
    from .hg19_portion import Hg19_Portion
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


    mode = sys.argv[1]
    if mode == "create":
        samfile = pysam.AlignmentFile(sys.argv[2], "rb")
        ad_seq_file = sys.argv[3]
        print("Create ad-data...", ad_seq_file,  file = sys.stderr)
        writer = AD_PersonDataWriter(ad_seq_file)
        tm0 = time.time()
        tm1 = None
        for fname in sys.argv[4:]:
            print("work", fname, file = sys.stderr)
            portion = Hg19_Portion.load(fname)
            ad_portion = AD_PortionMiner(portion, samfile)
            writer.addPortion(ad_portion)
            if tm1 is None:
                tm1 = time.time()
        writer.close()
        tm2 = time.time()
        print("Done at", tm1 - tm0, tm2 - tm0, file = sys.stderr)
    else:
        assert mode == "read"
        ad_seq_file = sys.argv[2]

    print("Reading", ad_seq_file, file = sys.stderr)
    reader = AD_PersonData(ad_seq_file)
    for chrom, pos, ref in sorted(test_data):
        print ("\t".join(["AD",str(chrom), str(pos), ref,
            str(reader.getAD(chrom, pos))]))
    reader.close()
