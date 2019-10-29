import pysam, json
from glob import glob
import numpy as np

from ..adlib.ad_lib import AD_LibReader
#========================================
class PysamList:
    def __init__(self, list_of_filenames):
        self.mSamFiles = []
        with open(list_of_filenames, "r") as inp:
            for line in inp:
                filename = line.partition('#')[0].strip()
                if not filename:
                    continue
                samfile = pysam.AlignmentFile(filename, "rb")
                print("Load pysam file:", filename, "\n",
                    samfile.check_index())
                self.mSamFiles.append(samfile)

    def mineAD(self, variant):
        if (len(variant.getRef()), len(variant.getAlt())) != (1, 1):
            return None, None
        ADfs, ADrs = [], []
        for samfile in self.mSamFiles:
            ADf, ADr = mineAD_ord(samfile, variant)
            ADfs.append(ADf)
            ADrs.append(ADr)
        return np.array(ADfs), np.array(ADrs)

#========================================
def _pileup(samfile, chrom, pos_from, pos_to):
    sam_chrom = (str(chrom) if 0 < chrom <= 22
        else {0: "M", 23: "X", 24: "Y"}[chrom])
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

#========================================
def mineAD_ord(samfile, variant):
    global MQ_thresh, BQ_thresh
    position = variant.getPos() - 1
    ADf, ADr = np.array([0.,0]), np.array([0., 0.])
    for pileupcolumn in _pileup(samfile, variant.getChromNum(),
            position, position + 1):
        if pileupcolumn.pos != position:
            continue
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue
            q_pos = pileupread.query_position
            MQ = pileupread.alignment.mapping_quality
            BQ = pileupread.alignment.query_qualities[q_pos]
            if MQ < MQ_thresh or BQ <BQ_thresh:
                continue
            if (variant.getRef().upper() ==
                    pileupread.alignment.query_sequence[q_pos].upper()):
                if pileupread.alignment.is_reverse:
                    ADr[0] += 1
                else:
                    ADf[0] += 1
            else:
                if pileupread.alignment.is_reverse:
                    ADr[1] += 1
                else:
                    ADf[1] += 1
    return ADf,ADr

#========================================
class AD_LibCollection:
    def __init__(self, lib_dir, dump_file):
        self.mLibSeq = []
        for fname in sorted(list(glob(lib_dir + "/*.ldx"))):
            self.mLibSeq.append(AD_LibReader(fname))
        self.mDumpFile = dump_file
        self.mDumpDict = dict()

    def mineAD(self, variant):
        if (len(variant.getRef()), len(variant.getAlt())) != (1, 1):
            if self.mDumpFile:
                key = "%d/%d" % (variant.getChromNum(), variant.getPos())
                if key not in self.mDumpDict:
                    self.mDumpDict[key] = [[[0., 0.]], [[0., 0.]]]
            return None, None
        ADfs, ADrs = [], []
        for lib in self.mLibSeq:
            seq = lib.getAD_seq(variant.getChromNum(), variant.getPos())
            if seq:
                for fam_vec in seq:
                    ADfs.append(fam_vec[0])
                    ADrs.append(fam_vec[1])
        if self.mDumpFile:
            key = "%d/%d" % (variant.getChromNum(), variant.getPos())
            if key not in self.mDumpDict:
                self.mDumpDict[key] = [
                    [[vec[0], vec[1]] for vec in ADfs],
                    [[vec[0], vec[1]] for vec in ADrs]]
        return np.array(ADfs), np.array(ADrs)

    def finishUp(self):
        if self.mDumpFile:
            with open(self.mDumpFile, "w") as outp:
                outp.write(json.dumps(self.mDumpDict,
                    indent = 4, sort_keys = True))

#========================================



