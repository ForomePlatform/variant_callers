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

import array
import numpy as np
import traceback

#========================================
# Portion of AD-data, contains one or two chunks
#========================================
class AD_Portion:
    PREFIX_LEN = 7
    BLOCK_PRE  = b"#Blk-0\n"
    BLOCK_BASE = b"#Block\n"

    def __init__(self, bin_input):
        self.mChunks = []
        title = bin_input.read(self.PREFIX_LEN)
        while title == self.BLOCK_PRE:
            self.mChunks.append(AD_PortionChunk(title, bin_input))
            title = bin_input.read(self.PREFIX_LEN)
        assert title == self.BLOCK_BASE
        self.mChunks.append(AD_PortionChunk(title, bin_input))
        self.mSize = sum([chunk.getSize() for chunk in self.mChunks])

    def getChrom(self):
        return self.mChunks[0].getChrom()

    def getShift(self):
        return self.mChunks[0].getShift()

    def getSize(self):
        return self.mSize

    def getInfo(self):
        return (self.mChunks[0].getChrom(), self.mChunks[0].getShift(),
            self.mChunks[0].getShift() + self.mSize)

    def isOf(self, chrom, pos):
        if self.getChrom() == chrom:
            for chunk in self.mChunks:
                if chunk.isOf(pos):
                    return True
        return False

    def isComplex(self):
        return len(self.mChunks) > 1

    def getAD(self, pos):
        for chunk in self.mChunks:
            if chunk.isOf(pos):
                return chunk.getAD(pos)
        assert False

    def toFile(self, bin_output):
        for chunk in self.mChunks:
            chunk.toFile(bin_output)

#========================================
class AD_PortionChunk:
    sZeroAD = np.array([[0., 0.], [0., 0.]])

    def __init__(self, prefix, bin_input):
        self.mPrefix = prefix
        self.mHead = array.array('L')
        self.mHead.fromfile(bin_input, 4)
        self.mChrom = self.mHead[0]
        self.mShift = self.mHead[1]
        self.mSize  = self.mHead[2]
        self.mPosRef = array.array('H')
        self.mPosRef.fromfile(bin_input, self.mSize)
        self.mAD_Tab = array.array('H')
        self.mAD_Tab.fromfile(bin_input, 4 * (self.mHead[3] - 1))

    def getChrom(self):
        return self.mChrom

    def getShift(self):
        return self.mShift

    def getSize(self):
        return self.mSize

    def isOf(self, pos):
        return self.mShift <= pos < (self.mShift + self.mSize)

    def getAD(self, pos):
        ref_idx = self.mPosRef[pos - self.mShift]
        if ref_idx == 0:
            return self.sZeroAD
        r_idx = (ref_idx - 1) << 2
        vv = [self.mAD_Tab[idx] for idx in range(r_idx, r_idx + 4)]
        return np.array([[vv[0], vv[1]], [vv[2], vv[3]]], float)

    def toFile(self, bin_output):
        bin_output.write(self.mPrefix)
        self.mHead.tofile(bin_output)
        self.mPosRef.tofile(bin_output)
        self.mAD_Tab.tofile(bin_output)

#========================================
# Reader of the person AD-data file
#========================================
class AD_PersonData:
    PREFIX = b"#SeqBlockAD\n"

    def __init__(self, fname, mine_mode = False):
        self.mFName = fname
        self.mInput = open(fname, 'rb')
        title =  self.mInput.read(12)
        assert title == self.PREFIX
        root_array = array.array('L')
        root_array.fromfile(self.mInput, 2)
        if mine_mode:
            if root_array[0] == 0:
                self.mTab = None
                return
        assert root_array[0] > 0
        self.mInput.seek(root_array[0])
        self.mTab = array.array('L')
        self.mTab.fromfile(self.mInput, 4 * root_array[1])
        self.mCurTabIdx = None
        self.mCurPortion = None

    def getFName(self):
        return self.mFName

    def isComplete(self):
        return self.mTab is not None

    def close(self):
        self.mInput.close()
        self.mInput = None

    def _setCurTab(self, idx):
        self.mCurTabIdx = idx
        self.mInput.seek(self.mTab[self.mCurTabIdx + 3])
        self.mCurPortion = AD_Portion(self.mInput)

    def directReadPortion(self):
        if self.mTab is not None:
            if self.mCurTabIdx is None:
                self.mCurTabIdx = 0
            else:
                self.mCurTabIdx += 4
            if self.mCurTabIdx >= len(self.mTab):
                return None
            try:
                self._setCurTab(self.mCurTabIdx)
                return self.mCurPortion
            except:
                traceback.print_exc()
                return False
        try:
            return AD_Portion(self.mInput)
        except:
            traceback.print_exc()
            return None

    def getAD(self, chrom, pos):
        assert self.mTab is not None
        if self.mCurPortion is not None and self.mCurPortion.isOf(chrom, pos):
            return self.mCurPortion.getAD(pos)
        self.mCurPortion, self.mCurTabIdx = None, None
        for idx0 in range(0, len(self.mTab), 4):
            if chrom != self.mTab[idx0]:
                continue
            if not (self.mTab[idx0 + 1] <= pos <
                    self.mTab[idx0 + 1] + self.mTab[idx0 + 2]):
                continue
            self._setCurTab(idx0)
            return self.mCurPortion.getAD(pos)
        return None

#========================================
if __name__=="__main__":
    import sys
    print("Check ad-data", sys.argv[1])
    rd = AD_PersonData(sys.argv[1], True)
    print("Complete ad-data:", rd.isComplete())
    cnt = 0
    cnt_bad = 0
    while True:
        portion = rd.directReadPortion()
        cnt += 1
        if portion is False:
            print("Failed portion:", cnt)
            cnt_bad += 1
            continue
        if portion is None:
            break
        if cnt % 100 == 0:
            print("Portions", cnt, "...", file = sys.stderr)
    print("Done", cnt, "bad=", cnt_bad)

