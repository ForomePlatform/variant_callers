#  Copyright (c) 2019. Partners HealthCare, Harvard Medical School’s
#  Department of Biomedical Informatics, Sergey Trifonov
#
#  Developed by Sergey Trifonov and Michael Bouzinier, based on contributions by:
#  Anwoy Kumar Mohanty, Andrew Bjonnes,
#  Ignat Leshchiner, Shamil Sunyaev and other members of Division of Genetics,
#  Brigham and Women's Hospital
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


import array, bz2
from io import BytesIO
from .ad_person import AD_Portion

#========================================
class AD_LibReader:
    PREFIX = b"#LibBlockAD\n"
    def __init__(self, fname):
        self.mInput = open(fname, 'rb')
        title =  self.mInput.read(12)
        assert title == self.PREFIX
        self.mSamples = []
        for line in self.mInput:
            sample_name = line.decode().rstrip()
            if sample_name:
                self.mSamples.append(sample_name)
            else:
                break
        root_array = array.array('Q')
        root_array.fromfile(self.mInput, 2)
        self.mInput.seek(root_array[0])
        self.mTab = array.array('Q')
        self.mTab.fromfile(self.mInput, root_array[1])
        self.mCurTabIdx = None
        self.mCurPortions = None

    def iterSampleNames(self):
        return iter(self.mSamples)

    def close(self):
        self.mInput.close()
        self.mInput = None

    def _nextPortions(self):
        if self.mCurTabIdx is None:
            self.mCurTabIdx = 0
        else:
            self.mCurTabIdx += 5
            if self.mCurTabIdx >= len(self.mTab):
                return None
        self._setupPortions()
        return self.mCurPortions[0].getInfo()

    def _getCurPortions(self):
        return self.mCurPortions

    def _setupPortions(self):
        offset = self.mTab[self.mCurTabIdx + 3]
        size = self.mTab[self.mCurTabIdx + 4]
        self.mInput.seek(offset)
        buffer = BytesIO(bz2.decompress(
            self.mInput.read(size)))
        self.mCurPortions = [AD_Portion(buffer)
            for sample_name in self.mSamples]

    def getAD_seq(self, chrom, pos):
        assert self.mTab is not None
        if (self.mCurPortions is not None and
                self.mCurPortions[0].isOf(chrom, pos)):
            return [portion.getAD(pos) for portion in self.mCurPortions]
        self.mCurPortions, self.mCurTabIdx = None, None
        for idx0 in range(0, len(self.mTab), 5):
            if chrom != self.mTab[idx0]:
                continue
            if not (self.mTab[idx0 + 1] <= pos <
                    self.mTab[idx0 + 1] + self.mTab[idx0 + 2]):
                continue
            self.mCurTabIdx = idx0
            self._setupPortions()
            return [portion.getAD(pos) for portion in self.mCurPortions]
        return None



