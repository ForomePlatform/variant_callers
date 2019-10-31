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


#========================================
# Reads file (hg19.fasta) and represents 
#  its content in serie of portions 
#========================================
class Hg19_Reader:
    def __init__(self, fname,
            block_size = 100000,
            chrom_is_int = False,
            upper_case = False):
        self.mInput = open(fname, "r")
        self.mBlockSize  = block_size
        self.mCurLine    = self.mInput.readline()
        self.mChromIsInt = chrom_is_int
        self.mUpperCase  = upper_case
        self.mCurChrom   = None
        self.mCurDiap    = None
        self.mCurLetters = None

    def close(self):
        self.mInput.close()
        self.mInput = None

    def getCurChrom(self):
        return self.mCurChrom

    def getCurDiap(self):
        return self.mCurDiap

    def getCurLetters(self):
        return self.mCurLetters

    def hasPosition(self, pos):
        return self.mCurDiap[0] <= pos < self.mCurDiap[1]

    def getLetter(self, pos):
        assert self.hasPosition(pos)
        return self.mCurLetters[pos - self.mCurDiap[0]]

    def read(self):
        if not self.mCurLine:
            self.mCurChrom = -1
            self.mCurLetters = None
            self.mCurDiap = None
            return False
        base_count = 0
        lines = []
        while base_count < self.mBlockSize and self.mCurLine:
            if self.mCurLine.startswith('>'):
                if base_count > 0:
                    break
                assert self.mCurLine.startswith('>chr')
                chrom = self.mCurLine.rstrip()[4:]
                if '_' in chrom:
                    # extra pseudo chromosome information, end of real work
                    self.mCurChrom = -1
                    self.mCurLetters = None
                    self.mCurDiap = None
                    return False
                if self.mChromIsInt:
                    if chrom.isdigit():
                        self.mCurChrom = int(chrom)
                    else:
                        self.mCurChrom = {"M": 0, "X": 23, "Y": 24}.get(chrom)
                else:
                    self.mCurChrom = chrom
                self.mCurDiap = [1, 1]
            elif self.mCurChrom is not None:
                letters = self.mCurLine.rstrip()
                base_count += len(letters)
                if self.mUpperCase:
                    lines.append(letters.upper())
                else:
                    lines.append(letters)
            self.mCurLine = self.mInput.readline()
        if self.mCurChrom is None:
            assert not self.mCurLine
            self.mCurChrom = -1
            self.mCurLetters = None
            self.mCurDiap = None
            return False
        assert all([len(letter_seq) == 50 for letter_seq in lines[:-1]])
        assert len(lines[-1]) == 50 or (len(lines[-1]) < 50 and (
                not self.mCurLine or self.mCurLine.startswith('>')))
        self.mCurLetters = ''.join(lines)
        first_pos = self.mCurDiap[1]
        self.mCurDiap = [first_pos, first_pos + base_count]
        return True

#========================================
# __main__ checks if reading is proper
#========================================
if __name__=="__main__":
    import sys
    test_chrom, test_pos, test_ref = "1", 142569315,  "t"

    rd = Hg19_Reader(sys.argv[1])

    cnt = 0
    while rd.read():
        cnt += 1
        if cnt % 100 == 0:
            print (cnt, "blocks...", file = sys.stderr)
        if rd.getCurChrom() == test_chrom and rd.hasPosition(test_pos):
            cur_ref = rd.getLetter(test_pos)
            print("Got it: chrom=", test_chrom, "pos=", test_pos,
                "ref=", cur_ref, "test=", cur_ref == test_ref)
            break
