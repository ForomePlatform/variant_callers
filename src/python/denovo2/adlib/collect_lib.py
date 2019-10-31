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

#========================================
class AD_LibBuilder:
    def __init__(self, fname, samples):
        self.mOutput = open(fname, 'wb')
        self.mOutput.write(AD_LibReader.PREFIX)
        self.mSampleCount = len(samples)
        for name in samples:
            self.mOutput.write(name.encode() + b"\n")
        self.mOutput.write(b"\n")
        self.mRootPos = self.mOutput.tell()
        array.array('Q', [0, 0]).tofile(self.mOutput)
        self.mTab = array.array('Q')

    def addPortions(self, portions):
        assert len(portions) == self.mSampleCount
        buffer = BytesIO()
        portion0 = portions[0]
        portion0.toFile(buffer)
        for portion in portions[1:]:
            assert portion.getChrom() == portion0.getChrom()
            assert portion.getShift() == portion0.getShift()
            assert portion.getSize() == portion0.getSize()
            portion.toFile(buffer)
        offset = self.mOutput.tell()
        self.mOutput.write(bz2.compress(buffer.getvalue()))
        size = self.mOutput.tell() - offset
        self.mTab.extend([portion0.getChrom(), portion0.getShift(),
            portion0.getSize(), offset, size])

    def close(self):
        root_array = array.array('Q', [self.mOutput.tell(), len(self.mTab)])
        self.mTab.tofile(self.mOutput)
        self.mOutput.seek(self.mRootPos)
        root_array.tofile(self.mOutput)
        self.mOutput.close()
        self.mOutput = None

#========================================
#========================================
if __name__=="__main__":
    import sys, os
    from .ad_person import AD_PersonData
    from .ad_lib import AD_LibReader

    if sys.argv[1] not in ("make", "info"):
        print("\n".join(["Available modes:",
            "make <out.ldx> [<in.ldx>] *idx",
            "info <.ldx>"]), file = sys.stderr)
        sys.exit()

    if sys.argv[1] == "info":
        inp_lib  = AD_LibReader(sys.argv[2])
        names = list(inp_lib.iterSampleNames())
        print("Info for library", sys.argv[2],
            "samples[%d]" % len(names))
        for name in names:
            print("\t", name)
        sys.exit()

    out_lib_name = sys.argv[2]
    inp_lib = None
    inp_readers = []
    sample_names = []
    if not out_lib_name.endswith('.ldx'):
        print("Output file extension must be .ldx:", out_lib_name,
            file = sys.stderr)
        sys.exit()
    assert sys.argv[1] == "make"
    files = sys.argv[3:]
    if files[0].endswith('.ldx'):
        inp_lib  = AD_LibReader(files[0])
        for nm in inp_lib.iterSampleNames():
            sample_names.append(nm)
        files = files[1:]
    for fname in files:
        if not fname.endswith('.idx'):
            print("Input file extension must be .idx:", fname,
                file = sys.stderr)
            sys.exit()
        inp_readers.append(AD_PersonData(fname, True))
        sample_names.append(os.path.basename(fname)[:-4])

    out_lib = AD_LibBuilder(out_lib_name, sample_names)
    cnt = 0
    while True:
        portions = []
        if inp_lib is not None:
            portion_info = inp_lib._nextPortions()
            if portion_info is None:
                for inp_rd in inp_readers:
                    if inp_rd.directReadPortion() is not None:
                        print ("Extra portion in", inp_rd.getFName())
                        sys.exit()
                break
            portions += inp_lib._getCurPortions()
        for inp_rd in inp_readers:
            portion = inp_rd.directReadPortion()
            if portion is None:
                if len(portions) > 0:
                    print ("Extra portion in", inp_rd.getFName())
                    sys.exit()
                for inp_rd1 in inp_readers[1:]:
                    if inp_rd1.directReadPortion() is not None:
                        print ("Extra portion in", inp_rd1.getFName())
                        sys.exit()
                break
            portions.append(portion)
        if len(portions) == 0:
            break
        out_lib.addPortions(portions)
        cnt += 1
        if cnt % 100 == 0:
            print(cnt, "portions...", file = sys.stderr)
            sys.stderr.flush()
    out_lib.close()
    print("Done:", cnt, file = sys.stderr)
