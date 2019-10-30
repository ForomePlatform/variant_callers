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

import array, bz2
from io import BytesIO

#========================================
# Builder of abstract model based on array.array
#========================================
class MDL_Builder:
    def __init__(self, fname, prefix, array_type, record_size,
            buffer_size = 50000):
        self.mOutput = open(fname, 'wb')
        self.mRecSize = record_size
        self.mArrayType = array_type
        self.mBufSize = buffer_size

        self.mOutput.write(prefix.encode())
        self.mOutput.write(("%s/%d\n" % (array_type, record_size)).encode())
        self.mRootPos = self.mOutput.tell()
        array.array('Q', [0, 0]).tofile(self.mOutput)
        self.mTab = array.array('Q')
        self.mCurChrom = None
        self.mBufPos = None
        self.mBufCur = None
        self.mBuffer = array.array(self.mArrayType)

    def flushBuffer(self):
        if len(self.mBuffer) == 0:
            return
        offset = self.mOutput.tell()
        io_buffer = BytesIO()
        self.mBuffer.tofile(io_buffer)
        offset = self.mOutput.tell()
        self.mOutput.write(bz2.compress(io_buffer.getvalue()))
        size = self.mOutput.tell() - offset
        self.mTab.extend([self.mCurChrom, self.mBufPos,
            self.mBufCur - self.mBufPos, offset, size])
        self.mBufPos = self.mBufCur
        self.mBuffer = array.array(self.mArrayType)

    def addRecord(self, chrom, pos, rec_data):
        assert len(rec_data) == self.mRecSize
        if chrom != self.mCurChrom:
            self.flushBuffer()
            self.mCurChrom = chrom
            self.mBufPos = self.mBufCur = pos
        else:
            assert self.mBufCur + 1 == pos
            self.mBufCur += 1
        self.mBuffer.extend(rec_data)
        if self.mBufCur >= self.mBufPos + self.mBufSize:
            self.flushBuffer()

    def close(self):
        self.flushBuffer()
        root_array = array.array('Q', [self.mOutput.tell(), len(self.mTab)])
        self.mTab.tofile(self.mOutput)
        self.mOutput.seek(self.mRootPos)
        root_array.tofile(self.mOutput)
        self.mOutput.close()
        self.mOutput = None

#========================================
class MDL_Reader:
    def __init__(self, fname, prefix, array_type, record_size):
        self.mInput = open(fname, 'rb')
        self.mRecSize = record_size
        self.mArrayType = array_type

        title =  self.mInput.read(len(prefix.encode())).decode()
        assert title == prefix, "Unexpected prefix in MDL file: " + title
        meta_line = self.mInput.readline().decode().rstrip()
        exp_meta = "%s/%d" % (self.mArrayType, self.mRecSize)
        assert meta_line == exp_meta, (
            "Wrong MDL model setup: " + meta_line + "/" + exp_meta)
        root_array = array.array('Q')
        root_array.fromfile(self.mInput, 2)
        self.mInput.seek(root_array[0])
        self.mTab = array.array('Q')
        self.mTab.fromfile(self.mInput, root_array[1])
        self.mCurTabIdx = None
        self.mCurBuffer = None

    def close(self):
        self.mInput.close()
        self.mInput = None

    def _directScanPortion(self):
        if self.mCurTabIdx is None:
            idx0 = 0
        else:
            idx0 = self.mCurTabIdx + 5
            if idx0 >= len(self.mTab):
                return None, None, None
        self._setupBuffer(idx0)
        return self.mTab[self.mCurTabIdx:self.mCurTabIdx + 3]

    def _setupBuffer(self, idx):
        self.mCurTabIdx = idx
        buf_chrom, buf_pos, buf_size, offset, byte_size = self.mTab[
            self.mCurTabIdx : self.mCurTabIdx + 5]
        self.mInput.seek(offset)
        io_buffer = BytesIO(bz2.decompress(
            self.mInput.read(byte_size)))
        self.mCurBuffer = array.array(self.mArrayType)
        self.mCurBuffer.fromfile(io_buffer, buf_size * self.mRecSize)

    def _getModel(self, chrom, pos):
        if self.mCurBuffer is None:
            return None
        buf_chrom, buf_pos, buf_size = self.mTab[
            self.mCurTabIdx : self.mCurTabIdx + 3]
        if buf_chrom == chrom and buf_pos <= pos < buf_pos + buf_size:
            shift = (pos - buf_pos) * self.mRecSize
            return self.mCurBuffer[shift: shift + self.mRecSize]
        return None

    def getModel(self, chrom, pos):
        assert self.mTab is not None
        ret = self._getModel(chrom, pos)
        if ret is not None:
            return ret
        self.mCurBuffer, self.mCurTabIdx = None, None
        for idx0 in range(0, len(self.mTab), 5):
            buf_chrom, buf_pos, buf_size = self.mTab[idx0: idx0 + 3]
            if buf_chrom == chrom and buf_pos <= pos < buf_pos + buf_size:
                self._setupBuffer(idx0)
                return self._getModel(chrom, pos)
        return None
