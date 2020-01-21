#!/bin/bash

## source /net/lib/python3.3/bin/activate
eval "$(/net/apps/x86_64/miniconda3/bin/conda shell.bash hook)"
# conda activate callers
export PYTHONPATH=/net/bgm/versions/develop/forome_callers/src/python
echo $PYTHONPATH
python $*

