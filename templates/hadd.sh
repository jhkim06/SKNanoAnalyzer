#!/bin/bash
export PATH=[MAMBA_BIN_PATH]:$PATH
export MAMBA_ROOT_PREFIX=[MAMBA_ROOT_PREFIX]
eval "$(micromamba shell hook -s bash)"
micromamba activate Nano

cd [WORKDIR]
hadd -f -j 8 [TARGET] output/hists_*.root
rm output/hists_*.root
exit $?
