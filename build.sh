#!/bin/tcsh
module purge
if ($SYS == "amd64_sles11") then
  module load intel/17.0
  mpdule load git
else if ($SYS == "amd64_sles12") then
  module load intel/17.0
  mpdule load git
else if ($SYS == "amd64_sles15") then
  module load intel/19.0.1
  module load git
endif
rm id
git rev-parse HEAD > id
make DEBUG=True IDA=True
make IDA=True
make DEBUG=True
make