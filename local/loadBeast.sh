#!/bin/bash -a
# -*- coding: utf-8 -*-
## @package runBeast.sh
# @author Sebastien Ravel

export OMP_NUM_THREADS=4
source /etc/profile.d/modules.sh
module load system/java/jre7 system/java/jdk6 compiler/gcc/4.9.2 bioinfo/beagle/4.0 bioinfo/beagle-lib/20150321 bioinfo/BEAST/1.8.1 mpi/openmpi/1.6.5
