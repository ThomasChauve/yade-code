#!/bin/bash

#OAR -n yade-layer5-00d-d-midY
#OAR -l /nodes=1/core=6,walltime=144:00:00

# Modules loading
source /soft/env.bash
module load yade/2018.02 

# Launch compute job
yade-2018-04-24.git-a8f9844 hydraulicInjection-layer.py > nohup.out
