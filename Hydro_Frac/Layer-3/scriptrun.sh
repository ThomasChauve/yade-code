#!/bin/bash

#OAR -n yade-layer-d-2Y-Y-2Y
#OAR -l /nodes=1/core=8,walltime=72:00:00

# Modules loading
source /soft/env.bash
module load yade/2018.02 

# Launch compute job
yade-2018-04-24.git-a8f9844 hydraulicInjection-layer.py > nohup.out
