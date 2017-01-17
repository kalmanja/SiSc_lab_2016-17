#!/usr/bin/env bash

#BSUB -J para-heat # job name
#BSUB -o run%J.qout # job output (use %J for job ID)
#BSUB -W 00:30 # limits in hours:minutes
#BSUB -M 512 # memory in MB
##BSUB -P lect0015 # (currently off) Use the job queue that has been created
##                 # for this particular lecture.
#BSUB -a openmpi -n 2 # Use MPI on n-cores
## BSUB -x # Run job exclusively NOT FOR DEBUGING, JUST FOR TIMING

# run the process
$MPIEXEC $FLAGS_MPI_BATCH ./heat>mk.txt
