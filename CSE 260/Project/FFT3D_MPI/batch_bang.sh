#!/bin/bash
# This script is intepreted by the Bourne Shell, sh
#
# Documentation for SGE is found in:
# http://docs.oracle.com/cd/E19279-01/820-3257-12/n1ge.html
#
# Tell SGE which shell to run the job script in rather than depending
# on SGE to try and figure it out.
#$ -S /bin/bash
#
# Export all my environment variables to the job
#$ -V
#
# Tun the job in the same directory from which you submitted it
#$ -cwd
#
#
# --- Don't change anything above this line ---
#
# Give a name to the job
#$ -N FFT
#
# Specify a time limit for the job, not more than 30 minutes
#$ -l h_rt=00:60:00
#
# Specify the parallel environment and number of cores,
# If not a multiple of 8, you'll get the whole node anyway
#$ -pe orte 8
#
# Join stdout and stderr so they are reported in job output file
#$ -j y
#
#
# Choose the queue to run the job
#
# Debug queue: only one node may be used at a time for up to 30 minutes
# Interactive or batch jobs, maximum of 1 job per user running at a time
#
# Normal queue: job may use all available compute nodes (256 cores)
# for up to 60 minutes
# Batch jobs, maximum of 2 jobs per user running at a time
# To use more than one node, specify the "normal" queue
#$ -q normal.q
# #$ -q debug.q
#
# Specifies the circumstances under which mail is to be sent to the job owner
# defined by -M option. For example, options "bea" cause mail to be sent at the 
# begining, end, and at abort time (if it happens) of the job.
# Option "n" means no mail will be sent.
#$ -m aeb
#
# *** Change to the address you want the notification sent to, and
# *** REMOVE the blank between the # and the $
# $ -M yqzhang@eng.ucsd.edu
#


echo
echo " *** Current working directory"
pwd
echo
echo " *** Compiler"
# Output which  compiler are we using and the environment
mpicc -v
echo
echo " *** Environment"
printenv

echo

echo ">>> Job Starts"
date

# -k     : No communication
# -v     : Verify the result
echo "-------------------------------------------------"
mpirun -np 4   ./driver -n 4   -x 2 -y 2
mpirun -np 4   ./driver -n 4   -x 2 -y 2 -k
mpirun -np 4   ./driver -n 4   -x 2 -y 2 -v
echo "-------------------------------------------------"
mpirun -np 4   ./driver -n 8   -x 2 -y 2
mpirun -np 4   ./driver -n 8   -x 2 -y 2 -k
mpirun -np 4   ./driver -n 8   -x 2 -y 2 -v
echo "-------------------------------------------------"
mpirun -np 4   ./driver -n 16  -x 2 -y 2
mpirun -np 4   ./driver -n 16  -x 2 -y 2 -k
mpirun -np 4   ./driver -n 16  -x 2 -y 2 -v
mpirun -np 16  ./driver -n 16  -x 4 -y 4
mpirun -np 16  ./driver -n 16  -x 4 -y 4 -k
mpirun -np 16  ./driver -n 16  -x 4 -y 4 -v
echo "-------------------------------------------------"
mpirun -np 4   ./driver -n 32  -x 2 -y 2
mpirun -np 4   ./driver -n 32  -x 2 -y 2 -k
mpirun -np 4   ./driver -n 32  -x 2 -y 2 -v
mpirun -np 16  ./driver -n 32  -x 4 -y 4
mpirun -np 16  ./driver -n 32  -x 4 -y 4 -k
mpirun -np 16  ./driver -n 32  -x 4 -y 4 -v
echo "-------------------------------------------------"
mpirun -np 4   ./driver -n 64  -x 2 -y 2
mpirun -np 4   ./driver -n 64  -x 2 -y 2 -k
mpirun -np 4   ./driver -n 64  -x 2 -y 2 -v
mpirun -np 16  ./driver -n 64  -x 4 -y 4
mpirun -np 16  ./driver -n 64  -x 4 -y 4 -k
mpirun -np 16  ./driver -n 64  -x 4 -y 4 -v
mpirun -np 64  ./driver -n 64  -x 8 -y 8
mpirun -np 64  ./driver -n 64  -x 8 -y 8 -k
mpirun -np 64  ./driver -n 64  -x 8 -y 8 -v
echo "-------------------------------------------------"
mpirun -np 4   ./driver -n 128 -x 2 -y 2
mpirun -np 4   ./driver -n 128 -x 2 -y 2 -k
mpirun -np 4   ./driver -n 128 -x 2 -y 2 -v
mpirun -np 16  ./driver -n 128 -x 4 -y 4
mpirun -np 16  ./driver -n 128 -x 4 -y 4 -k
mpirun -np 16  ./driver -n 128 -x 4 -y 4 -v
mpirun -np 64  ./driver -n 128 -x 8 -y 8
mpirun -np 64  ./driver -n 128 -x 8 -y 8 -k
mpirun -np 64  ./driver -n 128 -x 8 -y 8 -v
echo "-------------------------------------------------"
mpirun -np 4   ./driver -n 256 -x 2 -y 2
mpirun -np 4   ./driver -n 256 -x 2 -y 2 -k
mpirun -np 4   ./driver -n 256 -x 2 -y 2 -v
mpirun -np 16  ./driver -n 256 -x 4 -y 4
mpirun -np 16  ./driver -n 256 -x 4 -y 4 -k
mpirun -np 16  ./driver -n 256 -x 4 -y 4 -v
mpirun -np 64  ./driver -n 256 -x 8 -y 8
mpirun -np 64  ./driver -n 256 -x 8 -y 8 -k
mpirun -np 64  ./driver -n 256 -x 8 -y 8 -v
mpirun -np 256 ./driver -n 256 -x 16 -y 16
mpirun -np 256 ./driver -n 256 -x 16 -y 16 -k
mpirun -np 256 ./driver -n 256 -x 16 -y 16 -v
echo "-------------------------------------------------"
mpirun -np 4   ./driver -n 512 -x 2 -y 2
mpirun -np 4   ./driver -n 512 -x 2 -y 2 -k
mpirun -np 4   ./driver -n 512 -x 2 -y 2 -v
mpirun -np 16  ./driver -n 512 -x 4 -y 4
mpirun -np 16  ./driver -n 512 -x 4 -y 4 -k
mpirun -np 16  ./driver -n 512 -x 4 -y 4 -v
mpirun -np 64  ./driver -n 512 -x 8 -y 8
mpirun -np 64  ./driver -n 512 -x 8 -y 8 -k
mpirun -np 64  ./driver -n 512 -x 8 -y 8 -v
mpirun -np 256 ./driver -n 512 -x 16 -y 16
mpirun -np 256 ./driver -n 512 -x 16 -y 16 -k
mpirun -np 256 ./driver -n 512 -x 16 -y 16 -v
echo "-------------------------------------------------"

date
echo ">>> Job Ends"
