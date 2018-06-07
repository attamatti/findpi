#!/usr/bin/env python

# wrapper to run find_pi on lots of pdb files:
# outputs hits.txt with files that have stacks larger than 3
#filename.pdb n
#n = # hits

from glob import glob
import subprocess
from mpi4py import MPI

# MPI info
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

entfiles = glob('*.ent')


# get files and count
# decide how many files to assign to each MPI
filenum = len(entfiles)
groupsize = filenum//comm.size

# deal with case of # of MPIs > # of files
if groupsize == 0:
    groupsize = 1
    rangesize = filenum
else:
    rangesize = comm.size

# get the last few files that don't fit evenly into a group 
extras = filenum%comm.size

# run newstack on each MPI
for r in range(rangesize):
    if comm.rank == r:
        print(r*groupsize,(r+1)*groupsize)
        for i in entfiles[r*groupsize:(r+1)*groupsize]:
            print('{0} on rank {1}'.format(i,r))
            subprocess.Popen(['python find_pi.py {0} --h'.format(i)],  shell=True)
            
            
# finish the last files        
if comm.rank == 0:
    if extras >= 1:
        print('****')
        for i in entfiles[comm.size*groupsize:]:
            print('{0} finishing on rank 0'.format(i))
            subprocess.Popen(['python find_pi.py {0} --h'.format(i)],  shell=True)

    elif comm.size != len(entfiles): 
        print('Finishing {0} on rank 0'.format(entfiles[-1]))
        run_prog = subprocess.Popen(['python find_pi.py {0} --h'.format(i)],  shell=True, stdout=subprocess.PIPE)
        result = run_prog.stdout.readlines()
        if float(result[0].split()[1]) > 0:
            output.write('{0}\n'.format(i))

