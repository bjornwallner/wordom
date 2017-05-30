#!/usr/bin/env python2.7
import pickle
import os
import sys
from multiprocessing import Pool

def job((path,c,f)):
    cmd="{}/pymol_clusters.py -c {} -f {} output_wordom_50_800/avgpsn  {}/test/complex.pdb".format(path,c,f,path)
    print cmd
    os.system(cmd)
    return True

path_to_script=os.path.abspath(os.path.dirname(sys.argv[0]))

with open("cutoffs", 'rb') as f:
    cutoffs = pickle.load(f)

with open("freq", 'rb') as f:
    freq = pickle.load(f)
args=[]
for c in cutoffs:
    for f in freq:
        if c > 6:
            args.append((path_to_script,c,f))

print args
#sys.exit()
pool = Pool(processes=8)
pool.map(job,args)

#print cutoffs
#print freq
#for c in cutoffs:
#    for f in freq:
#        cmd="{}/pymol_clusters.py -c {} -f {} {}test/avgpsn {}test/complex.pdb".format(path_to_script,c,f,path_to_script,path_to_script)
#        print cmd
        #os.system(cmd)
