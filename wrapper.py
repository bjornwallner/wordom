#!/usr/bin/env python2.7
import pickle
import os
import sys
path_to_script=os.path.abspath(sys.arv[0])

with open("cutoffs", 'rb') as f:
    cutoffs = pickle.load(f)

with open("freq", 'rb') as f:
    freq = pickle.load(f)

print cutoffs
print freq
for c in cutoffs:
    for f in freq:
        cmd="{}/pymol_clusters.py -c {} -f {} {}test/avgpsn {}test/complex.pdb".format(path_to_script,c,f,path_to_script,path_to_script)
        print cmd
        #os.system(cmd)
