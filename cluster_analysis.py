#!/usr/bin/env python2.7

from wordom_parser import read_avg_strength, read_avg_clusters

import seaborn as sns
import matplotlib.pyplot as plt

import sys
import argparse
import numpy as np

#def get_largest_cluster(clusters):
#    sorted()

def Icrit_plot(cluster_all):
    vec={}
    for I in sorted(cluster_all.keys()):
        for frac in sorted(cluster_all[I].keys()):
            if not frac in vec:
                vec[frac]=[]

            s=sorted(cluster_all[I][frac], key=len,reverse=True)
            vec[frac].append((I,len(s[0])))
          #  print vec
            #print I, frac, len(cluster_all[I][frac]),len(s[0])
            #print s
    #print vec[70.0]
    #A=np.array(vec[70.0])
    #print A[:, 0]
    #print A[:, 1]
    ax = plt.axes()
    colors=sns.color_palette("Paired",n_colors=len(vec.keys()))
    n=0
    for frac in sorted(vec.keys()):
        A=np.array(vec[frac])
        plt.plot(A[:,0],A[:,1],label=str(frac),color=colors[n])
        n+=1
    handles, labels = ax.get_legend_handles_labels()
    print labels
    plt.legend()
    #ax.legend(handles, vec.keys())
    plt.savefig('test.png')

parser = argparse.ArgumentParser(description='Analysing clusters to determine Icrit')
parser.add_argument('-avg', default=None,type=str, nargs=1,required=True,
                    help='avg file from wordom PSN')

args = parser.parse_args()
#print args
infile=args.avg[0]

#sys.exit()
interactions = {}
clusters_all = {}
with open(infile, 'r') as infile:
    interactions = read_avg_strength(infile)
    #print(interactions)
    #sys.exit()
    infile.seek(0)
    clusters_all = read_avg_clusters(infile)


Icrit_plot(clusters_all)
