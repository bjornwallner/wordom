#!/usr/bin/env python2.7

from wordom_parser import read_avg_strength, read_avg_clusters, read_cross_corr

import seaborn as sns
import matplotlib.pyplot as plt

import sys
import argparse
import numpy as np
import os
import pandas as pd

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
   # print labels
    plt.legend()
    #ax.legend(handles, vec.keys())
    plt.savefig('Icrit.png')
    plt.clf()

parser = argparse.ArgumentParser(description='Analysing clusters to determine Icrit')
parser.add_argument('-avg', default=None,type=str, nargs=1,required=True,
                    help='avg file from wordom PSN')

args = parser.parse_args()
#print args
infile=args.avg[0]
cross_corr_file=os.path.dirname(infile) + "/cross-corr"
cross_corr_file2=os.path.dirname(infile) + "/cross-corr_lmi"
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
if os.path.exists(cross_corr_file):
    df=read_cross_corr(cross_corr_file)
    ax=plt.subplots(figsize=(10,10))
    g = sns.heatmap(df,cmap='YlOrRd')
    #plt.subplots_adjust(left=0.20,bottom=0.30)
    #plt.set_title('Pairwise correlation')
    plt.xticks(fontsize=5, rotation=90)
    plt.yticks(fontsize=5, rotation=0)
    #fig = ax.get_figure()
    #fig.tick_params(axis='x', which='major', labelsize=6)
    #plt.setp(ax.xaxis.get_major_ticks(), rotation=90,fontsize=6)
    #plt.setp(ax.yaxis.get_major_ticks(), rotation=0,fontsize=6)
    sns.plt.savefig('heatmap.eps')
    sns.plt.clf()
    df=read_cross_corr(cross_corr_file2)
    ax=plt.subplots(figsize=(10,10))
    g = sns.heatmap(df,cmap='YlOrRd')
    #plt.subplots_adjust(left=0.20,bottom=0.30)
    #plt.set_title('Pairwise correlation')
    plt.xticks(fontsize=5, rotation=90)
    plt.yticks(fontsize=5, rotation=0)
    #fig = ax.get_figure()
    #fig.tick_params(axis='x', which='major', labelsize=6)
    #plt.setp(ax.xaxis.get_major_ticks(), rotation=90,fontsize=6)
    #plt.setp(ax.yaxis.get_major_ticks(), rotation=0,fontsize=6)
    sns.plt.savefig('heatmap_lmi.eps')
    sns.plt.clf()

    sys.exit()
    g = sns.clustermap(df,method='average',cmap='YlOrRd')
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0,fontsize=6)
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90,fontsize=6)
#plt.subplots_adjust(left=0.04,right=0.83,bottom=0.17,top=0.95)
#plt.tight_layout()
    sns.plt.savefig('heatmap_cluster.eps')