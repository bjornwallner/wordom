#!/usr/bin/env python2.7
import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
from pymol import cmd
from wordom_parser import read_avg_strength, read_avg_clusters
import sys
import re
import seaborn as sns
#import pickle
import os
#from multiprocessing import Pool
import argparse
import pandas as pd
def visualize_path(path,pdb,outfile):
    nodes=path.split('=>')
    print outfile
    print nodes
    pymol.finish_launching()
    cmd.reinitialize()
    cmd.load(pdb)
    cmd.hide("everything")
#    cmd.show("ribbon")
    cmd.show("cartoon")

    #select the nodes.

def main():
    parser = argparse.ArgumentParser(description='visualising path from psnpath')
    parser.add_argument('-log', default=None, type=str, nargs=1, required=True,
                        help='log file from psnpath')
    parser.add_argument(
        "-pdb", nargs=1, metavar="PDBfile", type=str,help="PDB file to draw")
    args = parser.parse_args()
    logfile=args.log[0]
    pdb=args.pdb[0]
    print args
    df=pd.read_csv(logfile,sep='\s+',index_col=0)
    #print df.sort_values('HighFreqLen', ascending=False)
   # print df.HighFreqPath
    #column header gets shiften one to the right
    #Skip first column
    columns=df.columns[1:]
    #drop last empty column
    df.drop(df.columns[-1], axis=1, inplace=True)
    #shift the columns...
    df.columns=columns
    rows=df.NullFrames!='***SKIPPED***'
    for index,row in df[rows].sort_values('HighFreqLen', ascending=False).iterrows():
        #print row
        print row.HighFreqPath
        outfile="{}.{}.{}.len{}.png".format(logfile,row.Res1,row.Res2,row.HighFreqLen)
        visualize_path(row.HighFreqPath,pdb,outfile)
        #print row['HighFreqPath']
    #print df.HighFreqPath
    #print df
    #pymol.finish_launching()
    #cmd.reinitialize()
    #cmd.load(pdb)
    #cmd.hide("everything")
#    cmd.show("ribbon")
    #cmd.show("cartoon")


if __name__ == '__main__':
    main()