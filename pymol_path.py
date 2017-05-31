#!/usr/bin/env python2.7
import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
pymol.finish_launching()
from pymol import cmd
from wordom_parser import read_avg_strength, read_avg_clusters
import sys
import re
import seaborn as sns
#import pickle
import os
from multiprocessing import Pool
import argparse
import pandas as pd
def visualize_path((path,freq,pdb,outfile)):
    nodes=path.split('=>')
    #print outfile
    #print nodes
    node_color=sns.color_palette("RdBu",n_colors=len(nodes))
    #print node_color
    cmd.reinitialize()
    cmd.load(pdb)
    cmd.hide("everything")
#    cmd.show("ribbon")
    cmd.show("cartoon")
    resa = re.sub(r':\w', ':', nodes[0])
    a = "chain {} and resi {} and name CA".format(resa.split(':')[0],resa.split(':')[1])
    cmd.show("spheres", a)
    #cmd.select("sele",a)
    colorname="color"+str(0)
    cmd.set_color(colorname, node_color[0])
    cmd.color(colorname,a)
    cmd.label(a,'" %s:%s" % (resi, resn)')
    #cmd.label(a, '"%s" % (resi)')
    for i in range(1,len(nodes)):
        resb = re.sub(r':\w', ':', nodes[i])
        b = "chain {} and resi {} and name CA".format(resb.split(':')[0], resb.split(':')[1])
        #print a,b
        cmd.bond(a, b)
        cmd.set_bond("line_width",freq/5,a,b)
        cmd.set_bond("line_color","red",a,b)
        cmd.show("lines", a)
        cmd.show("spheres", a)
        cmd.show("lines", b)
        cmd.show("spheres", b)

        colorname="color"+str(i)
        cmd.set_color(colorname, node_color[i])
        cmd.color(colorname,b)
        cmd.label(b,'" %s:%s" % (resi, resn)')
        #cmd.label(b, '"%s" % (resi)')
        a=b
    #cmd.color(node_color[len)],b)
    #cmd.space("cmyk")
    #cmd.set("ray_shadow","off")
    #cmd.bg_color("white")
    #cmd.label_position([3,2,1])
    cmd.set("label_position", (2, 2, 2))
    sele = "chain A and resi 131 and not name H*"
    cmd.show("sticks",sele)
    cmd.util.cnc(sele)
    #cmd.center("chain A")
    cmd.deselect()
    ### cut below here and paste into script ###
    cmd.set_view ([\
    0.094197616,    0.035088558,   -0.994926095,\
    -0.611217141,    0.790879548,   -0.029978149,\
     0.785817504,    0.610941410,    0.095943026,\
    -0.002951451,   -0.004765198, -125.047027588,\
    30.400051117,   80.221611023,  147.147964478,\
    98.763908386,  151.776306152,  -20.000000000] )

### cut above here and paste into script ###

    #outfile="{}cluster{}-{}.png".format(outfolder,cut,freq)
    cmd.ray()
    cmd.sync(10)
    cmd.save(outfile)
    cmd.sync(10)
    outfile_pse=re.sub('.png','.pse',outfile)
    cmd.save(outfile_pse)
    cmd.sync(10)
    #
    #print outfile
    #os.system("sleep 2")
    #cmd.quit()

def main():
    parser = argparse.ArgumentParser(description='visualising path from psnpath')
    parser.add_argument('-log', default=None, type=str, nargs=1, required=True,
                        help='log file from psnpath')
    parser.add_argument(
        "-pdb", nargs=1, type=str, help="PDB file to draw")

    parser.add_argument(
        "-minlen",  type=int, default=1, help="minimum path length to consider")
    parser.add_argument(
        "-overwrite", action='store_true', help="overwrite if fig already exists")
    #pymol.finish_launching()
    args = parser.parse_args()
    logfile=args.log[0]
    pdb=args.pdb[0]
    #print args
    #print args.minlen
    #sys.exit()
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
    count=1
    jobs=[]
    for index,row in df[rows].sort_values('HighFreqLen', ascending=False).iterrows():
        #print row

        if row.HighFreqLen >= args.minlen:
            # print row.HighFreqPath
            dirname=os.path.dirname(logfile)
            basename=os.path.basename(logfile)
            outfile="{}/{}{}.{}.{}.freq{}.len{}.png".format(dirname,count,basename,row.Res1,row.Res2,row.HighFreqVal,row.HighFreqLen)
            if not os.path.exists(outfile) or args.overwrite:
                print outfile
                jobs.append((row.HighFreqPath,row.HighFreqVal,pdb,outfile))
                #visualize_path((row.HighFreqPath,row.HighFreqVal,pdb,outfile))
                visualize_path((row.HighFreqPath, 25, pdb, outfile))
            count+=1
            #sys.exit()
        #print row['HighFreqPath']

    #pool = Pool(processes=1)
    #pool.map(visualize_path,jobs)
#print df.HighFreqPath
    #print df
    #pymol.finish_launching()
    #cmd.reinitialize()
    #cmd.load(pdb)
    #cmd.hide("everything")
#    cmd.show("ribbon")
    #cmd.show("cartoon")
    cmd.quit()

if __name__ == '__main__':
    main()