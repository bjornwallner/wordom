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
import collections

def visualize_path((path,freq,pdb,outfile)):
    nodes=path.split('=>')
    save_pse=True

    #print outfile
    #print freq
    #sys.exit()
    #print nodes
    node_color=sns.color_palette("RdBu",n_colors=len(nodes))
    #print node_color
    cmd.reinitialize()
    cmd.load(pdb)
    cmd.hide("everything")
#    cmd.show("ribbon")
    cmd.show("cartoon")
    #resa = re.sub(r':\w', ':', nodes[0])
    (chain,resnum)=re.split(':\w',nodes[0])
    a = "chain {} and resi {} and name CA".format(chain,resnum)
    cmd.show("spheres", a)

    colorname="color"+str(0)
    cmd.set_color(colorname, node_color[0])
    cmd.color(colorname,a)
    cmd.label(a,'" %s:%s" % (resi, resn)')
    #cmd.label(a, '"%s" % (resi)')
    for i in range(1,len(nodes)):

        #resb = re.sub(r':\w', ':', nodes[i])
        (chain,resnum)=re.split(':\w',nodes[i])
        b = "chain {} and resi {} and name CA".format(chain,resnum)
        #print a,b
        cmd.bond(a, b)
        cmd.set_bond("line_width",5,a,b)
        cmd.set_bond("line_color","red",a,b)
        cmd.show("lines", a)
        cmd.show("spheres", a)

        cmd.show("lines", b)
        cmd.show("spheres", b)
        colorname="color"+str(i)
        cmd.set_color(colorname, node_color[i])
        cmd.color(colorname,b)
        #cmd.label(b,'" %s:%s" % (resi, resn)')
        cmd.label(b, '" %s%s" % (one_letter[resn],resi)')
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

    norm_factor=freq[nodes[0]]

    node_freqs={k: freq[k] for k in nodes}
    #print node_freqs
    #print freq[nodes]
    #print sorted(node_freqs.values())
    node_freq_min=min(node_freqs.values())
    #print node_freq_min
    #sys.exit()

    scale_cutoff=min(0.25,min(node_freqs.values())/norm_factor)

    for res in freq.keys():
        scale=freq[res]/norm_factor
        #print res,scale,scale_cutoff
        if scale>scale_cutoff:
            (chain,resnum)=re.split(':\w',res)
            sele = "chain {} and resi {} and name CA".format(chain,resnum)
            cmd.show("spheres", sele)
            cmd.set("sphere_scale", scale,selection=sele)
            if not res in nodes: #give residue not on highest freq path a different color
                last_color=len(nodes)-1
                colorname="color"+str(last_color)
                cmd.color(colorname,sele)
                cmd.label(sele, '" %s%s" % (one_letter[resn],resi)')

    cmd.util.cnc(sele)
    #cmd.center("chain A")
    cmd.deselect()
    ### cut below here and paste into script ###
    cmd.set_view([\
              0.354133159,   -0.915425777,    0.191264138,\
              -0.626513124,   -0.384070098,   -0.678214610,\
              0.694314659,    0.120348796,   -0.709537089,\
              -0.000188543,   -0.000061929, -119.472831726,\
              33.501632690,   81.432159424,  143.041992188,\
              87.118530273,  151.834289551,  -20.000000000\
              ])
    

### cut above here and paste into script ###

    #outfile="{}cluster{}-{}.png".format(outfolder,cut,freq)
   # if ray:
   #     cmd.ray()
   #     cmd.sync(10)
    cmd.save(outfile)
    cmd.sync(10)
    if save_pse:
        outfile_pse=re.sub('.png','.pse',outfile)
        cmd.save(outfile_pse)
        cmd.sync(10)
    #
    #print outfile
    #os.system("sleep 2")
    #cmd.quit()
def get_res_freq(frame_file):
    freq=collections.OrderedDict()
    freq['all']=collections.OrderedDict()
    with open(frame_file) as fp:
        for line in fp:
            line=line.rstrip()
            m=re.search('(\d+)\s+(.+=>.+$)',line)
            if m:
                #print frame_file,line
                frame=m.group(1)
                path=m.group(2)
                path_res=path.split('=>')
                if not frame in freq:
                    freq[frame]=collections.OrderedDict()

                for node in path_res:
                    if not node in freq[frame]:
                        freq[frame][node]=0
                    if not node in freq['all']:
                        freq['all'][node]=0
                    freq[frame][node] += 1
                    freq['all'][node] += 1



                #print re.split(':\w', path_res[0])[0]
   # for frame in freq.keys():
   #     print frame,freq[frame]

    total_count=sum(freq['all'].values())
    #print freq['all']
    for node in freq['all'].keys():
        freq['all'][node]=float(freq['all'][node])/float(total_count)
    #print freq['all']
  #  print sum(freq['all'].values())
    #sys.exit()
    return freq['all']

def main():
    parser = argparse.ArgumentParser(description='visualising path from psnpath')
    parser.add_argument('-log', default=None, type=str, nargs=1, required=True,
                        help='log file from psnpath')
    parser.add_argument(
        "-pdb", nargs=1, type=str, help="PDB file to draw")

    parser.add_argument(
        "-minlen",  type=int, default=1, help="minimum path length to consider")
    parser.add_argument(
        "-minfreq",  type=float, default=0.0, help="minimum path freq to consider")

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
    count={}
    jobs=[]
    for index,row in df[rows].sort_values('HighFreqLen', ascending=False).iterrows():
        #print row

        if row.HighFreqLen >= args.minlen and row.HighFreqVal >= args.minfreq:
            # print row.HighFreqPath
            #print row
            #continue
            dirname=os.path.dirname(logfile)
            basename=os.path.basename(logfile)
            basename_noext=re.sub('.log$','',basename)
            #print row
            name1=row.Res1.replace(':','_')
            name2=row.Res2.replace(':','_')
            if not name1 in count:
                count[name1]=1
            name_prefix = "{}/PSNPath-{}-{}-{}".format(dirname, basename_noext, name1, name2)
            name_prefix_count = "{}/{}PSNPath-{}-{}-{}".format(dirname, count[name1],basename_noext, name1, name2)
            frame_file=name_prefix +".frame"
            #stat_file=name_prefix + ".stat"
            res_freq=get_res_freq(frame_file)
            #continue
            #stat=get_stat(stat_file)
            outfile="{}.freq{}.len{}.png".format(name_prefix_count,row.HighFreqVal,row.HighFreqLen)
            if not os.path.exists(outfile) or args.overwrite:
                print outfile
                jobs.append((row.HighFreqPath,row.HighFreqVal,pdb,outfile))
                #visualize_path((row.HighFreqPath,row.HighFreqVal,pdb,outfile))
                visualize_path((row.HighFreqPath, res_freq, pdb, outfile))
                #sys.exit()
            count[name1]+=1
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
