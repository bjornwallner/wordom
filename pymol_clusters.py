#!/usr/bin/env python2.7
from __future__ import print_function
import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
from pymol import cmd
from wordom_parser import read_avg_strength, read_avg_clusters
#from visuals import color_selections
import sys
import re
import seaborn as sns
import pickle
import os
from multiprocessing import Pool

# Library functions
def select_clusters(clusters):
    clusternames = []
    for [cluster, cnum] in zip(clusters, range(len(clusters))):
        residues = []
        for resi in cluster:
            resi=re.sub(r':\w',':',resi)
            residues.append("(chain {} and resi {})".format(
                resi.split(':')[0], resi.split(':')[1]))
        residues = " or ".join(residues)
        clusternames.append("c{:>02d}".format(cnum))
        print(residues)
#        print(clusternames[-1])
        cmd.select(clusternames[-1], "({}) and name CA".format(residues))
    return clusternames


def color_selections(selections):
    cl = len(selections)
    colornames = []
    rgb_range = sns.color_palette("Set1",n_colors=cl)
    for [selection, color] in zip(selections, rgb_range):
       # color=rgb_range[0] #remove to not color everything with one color.
        cmd.set_color(selection, color)
        colornames.append(selection)
#        print(color)
        cmd.color(colornames[-1], selection)
        #print(selection + " " + str(color))
    return colornames




def get_color(resa,resb):
    shell1=["A:240"]
    shell2=["A:182","A:184","A:204","A:201","A:216"]
    shell3=["A:151","A:155","A:156","A:199","A:205","A:208","A:214","A:238","A:242"]
    shell4=["A:154","A:158","A:163","A:166","A:167","A:170","A:180","A:185","A:196","A:203","A:207","A:211","A:212","A:219","A:68"]
    shell5=["A:135","A:136","A:149","A:160","A:162","A:164","A:171","A:173","A:174","A:178","A:192","A:197","A:200","A:206","A:222","A:230","A:235","A:236","A:237","A:26","A:29","A:66","A:69"]
    shell6=["A:134","A:138","A:139","A:143","A:146","A:152","A:157","A:159","A:169","A:183","A:187","A:19","A:195","A:220","A:223","A:227","A:23","A:232","A:24","A:244","A:46","A:87","A:89","A:90","A:91"]


    
    colors=["red","orange","yellow","cyan","blue","magenta","hotpink","purple"] #"red","white","blue"]
    #for s in (90,70,50):
    #    colors.append("gray" + str(s))
    #    print(s)
    
    if resa in shell1 or resb in shell1:
        return colors[0]
    if resa in shell2 or resb in shell2:
        return colors[1]
    if resa in shell3 or resb in shell3:
        return colors[2]
    if resa in shell4 or resb in shell4:
        return colors[3]
    if resa in shell5 or resb in shell5:
        return colors[4]
    if resa in shell6 or resb in shell6:
        return colors[5]
    return colors[-1]
    

    

def bond_connections(clusters, interactions):
    cutoff=10
##    show_list=["A:240","A:182","A:184","A:204","A:201","A:216"]
#    shell1=["A:240"]
##    shell1=["A:204","A:240"]
#    shell2=["A:182","A:184","A:204","A:201","A:216"]
#    shell3=["A:151","A:155","A:156","A:199","A:205","A:208","A:214","A:238","A:242"]
#    shell4=["A:154","A:158","A:163","A:166","A:167","A:170","A:180","A:185","A:196","A:203","A:207","A:211","A:212","A:219","A:68"]
#    shell5=["A:135","A:136","A:149","A:160","A:162","A:164","A:171","A:173","A:174","A:178","A:192","A:197","A:200","A:206","A:222","A:230","A:235","A:236","A:237","A:26","A:29","A:66","A:69"]
#    shell6=["A:134","A:138","A:139","A:143","A:146","A:152","A:157","A:159","A:169","A:183","A:187","A:19","A:195","A:220","A:223","A:227","A:23","A:232","A:24","A:244","A:46","A:87","A:89","A:90","A:91"]
#     
#    show_list_wt=shell1 + shell2 + shell3 + shell4 + shell5 + shell6
#    show_list_S=shell1
#    show_list=shell1
#    #show_list=show_list_wt
#    #show_list.append(shell1)
#    #show_list.append(shell2)
#    #print(show_list)

    for cluster in clusters:
        #print(cluster)
        for i in range(len(cluster)):
            for j in range(i, len(cluster)):
                resa = cluster[i]
                resb = cluster[j]
                # Only draw bonds if interaction strength at all present
                if resa in interactions:
                    if resb in interactions[resa]:
                        resa_=re.sub(r':\w',':',resa)
                        resb_=re.sub(r':\w',':',resb)
                        a = "chain {} and resi {} and name CA".format(
                            resa_.split(':')[0], resa_.split(':')[1])
                        b = "chain {} and resi {} and name CA".format(
                            resb_.split(':')[0], resb_.split(':')[1])

                        
#                        if interactions[resa][resb] > cutoff and (resa in show_list or resb in show_list):
                        if interactions[resa][resb] > cutoff: #and (resa in show_list or resb in show_list):

                            #if resa not in show_list:
                            #    print(resa)
                            #if resb not in show_list:
                            #    print(resb)
#                            print(resa + " " + resb)
                            print(a + " " + b + " " + str(interactions[resa][resb]))
                            cmd.bond(a, b)
                            cmd.set_bond("line_width",
                                         (interactions[resa][resb]/3),
                                         a, b)
                            cmd.set_bond("line_color",get_color(resa,resb),a,b)



def show_cluster(clusters):
    for cluster in clusters:
        for resi in cluster:
            resi=re.sub(r':\w',':',resi)
            
            sele = "chain {} and resi {} and name CA".format(
                resi.split(':')[0], resi.split(':')[1])
            #print(resi,sele)
            cmd.show("lines", sele)
            cmd.show("spheres", sele)

#def visualize_in_pymol(input_args):
def visualize_in_pymol(args,outfolder,clusters,interactions,pdb,cut,freq):

    #(args,outfolder,clusters,interactions,pdb,cut,freq)=input_args
    #print(cut)
    #return
    pymol.finish_launching()
    cmd.reinitialize()
    cmd.load(pdb)
    cmd.hide("everything")
#    cmd.show("ribbon")
    cmd.show("cartoon")
    # Create bindings and selections, and color them
    bond_connections(clusters, interactions)
    selections = select_clusters(clusters)
    colors = color_selections(selections)
    print(colors)

#    colors = color_selections([sele])


    #cmd.color("gray70")
    #cmd.set("sphere_scale",0.75)
    cmd.space("cmyk")
    #cmd.set("ray_shadow","off")
    #cmd.bg_color("white")
    sele = "chain A and resi 131 and not name H*"
    cmd.show("sticks",sele)
    cmd.util.cnc(sele)
    #cmd.center("chain A")
    cmd.deselect()
    ### cut below here and paste into script ###
    cmd.set_view ([\
                   0.185403526,   -0.784912288,   -0.591205537,\
                   -0.726189077,    0.295874894,   -0.620551527,\
                   0.662004888,    0.544381559,   -0.515144050,\
                   -0.001714554,   -0.001769811, -125.023078918,\
                   30.330446243,   78.181671143,  146.038742065,\
                   98.763908386,  151.776306152,  -20.000000000] )
    outfile="{}cluster{}-{}.png".format(outfolder,cut,freq)
    if args.ray or args.generate_all:
        cmd.save(outfile)
### cut above here and paste into script ###

    # Show clusters
    if args.show[0] is None:
        show_cluster(clusters)
    else:
        shw = [int(c) for c in args.show[0].split(',')]
        for c in args.show[0]:
            show_cluster(clusters[c - 1])
    cmd.save(outfile)
    cmd.quit()



# Main; for callable scripts
def main():
    from argparse import ArgumentParser
    from sys import argv, stdin
    parser = ArgumentParser(
        description="Display Wordom PSN analysis clusters in PyMOL.")
    parser.add_argument(
        "-c", nargs=1, default=[None], metavar="float", 
        help="Cutoff to use (must be present in AVGfile)" +
        ", default=Use lowest found")
    parser.add_argument(
        "-f", nargs=1, default='50.0', metavar=float, 
        help="frequencey Cutoff to use (must be present in AVGfile)")
    
    parser.add_argument(
        "-show", nargs=1, default=[None], metavar="int[,int[...]]",
        help="Show specified clusters")
    parser.add_argument(
        "-ray", action='store_true', help="raytrace and save fig")
    parser.add_argument(
        "-generate_all", action='store_true', help="generate raytraced figs for all cutoffs and frequences")
    parser.add_argument(
        "avg", nargs=1, metavar="AVGfile", type=str,
        help="Wordom PSN avg file")
    parser.add_argument(
        "pdb", nargs=1, metavar="PDBfile", type=str,help="PDB file to draw")
    args = parser.parse_args()
    
    # Finish pymol launch
#    pymol.finish_launching()
#    cmd.feedback("disable","all","actions")
#    cmd.feedback("disable","all","results")
    # Set variables here
    pdb = args.pdb[0]
    avg = args.avg[0]
    cut = args.c[0]
    freq= float(args.f[0])
    shw = args.show[0]
    ray = args.ray or args.generate_all


    outfolder=os.path.abspath(avg)
    interactions = {}
    clusters_all = {}
    with open(avg, 'r') as infile:
        interactions = read_avg_strength(infile)
        #print(interactions)
        #sys.exit()
        infile.seek(0)
        clusters_all = read_avg_clusters(infile)
        #print(clusters)
        #sys.exit()
#    print("Hello")
    #sys.exit(1)
    # Select the cutoff
    #with open("cutoffs", 'wb') as f:
    #        pickle.dump(list(clusters.keys()), f)
    if cut is not None:
        # If provided
        cut = float(cut)
    else:
        # Default
        cutoffs = list(clusters_all.keys())
        cutoffs.sort()
        cut = cutoffs[0]

    cuts_to_sweep=[cut]
    freq_to_sweep=[freq]

    if args.generate_all:
        cuts_to_sweep=clusters_all.keys()
        freq_to_sweep=clusters_all[cut].keys()

    todo=[]
    for cut_ in cuts_to_sweep:
        for freq_ in freq_to_sweep:
            print("{} {}".format(cut_,freq_))
            todo.append((args,outfolder,clusters_all[cut_][freq_],interactions,pdb,cut_,freq_))
    #print(len(todo))
    #sys.exit()
    #pool = Pool(processes=8)
    #pool.map(visualize_in_pymol,todo)
    #print(todo)
    #visualize_in_pymol((todo))
    visualize_in_pymol(args,outfolder,clusters_all[cut][freq],interactions,pdb,cut,freq)
    #for cut_ in cuts_to_sweep:
    #    for freq_ in freq_to_sweep:


    if 1==0:
        # Select clusters
        print("Running for {} {}".format(cut_,freq_))
    #    clusters_freq = clusters_all[cut_]
    #    print(clusters_freq)
        clusters = clusters_all[cut_][freq_]
        #print(type(clusters_freq))
        #print(clusters_freq.keys())
        #with open("freq", 'wb') as f:
        #        pickle.dump(list(clusters_freq.keys()), f)
        #print(type(freq))
        #clusters = clusters_freq[80.0]
    #    print(clusters_freq)
    #    print(interactions)
         #ys.exit()
    #    print(clusters)
        cmd.reinitialize()
        cmd.load(pdb)
        cmd.hide("everything")
    #    cmd.show("ribbon")
        cmd.show("cartoon")

        # Create bindings and selections, and color them
        bond_connections(clusters, interactions)
        selections = select_clusters(clusters)
        colors = color_selections(selections)
        print(colors)

    #    colors = color_selections([sele])


        #cmd.color("gray70")
        #cmd.set("sphere_scale",0.75)
        cmd.space("cmyk")
        #cmd.set("ray_shadow","off")
        #cmd.bg_color("white")
        sele = "chain A and resi 131 and not name H*"
        cmd.show("sticks",sele)
        cmd.util.cnc(sele)
        #cmd.center("chain A")
        cmd.deselect()
        ### cut below here and paste into script ###
        cmd.set_view ([\
                       0.185403526,   -0.784912288,   -0.591205537,\
                       -0.726189077,    0.295874894,   -0.620551527,\
                       0.662004888,    0.544381559,   -0.515144050,\
                       -0.001714554,   -0.001769811, -125.023078918,\
                       30.330446243,   78.181671143,  146.038742065,\
                       98.763908386,  151.776306152,  -20.000000000] )
        outfile="{}cluster{}-{}.png".format(outfolder,cut,freq)
        if args.ray:
            cmd.save(outfile)
    ### cut above here and paste into script ###

        # Show clusters
        if shw is None:
            show_cluster(clusters)
        else:
            shw = [int(c) for c in shw.split(',')]
            for c in shw:
                show_cluster(clusters[c - 1])
        cmd.save(outfile)
        cmd.quit()

if __name__ == '__main__':
    main()
