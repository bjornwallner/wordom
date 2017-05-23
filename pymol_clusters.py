#!/usr/bin/env python2.7
from __future__ import print_function
import pymol
from pymol import cmd
from wordom_parser import read_avg_strength, read_avg_clusters
#from visuals import color_selections
import sys
import re
import seaborn as sns

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
        "avg", nargs=1, metavar="AVGfile",
        help="Wordom PSN avg file")
    parser.add_argument(
        "pdb", nargs=1, metavar="PDBfile", help="PDB file to draw")
    arguments = parser.parse_args(argv[1:])
    
    # Finish pymol launch
    pymol.finish_launching()
    
    # Set variables here
    pdb = arguments.pdb[0]
    avg = arguments.avg[0]
    cut = arguments.c[0]
    freq= float(arguments.f[0])
    shw = arguments.show[0]
    
    interactions = {}
    clusters = {}
    with open(avg, 'r') as infile:
        interactions = read_avg_strength(infile)
        #print(interactions)
        #sys.exit()
        infile.seek(0)
        clusters = read_avg_clusters(infile)
        print(clusters)
        #sys.exit()
#    print("Hello")
    #sys.exit(1)
    # Select the cutoff
    if cut is not None:
        # If provided
        cut = float(cut)
    else:
        # Default
        cutoffs = list(clusters.keys())
        cutoffs.sort()
        cut = cutoffs[0]
        
    # Select clusters
    print(cut)
    clusters_freq = clusters[cut]
#    print(clusters_freq)
    clusters = clusters[cut][freq]
    #print(type(clusters_freq))
    #print(clusters_freq.keys())
    #print(type(freq))
    #clusters = clusters_freq[80.0]
#    print(clusters_freq)
#    print(interactions)
     #ys.exit()
#    print(clusters)
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
    #cmd.space("cmyk")
    #cmd.set("ray_shadow","off")
    #cmd.bg_color("white")
    #sele = "chain A and resi 240 and not name H*"
    #cmd.show("sticks",sele)
    #cmd.util.cnc(sele)

    
    #cmd.select("thiopurine-binding_site", "resi 29+32+39+40+42+152+156+157+183+185+195+196+222+224+226+227+230+237")
    #cmd.select("active_site_loop", "resi 34-47")
    #cmd.select("sam_site", "resi 29+33+40+69+71+90+91+134+135+152+153+157")
#    cmd.deselect()
#    cmd.set_view ([\
#                   -0.305658966,   -0.713990152,    0.629898429,\
#                   -0.654916644,    0.637859344,    0.405215889,\
#                   -0.691106737,   -0.288665175,   -0.662581384,\
#                   0.001623333,    0.000937471,  -46.785041809,\
#                   46.941642761,   34.004577637,   31.637348175,\
#                   -13.265086174,  106.867988586,  -20.000000000] )
#    #S
#    if "S.pdb" in pdb:
#        print("PDB IS S.pdb")
#        cmd.set_view ([\
#                       -0.409751743,   -0.219807968,    0.885309339,\
#                       -0.638088524,    0.762622535,   -0.105983727,\
#                       -0.651861429,   -0.608328938,   -0.452747405,\
#                       0.001109511,    0.000494776,  -46.016716003,\
#                       47.185314178,   35.512237549,   32.690845490,\
#                       -14.064048767,  106.068992615,  -20.000000000] )
#
#    #C
#    if "C.pdb" in pdb:
#        print("PDB is C.pdb")
#        cmd.set_view ([\
#                        -0.258558095,   -0.645871639,    0.718325078,\
#                       -0.237979725,    0.763276696,    0.600630045,\
#                       -0.936205566,   -0.015643222,   -0.351061910,\
#                       0.001335760,   -0.000327773,  -39.497772217,\
#                       45.714607239,   35.260444641,   30.005773544,\
#                       -20.550725937,   99.582305908,  -20.000000000] )
#                      #                  -0.462704241,   -0.198452532,    0.864010990,\
#                  -0.555044293,    0.824802518,   -0.107798725,\
#                  -0.691241384,   -0.529438317,   -0.491791070,\
#                  0.000923563,    0.000070601,  -46.123336792,\
#                  46.376838684,   33.472648621,   30.712085724,\
#                  -13.914966583,  106.218093872,  -20.000000000] )
#
#   cmd.set_view ([-0.358574003,   -0.621906042,    0.696170390,\
#                  -0.551006794,    0.742994249,    0.379926175,\
#                  -0.753525376,   -0.247360140,   -0.609091520,\
#                  0.000761475,   -0.000055134,  -52.477943420,\
#                  45.892639160,   32.543830872,   30.253210068,\
#                  20.382640839,   84.580200195,  -20.000000000 ])
   # cmd.png(pdb+'.png')


#   cmd.set_view ([-0.313590229,    0.487615854,    0.814794183,\
#                 0.007288054,    0.859284937,   -0.511437297,\
#                 -0.949524581,   -0.154441506,   -0.273017645,\
#                 0.000000000,    0.000000000, -151.701156616,\
#                 33.086368561,   36.058006287,   31.844902039,\
#                 119.602340698,  183.799972534,  -20.000000000] )
#   
#

    print(shw)
    # Show clusters
    if shw is None:
        show_cluster(clusters)
    else:
        shw = [int(c) for c in shw.split(',')]
        for c in shw:
            show_cluster(clusters[c - 1])



if __name__ == '__main__':
    main()
