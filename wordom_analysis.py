#!/usr/bin/env python

import os
import sys
import glob
import argparse
import re
import shutil
import tempfile

parser = argparse.ArgumentParser(description='Running wordom analysis')
parser.add_argument('-traj', default=None,type=str, nargs='+', required=True,
                    help='trajectories for analysis')
parser.add_argument('-name', default='wordom',type=str,
                    help='name of analysis, will used in some output files')
parser.add_argument('-prefix', default='',type=str,
                    help='the prefix will be added to outfiles')
parser.add_argument('-start',  default=None,type=float, required=True,
                    help='start time')
parser.add_argument('-end', default=None,type=float,  required=True,
                    help='end frame')
parser.add_argument('-ref',  metavar= 'REF.pdb', default=None,type=str,
                    help='Reference pdb for output visualisation')
parser.add_argument('-psn',  action='store_true',
                    help='Do protein structure network analysis, requires -ref REF.pdb')
parser.add_argument('-psnpath',  action='store_true',
                    help='Do PSNPATH analysis, requires -ref REF.pdb')
parser.add_argument('-corr',  action='store_true',
                    help='Do cross correlation analysis, requires -ref REF.pdb')
parser.add_argument('-align',  action='store_true',
                    help='align trajectories, requires -ref REF.pdb')
parser.add_argument('-timestep', default=0.2 ,type=float, 
                    help='time step per frame in nanoseconds, works for the skip100')

args = parser.parse_args()
if not args.ref and (args.psn or args.psnpath or args.corr or args.align):
    sys.stderr.write("-ref REF.pdb need for the choosen analysis\n")
    sys.exit()

if args.ref and not os.path.exists(args.ref):
    sys.stderr.write("Can not find -ref {}\n".format(args.ref))
    sys.exit()
#print args.traj
print args
#for infile in glob.glob("test/*.skip100.xtc")



start_frame=int(args.start/args.timestep)
end_frame=int(args.end/args.timestep)
#print args.start,start_frame
#s=args.start
#print "start time: {} corresponds to frame {}".format("hej","hej")
print "start time: {} ns corresponds to frame {} assuming {} ns/frame".format(args.start,start_frame,args.timestep)
print "end time: {} ns corresponds to frame {} assuming {} ns/frame".format(args.end,end_frame,args.timestep)




trjs=[]
print "Converting individual trajectories..."
for intraj in args.traj:
    if os.path.exists(intraj):
        outtraj=re.sub(r'.xtc$','.'+str(start_frame)+'-'+str(end_frame)+'.dcd',intraj)
        print "{} -> {}".format(intraj,outtraj)
        cmd = "wordom -F range -beg {} -end {}  -itrj {} -otrj {}".format(start_frame,end_frame,intraj,outtraj)
        if not os.path.exists(outtraj):
            print cmd
            os.system(cmd)
        trjs.append(outtraj)
    else:
        sys.stderr.write("Could not find {}, skipping...\n".format(intraj))

print "{} trajectories converted".format(len(trjs))
if len(trjs) > 0:    
    print "Merging indiviual trajectories to one large trajectory file..."
    outfile="{}{}.{}-{}.dcd".format(args.prefix,"merged",start_frame,end_frame)
    print outfile
    shutil.copy(trjs[0],outfile)
    for i in range(1,len(trjs)):
        cmd = "wordom -atrj {} -otrj {}".format(trjs[i],outfile)
        print cmd
        os.system(cmd)

    cmd="wordom -F range -beg {} -end {}  -itrj {} -otrj {}".format(start_frame,end_frame,intraj,outtraj)
    outdir="output_{}_{}_{}/".format(args.name,start_frame,end_frame)
    psn_raw=outdir +"rawpsn"
    if args.psn or args.psnpath:
        print "Do PSN"
        outdir="output_{}_{}_{}/".format(args.name,start_frame,end_frame)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        dirpath = tempfile.mkdtemp()
        print dirpath
        print outdir
        psn_input=args.prefix+'psn.inp'
        f=open(psn_input,'w')
        f.write("""BEGIN psn
--TITLE psn
--SELE /*/*/* 
--INTMIN 0.0:10.0:0.5 
--DISTCUTOFF 4.5 
--STABLECUTOFF 0.5 
--HUBCONTCUTOFF 3 
--TERMINI 0 
--VERBOSE 1 
END
""")
        f.close()
        psn_input_auto=args.prefix+'psn.auto.inp'
        f=open(psn_input_auto,'w')
        f.write("""BEGIN psn
--TITLE psn
--SELE /*/*/* 
--INTMIN AUTO 
--DISTCUTOFF 4.5 
--STABLECUTOFF 0.5 
--HUBCONTCUTOFF 3 
--TERMINI 0 
--VERBOSE 1 
END
""")
        f.close()
        
        cmd = "cd {};wordom -iA {} -itrj {} -imol {} -nopbc".format(dirpath,os.path.abspath(psn_input),os.path.abspath(outfile),os.path.abspath(args.ref))
        print cmd
        #sys.exit()
        os.system(cmd)
        os.system("mv {}/* {}".format(dirpath,outdir))
        os.rmdir(dirpath)
        print "Output from PSN written to {}".format(outdir)

    if args.align:
        outfile_aligned=re.sub(r'.dcd$','.alignd.dcd',outfile)
        cmd = 'wordom -ia rmsd --TITLE rmsd1 --SELE "/*/*/CA" --TRJOUT {} -imol {} -itrj {}'.format(outfile_aligned,args.ref, outfile)
        print cmd
        os.system(cmd)
    if args.corr or args.psnpath:
        outfile_aligned=re.sub(r'.dcd$','.alignd.dcd',outfile)
        cmd = 'wordom -ia rmsd --TITLE rmsd1 --SELE "/*/*/CA" --TRJOUT {} -imol {} -itrj {}'.format(outfile_aligned,args.ref, outfile)
        print cmd
        os.system(cmd)
        cmd = 'wordom -ia corr --TITLE cross-corr --SELE "/*/*/*" --TYPE DCC --LEVEL RES -imol {} -itrj {} -nopbc'.format(args.ref,outfile_aligned)
        print cmd
        os.system(cmd)
        os.system("mv cross-corr {}".format(outdir))

#    if args.psnpath:
        
        
        
    #wordom -ia corr --TITLE cross-corr --SELE "/*/*/*" --TYPE DCC --LEVEL RES -imol ref.pdb -itrj whole_traj.127.aligned.dcd -nopbc

else:
    print "No trajectories found..."
    




