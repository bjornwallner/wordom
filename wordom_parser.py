import re
import sys
import pandas as pd
import collections

def read_cross_corr(infile):
    dict=collections.OrderedDict()
    for line in open(infile,'r'):
        if line.startswith("#"):
            continue
        line = line.rstrip().lstrip()
        #print line
        #line = line.rstrip()
        (i,j,resi,resj,corr)=line.split()
        #remove chain
        resi = re.sub('^\w:\w', '', resi)
        resj = re.sub('^\w:\w', '', resj)
        #print resi
        if len(resi) == 2:
            resi = "0" + resi
        if len(resj) == 2:
            resj = "0" + resj
        if not resi in dict:
            dict[resi]=collections.OrderedDict()
        #corr_filtered=corr
        dict[resi][resj] = 0
        if float(corr) > 0.5:
            dict[resi][resj] = float(corr)

        #t=line.split()
        #print t
        #print line

    #dict2={} #collections.OrderedDict()
    ##for k1 in dict.keys():
    #    dict2[k1]={} #collections.OrderedDict()
    #    for k2 in dict.keys():
    #        dict2[k1][k2]=dict[k1][k2]
    df=pd.DataFrame.from_dict(dict,orient='index') #columns')
    ##print df
    #sys.exit()
    return(df)


# Read residue interaction strength from PSN analysis avg output files
def read_avg_strength(infile):
    m_start = re.compile("^\*\*\* Averaged Interaction Strength \*\*\*")
    m_end = re.compile("^========================")
    m_entry = re.compile("^\s*.:\w*\d+\s+.:\w*\d+\s+\d+\.\d+\s+\d+\.\d+\s*$")
    interactions = {}
    reading = False
    for line in infile:
        #print "B" + line
        if reading:
            # Stop reading if end of interaction strength section
            if m_end.search(line):
                break
            else:
                #print "IN search " + line
                if m_entry.search(line):
                 #   print "C" + line
                    [a, b, strength,seq_something] = line.split()
                  #  print a_tmp,b_tmp
#                    a=re.sub(r':\w',':',a_tmp)
#                    b=re.sub(r':\w',':',b_tmp)
                   # print a,b
                    if not a in interactions:
                        interactions[a] = {}
                    if not b in interactions:
                        interactions[b] = {}
                    # Assign symmetrically
                    interactions[a][b] = interactions[b][a] = float(strength)
        # Start reading when header found
        elif m_start.search(line):
            reading = True

    return interactions 



# Read clusters from PSN analysis avg output files
def read_avg_clusters(infile):
    m_start = re.compile("^\*\*\* Stable Cluster Compositions \*\*\*")
    m_new = re.compile("^Imin:")
    m_new_freq = re.compile("^Freq:")
    m_end = re.compile("^========================")
    m_entry = re.compile("^C\s*\d+:")
    reading = False
    clusters = {}
    current = None
    freq = None
    for line in infile:
#        print "B" + line
        if reading:
            if m_end.search(line):
                return clusters
            else:
                if m_entry.search(line) and freq:
                    entries = ":".join(line.split(':')[1:])
                    #cluster = int(line.split(':')[0][1:])
#                    if freq == 80:
#                        print "B " + line
                    clusters[current][freq].append(entries.split())
                        #   print "D " + line
                   # print current + " " + entries
                elif m_new.search(line):
                    current = round(float(line.split()[1]),1)
                 #   print(current)
                    clusters[current] = {}
                elif m_new_freq.search(line):
                    freq = round(float(line.split()[1]),1)
                    clusters[current][freq]=[]
        

        elif m_start.search(line):
            reading = True
    #sys.exit(1)
    return clusters
