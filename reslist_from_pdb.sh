grep CA wt_local/ref.pdb  | grep ATOM | awk '{print :}' > reslist.txt
