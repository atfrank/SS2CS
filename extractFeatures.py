import pandas as pd
import numpy as np

def get_resname_int(resname):
    if resname == 'A' or resname == 'ADE':
        return 1
    elif resname == 'G' or resname == 'GUA':
        return 2
    elif resname == 'C' or resname == 'CYT':
        return 3
    elif resname == 'U' or resname == 'URA':
        return 4

def get_resname_char(resname):
    if resname == 1:
        return 'ADE'
    elif resname == 2:
        return 'GUA'
    elif resname == 3:
        return 'CYT'
    elif resname == 4:
        return 'URA'

def extractCT(ct_path,rna):
    # function to extract features from DSSR annotated secondary structure (.ct)
    df = pd.read_csv(ct_path, delim_whitespace=True,header=None,skiprows=1)
    length = len(df)
    df.columns = ['resid','resname','i_minus_1','i_plus_1','i_bp','resid2']

    #residue name of residue i
    i_resname = df['resname'].apply(lambda x: get_resname_int(x)) # Series
    i_resname_char = df['resname'].apply(lambda x: get_resname_char(get_resname_int(x))) 

    #residue that base paired to i is (resname instead of resid):
    i_bp = df['i_bp']
    i_bp_resname = []
    for i in range(length):
      if i_bp[i] == 0:
        i_bp_resname.append(0)
      else:
        i_bp_resname.append(get_resname_int(df[df['resid']==i_bp[i]]['resname'].values[0]))
    i_bp_resname = pd.Series(i_bp_resname)

    #residue that base paired to i-1 is:
    i_minus_1_bp = df['i_bp'].values.tolist()
    i_minus_1_bp = [0]+i_minus_1_bp # add 0 at the begining
    i_minus_1_bp.pop() # delete last element
    i_minus_1_bp_resname = []
    for i in range(length):
      if i_minus_1_bp[i] == 0:
        i_minus_1_bp_resname.append(0)
      else:
        i_minus_1_bp_resname.append(get_resname_int(df[df['resid']==i_minus_1_bp[i]]['resname'].values[0]))
    i_minus_1_bp_resname = pd.Series(i_minus_1_bp_resname)

    #residue that base paired to i+1 is:
    i_plus_1_bp = df['i_bp'].values.tolist()
    i_plus_1_bp = i_plus_1_bp+[0] # add 0 at the end
    i_plus_1_bp.pop(0) # delete last element
    i_plus_1_bp_resname = []
    for i in range(length):
      if i_plus_1_bp[i] == 0:
        i_plus_1_bp_resname.append(0)
      else:
        i_plus_1_bp_resname.append(get_resname_int(df[df['resid']==i_plus_1_bp[i]]['resname'].values[0]))
    i_plus_1_bp_resname = pd.Series(i_plus_1_bp_resname)

    #residue name of i-1
    i_minus_1_resname = i_resname.values.tolist()
    i_minus_1_resname = [0] + i_minus_1_resname
    i_minus_1_resname.pop()
    i_minus_1_resname = pd.Series(i_minus_1_resname)

    #residue name of i+1
    i_plus_1_resname = i_resname.values.tolist()
    i_plus_1_resname = i_plus_1_resname+[0]
    i_plus_1_resname.pop(0)
    i_plus_1_resname = pd.Series(i_plus_1_resname) 

    #the previous residue of the residue base paired to i: i_bp_minus_1
    i_bp_minus_1 = list(map(lambda x: x-1, i_bp.values.tolist()))
    i_bp_minus_1_resname = []
    for i in range(length):
      if i_bp_minus_1[i]==-1 or i_bp_minus_1[i]==0:
        i_bp_minus_1_resname.append(0)
      else:
        i_bp_minus_1_resname.append(get_resname_int(df[df['resid']==i_bp_minus_1[i]]['resname'].values[0]))
    i_bp_minus_1_resname = pd.Series(i_bp_minus_1_resname)

    #the next residue of the residue base paired to i: i_bp_plus_1
    i_bp_plus_1 = list(map(lambda x: x+1, i_bp.values.tolist()))
    i_bp_plus_1_resname = []
    for i in range(length):
      if i_bp_plus_1[i]==1 or i_bp_plus_1[i]==(length+1):
        i_bp_plus_1_resname.append(0)
      else:
        i_bp_plus_1_resname.append(get_resname_int(df[df['resid']==i_bp_plus_1[i]]['resname'].values[0]))  
    i_bp_plus_1_resname = pd.Series(i_bp_plus_1_resname)

    #i_bp_prev_bp
    i_bp_minus_1_bp_resname = []
    for i in range(length):
      if i_bp_minus_1[i]==-1 or i_bp_minus_1[i]==0:
        i_bp_minus_1_bp_resname.append(0)
      else:
        i_bp_minus_1_bp = df[df['resid']==i_bp_minus_1[i]]['i_bp'].values[0]
        if i_bp_minus_1_bp == 0:
          i_bp_minus_1_bp_resname.append(0)
        else:
          i_bp_minus_1_bp_resname.append(get_resname_int(df[df['resid']==i_bp_minus_1_bp]['resname'].values[0]))
    i_bp_minus_1_bp_resname = pd.Series(i_bp_minus_1_bp_resname)

    #i_bp_next_bp
    i_bp_plus_1_bp_resname = []
    for i in range(length):
      if i_bp_plus_1[i]==1 or i_bp_plus_1[i]==(length+1):
        i_bp_plus_1_bp_resname.append(0)
      else:
        i_bp_plus_1_bp = df[df['resid']==i_bp_plus_1[i]]['i_bp'].values[0]
        if i_bp_plus_1_bp == 0:
          i_bp_plus_1_bp_resname.append(0)
        else:
          i_bp_plus_1_bp_resname.append(get_resname_int(df[df['resid']==i_bp_plus_1_bp]['resname'].values[0]))  
    i_bp_plus_1_bp_resname = pd.Series(i_bp_plus_1_bp_resname) 

    # create other information as Series
    rnaid = pd.Series([rna]*length)
    total_length = pd.Series([length]*length)
    resid = df['resid']
    resname = df['resname']

    #print(features)
    features = pd.concat([rnaid, total_length, resid, resname, i_resname, i_minus_1_resname,
    i_plus_1_resname, i_bp_resname, i_minus_1_bp_resname, i_plus_1_bp_resname, i_bp_minus_1_resname,
    i_bp_plus_1_resname, i_bp_minus_1_bp_resname, i_bp_plus_1_bp_resname], axis=1)
    features.columns = ["id","length","resid","i_resname_char","i_resname","i_minus_1_resname","i_plus_1_resname",
                        "i_bp_resname", "i_minus_1_bp_resname", "i_plus_1_bp_resname", "i_bp_minus_1_resname",
                        "i_bp_plus_1_resname", "i_bp_minus_1_bp_resname", "i_bp_plus_1_bp_resname"]
    return features
  