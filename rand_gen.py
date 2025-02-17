import random
import pandas as pd
import os
import math
import numpy as np

os.getcwd()


def gen_meta_seq(KR, AILV, FWY, Length):
    # Define groups
    group_KR = '(KR)'
    group_AILV = '(ILV)'
    group_FWY = '(FW)'
    group_NQSTG = '(QTS)'

    while True:
        # Randomly shuffle positions
        positions = list(range(Length))
        random.shuffle(positions)

        # Assign positions for each group
        kr_positions = positions[:KR]
        positions = positions[KR:]

        ailv_positions = positions[:AILV]
        positions = positions[AILV:]

        fwy_positions = positions[:FWY]
        positions = positions[FWY:]

        nqstg_positions = positions

        # Create a sequence list with placeholders
        sequence = [''] * Length

        # Fill in the placeholders with the corresponding group
        for pos in kr_positions:
            sequence[pos] = group_KR
        for pos in ailv_positions:
            sequence[pos] = group_AILV
        for pos in fwy_positions:
            sequence[pos] = group_FWY
        for pos in nqstg_positions:
            sequence[pos] = group_NQSTG

        # Join the sequence into a string
        metaseq = ''.join(sequence)
        
        if KR>3:
        # Check for consecutive or interspersed AILV and FWY groups
            if not has_consecutive_or_interspersed(metaseq, group_AILV, group_FWY, 4):
                return metaseq
        else:
            if not has_consecutive_or_interspersed(metaseq, group_AILV, group_FWY, 5):
                return metaseq            

def has_consecutive_or_interspersed(sequence, group1, group2, limit):
    count = 0
    for elem in sequence.split(')'):
        elem += ')'
        if elem in [group1, group2]:
            count += 1
            if count >= limit:
                return True
        else:
            count = 0
    return False

def genmetaseqs(KR, AILV, FWY,Length,number):
    i=0
    set=[]
    for i in range(int(number)):
        #阳离子 脂肪族 芳香族 
        set.append(gen_meta_seq(KR, AILV, FWY, Length))
        i+=1
    return set

def sixmeta(net,num):
    set=[]                
    if net==2:
        for AILV in range(1,5):
            for FWY in range(0,2):
                set.append(genmetaseqs(2,AILV,FWY,6,num))
                        
    if net==3:
        for AILV in range(1,4):
            for FWY in range(0,2):
                set.append(genmetaseqs(3,AILV,FWY,6,num))
                
    if net==4:
        for AILV in range(1,3):
            for FWY in range(0,2):
                set.append(genmetaseqs(4,AILV,FWY,6,num))
    
    seqs = sum(set,[]) 
#     print(len(seqs)) 
    return seqs

def eightmeta(net,num):
    set=[]                
    if net==2:
        for AILV in range(2,7):
            for FWY in range(0,2):
                set.append(genmetaseqs(2,AILV,FWY,8,num))
                        
    if net==3:
        for AILV in range(1,6):
            for FWY in range(0,2):
                set.append(genmetaseqs(3,AILV,FWY,8,num))
                
    if net==4:
        for AILV in range(1,5):
            for FWY in range(0,2):
                set.append(genmetaseqs(4,AILV,FWY,8,num))
    
    seqs = sum(set,[]) 
#     print(len(seqs)) 
    return seqs

def tenmeta(net,num):
    set=[]                
    if net==2:
        for AILV in range(4,9):
            for FWY in range(0,3):
                set.append(genmetaseqs(2,AILV,FWY,10,num))
                        
    if net==3:
        for AILV in range(3,8):
            for FWY in range(0,3):
                set.append(genmetaseqs(3,AILV,FWY,10,num))
                
    if net==4:
        for AILV in range(2,7):
            for FWY in range(0,3):
                set.append(genmetaseqs(4,AILV,FWY,10,num))
                
    if net==5:
        for AILV in range(1,6):
            for FWY in range(0,3):
                set.append(genmetaseqs(5,AILV,FWY,10,num))  
                
    seqs = sum(set,[]) 
#     print(len(seqs)) 
    return seqs

def thirteenmeta(net,num):
    set=[]

    if net==2:
        for AILV in range(7,11):
            for FWY in range(0,2):
                set.append(genmetaseqs(2,AILV,FWY,13,num))
                        
    if net==3:
        for AILV in range(6,10):   
            for FWY in range(0,2):
                set.append(genmetaseqs(3,AILV,FWY,13,num))
                
    if net==4:
        for AILV in range(5,9):
            for FWY in range(0,2):
                set.append(genmetaseqs(4,AILV,FWY,13,num))
                        
    if net==5:
        for AILV in range(4,8):
            for FWY in range(0,2):
                set.append(genmetaseqs(5,AILV,FWY,13,num))
            
    seqs = sum(set,[]) 
#     print(len(seqs)) 
    return seqs

def merge_meta(num,net,l):
    df6_2 = None
    df6_3 = None
    df8_2 = None
    df8_3 = None
    df8_4 = None
    df10_2 = None
    df10_3 = None
    df10_4 = None
    df10_5 = None
    df13_2 = None
    df13_3 = None
    df13_4 = None
    df13_5 = None

    if l == 6:
        if net ==2:
            df6_2 = sixmeta(2, num)
        elif net==3:
            df6_3 = sixmeta(3, num)
    elif l == 8:
        if net ==2:
            df8_2 = eightmeta(2, num)
        elif net==3:
            df8_3 = eightmeta(3, num)
        elif net ==4:
            df8_4 = eightmeta(4, num)
    elif l == 10:
        if net ==2:
            df10_2 = tenmeta(2, num)
        elif net==3:
            df10_3 = tenmeta(3, num)
        elif net ==4:
            df10_4 = tenmeta(4, num)
        elif net ==5:
            df10_5 = tenmeta(5, num)
    
    elif l == 13:
        if net ==2:
            df13_2 = thirteenmeta(2, num)
        elif net==3:
            df13_3 = thirteenmeta(3, num)
        elif net ==4:
            df13_4 = thirteenmeta(4, num)
        elif net ==5:
            df13_5 = thirteenmeta(5, num)
    
    seqall = []
    if df6_2 is not None:
        seqall += df6_2
    if df6_3 is not None:
        seqall += df6_3
    if df8_2 is not None:
        seqall += df8_2
    if df8_3 is not None:
        seqall += df8_3
    if df8_4 is not None:
        seqall += df8_4
    if df10_2 is not None:
        seqall += df10_2
    if df10_3 is not None:
        seqall += df10_3
    if df10_4 is not None:
        seqall += df10_4
    if df10_5 is not None:
        seqall += df10_5

    if df13_2 is not None:
        seqall += df13_2
    if df13_3 is not None:
        seqall += df13_3
    if df13_4 is not None:
        seqall += df13_4
    if df13_5 is not None:
        seqall += df13_5

    seqall = list(set(seqall))
#     print(len(seqall))
    seqall = [seq for seq in set(seqall) if (seq[0] not in "PST") and (seq[-1] not in "FWY")]
    df_all = pd.DataFrame(seqall,columns=['Sequence'])
    
#     df_all.to_csv('Rand_gen_seq/rand_gen_seq.csv',index=False,header=True)
    return df_all



delta = 100
def acquireH(amino):
    if amino == 'A':
        H = 0.310
        return H
    if amino == 'C':
        H = 1.540
        return H
    if amino == 'D':
        H = -0.770
        return H
    if amino == 'E':
        H = -0.640
        return H
    if amino == 'F':
        H = 1.790
        return H
    if amino == 'G':
        H = 0.000
        return H
    if amino == 'H':
        H = 0.130
        return H
    if amino == 'I':
        H = 1.800
        return H
    if amino == 'K':
        H = -0.990
        return H
    if amino == 'L':
        H = 1.700
        return H
    if amino == 'M':
        H = 1.230
        return H
    if amino == 'N':
        H = -0.600
        return H
    if amino == 'P':
        H = 0.720
        return H
    if amino == 'Q':
        H = -0.220
        return H
    if amino == 'R':
        H = -1.010
        return H
    if amino == 'S':
        H = -0.040
        return H
    if amino == 'T':
        H = 0.260
        return H
    if amino == 'V':
        H = 1.220
        return H
    if amino == 'W':
        H = 2.250
        return H
    if amino == 'Y':
        H = 0.960
        return H

    return -1

def computeMuH(sequence):
    miuH={}
    H_sin_sum = 0
    H_cos_sum = 0
    N = len(sequence)
    space = 0
    amino_error = False
    for i in range(N):
        amino = sequence[i]
        if amino == ' ':
            space += 1
        H = acquireH(amino)
        if H == -1:
            amino_error = True
            break
        ndelta = math.pi * ((i + 1) * delta / 180)
        H_sin = math.sin(ndelta) * H
        H_cos = math.cos(ndelta) * H
        H_sin_sum += H_sin
        H_cos_sum += H_cos
    if amino_error == True:
        return ''
    H_sum = math.pow(H_sin_sum, 2) + math.pow(H_cos_sum, 2)
    muH = math.pow(H_sum, 0.5) * (1.0 / (N - space))
#     miuH['Sequence'] = sequence
#     miuH['miuH'] = muH

    return muH

def computeH(sequence):
        H_sin_sum = 0
        H_cos_sum = 0
        N = len(sequence)
        space = 0
        amino_error = False
        sum_ = 0
        for i in range(N):
            amino = sequence[i]
            if amino == ' ':
                space += 1
            H = acquireH(amino)
            if H == -1:
                amino_error = True
                break
            sum_ += H
        H_result = sum_ / N
        
        return H_result
    
def miuHseq(seq):        
    return computeMuH(seq)
def Hseq(seq):        
    return computeH(seq)

def miuH(seqs):
    H=[]
    b=[]
    for i in seqs:
        d=computeMuH(i)
        c=d.copy()
        b.append(c)
        a=pd.DataFrame(b)
    return a

def select_subseq_H(data=None,H=None):
    data['H']=data['Sequence'].map(Hseq)
    subdf = data[data.H>=H]
    return subdf

def Length(x):
    return len(x)

def RKHnum(x):
    num = 0
    for i in x:
        if i == 'K' or i == 'R'or i == 'H':
            num += 1
    return num

def select_seq(data=None,seq_length=None,RKH=None,limit=None,H=None):
    data = select_subseq_H(data=data,H=H)
    print('after select subseq H ',len(data))
    data.loc[:,'Length']=data['Sequence'].map(Length)
    print('After calculating Length:', len(data))
    data.loc[:,'miuH']=data['Sequence'].map(miuHseq)
    print('After calculating miuH:', len(data))
    data.loc[:,'RKH']=data['Sequence'].map(RKHnum)
    print('After RKH calculation:', len(data))
    miuH_max = data[(data.Length==seq_length)&(data.RKH==RKH)]['miuH'].max(axis=0)
    print('Max miuH:', miuH_max)
    miuH_min = miuH_max - limit
    print('Min miuH threshold:', miuH_min)

    subdf = data[(data.Length==seq_length)&(data.RKH==RKH)]
    print('After final selection:', len(subdf))
  
    return subdf



