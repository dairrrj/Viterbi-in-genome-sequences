# -*- coding: utf-8 -*-
import numpy as np
import sys
import re
import math

log=math.log
# read student-ID.py <input HMM file> <input fasta file>
# read HMM file 
hmm_name=sys.argv[1]
hmm_name=hmm_name.strip('<')
hmm_name=hmm_name.strip('>')
fid = open(hmm_name)
lines = fid.readlines()
line1 = lines[0].split(' ')
n_state = int(line1[0])
n_symbols = int(line1[1])
symbols = [] # four symbols ACGT
inip=[] # initial probability
em = [[0] * n_symbols for i in range(n_state)]# emission matrix
tm = [[0] * n_state for i in range(n_state)] # transition matrix  
for s in line1[2].strip('\n'):
    symbols.append(s)   
line2 = lines[1].split(' ')
for s in line2:
    inip.append(float(s))
line3 = lines[2].split(' ')
while '' in line3:
    line3.remove('')
    line4 = lines[3].split(' ')
while '' in line4:
    line4.remove('')
    
for i in range(n_state):
    tm[0][i] = float(line3[i])
    tm[1][i] = float(line4[i])
for j in range(n_symbols):
    em[0][j] = float(line3[j+n_state])
    em[1][j] = float(line4[j+n_state])
 
# read fasta file
fasta_name=sys.argv[2]
fasta_name=fasta_name.strip('<')
fasta_name=fasta_name.strip('>')
fid = open(fasta_name)
lines2 = fid.readlines()
observed=''
observe_seq = []
for ii,line in enumerate(lines2):
    line = line.strip('\n')
    g = re.match('^>', line)
    if g is None:
        observed = observed+lines2[ii].strip()
for o in observed:
    if o == "a" or o == "A":
        observe_seq.append(0)
    elif o == "c" or o == "C":
        observe_seq.append(1)
    elif o == "g" or o == "G":
        observe_seq.append(2)
    elif o == "t" or o == "T":
        observe_seq.append(3)
    
tm = np.array(tm)
em = np.array(em)
inip = np.array(inip)
L = len(observed)

# table initialization
max_p = np.zeros(shape=(L,n_state))
back_path = np.zeros(shape=(L,n_state))

for j in range(n_state):
    max_p[0][j] = log(inip[j])+log(em[j][observe_seq[0]])
    back_path[0][j] = j

# table filling
for k in range(1,L):
    for j1 in range(n_state):
        best = float("-inf")
        for j2 in range(n_state):
            nprob = max_p[k-1][j2] + log(tm[j2][j1]) + log(em[j1][observe_seq[k]])
            if nprob > best:
                best = nprob
                state = j2
        max_p[k][j1] = best
        back_path[k][j1] = state

# traceback
output = np.zeros(L)
output[L-1] = max_p[L-1,:].argmax()   
for x in range(L-2,-1,-1):

    j = back_path[x+1][int(output[x+1])]
    output[x]=j


# output the state information and number of state B segement
nB=0              
i=0 
file=r"result.position"
pos = []

while i<L:
    if output[i] == 0:
        j = i+1
        while j<L:
            if output[j]==1:
                pos.append([i+1,j,0])
                break 
            elif output[j]==0:
                j+=1   
        if j==L:
            pos.append([i+1,j,0])
        i=j     
            
    elif output[i] == 1:
        j = i+1
        while j<L:
            if output[j]==0:
                pos.append([i+1,j,1])
                break
            elif output[j]==1:
                j+=1 
        if j==L:
            pos.append([i+1,j,0])
        i=j
        
with open(file,'w') as f:
    for i in pos:
        if i[2]==0:
            f.write(str(i[0])+" " + str(i[1])+ " state A"+ "\n")
            print(str(i[0])+" " + str(i[1])+ " state A")
        else:
            f.write(str(i[0])+" " + str(i[1])+ " state B"+ "\n")
            print(str(i[0])+" " + str(i[1])+ " state A")
            nB+=1
    f.write("number of state B genome segement is %s" %nB)
            
print("number of state B genome segement is %s" %nB)


