#*****************************************************************************
import sys
import numpy as np

# dependency on pytools

if __name__ == '__main__':

    inputfile = sys.argv[1]

    fmt = '%14.14e'

    with open(inputfile,'r') as finp:
        lines = finp.readlines()

    fout=open('pickup.vec','w')
    for line in lines: 
        picked=line.strip().split('_')[1]
        #print(picked)
        fout.write(str(picked) + '\n')
    fout.close()

#    for i in range(L):
#        for j in range(L):
#            fout.write(str(i)+' '+str(j)+' '+str(fmt%mat[i,j])+'\n')
#    fout.close()


