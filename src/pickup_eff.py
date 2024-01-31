#*****************************************************************************
import sys
import numpy as np

# dependency on pytools

if __name__ == '__main__':

    inputfile = sys.argv[1]

    fmt = '%14.14e'

    with open(inputfile,'r') as finp:
        lines = finp.readlines()

    fout=open('pickup.pos','w')
    for line in lines: 
        picked=line.strip().split('-')
        #print(picked)
        fout.write(str(picked[0])+' '+str(picked[1])+' '+str(picked[2])+'\n')
    fout.close()

#    for i in range(L):
#        for j in range(L):
#            fout.write(str(i)+' '+str(j)+' '+str(fmt%mat[i,j])+'\n')
#    fout.close()


