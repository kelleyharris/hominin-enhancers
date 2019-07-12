import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from math import sqrt
from math import exp
import scipy.stats 
import numpy as np

neand_rare, control_rare, neand_common, control_common=dict({}),dict({}), dict({}), dict({})

infile=open('rare_vars_v_celltype_pleio.txt')
lines=infile.readlines()
infile.close()

pleio_bounds=[1,2,4,6,9,15,25]
tissues=[]
for line in lines[1:]:
    s=line.strip('\n').split(',')
    if int(s[0]) in pleio_bounds:
        tissues.append(s[0])
        t=tissues[-1]
        neand_rare[t]=int(s[1])
        neand_common[t]=int(s[3])+int(s[2])
        control_rare[t]=int(s[4])
        control_common[t]=int(s[6])+int(s[5])
    else:
        t=tissues[-1]
        neand_rare[t]+=int(s[1])
        neand_common[t]+=int(s[2])+int(s[3])
        control_rare[t]+=int(s[4])
        control_common[t]+=int(s[5])+int(s[6])


odds_info=[]
for tissue in tissues:
    print tissue
    chi2_info=scipy.stats.chi2_contingency(np.array([[neand_rare[tissue],neand_common[tissue]],[control_rare[tissue],control_common[tissue]]]))
    odds_info.append([1.0*(control_common[tissue])*neand_rare[tissue]/(control_rare[tissue]*(neand_common[tissue])),sqrt(1.0/neand_rare[tissue]+1.0/control_rare[tissue]+1.0/(control_common[tissue])+1.0/(neand_common[tissue])),tissue, chi2_info])

#odds_info.sort()
x,y,ymin,ymax=[],[],[],[]
xtotal,xlabel=[],[]

counter=0
for oi in odds_info:
    print oi[0], oi[0]*exp(-2*oi[1]), oi[0]*exp(2*oi[1]),oi[2], oi[-1]
    if xlabel==[]:
        xlabel.append('1')
    else:
        xlabel.append(str(oi[2])+'-'+str(pleio_bounds[counter+1]-1))
    y.append(oi[0])
    ymin.append(oi[0]-oi[0]*exp(-2*oi[1]))
    ymax.append(oi[0]*exp(2*oi[1])-oi[0])
    counter+=1
    x.append(counter)

plt.ylim(ymin=1.03)
plt.ylim(ymax=1.11)
plt.xlim(xmax=len(xlabel)+1)
xtotal=np.arange(1,len(xlabel)+1)
#plt.axhline(y=1,linestyle='dotted',color='grey')
plt.ylabel('Singleton enrichment')
plt.xticks(xtotal,xlabel,rotation='vertical')
plt.xlabel('Enhancer pleiotropy')
plt.errorbar(x,y,yerr=[ymin,ymax],color='blue',fmt='o')
fig=plt.gcf()
fig.subplots_adjust(bottom=0.55)
#fig.set_size_inches((3.5,4.5))
plt.savefig('oddssingle_v_celltype_pleio.png')
