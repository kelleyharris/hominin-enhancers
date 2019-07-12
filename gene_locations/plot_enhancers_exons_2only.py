import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from math import sqrt
from math import exp
import scipy.stats
import numpy as np

param_sets=[]
populs=['oceania','westeurasia','eastasia','southasia','centralasia','america']
for popul in populs:
    for oneortwo in ['2','1minus2']:
        for species in ['neandertal','denisova']:
            param_sets.append((oneortwo,species,popul))

pop_lw, species_color=dict({}),dict({})
popmarker=dict({})
for p in zip(populs,['purple','magenta','orange','red','blue','black']):
    species_color[p[0]]=p[1]
for p in zip(['1','2'],['orange','red']):
    species_color[p[0]]=p[1]
for p in zip(['oceania','westeurasia','southasia','eastasia','centralasia','america'],['solid','solid','dashdot','dashed','dotted','dotted']):
    pop_lw[p[0]]=p[1]
for p in zip(['oceania','westeurasia','southasia','eastasia','centralasia','america'],['o','<','v','>','x','^']):
    popmarker[p[0]]=p[1]


offset=dict({})
populs+=['denisova']
for i in range(len(populs)):
    offset[populs[i]]=-0.25+0.1*i
populs=['oceania','westeurasia','southasia','eastasia','centralasia','america']

"""
fncolor,spmarker=dict({}),dict({})
for fnpair in zip([('exon','2'),('exon','1minus2'),('enhancer','2'),('enhancer','1minus2')],['red','magenta','purple','blue']):
    fncolor[fnpair[0]]=fnpair[1]
for sppair in zip(['neandertal','denisova'],['o','x']):
    spmarker[sppair[0]]=sppair[1]
"""

x, y, xmin, xmax, ymin, ymax=dict({}), dict({}), dict({}), dict({}), dict({}),dict({})
#for t1 in ['exon','enhancer']:
#    for t2 in ['2','1minus2']:
for t1 in populs:
    for t2 in ['neandertal','denisova']:
        t=(t1,t2)
        x[t],y[t],xmin[t],xmax[t],ymin[t],ymax[t]=[],[],[],[],[],[]

xtotal,xlabel=[],[]
count_neand, count_control=dict({}),dict({})

marker=0
textheight=1.6
oneortwo='2'
for popul in populs:
    infile=open(popul+'_enhancer_exon_summary_counts.txt')
    lines=infile.readlines()
    infile.close()
    for line in lines[1:]:
        s=line.strip('\n').split(' ')
        count_neand[(s[0],s[1],s[2])]=int(s[3])
        count_control[(s[0],s[1],s[2])]=int(s[4])
    for fn in ['exon','enhancer']:
        for species in ['neandertal','denisova']:
            t=(popul,species)
            odds_info=(1.0*(count_control[(species,'total',oneortwo)]-count_control[(species,fn,oneortwo)])*count_neand[(species,fn,oneortwo)]/(count_control[(species,fn,oneortwo)]*(count_neand[(species,'total',oneortwo)]-count_neand[(species,fn,oneortwo)])),sqrt(1.0/count_neand[(species,fn,oneortwo)]+1.0/count_control[(species,fn,oneortwo)]+1.0/(count_neand[(species,'total',oneortwo)]-count_neand[(species,fn,oneortwo)])+1.0/(count_control[(species,'total',oneortwo)]-count_control[(species,fn,oneortwo)])))
            xtotal.append(len(xtotal))
#            print y.keys()
#            print y[t]
#            print offset[t[0]]
            y[t].append(len(y[t])+offset[t[0]])
            xlabel.append(fn+' ('+species.capitalize()+')')
            x[t].append(odds_info[0])
            print t, len(x[t])
            xmin[t].append(odds_info[0]*(1-exp(-2*odds_info[1])))
            xmax[t].append(odds_info[0]*(exp(2*odds_info[1])-1))
print x.keys()

"""
plt.axhline(y=1,linestyle='dotted',color='grey')
plt.xticks(xtotal,xlabel,rotation='vertical')
plt.xlim(xmin=-0.5)
plt.xlim(xmax=len(xtotal)-0.5)
for pos in [8,24]:
    plt.axvspan(pos-0.5,pos+7.5,facecolor='blue',alpha=0.25)
"""
#    plt.axvline(x=pos-0.5,color='grey')
#for pos in [0,4,8,12]:
#    plt.text(


#plt.text(0.5,textheight,'Oceania')
#plt.text(6.5,textheight,'West Eurasia')
#plt.text(12.5,textheight,'East Asia')
#plt.text(18.5,textheight,'South Asia')

poplabel=dict({})
for popul in ['oceania','westeurasia','southasia','eastasia','centralasia','america']:
#    param_sets.append((species,popul))
    if popul in ['oceania','america']:
        poplabel[popul]=popul.capitalize()
    elif popul=='westeurasia':
        poplabel[popul]='West Eurasia'
    else:
        poplabel[popul]=popul[:-4].capitalize()+' Asia'

for popul in populs:
    t=(popul,'neandertal')
    plt.errorbar(y[t],x[t],yerr=[xmin[t],xmax[t]],color=species_color[popul],linestyle='',label=poplabel[popul],marker=popmarker[popul],fmt='o')

t=('oceania','denisova')
y['denisova']=[offset['denisova'],offset['denisova']+1]
plt.errorbar(y['denisova'],x[t],yerr=[xmin[t],xmax[t]],color='grey',linestyle='',label=poplabel['oceania']+'(Denisova)',marker=popmarker[popul],fmt='o')

plt.xticks([0,1],['Exons','Enhancers'])
#plt.legend(loc='lower right',ncol=2)

#for t in x.keys():
#    plt.errorbar(y[t],x[t],yerr=[xmin[t],xmax[t]],fmt='o',color=fncolor[t])
#    print len(x[t])

plt.ylabel('Archaic variant depletion')
fig=plt.gcf()
fig.subplots_adjust(bottom=0.5)
plt.savefig('Exon_v_enhancer_depletion_only2.png')





