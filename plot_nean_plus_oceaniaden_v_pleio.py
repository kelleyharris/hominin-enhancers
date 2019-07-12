import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from math import sqrt
from math import exp
import scipy.stats 
import numpy as np
import sys

param_sets=[]
poplabel=dict({})

poplabel['denisova']=' (D)'
poplabel['neandertal']=' (N)'

poplabel['1']=' (Distal)'
poplabel['2']=' (Proximal)'

offset=dict({})
populs=['oceania','westeurasia','southasia','eastasia','centralasia','america']
for i in range(len(populs)):
    offset[populs[i]]=-0.25+0.1*i


for species in ['denisova','neandertal']:
    for oneortwo in ['1','2']:
        for popul in ['oceania','westeurasia','southasia','eastasia','centralasia','america']:
            param_sets.append((str(oneortwo),species,popul))
            if popul in ['oceania','america']:
                poplabel[(str(oneortwo),popul)]=popul.capitalize()
            elif popul=='westeurasia':
                poplabel[(oneortwo,popul)]='West Eurasia'
            else:
                poplabel[(oneortwo,popul)]=popul[:-4].capitalize()+' Asia'

binstarts=[1,2,4,6,9,14]
pop_lw, species_color=dict({}),dict({})
popmarker=dict({})
species_color[('neandertal','2')]='blue'
species_color[('neandertal','1')]='cyan'
species_color[('denisova','2')]='red'
species_color[('denisova','1')]='magenta'

for p in zip(populs,['purple','magenta','orange','red','blue','black']):
    species_color[p[0]]=p[1]
for p in zip(['1','2'],['orange','red']):
    species_color[p[0]]=p[1]
for p in zip(['oceania','westeurasia','southasia','eastasia','centralasia','america'],['solid','solid','dashdot','dashed','dotted','dotted']):
    pop_lw[p[0]]=p[1]
for p in zip(['oceania','westeurasia','southasia','eastasia','centralasia','america'],['o','<','v','>','x','^']):
    popmarker[p[0]]=p[1]

for (oneortwo,species,popul) in param_sets:
    if oneortwo=='2':
        infile=open(species+'_'+popul+oneortwo+'_snps.txt')
    else:
        infile=open(popul+'_'+species+'_1minusboth2_snps.txt')
    lines=infile.readlines()
    num_neand=len(lines)
    infile.close()

    if oneortwo=='2':
        infile=open(species+'_'+popul+oneortwo+'_controls.txt')
    else:
        infile=open(popul+'_'+species+'_1minus2_controls.txt')
    lines=infile.readlines()
    num_control=len(lines)
    infile.close()

    neand_nonenhancer, control_nonenhancer, neand_enhancer, control_enhancer=dict({}),dict({}), dict({}), dict({})

    if oneortwo=='2':
        infile=open(popul+'_'+species+oneortwo+'_analysis/'+popul+'_'+species+oneortwo+'_counts_v_celltype_pleio_number.txt')
    else:
        infile=open('bugtest_1minus2/'+popul+'_'+species+'_1minus2_counts_v_celltype_pleio_number.txt')
    lines=infile.readlines()
    infile.close()

    pleio_vals=[]
    neand_count=[]
    control_count=[]
    for i in range(23):
        neand_count.append(0)
        control_count.append(0)

    for line in lines[1:]:
        s=line.strip('\n').split(',')
        if not int(s[0]) in pleio_vals:
            pleio_vals.append(int(s[0]))
            p=pleio_vals[-1]
            neand_enhancer[p]=int(s[2])
            control_enhancer[p]=int(s[3])
        else:
            p=int(s[0])
            neand_enhancer[p]+=int(s[2])
            control_enhancer[p]+=int(s[3])
        control_count[int(s[0])-1]+=int(s[3])
        neand_count[int(s[0])-1]+=int(s[2])

    odds_info=[]
    ind0,ind1=0,1
    while ind1<len(pleio_vals):
        if pleio_vals[ind1] in binstarts:
            p=pleio_vals[ind0]
            chi2_info=scipy.stats.fisher_exact(np.array([[neand_enhancer[p],num_neand-neand_enhancer[p]],[control_enhancer[p],num_control-control_enhancer[p]]]))
            odds_info.append([p,1.0*(num_control-control_enhancer[p])*neand_enhancer[p]/(control_enhancer[p]*(num_neand-neand_enhancer[p])),sqrt(1.0/neand_enhancer[p]+1.0/control_enhancer[p]+1.0/(num_control-control_enhancer[p])+1.0/(num_neand-neand_enhancer[p])), chi2_info])
            ind0=ind1
            ind1=ind0+1
        else:
            neand_enhancer[pleio_vals[ind0]]+=neand_enhancer[pleio_vals[ind1]]
            control_enhancer[pleio_vals[ind0]]+=control_enhancer[pleio_vals[ind1]]
            ind1+=1
    p=pleio_vals[ind0]
    chi2_info=scipy.stats.fisher_exact(np.array([[neand_enhancer[p],num_neand-neand_enhancer[p]],[control_enhancer[p],num_control-control_enhancer[p]]]))
    odds_info.append([p,1.0*(num_control-control_enhancer[p])*neand_enhancer[p]/(control_enhancer[p]*(num_neand-neand_enhancer[p])),sqrt(1.0/neand_enhancer[p]+1.0/control_enhancer[p]+1.0/(num_control-control_enhancer[p])+1.0/(num_neand-neand_enhancer[p])), chi2_info])

    odds_info.sort()
    x, y, ymin, ymax=[],[],[],[]
    xtotal,xlabel=[],[]
    odds_info.append([pleio_vals[-1]])

    for ind in range(len(odds_info)-1):
        oi=odds_info[ind]
#for oi in odds_info:
        x.append(len(x)+offset[popul])
        if odds_info[ind+1][0]>oi[0]+1:
            xlabel.append(str(oi[0])+'-'+str(odds_info[ind+1][0]-1))
        else:
            xlabel.append(str(oi[0]))
        y.append(oi[1])
        ymin.append(oi[1]-oi[1]*exp(-2*oi[2]))
        ymax.append(oi[1]*exp(2*oi[2])-oi[1])
    if species=='denisova' and popul=='oceania' and oneortwo=='2':
        plt.errorbar(x,y,yerr=[ymin,ymax],color='grey',fmt='o',linestyle=pop_lw[popul],label=poplabel[(oneortwo,popul)]+' (Denisova)',marker=popmarker[popul])
    elif species=='neandertal' and oneortwo=='2':
        if popul=='oceania':
            poplabel[(oneortwo,popul)]+=' (Neandertal)'
        plt.errorbar(x,y,yerr=[ymin,ymax],color=species_color[popul],fmt='o',linestyle=pop_lw[popul],label=poplabel[(oneortwo,popul)],marker=popmarker[popul])

plt.legend(loc='upper right',fontsize=10)

plt.axhline(y=1,linestyle='dotted',color='grey')
plt.ylabel(species.capitalize()+'-to-control variant odds ratio')
plt.xticks(x,xlabel)
plt.xlabel('Enhancer pleiotropy')
fig=plt.gcf()
fig.subplots_adjust(bottom=0.25)
fig.set_size_inches((5,5))
plt.tight_layout()
plt.savefig('archaic_v_pleio_neand_plus_oceaniaden.png')
#plt.savefig(popul+'_'+species+oneortwo+'_analysis/'+popul+'_odds_'+species+oneortwo+'_v_pleio_allcelltypes_binned.png')
