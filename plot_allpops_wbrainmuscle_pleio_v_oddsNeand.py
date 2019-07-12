import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from math import sqrt
from math import exp
import scipy.stats 
import numpy as np
import sys

binstarts=[1,2,4,6,9,14]
colors=['blue','cyan','red','magenta']

offset=0
for species in ['neandertal','denisova']:
    infile=open(species+'_allpops_2_snps_plus_bstat.txt')
    lines=infile.readlines()
    num_neand=len(lines)
    infile.close()

    infile=open('allpops_'+species+'_2_controls.txt')
    lines=infile.readlines()
    num_control=len(lines)
    infile.close()

    for insert in ['','_withfbrain_and_fmuscle']:
        neand_nonenhancer, control_nonenhancer, neand_enhancer, control_enhancer=dict({}),dict({}), dict({}), dict({})
    
        infile=open('allpop_analysis/allpops_'+species+'_2_counts'+insert+'_v_celltype_pleio_number.txt')
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
            x.append(len(x)+offset*0.125)
            if odds_info[ind+1][0]>oi[0]+1:
                xlabel.append(str(oi[0])+'-'+str(odds_info[ind+1][0]-1))
            else:
                xlabel.append(str(oi[0]))
            y.append(oi[1])
            ymin.append(oi[1]-oi[1]*exp(-2*oi[2]))
            ymax.append(oi[1]*exp(2*oi[2])-oi[1])
        if insert=='':
            plt.errorbar(x,y,yerr=[ymin,ymax],fmt='o',label=species.capitalize()+' (All enhancers)',color=colors[offset])
        else:
            plt.errorbar(x,y,yerr=[ymin,ymax],fmt='o',label=species.capitalize()+' (Active in fetal muscle or fetal brain)',color=colors[offset])
        offset+=1
            
plt.legend(loc='upper right',fontsize=8)
plt.xlim(xmin=-1)
plt.xlim(xmax=len(xlabel)+1)
plt.axhline(y=1,linestyle='dotted',color='grey')
plt.ylabel('Archaic-to-control variant odds ratio')
plt.xticks(x,xlabel)
plt.xlabel('Enhancer pleiotropy')
fig=plt.gcf()
fig.subplots_adjust(bottom=0.25)
fig.set_size_inches((5,5))
plt.savefig('archaic_v_pleio_allspecies_allpops_wbrainmuscle.png')
#plt.savefig(popul+'_'+species+oneortwo+'_analysis/'+popul+'_odds_'+species+oneortwo+'_v_pleio_allcelltypes_binned.png')
