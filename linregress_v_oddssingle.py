import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from math import sqrt
from math import exp
import scipy.stats 
import numpy as np
from scipy.stats import linregress
param_sets=[]
for oneortwo in ['1','2']:
    for species in ['neandertal','denisova']:
        for popul in ['oceania','westeurasia','eastasia','southasia']:
            param_sets.append((oneortwo,species,popul))

tcolor=dict({})
for cpair in [('Fetal','red'),('Adult','blue'),('Embryonic','orange'),('Exons','green'),('HARs','black')]:
    tcolor[cpair[0]]=cpair[1]

for (oneortwo,species,popul) in param_sets:
    neand_nonenhancer, control_nonenhancer, neand_enhancer, control_enhancer=dict({}),dict({}), dict({}), dict({})
    infile=open(popul+'_'+species+oneortwo+'_analysis/exon_har_variant_counts.txt')
    lines=infile.readlines()
    infile.close()
    s=lines[1].split(' ')
    num_neand,num_control=int(s[1]),int(s[2])
    s=lines[2].split(' ')
    neand_enhancer[('Exons','Exons')],control_enhancer[('Exons','Exons')]=int(s[1]),int(s[2])
    s=lines[3].split(' ')
    neand_enhancer[('HARs','HARs')],control_enhancer[('HARs','HARs')]=int(s[1]),int(s[2])

    infile=open(popul+'_'+species+oneortwo+'_analysis/'+popul+'_'+species+oneortwo+'_counts_alt_cell_types.txt')
    lines=infile.readlines()
    infile.close()

    tissues=[]
    for line in lines[2:]:
        s=line.strip('\n').split(',')
        if s[2] in ['ESC','iPSC']:
            s[1]='Embryonic'
        if not (s[1],s[2]) in tissues and not s[2][:2]=='Un':
            tissues.append((s[1],s[2]))
            t=tissues[-1]
            neand_enhancer[t]=int(s[4])
            control_enhancer[t]=int(s[5])
        else:
            t=(s[1],s[2])
            neand_enhancer[t]+=int(s[4])
            control_enhancer[t]+=int(s[5])


    odds_info_neand=dict({})
    for tissue in tissues:
        print tissue
        odds_info_neand[tissue]=[1.0*(num_control-control_enhancer[tissue])*neand_enhancer[tissue]/(control_enhancer[tissue]*(num_neand-neand_enhancer[tissue])),sqrt(1.0/neand_enhancer[tissue]+1.0/control_enhancer[tissue]+1.0/(num_control-control_enhancer[tissue])+1.0/(num_neand-neand_enhancer[tissue]))]

    enhancer_common, control_common, enhancer_rare, control_rare=dict({}),dict({}), dict({}), dict({})

    for t in tissues:
        enhancer_common[t],control_common[t],enhancer_rare[t],control_rare[t]=0,0,0,0
    enhancer_common[('HARs','HARs')],control_common[('HARs','HARs')],enhancer_rare[('HARs','HARs')],control_rare[('HARs','HARs')]=6130,6292,3745,3655
    enhancer_common[('Exons','Exons')],control_common[('Exons','Exons')],enhancer_rare[('Exons','Exons')],control_rare[('Exons','Exons')]=463826, 497059, 378496, 320136

    infile=open('../rare_common_vars_alt_cell_types_enhancer_control.txt')
    lines=infile.readlines()
    infile.close()

    tissues=[]
    for line in lines[2:]:
        s=line.strip('\n').split(',')
        if s[2] in ['ESC','iPSC']:
            s[1]='Embryonic'
        if not (s[1],s[2]) in tissues and not s[2][:3]==' Un':
            tissues.append((s[1],s[2]))
            t=tissues[-1]
            enhancer_rare[t]=int(s[3])
            enhancer_common[t]=int(s[4])+int(s[5])
            control_rare[t]=int(s[6])
            control_common[t]=int(s[7])+int(s[8])
        elif not s[2][:3]==' Un':
            t=(s[1],s[2])
            enhancer_rare[t]+=int(s[3])
            enhancer_common[t]+=int(s[4])+int(s[5])
            control_rare[t]+=int(s[6])
            control_common[t]+=int(s[7])+int(s[8])

    x, y, xmin, xmax, ymin, ymax=dict({}), dict({}), dict({}), dict({}), dict({}),dict({})
    x_nobrain, y_nobrain=dict({}),dict({})
    for t in ['Adult','Fetal','Embryonic','HARs','Exons']:
        x[t],y[t],xmin[t],xmax[t],ymin[t],ymax[t]=[],[],[],[],[],[]
        x_nobrain[t],y_nobrain[t]=[],[]

    print tissues
    allneand,allsfs=[],[]
    for t in tissues:
        x[t[0]].append(odds_info_neand[t][0])
        allneand.append(odds_info_neand[t][0])
        xmin[t[0]].append(odds_info_neand[t][0]*(1-exp(-2*odds_info_neand[t][1])))
        xmax[t[0]].append(odds_info_neand[t][0]*(exp(2*odds_info_neand[t][1])-1))
        omean, ovar=1.0*(enhancer_rare[t]*control_common[t])/(control_rare[t]*enhancer_common[t]),sqrt(1.0/(enhancer_rare[t])+1.0/(control_common[t])+1.0/(control_rare[t])+1.0/(enhancer_common[t]))
        y[t[0]].append(omean)
        allsfs.append(omean)
        ymin[t[0]].append(omean*(1-exp(-2*ovar)))
        ymax[t[0]].append(omean*(exp(2*ovar)-1))
        print t[0],t[1],y[t[0]][-1],y[t[0]][-1]-ymin[t[0]][-1],y[t[0]][-1]+ymax[t[0]][-1]
        if not (t[0]=='Fetal' and t[1] in ['BRAIN','MUSCLE']) and not (t[0]=='Adult' and t[1]=='THYMUS'):
            x_nobrain[t[0]].append(x[t[0]][-1])
            y_nobrain[t[0]].append(y[t[0]][-1])
#    if x[t[0]][-1]<0.87 or x[t[0]][-1]>0.97 or y[t[0]][-1]<1.075 or y[t[0]][-1]>1.09:
        if t in [('Fetal', 'Blood & T-cell'),('Adult','Thymus'),('Fetal','Muscle'),('Fetal','Brain'),('Adult','Neurosph'),('Adult','Brain'),('Adult','Mesench'),('Adult','Blood & T-cell'),('Adult','Adipose'),('Fetal','Heart')]:
            plt.text(y[t[0]][-1]+0.001,x[t[0]][-1]+0.005,str(t[1]))
#print linregress(allneand,allsfs)        

    matplotlib.rcParams.update({'font.size':18})
        
    slope,yint,rval,pval,sterror= linregress(allsfs,allneand)
    plt.plot([1.04,1.11],[yint+slope*1.04,yint+slope*1.11],color='grey',linestyle='dashed')

#    plt.ylim(ymin=0.75)
#    plt.ylim(ymax=1.15)
    plt.text(1.05,0.7,'r='+str(0.001*int(1000*rval))+', p<'+str(0.00001*int(100000*pval)),color='grey')

    for t in ['Fetal','Adult','Embryonic']: #,'Exons']: #,'HARs']:
        plt.errorbar(y[t],x[t],yerr=[xmin[t],xmax[t]],xerr=[ymin[t],ymax[t]],fmt='o',color=tcolor[t],label=t)
    plt.legend(loc='upper right',ncol=2,fontsize=14)

    plt.axhline(y=1,linestyle='dotted',color='grey')

    plt.xlabel('Singleton enrichment')
    plt.ylabel(species.capitalize()+' variant depletion')
    plt.savefig(popul+'_'+species+oneortwo+'_analysis/'+popul+'_odds'+species+oneortwo+'_v_oddssingle.png')
    plt.clf()
