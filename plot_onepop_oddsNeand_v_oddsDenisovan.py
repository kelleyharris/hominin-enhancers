import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from math import sqrt
from math import exp
import scipy.stats 
import numpy as np

import sys
popul=sys.argv[1]

def get_min_mean_max(neand,control,tot_neand,tot_control):
    odds_ratio=1.0*(tot_control-control)*neand/(control*(tot_neand-neand))
    std_error=sqrt(1.0/neand+1.0/control+1.0/tot_neand+1.0/tot_control)
    return odds_ratio*(1-exp(-2*std_error)),odds_ratio,odds_ratio*(exp(2*std_error)-1)
                   
neand_nonenhancer, control_nonenhancer, neand_enhancer, control_enhancer=dict({}),dict({}), dict({}), dict({})

infile=open('denisova_'+popul+'2_snps.txt')
lines=infile.readlines()
num_neand=len(lines)
infile.close()

infile=open('denisova_'+popul+'2_controls.txt')
lines=infile.readlines()
num_control=len(lines)
infile.close()

infile=open(popul+'_enhancer_exon_summary_counts.txt')
lines=infile.readlines()
infile.close()

ind=0
while not lines[ind].startswith('neandertal exon 2'):
    ind+=1
s=lines[ind].split(' ')
xmin_exon, neand_exon, xmax_exon = get_min_mean_max(int(s[3]),int(s[4]),num_neand,num_control)
while not lines[ind].startswith('denisova exon 2'):
    ind+=1
s=lines[ind].split(' ')
ymin_exon, denisova_exon, ymax_exon = get_min_mean_max(int(s[3]),int(s[4]),num_neand,num_control)
plt.errorbar([neand_exon],[denisova_exon],xerr=[[xmin_exon],[xmax_exon]],yerr=[[ymin_exon],[ymax_exon]],fmt='o',color='cyan',label='Exon')
#plt.errorbar

infile=open(popul+'_denisova2_analysis/'+popul+'_denisova2_counts_alt_cell_types.txt')
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

odds_info_deni=dict({})
for tissue in tissues:
    print tissue
    odds_info_deni[tissue]=[1.0*(num_control-control_enhancer[tissue])*neand_enhancer[tissue]/(control_enhancer[tissue]*(num_neand-neand_enhancer[tissue])),sqrt(1.0/neand_enhancer[tissue]+1.0/control_enhancer[tissue]+1.0/(num_control-control_enhancer[tissue])+1.0/(num_neand-neand_enhancer[tissue]))]


infile=open('neandertal_'+popul+'2_snps.txt')
lines=infile.readlines()
num_neand=len(lines)
infile.close()

infile=open('neandertal_'+popul+'2_controls.txt')
lines=infile.readlines()
num_control=len(lines)
infile.close()

infile=open(popul+'_neandertal2_analysis/'+popul+'_neandertal2_counts_alt_cell_types.txt')
lines=infile.readlines()
infile.close()

neand_enhancer,control_enhancer=dict({}),dict({})
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


x, y, xmin, xmax, ymin, ymax=dict({}), dict({}), dict({}), dict({}), dict({}),dict({})
x_nobrain, y_nobrain=dict({}),dict({})
for t in ['Adult','Fetal','Embryonic','HARs','Exons']:
    x[t],y[t],xmin[t],xmax[t],ymin[t],ymax[t]=[],[],[],[],[],[]
    x_nobrain[t],y_nobrain[t]=[],[]

parity=dict({})
for pkey in [('Fetal', 'Blood & T-cell',1),('Adult','Thymus',-1),('Fetal','Muscle',1),('Fetal','Brain',-1),('Fetal','Neurosph',1),('Adult','Mesench',-1),('Adult','Epithelial',1),('Fetal','Heart',1),('Fetal','Thymus',-1),('Embryonic','ESC',1),('Adult','Ovary',1),('Adult','Muscle',-1),('Adult','Digestive',1),('Fetal','HSC & B-cell',1),('Adult','HSC & B-cell',1),('Adult','Sm. Muscle',1),('Adult','Digestive',1),('Fetal','Digestive',1),('Adult','Adipose',-2)]:
    parity[tuple(pkey[:2])]=pkey[2]

tcolor=dict({})
for cpair in [('Fetal','red'),('Adult','blue'),('Embryonic','orange'),('Exons','green')]: #,('HARs','black')]:
    tcolor[cpair[0]]=cpair[1]

print tissues
for t in tissues:
    x[t[0]].append(odds_info_deni[t][0])
    xmin[t[0]].append(odds_info_deni[t][0]*(1-exp(-2*odds_info_deni[t][1])))
    xmax[t[0]].append(odds_info_deni[t][0]*(exp(2*odds_info_deni[t][1])-1))
    omean, ovar = 1.0*(num_control-control_enhancer[t])*neand_enhancer[t]/(control_enhancer[t]*(num_neand-neand_enhancer[t])),sqrt(1.0/neand_enhancer[t]+1.0/control_enhancer[t]+1.0/(num_control-control_enhancer[t])+1.0/(num_neand-neand_enhancer[t]))
    y[t[0]].append(omean)
    ymin[t[0]].append(omean*(1-exp(-2*ovar)))
    ymax[t[0]].append(omean*(exp(2*ovar)-1))
    print t[0],t[1],y[t[0]][-1],y[t[0]][-1]-ymin[t[0]][-1],y[t[0]][-1]+ymax[t[0]][-1]
    if not (t[0]=='Fetal' and t[1] in ['BRAIN','MUSCLE']) and not (t[0]=='Adult' and t[1]=='THYMUS'):
        x_nobrain[t[0]].append(x[t[0]][-1])
        y_nobrain[t[0]].append(y[t[0]][-1])
#    if x[t[0]][-1]<0.87 or x[t[0]][-1]>0.97 or y[t[0]][-1]<1.075 or y[t[0]][-1]>1.09:
    if t in parity.keys():
        if parity[t]==1:
            plt.text(y[t[0]][-1]+0.001,x[t[0]][-1]+0.005,str(t[1]),fontsize='10',color=tcolor[t[0]])
        elif parity[t]==-2:
            plt.text(y[t[0]][-1]-0.017,x[t[0]][-1]+0.01,str(t[1]),fontsize='10',color=tcolor[t[0]])
        else:
            plt.text(y[t[0]][-1]+0.001,x[t[0]][-1]-0.015,str(t[1]),fontsize='10',color=tcolor[t[0]])
        

for t in ['Fetal','Adult','Embryonic']: #,'Exons']: #,'HARs']:
    plt.errorbar(y[t],x[t],yerr=[xmin[t],xmax[t]],xerr=[ymin[t],ymax[t]],fmt='o',color=tcolor[t],label=t)
plt.legend(loc='lower right')

plt.plot([0.4,1.1],[0.4,1.1],linestyle='dotted',color='grey')

plt.axhline(y=1,linestyle='dotted',color='grey')
plt.axvline(x=1,linestyle='dotted',color='grey')
plt.xlim(xmax=1)
#plt.ylim(ymax=1.05)
plt.xlim(xmin=0.7)
#plt.ylim(ymin=0.45)

plt.title(popul.capitalize())
plt.xlabel('Neandertal variant depletion')
plt.ylabel('Denisova variant depletion')
plt.savefig(popul+'_allcelltypes_fetal_adult_oddsneand_v_oddsdenisovan.png')
