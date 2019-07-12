import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from math import sqrt
from math import exp
import scipy.stats 
import numpy as np
from scipy.stats import linregress

import sys
species=sys.argv[1]

neand_nonenhancer, control_nonenhancer, neand_enhancer, control_enhancer=dict({}),dict({}), dict({}), dict({})

def get_min_mean_max(neand,control,tot_neand,tot_control):
    odds_ratio=1.0*(tot_control-control)*neand/(control*(tot_neand-neand))
    std_error=sqrt(1.0/neand+1.0/control+1.0/tot_neand+1.0/tot_control)
    return odds_ratio*(1-exp(-2*std_error)),odds_ratio,odds_ratio*(exp(2*std_error)-1)

infile=open(species+'_allpops_2_snps_plus_bstat.txt')
lines=infile.readlines()
num_neand=len(lines)
infile.close()

infile=open('allpops_'+species+'_2_controls.txt')
lines=infile.readlines()
num_control=len(lines)
infile.close()

enhancer_rare, enhancer_common, control_rare, control_common = dict({}), dict({}), dict({}), dict({})

enhancer_common[('HARs','HARs')],control_common[('HARs','HARs')],enhancer_rare[('HARs','HARs')],control_rare[('HARs','HARs')]=6130,6292,3745,3655
enhancer_common[('Exons','Exons')],control_common[('Exons','Exons')],enhancer_rare[('Exons','Exons')],control_rare[('Exons','Exons')]=463826, 497059, 378496, 320136

ymin_rare, rare_exon, ymax_rare = get_min_mean_max(enhancer_rare[('Exons','Exons')],enhancer_common[('Exons','Exons')],enhancer_rare[('Exons','Exons')]+control_rare[('Exons','Exons')],enhancer_common[('Exons','Exons')]+control_common[('Exons','Exons')])

infile=open('allpop_analysis/'+species+'_exon_har_variant_counts.txt')
lines=infile.readlines()
infile.close()

ind=0
while not lines[ind].startswith('Exons'):
    ind+=1
s=lines[ind].split(' ')
xmin_exon, neand_exon, xmax_exon = get_min_mean_max(int(s[1]),int(s[2]),num_neand,num_control)
ind+=1
s=lines[ind].split(' ')
xmin_har, neand_har, xmax_har = get_min_mean_max(int(s[1]),int(s[2]),num_neand,num_control)

infile=open('allpop_analysis/allpops_'+species+'_2_counts_alt_cell_types.txt')
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

x, y, xmin, xmax, ymin, ymax=dict({}), dict({}), dict({}), dict({}), dict({}),dict({})
for t in ['Adult','Fetal','Embryonic','HARs','Exons']:
    x[t],y[t],xmin[t],xmax[t],ymin[t],ymax[t]=[],[],[],[],[],[]

parity=dict({})
for pkey in [('Fetal', 'Blood & T-cell',1),('Adult','Thymus',-1),('Fetal','Muscle',-2),('Fetal','Brain',1),('Fetal','Neurosph',1),('Adult','Mesench',1),('Embryonic','iPSC',1),('Embryonic','ESC',1),('Adult','Brain',1)]:
    parity[tuple(pkey[:2])]=pkey[2]

tcolor=dict({})
for cpair in [('Fetal','red'),('Adult','blue'),('Embryonic','orange'),('Exons','green')]: #,('HARs','black')]:
    tcolor[cpair[0]]=cpair[1]

infile=open('../rare_common_vars_alt_cell_types_enhancer_control.txt')
lines=infile.readlines()
infile.close()

tissues=[]
for line in lines[2:]:
    s=line.strip('\n').split(',')
    if s[2] in ['ESC','iPSC']:
        s[1]='Embryonic'
    elif s[2]=='Neurosph':
        s[1]='Fetal'
    if not (s[1],s[2]) in tissues and not s[2][:3]==' Un':
        tissues.append((s[1],s[2]))
        t=tissues[-1]
        enhancer_rare[t]=int(s[3])
        enhancer_common[t]=int(s[5])+int(s[4])
        control_rare[t]=int(s[6])
        control_common[t]=int(s[7])+int(s[8])
    elif not s[2][:3]==' Un':
        t=(s[1],s[2])
        enhancer_rare[t]+=int(s[3])
        enhancer_common[t]+=int(s[4])+int(s[5])
        control_rare[t]+=int(s[6])
        control_common[t]+=int(s[7])+int(s[8])
    if s[2] in ['Brain','Neurosph']:
        tissues.pop()

for t in tissues:
    x[t[0]].append(odds_info_deni[t][0])
    xmin[t[0]].append(odds_info_deni[t][0]*(1-exp(-2*odds_info_deni[t][1])))
    xmax[t[0]].append(odds_info_deni[t][0]*(exp(2*odds_info_deni[t][1])-1))
    raremin, rare, raremax = get_min_mean_max(enhancer_rare[t],enhancer_common[t],enhancer_rare[t]+control_rare[t],enhancer_common[t]+control_common[t])
#    omean, ovar = 1.0*(num_control-control_enhancer[t])*neand_enhancer[t]/(control_enhancer[t]*(num_neand-neand_enhancer[t])),sqrt(1.0/neand_enhancer[t]+1.0/control_enhancer[t]+1.0/(num_control-control_enhancer[t])+1.0/(num_neand-neand_enhancer[t]))
    y[t[0]].append(rare)
    ymin[t[0]].append(raremin)
    ymax[t[0]].append(raremax)
    if t in parity.keys():
        if parity[t]==1:
            plt.text(x[t[0]][-1]+0.001,y[t[0]][-1]+0.005,str(t[1]),fontsize='10',color=tcolor[t[0]])
        elif parity[t]==-2:
            plt.text(x[t[0]][-1]-0.017,y[t[0]][-1]+0.01,str(t[1]),fontsize='10',color=tcolor[t[0]])
        else:
            plt.text(x[t[0]][-1]+0.001,y[t[0]][-1]-0.015,str(t[1]),fontsize='10',color=tcolor[t[0]])
        

for t in ['Fetal','Adult','Embryonic']: #,'Exons']: #,'HARs']:
    plt.errorbar(x[t],y[t],xerr=[xmin[t],xmax[t]],yerr=[ymin[t],ymax[t]],fmt='o',color=tcolor[t],label=t)
plt.legend(loc='upper right')

plt.axvline(x=1,linestyle='dotted',color='grey')

slope,yint,rval,pval,sterror= linregress(x['Fetal']+x['Adult']+x['Embryonic'],y['Fetal']+y['Adult']+y['Embryonic'])
plt.plot([0.6,1.05],[yint+slope*0.6,yint+slope*1.05],color='grey',linestyle='dashed')

if species=='denisova':
    plt.xlim(xmin=0.6)
    plt.xlim(xmax=1.05)
    plt.text(0.7,1.05,r'$r^2=$'+str(0.001*int(1000*rval**2))+',\n p<'+str(0.00001*int(100000*pval)),color='grey')
else:
    plt.xlim(xmin=0.75)
    plt.xlim(xmax=0.95)
    plt.text(0.8,1.055,r'$r^2=$'+str(0.001*int(1000*rval**2))+',\n p<'+str(0.00001*int(100000*pval)),color='grey')
    

plt.title('All tissues except fetal brain, neurospheres, and adult brain')
plt.xlabel(species.capitalize()+' variant depletion')
plt.ylabel('Rare variant enrichment')
plt.savefig('allpop_analysis/allcelltypes_fetal_adult_odds'+species+'_v_oddsrare_nobrain.png')
