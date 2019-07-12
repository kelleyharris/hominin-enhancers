import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from math import sqrt
from math import exp
import scipy.stats 
import numpy as np
import subprocess
import sys
from subprocess import PIPE

if sys.argv[1]=='allpops':
    [oneortwo, species]=sys.argv[2:4]
    popul='all'
elif len(sys.argv)>3:
    [oneortwo,species,popul]=sys.argv[1:4]
elif len(sys.argv)>2:
    [popul,species]=sys.argv[1:3]
    oneortwo='minus'
else:
    oneortwo='both'
    popul=sys.argv[1]

if popul=='all':
    infile=open(species+'_allpops_'+oneortwo+'_snps_plus_bstat.txt')
elif oneortwo=='both':
    infile=open(popul+'_combined_archaic_1minusboth2_snps.txt')
elif oneortwo=='minus':
    infile=open(popul+'_'+species+'_1minus2_snps.txt')
else:
    infile=open(species+'_'+popul+oneortwo+'_snps.txt')
lines=infile.readlines()
num_neand=len(lines)
infile.close()

if popul=='all':
    infile=open('allpops_'+species+'_'+oneortwo+'_controls.txt')
elif oneortwo=='both':
    infile=open(popul+'_1minus2_controls.txt')
elif oneortwo=='minus':
    infile=open(popul+'_'+species+'_1minus2_controls.txt')
else:
    infile=open(species+'_'+popul+oneortwo+'_controls.txt')
lines=infile.readlines()
num_control=len(lines)
infile.close()

print num_neand,num_control

neand_nonenhancer, control_nonenhancer, neand_enhancer, control_enhancer=dict({}),dict({}), dict({}), dict({})
#num_neand,num_control=0,0

#infile=open('eur_neand_control_counts_alt_cell_types.txt')
if popul=='all':
    infile=open('allpop_analysis/allpops_'+species+'_'+oneortwo+'_counts_alt_cell_types.txt')
elif oneortwo=='both':
    infile=open(popul+'_1minus2_analysis/'+popul+'_1minus2_counts_alt_cell_types.txt')
elif oneortwo=='minus':
    infile=open('bugtest_1minus2/'+popul+'_'+species+'_1minus2_counts_alt_cell_types.txt')
else:
    infile=open(popul+'_'+species+oneortwo+'_analysis/'+popul+'_'+species+oneortwo+'_counts_alt_cell_types.txt')
lines=infile.readlines()
infile.close()

tissues=[]
for line in lines[2:]:
    s=line.strip('\n').split(',')
    if not (s[1],s[2]) in tissues:
        tissues.append((s[1],s[2]))
        t=tissues[-1]
        neand_enhancer[t]=int(s[4])
        control_enhancer[t]=int(s[5])
    else:
        t=(s[1],s[2])
        neand_enhancer[t]+=int(s[4])
        control_enhancer[t]+=int(s[5])
print neand_enhancer
print control_enhancer
odds_info=[]
for tissue in tissues:
    if neand_enhancer[tissue]*control_enhancer[tissue]==0:
        del neand_enhancer[tissue]
        del control_enhancer[tissue]
    else:
        chi2_info=scipy.stats.fisher_exact(np.array([[neand_enhancer[tissue],num_neand-neand_enhancer[tissue]],[control_enhancer[tissue],num_control-control_enhancer[tissue]]]))
        odds_info.append([1.0*(num_control-control_enhancer[tissue])*neand_enhancer[tissue]/(control_enhancer[tissue]*(num_neand-neand_enhancer[tissue])),sqrt(1.0/neand_enhancer[tissue]+1.0/control_enhancer[tissue]+1.0/(num_control-control_enhancer[tissue])+1.0/(num_neand-neand_enhancer[tissue])),tissue, chi2_info])

odds_info.sort()
x, y, ymin, ymax=dict({}),dict({}),dict({}),dict({})
xtotal,xlabel=[],[]
for t in ['Adult','Fetal','Embryonic']:
    x[t],y[t],ymin[t],ymax[t]=[],[],[],[]
for oi in odds_info:
    print oi[0], oi[0]*exp(-2*oi[1]), oi[0]*exp(2*oi[1]),oi[2], oi[-1]
    if oi[2][1].startswith('ESC') or oi[2][1]=='iPSC':
        oi[2]=('Embryonic',oi[2][1])
#        oi[2][0]='Embryonic'
    x[oi[2][0]].append(len(x['Fetal'])+len(x['Adult'])+len(x['Embryonic'])+1)
    xlabel.append(oi[2][1])
    y[oi[2][0]].append(oi[0])
    ymin[oi[2][0]].append(oi[0]-oi[0]*exp(-2*oi[1]))
    ymax[oi[2][0]].append(oi[0]*exp(2*oi[1])-oi[0])
#    print t, oi[2][1], ymin[oi[2][0]][-1],y[oi[2][0]][-1],ymax[oi[2][0]][-1]
#    if oi[-1][-1]<0.05:
#        if len(x)%2==0:
#            plt.text(x[oi[2][0]][-1],1.2*ymax[oi[2][0]][-1]+y[oi[2][0]][-1],'p<'+('%0.3g' % oi[-1][-1]),fontsize=6)
#        else:
#            plt.text(x[oi[2][0]][-1],-1.2*ymin[oi[2][0]][-1]+y[oi[2][0]][-1],'p<'+('%0.3g' % oi[-1][-1]),fontsize=6)
#errorbars=[]
#for i in range(len(x)):
#    errorbars.append([ymin[i],ymax[i]])
plt.xlim(xmax=len(xlabel)+1)
xtotal=np.arange(1,len(xlabel)+1)
#plt.ylim(ymax=1.25)
#plt.ylim(ymin=0.8)
plt.axhline(y=1,linestyle='dotted',color='grey')
if oneortwo=='both':
    plt.ylabel('archaic-to-control variant odds ratio')
else:
    plt.ylabel(species+'-to-control variant odds ratio')
plt.xticks(xtotal,xlabel,rotation='vertical')
plt.errorbar(x['Fetal'],y['Fetal'],yerr=[ymin['Fetal'],ymax['Fetal']],color='red',fmt='o')
plt.errorbar(x['Adult'],y['Adult'],yerr=[ymin['Adult'],ymax['Adult']],color='blue',fmt='o')
plt.errorbar(x['Embryonic'],y['Embryonic'],yerr=[ymin['Embryonic'],ymax['Embryonic']],color='orange',fmt='o')
fig=plt.gcf()
fig.subplots_adjust(bottom=0.35)

plt.ylim(ymin=0.6)
plt.ylim(ymax=1.05)

if popul=='all':
    plt.savefig('allpop_analysis/allpops_'+species+'_'+oneortwo+'_allpleio_allcelltypes_fetal_adult.png')
elif oneortwo=='both':
    plt.savefig(popul+'_1minus2_analysis/'+popul+'_1minus2_allpleio_allcelltypes_fetal_adult.png')
elif oneortwo=='minus':
    plt.savefig('bugtest_1minus2/'+popul+'_'+species+'_1minusboth2_allpleio_allcelltypes_fetal_adult.png')
else:
    plt.savefig(popul+'_'+species+oneortwo+'_analysis/'+popul+'_'+species+oneortwo+'_allpleio_allcelltypes_fetal_adult.png')
