import subprocess
import sys
from subprocess import PIPE

focal_celltypes=[('Fetal','Brain'),('Fetal','Muscle')]

if sys.argv[1]=='allpops':
    [oneortwo, species]=sys.argv[2:4]
    popul='all'
elif len(sys.argv)>3:
    [oneortwo,species,popul]=sys.argv[1:4]
elif len(sys.argv)>2:
    [popul,species]=sys.argv[1:3]
    oneortwo='minus'
else:
    popul=sys.argv[1]
    oneortwo='both'

infile=open('Consolidated_EpigenomeIDs_summary_Table.csv')
cell_lines=infile.readlines()
infile.close()

infile=open('enhancer_states.txt')
enhancer_lines=infile.readlines()
infile.close()
 
celltypes=[]
celltype=dict({})
for line in cell_lines:
    s=line.strip('\n').split(',')
    if s[15].startswith('Fetal') or s[19].startswith('Fetus'):
        celltype[s[1]]=('Fetal',s[3])
    else:
        celltype[s[1]]=('Adult',s[3])
    if s[3] in ['IMR90','GROUP','ENCODE2012','','Other']:
        celltype[s[1]]='Other'
    elif not s[3] in celltypes:
        celltypes.append(s[3])
s=enhancer_lines[0].strip('\n').split(',')
max_pleio=len(s)-3

for i in range(3,len(s)):
    celltype[i]=celltype[s[i]]

num_enhancers, num_controls, num_neand = dict({}),dict({}),dict({})
for i in range(1,max_pleio):
    num_enhancers[i]=0
    num_controls[i]=0
    num_neand[i]=0

if popul=='all':
    infile=open(species+'_allpops_'+oneortwo+'_snps_plus_bstat.txt')
elif oneortwo=='both':
    infile=open(popul+'_combined_archaic_1minusboth2_snps.txt')
elif oneortwo=='minus':
    infile=open(popul+'_'+species+'_1minusboth2_snps.txt')
#    infile=open(popul+'_'+species+'_1minus2_snps.txt')
else:
    infile=open(species+'_'+popul+oneortwo+'_snps.txt')
#infile=open('neanderthal_eastasian_snps.txt')
#infile=open('neanderthal_westeurasian_snps.txt')
neand_vars=[]
infile.readline()
for line in infile:
    s=line.strip('\n').split('\t')
    neand_vars.append((int(s[1]),int(s[3])))
infile.close()

#infile=open('neanderthal_eastasian_controls.txt')
#infile=open('neanderthal_westeurasian_controls.txt')
if popul=='all':
    infile=open('allpops_'+species+'_'+oneortwo+'_controls.txt')
elif oneortwo=='both':
    infile=open(popul+'_1minus2_controls.txt')
elif oneortwo=='minus':
    infile=open(popul+'_'+species+'_1minus2_controls.txt')
else:
    infile=open(species+'_'+popul+oneortwo+'_controls.txt')
control_vars=[]
infile.readline()
for line in infile:
    s=line.strip('\n').split('\t')
    control_vars.append((int(s[1]),int(s[3])))
infile.close()

neand_vars.sort()
control_vars.sort()

print 'Neand', len(neand_vars)
print 'control',len(control_vars)

neand_ind,control_ind=0,0
highest_pleio=1

start_X=0
while not enhancer_lines[start_X][0]=='X':
    start_X+=1

total_neand,total_control=0,0
for line in enhancer_lines[1:start_X]:
    active_types=[]
    s=line.strip('\n').split(',')
    for i in range(3,len(s)):
        if not s[i]=='0' and not celltype[i]=='Other' and not celltype[i] in active_types:
            active_types.append(celltype[i])
    chrom, start, end=int(s[0]),int(s[1]),int(s[2])
    while neand_ind<len(neand_vars) and neand_vars[neand_ind][0]<chrom:
        neand_ind+=1
    while control_ind<len(control_vars) and control_vars[control_ind][0]<chrom:
        control_ind+=1
    while neand_ind<len(neand_vars) and neand_vars[neand_ind][0]==chrom and neand_vars[neand_ind][1]<start:
        neand_ind+=1
    while control_ind<len(control_vars) and control_vars[control_ind][0]==chrom and control_vars[control_ind][1]<start:
        control_ind+=1
    neanders,controls=0,0
    while neand_ind<len(neand_vars) and neand_vars[neand_ind][0]==chrom and neand_vars[neand_ind][1]<=end:
        neanders+=1
        neand_ind+=1
    while control_ind<len(control_vars) and control_vars[control_ind][0]==chrom and control_vars[control_ind][1]<=end:
        controls+=1
        control_ind+=1
    if len(active_types)>highest_pleio:
        highest_pleio=len(active_types)
    total_neand+=neanders
    total_control+=controls
    if len(active_types)>0 and (focal_celltypes[0] in active_types or focal_celltypes[1] in active_types):
        num_enhancers[len(active_types)]+=1
        num_neand[len(active_types)]+=neanders
        num_controls[len(active_types)]+=controls

print total_neand,total_control
    
output='Pleiotropy_number Num_enhancers Neand_vars Control_vars\n'
for pleio in range(1,highest_pleio+1):
    output+=str(pleio)+','+str(num_enhancers[pleio])+','+str(num_neand[pleio])+','+str(num_controls[pleio])+'\n'

#outfile=open('eur_neand_control_counts_v_celltype_pleio_number.txt','w')
if popul=='all':
    outfile=open('allpop_analysis/allpops_'+species+'_'+oneortwo+'_counts_withfbrain_and_fmuscle_v_celltype_pleio_number.txt','w')
elif oneortwo=='both':
    outfile=open(popul+'_1minus2_analysis/'+popul+'_1minus2_counts_withfbrain_and_fmuscle_v_celltype_pleio_number.txt','w')
elif oneortwo=='minus':
    outfile=open('bugtest_1minus2/'+popul+'_'+species+'_1minus2_counts_withfbrain_and_fmuscle_v_celltype_pleio_number.txt','w')
else:
    outfile=open(popul+'_'+species+oneortwo+'_analysis/'+popul+'_'+species+oneortwo+'_counts_withfbrain_and_fmuscle_v_celltype_pleio_number.txt','w')
outfile.write(output)
outfile.close()
            
