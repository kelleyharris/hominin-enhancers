import subprocess
import sys
from subprocess import PIPE

if sys.argv[1]=='allpops':
    [oneortwo,species]=sys.argv[2:4]
    popul='all'
elif len(sys.argv)>3:
    [oneortwo,species,popul]=sys.argv[1:4]
elif len(sys.argv)>2:
    [popul,species]=sys.argv[1:3]
    oneortwo='minus'
else:
    oneortwo='both'
    popul=sys.argv[1]

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
    elif s[3] in ['ESC','iPSC']:
        celltype[s[1]]=('Embryonic',s[3])
    else:
        celltype[s[1]]=('Adult',s[3])
#        print 'Adult', line
    if s[3]=='Other' and s[5]=='OVRY':
        celltype[s[1]]=('Adult','Ovary')
    elif s[3] in ['IMR90','GROUP','ENCODE2012','','Other']:
        celltype[s[1]]='Other'
    if not celltype[s[1]] in celltypes and not celltype[s[1]]=='Other':
        celltypes.append(celltype[s[1]])
s=enhancer_lines[0].strip('\n').split(',')
max_pleio=len(s)-3

for i in range(3,len(s)):
    celltype[i]=celltype[s[i]]

#print celltypes
num_enhancers, num_controls, num_neand = dict({}),dict({}),dict({})
for i in range(1,max_pleio):
#    for stage in ['Fetal','Adult']:
    for c in celltypes:
        thiskey=(i,c[0],c[1])
        num_enhancers[thiskey]=0
        num_controls[thiskey]=0
        num_neand[thiskey]=0
#print num_enhancers.keys()

#infile=open('neanderthal_eastasian_snps.txt')
#infile=open('neanderthal_westeurasian_snps.txt')
if popul=='all':
    infile=open(species+'_allpops_'+oneortwo+'_snps_plus_bstat.txt')
elif oneortwo=='both':
    infile=open(popul+'_combined_archaic_1minusboth2_snps.txt')
elif oneortwo=='minus':
    infile=open(popul+'_'+species+'_1minus2_snps.txt')
else:
    infile=open(species+'_'+popul+oneortwo+'_snps.txt')
neand_vars=[]
infile.readline()
for line in infile:
    s=line.strip('\n').split('\t')
    neand_vars.append((int(s[1]),int(s[3])))
infile.close()

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
    for atype in active_types:
        thiskey=(len(active_types),atype[0],atype[1])
        num_enhancers[thiskey]+=1
        num_neand[thiskey]+=neanders
        num_controls[thiskey]+=controls

print total_neand,total_control
print num_neand
print num_controls
#alias=dict({})
#for c in celltypes:
#    if c=='Blood & T-cell':
#        alias
    
output='Pleiotropy_number Tissue Num_enhancers Neand_vars Control_vars\n'
for pleio in range(1,highest_pleio+1):
#    for stage in ['Adult','Fetal']:
    for c in celltypes:
        if c[1]=='Neurosph':
            output+=str(pleio)+',Fetal,'+c[1]+','+str(num_enhancers[(pleio,c[0],c[1])])+','+str(num_neand[(pleio,c[0],c[1])])+','+str(num_controls[(pleio,c[0],c[1])])+'\n'
        else:
            output+=str(pleio)+','+c[0]+','+c[1]+','+str(num_enhancers[(pleio,c[0],c[1])])+','+str(num_neand[(pleio,c[0],c[1])])+','+str(num_controls[(pleio,c[0],c[1])])+'\n'

if popul=='all':
    outfile=open('allpop_analysis/allpops_'+species+'_'+oneortwo+'_counts_alt_cell_types.txt','w')
elif oneortwo=='both':
    outfile=open(popul+'_1minus2_analysis/'+popul+'_1minus2_counts_alt_cell_types.txt','w')
elif oneortwo=='minus':
    outfile=open('bugtest_1minus2/'+popul+'_'+species+'_1minus2_counts_alt_cell_types.txt','w')
else:
    outfile=open(popul+'_'+species+oneortwo+'_analysis/'+popul+'_'+species+oneortwo+'_counts_alt_cell_types.txt','w')

outfile.write(output)
outfile.close()
            
