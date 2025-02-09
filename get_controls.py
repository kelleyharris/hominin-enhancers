import random
from numpy import arange
import sys

[oneortwo,species,popul]=sys.argv[1:4]

if len(sys.argv)>4:
    infile=open('allpop_component_category_counts_'+species+'_'+popul+oneortwo+'_snps.txt')
else:
    infile=open('category_counts_'+species+'_'+popul+oneortwo+'_snps.txt')
lines=infile.readlines()
infile.close()

selected_indices=dict({})
observed_vars=dict({})
next_selected=dict({})
categories=[(1,0)]

residual=0
for line in lines[1:]:
    s=line.strip('\n').split(' ')
    if not (int(s[0]),int(s[1]))==categories[-1]:
        categories.append((int(s[0]),int(s[1])))
    if 2*int(s[2])>int(s[3])+residual:
        print line
        selected_indices[(int(s[0]),int(s[1]))]=arange(1,int(s[3])+1)
        residual+=2*int(s[2])-int(s[3])-residual
    else:
        selected_indices[(int(s[0]),int(s[1]))]=random.sample(range(1,int(s[3])+1),2*int(s[2])+residual)
        residual=0
    selected_indices[(int(s[0]),int(s[1]))].sort()
    observed_vars[(int(s[0]),int(s[1]))]=0
    next_selected[(int(s[0]),int(s[1]))]=0

#print selected_indices

bstat_range=dict({})
for i in range(1,categories[-1][0]+1):
    bstat_range[i]=[]
for i in range(len(categories)):
    bstat_range[categories[i][0]].append(categories[i][1])

output=''

for chrom in range(1,23):
    infile=open('bkgd/hg19_Bstats_bin50_chrom'+str(chrom)+'.bed')
    blines=infile.readlines()
    infile.close()
    bs=blines[0].strip('\n').split('\t')
    bind=0

    infile=open('summaries/'+oneortwo+'/'+species+'/'+popul+'/summaries/pred.'+str(chrom)+'.thresh-50.length-0.00')
    control_lines=infile.readlines()
    infile.close()
    cind=0
    while control_lines[cind][0]=='#':
        cind+=1

    for line in control_lines[cind:]:
        s=line.strip('\n').split('\t')
        if float(s[-4])==0 and int(s[6])<=categories[-1][0] and int(s[6])>=1 and len(bstat_range[int(s[6])])>0:            
            while bind<len(blines)-1 and int(bs[1])<int(s[3]):
                bind+=1
                bs=blines[bind].strip('\n').split('\t')
            this_bstat=int(bs[2])
            if len(bstat_range[int(s[6])])==20 or (this_bstat>=bstat_range[int(s[6])][0] and this_bstat<=bstat_range[int(s[6])][-1] and this_bstat in bstat_range[int(s[6])]):
                observed_vars[(int(s[6]),this_bstat)]+=1
                if next_selected[(int(s[6]),this_bstat)]<len(selected_indices[(int(s[6]),this_bstat)]) and observed_vars[(int(s[6]),this_bstat)]==selected_indices[(int(s[6]),this_bstat)][next_selected[(int(s[6]),this_bstat)]]:
                    output+=line.strip('\n')+'\t'+str(this_bstat)+'\n'
                    next_selected[(int(s[6]),this_bstat)]+=1
                    
    print chrom

if len(sys.argv)>4:
    outfile=open('allpop_component_'+species+'_'+popul+oneortwo+'_controls.txt','w')
else:
    outfile=open(species+'_'+popul+oneortwo+'_controls.txt','w')
outfile.write(output)
outfile.close()
    
