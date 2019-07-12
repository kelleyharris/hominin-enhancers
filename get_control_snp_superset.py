import sys

[oneortwo,species,popul]=sys.argv[1:4]
if len(sys.argv)>4:
    infile=open('bfsorted_allpop_component_'+popul+'_'+species+'_'+oneortwo+'_snps.txt')
else:
    infile=open('bstat_freq_sorted_'+species+'_'+popul+oneortwo+'_snps.txt')
lines=infile.readlines()
infile.close()

den_count,nonden_count=dict({}),dict({})
categories=[(1,0)]
den_count[(1,0)],nonden_count[(1,0)]=0,0
for line in lines:
    s=line.strip('\n').split('\t')
    print s, categories[-1]
    if not (int(s[6]),int(s[17]))==categories[-1]:
        categories.append((int(s[6]),int(s[17])))
        den_count[categories[-1]]=0
        nonden_count[categories[-1]]=0
    den_count[categories[-1]]+=1                          

bstat_range=dict({})
for i in range(1,categories[-1][0]+1):
    bstat_range[i]=[]
    for j in range(20):
        if (i,j*50) in categories:
            bstat_range[i].append(j*50)
    print i, bstat_range[i]

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
#        print s
        if float(s[13])==0 and int(s[6])<=categories[-1][0] and int(s[6])>=1 and len(bstat_range[int(s[6])])>0:            
            while bind<len(blines)-1 and int(bs[1])<int(s[3]):
                bind+=1
                bs=blines[bind].strip('\n').split('\t')
            this_bstat=int(bs[2])
            if len(bstat_range[int(s[6])])==20 or (this_bstat>=bstat_range[int(s[6])][0] and this_bstat<=bstat_range[int(s[6])][-1] and this_bstat in bstat_range[int(s[6])]):
                nonden_count[(int(s[6]),this_bstat)]+=1
    print chrom

if len(sys.argv)>4:
    outfile=open('allpop_component_category_counts_'+species+'_'+popul+oneortwo+'_snps.txt','w')
else:
    outfile=open('category_counts_'+species+'_'+popul+oneortwo+'_snps.txt','w')
output='Allele_count Bstat Neanderthal non-Neanderthal\n'
for cat in categories:
    output+=str(cat[0])+' '+str(cat[1])+' '+str(den_count[cat])+' '+str(nonden_count[cat])+'\n'
outfile.write(output)
outfile.close()
