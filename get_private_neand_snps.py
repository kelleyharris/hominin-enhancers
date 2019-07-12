import sys

[oneortwo,species,popul]=sys.argv[1:4]

output=''

for chrom in range(1,23):
    variants=0
    infile=open('summaries/'+oneortwo+'/'+species+'/'+popul+'/summaries/pred.'+str(chrom)+'.thresh-50.length-0.00')
    lines=infile.readlines()
    infile.close()

    ind=0
    while lines[ind].startswith('#'):
        ind+=1
    for line in lines[ind:]:
        s=line.strip('\n').split('\t')
        if s[7]=='1':
            output+=line
            variants+=1
    print chrom, variants

outfile=open(species+'_'+popul+oneortwo+'_snps.txt','w')
outfile.write(output)
outfile.close()
