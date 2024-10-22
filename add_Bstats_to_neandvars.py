import sys
[oneortwo,species,popul]=sys.argv[1:4]

infile=open(species+'_'+popul+oneortwo+'_snps.txt')
lines=infile.readlines()
infile.close()

output=''

chrom='0'

for line in lines:
    s=line.strip('\n').split('\t')
    if not s[1]==chrom:
        chrom=s[1]
        infile=open('bkgd/hg19_Bstats_bin50_chrom'+chrom+'.bed')
        blines=infile.readlines()
        infile.close()
        bs=blines[0].strip('\n').split('\t')
        bind=0
#    print bs
#    print s
    while bind<len(blines)-1 and int(bs[1])<int(s[3]):
        bind+=1
        bs=blines[bind].strip('\n').split('\t')
    output+=line.strip('\n')+'\t'+bs[2]+'\n'
    if not (int(bs[0])<=int(s[3]) and int(bs[1])>= int(s[3])):
        print chrom, bs

outfile=open(species+'_'+popul+oneortwo+'_snps_plus_bstat.txt','w')
outfile.write(output)
outfile.close()
