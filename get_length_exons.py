import sys
chrom=sys.argv[1]

infile=open('exons_chr'+chrom+'.txt')
lines=infile.readlines()
infile.close()

length=0
for line in lines:
    s=line.split('\t')
    length+=int(s[1])-int(s[0])+2
print length
