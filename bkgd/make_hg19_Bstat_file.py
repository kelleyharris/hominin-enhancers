blist=[]
for i in range(20):
    infile=open('hg19_bin50_B'+str(50*i)+'.bed')
    lines=infile.readlines()
    infile.close()
    for line in lines:
        s=line.strip('\n').split('\t')
        blist.append((s[0],int(s[1]),int(s[2]),50*i))
blist.sort()

this_chrom=blist[0][0]
output=''
for i in range(len(blist)):
    if not blist[i][0]==this_chrom:
        outfile=open('hg19_Bstats_bin50_chrom'+this_chrom[3:]+'.bed','w')
        outfile.write(output)
        outfile.close()
        output=''
        this_chrom=blist[i][0]
    output+=str(blist[i][1])+'\t'+str(blist[i][2])+'\t'+str(blist[i][3])+'\n'

outfile=open('hg19_Bstats_bin50_chrom'+blist[-1][0][3:]+'.bed','w')
outfile.write(output)
outfile.close()
