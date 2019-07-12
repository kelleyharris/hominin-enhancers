import sys
chrom=sys.argv[1]

infile=open('refGene.txt')
lines=infile.readlines()
infile.close()

starts=[]
end=dict({})
direction=dict({})

for line in lines:
    s=line.split('\t')
    if s[2]=='chr'+chrom:
        starts.append(int(s[4]))
        end[starts[-1]]=int(s[5])
        direction[starts[-1]]=s[3]
starts.sort()

i=1
while i<len(starts):
    if starts[i]<=end[starts[i-1]]:
        end[starts[i-1]]=max(end[starts[i-1]],end[starts[i]])
        starts.pop(i)
    else:
        i+=1

output=''
for i in range(len(starts)):
    output+=str(starts[i])+'\t'+str(end[starts[i]])+'\t'+direction[starts[i]]+'\n'
outfile=open('genes_chr'+chrom+'.txt','w')
outfile.write(output)
outfile.close()
