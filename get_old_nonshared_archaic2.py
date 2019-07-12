import sys
species=sys.argv[1]

for popul in ['westeurasia','oceania','eastasia','southasia','america','centralasia']:
    snps1, snps2=dict({}),dict({})
    infile=open(species+'_'+popul+'1_snps_plus_bstat.txt')
    snps1[species]=infile.readlines()
    infile.close()

    for sp in ['denisova','neandertal']:
        infile=open(species+'_'+popul+'2_snps_plus_bstat.txt')
        snps2[sp]=infile.readlines()
        infile.close()

    for (sp1, sp2) in [(species,'neandertal'),(species,'denisova')]: #[('neandertal','neandertal'),('neandertal','denisova'),('denisova','neandertal'),('denisova','denisova')]:
        ind1, ind2=0,0
        while ind1<len(snps1[sp1]) and ind2<len(snps2[sp2]):
            s1=snps1[sp1][ind1].split('\t')
            s2=snps2[sp2][ind2].split('\t')
            if int(s1[1])<int(s2[1]):
                ind1+=1
            elif int(s2[1])<int(s1[1]):
                ind2+=1
            elif int(s1[3])<int(s2[3]):
                ind1+=1
            elif int(s2[3])<int(s1[3]):
                ind2+=1
            else:
                snps1[sp1].pop(ind1)
                ind2+=1
    outfile=open(popul+'_'+species+'_1minusboth2_snps.txt','w')
    for line in snps1[species]:
        outfile.write(line)
    outfile.close()
    
