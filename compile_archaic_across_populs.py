import sys
species=sys.argv[1] #valid arguments are 'neandertal' or 'denisova'
oneortwo=sys.argv[2] #valid arguments are '2' or '1minus2'

populs=['westeurasia','eastasia','southasia','oceania','america','centralasia']

if oneortwo=='2':
    infile=open(species+'_'+populs[0]+'2_snps_plus_bstat.txt')
else:
    infile=open(populs[0]+'_'+species+'_1minusboth2_snps.txt')

archaic_snps=infile.readlines()
infile.close()

for i in range(len(archaic_snps)):
    archaic_snps[i]=archaic_snps[i].strip('\n')+'\t'+populs[0]+'\n'

for i in range(1,len(populs)):
    master_ind=0
    popul_ind=0

    if oneortwo=='2':
        infile=open(species+'_'+populs[i]+'2_snps_plus_bstat.txt')
    else:
        infile=open(populs[i]+'_'+species+'_1minusboth2_snps.txt')
    popul_lines=infile.readlines()
    infile.close()

    while master_ind<len(archaic_snps) and popul_ind<len(popul_lines):
        master_s=archaic_snps[master_ind].split('\t')
        popul_s=popul_lines[popul_ind].split('\t')
        if int(master_s[1])<int(popul_s[1]):
            master_ind+=1
        elif int(master_s[1])>int(popul_s[1]):
            archaic_snps.insert(master_ind,popul_lines[popul_ind].strip('\n')+'\t'+populs[i]+'\n')
            master_ind+=1
            popul_ind+=1
        elif int(master_s[3])<int(popul_s[3]):
            master_ind+=1
        elif int(master_s[3])>int(popul_s[3]):
            archaic_snps.insert(master_ind,popul_lines[popul_ind].strip('\n')+'\t'+populs[i]+'\n')
            master_ind+=1
            popul_ind+=1
        else:
            popul_ind+=1

outfile=open(species+'_allpops_'+oneortwo+'_snps_plus_bstat.txt','w')
for line in archaic_snps:
    outfile.write(line)
outfile.close()



