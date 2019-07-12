import sys

[oneortwo,species,popul]=sys.argv[1:4]

neand_tot_1chrom, control_tot_1chrom=dict({}),dict({})
neand_common, control_common=dict({}), dict({})
for chrom in range(1,23):
    neand_tot_1chrom[chrom], control_tot_1chrom[chrom]=0,0
    neand_common[chrom],control_common[chrom]=0,0


if popul=='allpops':
    infile=open(species+'_allpops_'+oneortwo+'_snps_plus_bstat.txt')
else:
    infile=open(species+'_'+popul+oneortwo+'_snps.txt')
#infile=open('neanderthal_eastasian_snps.txt')
neand_vars=[]
infile.readline()
for line in infile:
    s=line.strip('\n').split('\t')
    neand_vars.append((int(s[1]),int(s[3]),int(s[-3])))
    neand_tot_1chrom[int(s[1])]+=1
    if int(s[-3])>5:
        neand_common[int(s[1])]+=1
infile.close()

#infile=open('neanderthal_eastasian_controls.txt')
if popul=='allpops':
    infile=open('allpops_'+species+'_'+oneortwo+'_controls.txt')
else:
    infile=open(popul+'_'+species+'_'+oneortwo+'_controls.txt')
control_vars=[]
infile.readline()
for line in infile:
    s=line.strip('\n').split('\t')
#    print s
    control_vars.append((int(s[1]),int(s[3]),int(s[6])))
    control_tot_1chrom[int(s[1])]+=1
    if int(s[6])>5:
        control_common[int(s[1])]+=1
infile.close()

neand_vars.sort()
control_vars.sort()

output='Variant category Archaic Control\n'
output+='Total '+str(len(neand_vars))+' '+str(len(control_vars))+'\n'

neand_exons, control_exons=0,0
neand_ind=0

neand_1chrom, control_1chrom=dict({}),dict({})

for chrom in range(1,23):
    neand_1chrom[chrom],control_1chrom[chrom]=0,0

#for neand_ind in range(1,len(neand_vars)):
#    if not neand_vars[neand_ind][0]==neand_vars[neand_ind-1][0]:
#        print neand_vars[neand_ind]

infile=open('promoter_states.txt')
plines=infile.readlines()
infile.close()

neand_promoters, control_promoters=0,0
neand_ind,control_ind=0,0
for line in plines[1:]:
    s=line.strip('\n').split(',')
    while neand_ind<len(neand_vars) and neand_vars[neand_ind][0]<int(s[0]):
        neand_ind+=1
    while neand_ind<len(neand_vars) and neand_vars[neand_ind][0]==int(s[0]) and neand_vars[neand_ind][1]<int(s[1]):
        neand_ind+=1
    while neand_ind<len(neand_vars) and neand_vars[neand_ind][0]==int(s[0]) and neand_vars[neand_ind][1]>=int(s[1]) and neand_vars[neand_ind][1]<=int(s[2]):
        neand_promoters+=1
        neand_ind+=1

    while control_ind<len(control_vars) and control_vars[control_ind][0]<int(s[0]):
        control_ind+=1
    while control_ind<len(control_vars) and  control_vars[control_ind][0]==int(s[0]) and control_vars[control_ind][1]<int(s[1]):
        control_ind+=1
    while control_ind<len(control_vars) and  control_vars[control_ind][0]==int(s[0]) and control_vars[control_ind][1]>=int(s[1]) and control_vars[control_ind][1]<=int(s[2]):
        control_promoters+=1
        control_ind+=1

neand_ind=0
for chrom in range(1,23):
#    print 'chrom, neand_ind:', chrom, neand_ind,len(neand_vars)
    infile=open('gene_locations/exons_chr'+str(chrom)+'.txt')
    genelines=infile.readlines()
    infile.close()

    gene_ind=0
    while neand_ind<len(neand_vars)-1 and neand_vars[neand_ind][0]<chrom:
        neand_ind+=1
    while gene_ind<len(genelines) and neand_ind<len(neand_vars) and neand_vars[neand_ind][0]==chrom:
        s=genelines[gene_ind].split('\t')
        if int(s[0])>neand_vars[neand_ind][1]:
#            print s[0], neand_vars[neand_ind][1]
            neand_ind+=1
        elif int(s[1])<neand_vars[neand_ind][1]:
#            print chrom, s, neand_vars[neand_ind]
            gene_ind+=1
        else:
#            print neand_vars[neand_ind], s
            neand_exons+=1
#            if neand_vars[neand_ind][-1]>1:
#                print neand_vars[neand_ind]
            neand_1chrom[chrom]+=1
            neand_ind+=1

print 'neanderthal vars in exons vs not: ',neand_exons,len(neand_vars)-neand_exons
neand_ind=0

for chrom in range(1,23):
    infile=open('gene_locations/exons_chr'+str(chrom)+'.txt')
    genelines=infile.readlines()
    infile.close()

    gene_ind=0
    while control_vars[neand_ind][0]<chrom and neand_ind<len(control_vars)-1:
        neand_ind+=1
    while gene_ind<len(genelines) and neand_ind<len(control_vars) and control_vars[neand_ind][0]==chrom:
        s=genelines[gene_ind].split('\t')
        if int(s[0])>control_vars[neand_ind][1]:
            neand_ind+=1
        elif int(s[1])<control_vars[neand_ind][1]:
#            print chrom, s, neand_vars[neand_ind]
            gene_ind+=1
        else:
#            print control_vars[neand_ind],s
#            if control_vars[neand_ind][-1]>1:
#                print control_vars[neand_ind]
            control_1chrom[chrom]+=1
            control_exons+=1
            neand_ind+=1
print 'Control vars in exons vs not: ',control_exons,len(control_vars)-control_exons
output+='Exons '+str(neand_exons)+' '+str(control_exons)+'\n'
for chrom in range(1,23):
    print chrom, neand_1chrom[chrom],control_1chrom[chrom], neand_tot_1chrom[chrom],control_tot_1chrom[chrom], neand_common[chrom],control_common[chrom], 1.0*neand_1chrom[chrom]/control_1chrom[chrom], 1.0*neand_tot_1chrom[chrom]/control_tot_1chrom[chrom], 1.0*neand_common[chrom]/control_common[chrom]

neand_hars, control_hars=0,0

infile=open('human_accelerated_enhancers_sorted.txt')
harlines=infile.readlines()
infile.close()

neand_ind=0
har_ind=0
while not harlines[har_ind][0]=='1':
    har_ind+=1
print har_ind,len(harlines),harlines[har_ind]
for line in harlines[har_ind:]:
    s=line.strip('\n').split('\t')
    while neand_ind<len(neand_vars)-1 and int(s[0])>neand_vars[neand_ind][0]:
        neand_ind+=1
    while neand_vars[neand_ind][0]==int(s[0]) and neand_ind<len(neand_vars)-1 and int(s[1])>neand_vars[neand_ind][1]:
        neand_ind+=1
    while neand_ind<len(neand_vars)-1 and int(s[0])==neand_vars[neand_ind][0] and int(s[1])<=neand_vars[neand_ind][1] and int(s[2])>=neand_vars[neand_ind][1]:
        neand_hars+=1
#        print'N', s[0],s[1]
        neand_ind+=1

print 'neanderthal vars in HARs vs not', neand_hars, len(neand_vars)-neand_hars

neand_ind=0

for line in harlines[har_ind:]:
    s=line.strip('\n').split('\t')
    while neand_ind<len(control_vars)-1 and int(s[0])>control_vars[neand_ind][0]:
        neand_ind+=1
    while control_vars[neand_ind][0]==int(s[0]) and neand_ind<len(control_vars)-1 and int(s[1])>control_vars[neand_ind][1]:
        neand_ind+=1
    while neand_ind<len(control_vars)-1 and int(s[0])==control_vars[neand_ind][0] and int(s[1])<=control_vars[neand_ind][1] and int(s[2])>=control_vars[neand_ind][1]:
        control_hars+=1
#        print'C', s[0],s[1]
        neand_ind+=1

output+='HARs '+str(neand_hars)+' '+str(control_hars)+'\n'
print 'control vars in HARs vs not', control_hars, len(control_vars)-control_hars

output+='Promoters '+str(neand_promoters)+' '+str(control_promoters)+'\n'
print 'neand v control in promoters', neand_promoters, control_promoters

if popul=='allpops':
    outfile=open('allpop_analysis/'+species+'_exon_har_variant_counts.txt','w')
else:
    outfile=open(popul+'_'+species+oneortwo+'_analysis/exon_har_variant_counts.txt','w')
outfile.write(output)
outfile.close()
