import sys

popul=sys.argv[1]

neand_count, control_count=dict({}),dict({})
neand_snps, control_snps=dict({}),dict({})
for sp in ['denisova2','neandertal2','denisova1minus2','neandertal1minus2']:
    for region in ['enhancer','exon']:
        neand_count[(sp,region)], control_count[(sp,region)]=0,0
    neand_snps[sp], control_snps[sp]=[],[]

for sp in ['neandertal','denisova']:
    infile=open(popul+'_'+sp+'_1minusboth2_snps.txt')
    lines=infile.readlines()
    infile.close()
    neand_count[(sp+'1minus2','total')]=len(lines)

    for line in lines:
        s=line.strip('\n').split('\t')
        neand_snps[sp+'1minus2'].append((int(s[1]),int(s[3])))

    infile=open(sp+'_'+popul+'2_snps.txt')
    lines=infile.readlines()
    infile.close()
    neand_count[(sp+'2','total')]=len(lines)

    for line in lines:
        s=line.strip('\n').split('\t')
        neand_snps[sp+'2'].append((int(s[1]),int(s[3])))

    infile=open(popul+'_'+sp+'_1minus2_controls.txt')
    lines=infile.readlines()
    infile.close()
    control_count[(sp+'1minus2','total')]=len(lines)

    for line in lines:
        s=line.strip('\n').split('\t')
        control_snps[sp+'1minus2'].append((int(s[1]),int(s[3])))

    infile=open(sp+'_'+popul+'2_controls.txt')
    lines=infile.readlines()
    infile.close()
    control_count[(sp+'2','total')]=len(lines)

    for line in lines:
        s=line.strip('\n').split('\t')
        control_snps[sp+'2'].append((int(s[1]),int(s[3])))

for sp in neand_snps.keys():    
    neand_ind=0
    for chrom in range(1,23):
        infile=open('gene_locations/exons_chr'+str(chrom)+'.txt')
        genelines=infile.readlines()
        infile.close()

        gene_ind=0
        while neand_ind<len(neand_snps[sp])-1 and neand_snps[sp][neand_ind][0]<chrom:
            neand_ind+=1
        while gene_ind<len(genelines) and neand_ind<len(neand_snps[sp]) and neand_snps[sp][neand_ind][0]==chrom:
            s=genelines[gene_ind].split('\t')
            if int(s[0])>neand_snps[sp][neand_ind][1]:
                neand_ind+=1
            elif int(s[1])<neand_snps[sp][neand_ind][1]:
                gene_ind+=1
            else:
                neand_count[(sp,'exon')]+=1
                neand_ind+=1
    neand_ind=0
    for chrom in range(1,23):
        infile=open('gene_locations/exons_chr'+str(chrom)+'.txt')
        genelines=infile.readlines()
        infile.close()

        gene_ind=0
        while neand_ind<len(control_snps[sp])-1 and control_snps[sp][neand_ind][0]<chrom:
            neand_ind+=1
        while gene_ind<len(genelines) and neand_ind<len(control_snps[sp]) and control_snps[sp][neand_ind][0]==chrom:
            s=genelines[gene_ind].split('\t')
            if int(s[0])>control_snps[sp][neand_ind][1]:
                neand_ind+=1
            elif int(s[1])<control_snps[sp][neand_ind][1]:
                gene_ind+=1
            else:
                control_count[(sp,'exon')]+=1
                neand_ind+=1

infile=open('enhancer_states.txt')
plines=infile.readlines()
infile.close()


for sp in neand_snps.keys():    
    neand_ind=0
    for line in plines[1:]:
        s=line.strip('\n').split(',')
        while neand_ind<len(neand_snps[sp]) and neand_snps[sp][neand_ind][0]<int(s[0]):
            neand_ind+=1
        while neand_ind<len(neand_snps[sp]) and neand_snps[sp][neand_ind][0]==int(s[0]) and neand_snps[sp][neand_ind][1]<int(s[1]):
            neand_ind+=1
        while neand_ind<len(neand_snps[sp]) and neand_snps[sp][neand_ind][0]==int(s[0]) and neand_snps[sp][neand_ind][1]>=int(s[1]) and neand_snps[sp][neand_ind][1]<=int(s[2]):
            neand_count[(sp,'enhancer')]+=1
            neand_ind+=1

    neand_ind=0
    for line in plines[1:]:
        s=line.strip('\n').split(',')
        while neand_ind<len(control_snps[sp]) and control_snps[sp][neand_ind][0]<int(s[0]):
            neand_ind+=1
        while neand_ind<len(control_snps[sp]) and control_snps[sp][neand_ind][0]==int(s[0]) and control_snps[sp][neand_ind][1]<int(s[1]):
            neand_ind+=1
        while neand_ind<len(control_snps[sp]) and control_snps[sp][neand_ind][0]==int(s[0]) and control_snps[sp][neand_ind][1]>=int(s[1]) and control_snps[sp][neand_ind][1]<=int(s[2]):
            control_count[(sp,'enhancer')]+=1
            neand_ind+=1
    
output='Species Function inference_set Count_archaic Count_control\n'

for sp in ['neandertal','denisova']:
    for function in ['exon','enhancer','total']:
        for oneortwo in ['2','1minus2']:
            output+=sp+' '+function+' '+oneortwo+' '+str(neand_count[(sp+oneortwo,function)])+' '+str(control_count[(sp+oneortwo,function)])+'\n'
        
outfile=open(popul+'_enhancer_exon_summary_counts.txt','w')
outfile.write(output)
outfile.close()
