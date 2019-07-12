import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

populs=['oceania','southasia','westeurasia','eastasia','centralasia','america']
poplabels=['Oceania','South\nAsia','West\nEurasia','East\nAsia','Central\nAsia','America']

count=dict({})
y=dict({})
x=dict({})
for species in ['neandertal','denisova']:
    x[species]=[]
    for oneortwo in ['2','1minus2']:
        y[(species,oneortwo)]=[]

for species in ['neandertal','denisova']:
    for popul in populs:
        infile=open(species+'_'+popul+'2_snps.txt')
        lines=infile.readlines()
        infile.close()
        y[(species,'2')].append(len(lines))

        infile=open(popul+'_'+species+'_1minusboth2_snps.txt')
        lines=infile.readlines()
        infile.close()
        y[(species,'1minus2')].append(len(lines))

for i in range(len(populs)):
    x['denisova'].append(i+0.33)
    x['neandertal'].append(i)

width=0.33
colors=['blue','cyan','red','magenta']
colind=0
for species in ['neandertal','denisova']:
    plt.bar(x[species],y[(species,'2')],width,label='Proximal '+species.capitalize(),color=colors[colind])
    colind+=1
    plt.bar(x[species],y[(species,'1minus2')],width,bottom=y[(species,'2')],label='Distal '+species.capitalize(),color=colors[colind])
    colind+=1

plt.ylabel('SNP count',fontsize=15)
#plt.xticks(x['neandertal'],poplabels)
plt.xticks(x['denisova'],poplabels,fontsize=15)
plt.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
plt.legend(loc='upper right')
matplotlib.rcParams.update({'font.size':13})
#plt.tight_layout()
plt.savefig('snpcall_counts_2_1minus2.pdf',format='pdf')
