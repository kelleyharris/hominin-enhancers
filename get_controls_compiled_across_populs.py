import subprocess
import sys
from subprocess import PIPE

[oneortwo,species]=sys.argv[1:3]

infile=open(species+'_allpops_'+oneortwo+'_snps_plus_bstat.txt')
lines=infile.readlines()
infile.close()

populs=[]
output=dict({})
for line in lines:
    s=line.strip('\n').split('\t')
    if len(s)>1:
        if not s[-1] in populs:
            populs.append(s[-1])
            output[s[-1]]=''
        output[s[-1]]+=line

print populs

for popul in populs:
    outfile=open('allpop_component_'+species+'_'+popul+oneortwo+'_snps_plus_bstat.txt','w')
    outfile.write(output[popul])
    outfile.close()

    scratch3=subprocess.Popen(['sort','-k7,7n','-k18,18n','allpop_component_'+species+'_'+popul+oneortwo+'_snps_plus_bstat.txt'],stdout=PIPE)
    outfile=open('bfsorted_allpop_component_'+popul+'_'+species+'_'+oneortwo+'_snps.txt','w')
    outfile.write(scratch3.communicate()[0])
    outfile.close()

    if oneortwo=='1minus2':
        print 'get_combined1_snp_superset'
        scratch4=subprocess.Popen(['python','get_1species_1minus2_snp_superset.py',popul,species,'allpops'],stdout=PIPE)
        dummy=scratch4.communicate()[0]
        print 'get_combined1_controls'
        scratch5=subprocess.Popen(['python','get_1species_combined1_controls.py',popul, species,'allpops'],stdout=PIPE)
        dummy=scratch5.communicate()[0]
    else:
        scratch4=subprocess.Popen(['python','get_control_snp_superset.py',oneortwo,species,popul,'allpops'],stdout=PIPE)
        dummy=scratch4.communicate()[0]
        scratch5=subprocess.Popen(['python','get_controls.py',oneortwo,species,popul,'allpops'],stdout=PIPE)
        dummy=scratch5.communicate()[0]


catline=['cat']
for popul in populs:
    if oneortwo=='1minus2':
        catline.append('allpop_component_'+popul+'_'+species+'_1minus2_controls.txt')
    else:
        catline.append('allpop_component_'+species+'_'+popul+oneortwo+'_controls.txt')

scratch6=subprocess.Popen(catline,stdout=PIPE)
outfile=open('temp_'+oneortwo+'_'+species+'.txt','w')
outfile.write(scratch6.communicate()[0])
outfile.close()
scratch7=subprocess.Popen(['sort','-k2,2n','-k4,4n','temp_'+oneortwo+'_'+species+'.txt'],stdout=subprocess.PIPE)
outfile=open('allpops_'+species+'_'+oneortwo+'_controls.txt','w')
outfile.write(scratch7.communicate()[0])
outfile.close()    
