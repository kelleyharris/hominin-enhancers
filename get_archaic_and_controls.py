import subprocess
import sys
from subprocess import PIPE

[oneortwo,species,popul]=sys.argv[1:4]

scratch1=subprocess.Popen(['python','get_private_neand_snps.py',oneortwo,species,popul],stdout=PIPE)
dummy=scratch1.communicate()[0]
scratch2=subprocess.Popen(['python','add_Bstats_to_neandvars.py',oneortwo,species,popul],stdout=PIPE)
dummy=scratch2.communicate()[0]
scratch3=subprocess.Popen(['sort','-k7,7n','-k18,18n',species+'_'+popul+oneortwo+'_snps_plus_bstat.txt'],stdout=PIPE)
outfile=open('bstat_freq_sorted_'+species+'_'+popul+oneortwo+'_snps.txt','w')
outfile.write(scratch3.communicate()[0])
outfile.close()
scratch4=subprocess.Popen(['python','get_control_snp_superset.py',oneortwo,species,popul],stdout=PIPE)
dummy=scratch4.communicate()[0]
scratch5=subprocess.Popen(['python','get_controls.py',oneortwo,species,popul],stdout=PIPE)
dummy=scratch5.communicate()[0]
scratch6=subprocess.Popen(['mkdir',popul+'_'+species+oneortwo+'_analysis'])
dummy=scratch6.communicate()[0]
scratch7=subprocess.Popen(['python2.7','count_enhancers_v_celltype_pleio_category.py', oneortwo, species,popul])
dummy=scratch7.communicate()[0]
scratch8=subprocess.Popen(['python2.7','count_neand_control_celltypes.py', oneortwo, species,popul])
dummy=scratch8.communicate()[0]
#scratch9=subprocess.Popen(['python2.7','plot_allpleio_allcelltypes_oddsNeand.py', oneortwo, species,popul])
#dummy=scratch9.communicate()[0]



