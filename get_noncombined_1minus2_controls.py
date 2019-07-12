import subprocess
import sys
from subprocess import PIPE

popul=sys.argv[1]
species=sys.argv[2]

scratch3=subprocess.Popen(['sort','-k7,7n','-k18,18n',popul+'_'+species+'_1minusboth2_snps.txt'],stdout=PIPE)
outfile=open('bstat_freq_sorted_'+popul+'_'+species+'_1minus2_snps.txt','w')
outfile.write(scratch3.communicate()[0])
outfile.close()
print 'get_combined1_snp_superset'
scratch4=subprocess.Popen(['python','get_1species_1minus2_snp_superset.py',popul,species],stdout=PIPE)
dummy=scratch4.communicate()[0]
print 'get_combined1_controls'
scratch5=subprocess.Popen(['python','get_1species_combined1_controls.py',popul, species],stdout=PIPE)
dummy=scratch5.communicate()[0]
#scratch6=subprocess.Popen(['mkdir',popul+'_1minus2_analysis'])
#dummy=scratch6.communicate()[0]
scratch7=subprocess.Popen(['python2.7','count_enhancers_v_celltype_pleio_category.py', popul,species])
dummy=scratch7.communicate()[0]
scratch8=subprocess.Popen(['python2.7','count_neand_control_celltypes.py', popul,species])
dummy=scratch8.communicate()[0]
scratch9=subprocess.Popen(['python2.7','plot_allpleio_allcelltypes_oddsNeand.py', popul,species])
dummy=scratch9.communicate()[0]



