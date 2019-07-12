import subprocess
import sys
from subprocess import PIPE

popul=sys.argv[1]

#scratch2=subprocess.Popen(['sort','-k2,2n','-k4,4n',popul+'_combined_archaic_1minus2.txt'],stdout=PIPE)
#outfile=open(popul+'_combined_archaic_1minus2_snps.txt','w')
#outfile.write(scratch2.communicate()[0])
#outfile.close()
scratch3=subprocess.Popen(['sort','-k7,7n','-k18,18n',popul+'_combined_archaic_1minus2_snps.txt'],stdout=PIPE)
outfile=open('bstat_freq_sorted_'+popul+'_combined_archaic_1minus2.txt','w')
outfile.write(scratch3.communicate()[0])
outfile.close()
print 'get_combined1_snp_superset'
scratch4=subprocess.Popen(['python','get_combined1_snp_superset.py',popul],stdout=PIPE)
dummy=scratch4.communicate()[0]
print 'get_combined1_controls'
scratch5=subprocess.Popen(['python','get_combined1_controls.py',popul],stdout=PIPE)
dummy=scratch5.communicate()[0]
scratch6=subprocess.Popen(['mkdir',popul+'_1minus2_analysis'])
dummy=scratch6.communicate()[0]
scratch7=subprocess.Popen(['python2.7','count_enhancers_v_celltype_pleio_category.py', popul])
dummy=scratch7.communicate()[0]
scratch8=subprocess.Popen(['python2.7','count_neand_control_celltypes.py', popul])
dummy=scratch8.communicate()[0]
scratch9=subprocess.Popen(['python2.7','plot_allpleio_allcelltypes_oddsNeand.py', popul])
dummy=scratch9.communicate()[0]



