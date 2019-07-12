# hominin-enhancers
A code repository for "Selection against archaic DNA in human regulatory regions" by Telis, Aguilar, &amp; Harris 2019

This repository contains everything you'll need to generate the plots in the main text of Telis, et al. 2019. If you need additional info, please contact Kelley Harris at harriske@uw.edu. Scripts are written to run with python 2.7.

Many of these scripts depend on the file "enhancer_states.txt" that was downloaded from the RoadMap annotation database. Due to GitHub file size constraints, this file is provided in gzipped form. You will need to decompress it before running the upstream parts of this pipeline.

All figures in the paper that summarize archaic variation (1B, 2, 3A, 3B, 4A, 4B, 5C, 5D) require you to first compile counts of Neanderthal and Denisovan variation in enhancers. Figures 1B, 2, 4A, and 4B utilize counts from specific to the Simons Genome Diversity Project (SGDP) populations, which are Westeurasia, Southasia, Eastasia, Centralasia, America, and Oceania. The summary counts you need for Figures 1, 2, and 3 are available in the repository in the following files:

  SPECIES_POPUL2_snps.txt (1)
  
  SPECIES_POPUL2_controls.txt (2)
  
where SPECIES equals neandertal or denisova and popul is one of the SGDP population names (all lowercase). Throughout the README, line numbers are provided for reference but are not intended to be typed into the commands. For Figure 4, you also need:

  POPUL_SPECIES_1minusboth2_snps.txt (3)
  
  POPUL_SPECIES_1minus2_controls.txt (4)

If you want to generate these count summaries from scratch, you will need to download the archaic call repository available at https://sriramlab.cass.idre.ucla.edu/public/sankararaman.curbio.2016/summaries.tgz and copy the resulting "summaries" directory into your clone of this "hominin-enhancers" directory. 

The following command generates files (1) and (2):

  python get_archaic_and_controls.py 2 SPECIES POPUL

The following command will generate the count summaries (3) and (4), the first script looping through all SGDP populations automatically:

  python get_old_nonshared_archaic2.py SPECIES (5)
  
  python get_noncombined_1minus2_controls.py POPUL SPECIES (6)
  
To generated Figure 5, you'll need to generate the following combined sets of variation compiled across all populations:
  SPECIES_allpops_2_snps_plus_bstat.txt
  allpops_SPECIES_2_controls.txt
  
To generate these files from the raw calls, the first step is to generate the population-specific counts that we already discussed. The second step is to run the commands

  python compile_archaic_across_populs.py SPECIES 2
  
  python get_controls_compiled_across_populs.py 2 SPECIES

Figures 5C, 5D, and 5E also utilize some site frequency spectrum data from the 1000 Genomes African individuals. The file rare_common_vars_alt_cell_types_enhancer_control.txt contains summary data utlized for Figures 5C and 5D, while rare_vars_v_celltype_pleio.txt contains summary data that go into Figure 5E. 

The plotting commands for generating the figures from the paper are as follows:

Figure 1B:

  python count_1and2_neand_control_enhancers_exons.py america
  
  python count_1and2_neand_control_enhancers_exons.py centralasia
  
  python count_1and2_neand_control_enhancers_exons.py eastasia
  
  python count_1and2_neand_control_enhancers_exons.py southasia  
  
  python count_1and2_neand_control_enhancers_exons.py westeurasia
  
  python count_1and2_neand_control_enhancers_exons.py oceania
  
  python plot_enhancers_exons_2only.py

Figure 2:

  python plot_nean_plus_oceanianden_v_pleio.py 

Figure 3A:

  python count_neand_control_celltypes.py allpops 2 neandertal
  
  python count_neand_control_celltypes.py allpops 2 denisova
  
  python count_neand_control_exons_HARs.py 2 neandertal allpops
  
  python count_neand_control_exons_HARs.py 2 denisova allpops
  
  python plot_allpops_oddsNeand_v_oddsDenisovan.py
  
Figure 3B:

  python count_enhancers_v_celltype_pleio_category.py allpops 2 neandertal
  
  python count_enhancers_v_celltype_pleio_category.py allpops 2 denisova
  
  python count_enhancers_onecelltype_plus_pleio_category.py
  
  python plot_allpops_wbrainmuscle_pleio_v_oddsNeand.py

Figure 4A:

  python count_1and2_snpcalls.py

Figure 5C:

  python plot_allpops_archaic_v_oddsrare_wExon.py denisova

Figure 5D:

  python plot_allpops_archaic_v_oddsrare_wExon.py neandertal 
  
Figure 5E:

  python plot_oddssingle_v_celltype_pleio.py
  
