# use PLINK to detect possible linkage between selected markers

# default settings: 1000 kbp window, omitted r2 < 0.1 from ld report
# note markers listed as potentially linked below should not be used in tandem
plink --vcf ../LDtest_NAMED.vcf.gz --make-founders -const-fid --allow-extra-chr --r2 --ld-window-r2 0.1 --out candidate_PLINK_LD
# CHR_A         BP_A                   SNP_A  CHR_B         BP_B                   SNP_B           R2 
# scaffold_1      3515865  scaffold_1_3515865_A_G scaffold_1      3626953  scaffold_1_3626953_G_A     0.114403 
# scaffold_1     25311758 scaffold_1_25311758_G_A scaffold_1     26305924 scaffold_1_26305924_A_G     0.482621 
# scaffold_2     13443561 scaffold_2_13443561_G_A scaffold_2     14360108 scaffold_2_14360108_C_T     0.177941 
# scaffold_3     14046673 scaffold_3_14046673_G_T scaffold_3     14046676 scaffold_3_14046676_A_C            1 
# scaffold_5     16132312 scaffold_5_16132312_A_G scaffold_5     16132319 scaffold_5_16132319_T_G     0.887139 
# scaffold_5     16132312 scaffold_5_16132312_A_G scaffold_5     16909120 scaffold_5_16909120_T_G     0.177941 
# scaffold_5     16132319 scaffold_5_16132319_T_G scaffold_5     16909120 scaffold_5_16909120_T_G        0.225 
# scaffold_7      6810378  scaffold_7_6810378_T_G scaffold_7      7632091  scaffold_7_7632091_T_C     0.106838 
# scaffold_8     14118570 scaffold_8_14118570_C_T scaffold_8     14399989 scaffold_8_14399989_C_A     0.185357 
# scaffold_8     14118570 scaffold_8_14118570_C_T scaffold_8     14466046 scaffold_8_14466046_C_T     0.141651 
# scaffold_8     14218346 scaffold_8_14218346_G_A scaffold_8     14532835 scaffold_8_14532835_A_G     0.100833 
# scaffold_8     14466046 scaffold_8_14466046_C_T scaffold_8     14532835 scaffold_8_14532835_A_G     0.169183 
