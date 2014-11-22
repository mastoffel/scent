
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Users/Martin/Downloads/WinArl35_extracted/WinArl35/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (arp_genotypes.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 13/11/14 at 14:14:42", "arp_genotypes.htm#13_11_14at14_14_42"))
	insDoc(aux1, gLnk("R", "Settings", "arp_genotypes.htm#13_11_14at14_14_42_run_information"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "Beach_1", "arp_genotypes.htm#13_11_14at14_14_42_group0"))
		insDoc(aux2, gLnk("R", "Beach_2", "arp_genotypes.htm#13_11_14at14_14_42_group1"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "arp_genotypes.htm#13_11_14at14_14_42_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "arp_genotypes.htm#13_11_14at14_14_42_comp_sum_het"))
		insDoc(aux2, gLnk("R", "Theta(H)", "arp_genotypes.htm#13_11_14at14_14_42_comp_sum_thetaH"))
		insDoc(aux2, gLnk("R", "No. of alleles", "arp_genotypes.htm#13_11_14at14_14_42_comp_sum_numAll"))
		insDoc(aux2, gLnk("R", "Allelic range", "arp_genotypes.htm#13_11_14at14_14_42_comp_sum_allRange"))
		insDoc(aux2, gLnk("R", "Garza-Williamson index", "arp_genotypes.htm#13_11_14at14_14_42_comp_sum_GW"))
		insDoc(aux2, gLnk("R", "Garza-Williamson modified index", "arp_genotypes.htm#13_11_14at14_14_42_comp_sum_GWN"))
		aux2 = insFld(aux1, gFld("Genetic structure (samp=pop)", "arp_genotypes.htm#13_11_14at14_14_42_pop_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "arp_genotypes.htm#13_11_14at14_14_42_pop_pairw_diff"))
	aux1 = insFld(foldersTree, gFld("Run of 13/11/14 at 14:21:02", "arp_genotypes.htm#13_11_14at14_21_02"))
	insDoc(aux1, gLnk("R", "Settings", "arp_genotypes.htm#13_11_14at14_21_02_run_information"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "Beach_1", "arp_genotypes.htm#13_11_14at14_21_02_group0"))
		insDoc(aux2, gLnk("R", "Beach_2", "arp_genotypes.htm#13_11_14at14_21_02_group1"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "arp_genotypes.htm#13_11_14at14_21_02_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "arp_genotypes.htm#13_11_14at14_21_02_comp_sum_het"))
		insDoc(aux2, gLnk("R", "Theta(H)", "arp_genotypes.htm#13_11_14at14_21_02_comp_sum_thetaH"))
		insDoc(aux2, gLnk("R", "No. of alleles", "arp_genotypes.htm#13_11_14at14_21_02_comp_sum_numAll"))
		insDoc(aux2, gLnk("R", "Allelic range", "arp_genotypes.htm#13_11_14at14_21_02_comp_sum_allRange"))
		insDoc(aux2, gLnk("R", "Garza-Williamson index", "arp_genotypes.htm#13_11_14at14_21_02_comp_sum_GW"))
		insDoc(aux2, gLnk("R", "Garza-Williamson modified index", "arp_genotypes.htm#13_11_14at14_21_02_comp_sum_GWN"))
		aux2 = insFld(aux1, gFld("Genetic structure (samp=pop)", "arp_genotypes.htm#13_11_14at14_21_02_pop_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "arp_genotypes.htm#13_11_14at14_21_02_pop_pairw_diff"))
	aux1 = insFld(foldersTree, gFld("Run of 13/11/14 at 14:22:49", "arp_genotypes.htm#13_11_14at14_22_49"))
	insDoc(aux1, gLnk("R", "Settings", "arp_genotypes.htm#13_11_14at14_22_49_run_information"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "Beach_1", "arp_genotypes.htm#13_11_14at14_22_49_group0"))
		insDoc(aux2, gLnk("R", "Beach_2", "arp_genotypes.htm#13_11_14at14_22_49_group1"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "arp_genotypes.htm#13_11_14at14_22_49_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "arp_genotypes.htm#13_11_14at14_22_49_comp_sum_het"))
		insDoc(aux2, gLnk("R", "Theta(H)", "arp_genotypes.htm#13_11_14at14_22_49_comp_sum_thetaH"))
		insDoc(aux2, gLnk("R", "No. of alleles", "arp_genotypes.htm#13_11_14at14_22_49_comp_sum_numAll"))
		insDoc(aux2, gLnk("R", "Allelic range", "arp_genotypes.htm#13_11_14at14_22_49_comp_sum_allRange"))
		insDoc(aux2, gLnk("R", "Garza-Williamson index", "arp_genotypes.htm#13_11_14at14_22_49_comp_sum_GW"))
		insDoc(aux2, gLnk("R", "Garza-Williamson modified index", "arp_genotypes.htm#13_11_14at14_22_49_comp_sum_GWN"))
		aux2 = insFld(aux1, gFld("Genetic structure (samp=pop)", "arp_genotypes.htm#13_11_14at14_22_49_pop_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "arp_genotypes.htm#13_11_14at14_22_49_pop_pairw_diff"))