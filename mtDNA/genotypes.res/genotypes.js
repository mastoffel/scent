
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Program%20Files/Arlequin/WinArl35/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (genotypes.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 16/11/14 at 14:37:16", "genotypes.htm#16_11_14at14_37_16"))
	insDoc(aux1, gLnk("R", "Settings", "genotypes.htm#16_11_14at14_37_16_run_information"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "Beach_1", "genotypes.htm#16_11_14at14_37_16_group0"))
		insDoc(aux2, gLnk("R", "Beach_2", "genotypes.htm#16_11_14at14_37_16_group1"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "genotypes.htm#16_11_14at14_37_16_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "genotypes.htm#16_11_14at14_37_16_comp_sum_het"))
		insDoc(aux2, gLnk("R", "Theta(H)", "genotypes.htm#16_11_14at14_37_16_comp_sum_thetaH"))
		insDoc(aux2, gLnk("R", "No. of alleles", "genotypes.htm#16_11_14at14_37_16_comp_sum_numAll"))
		aux2 = insFld(aux1, gFld("Genetic structure (samp=pop)", "genotypes.htm#16_11_14at14_37_16_pop_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "genotypes.htm#16_11_14at14_37_16_pop_pairw_diff"))
	aux1 = insFld(foldersTree, gFld("Run of 16/11/14 at 14:42:45", "genotypes.htm#16_11_14at14_42_45"))
	insDoc(aux1, gLnk("R", "Settings", "genotypes.htm#16_11_14at14_42_45_run_information"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "Beach_1", "genotypes.htm#16_11_14at14_42_45_group0"))
		insDoc(aux2, gLnk("R", "Beach_2", "genotypes.htm#16_11_14at14_42_45_group1"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "genotypes.htm#16_11_14at14_42_45_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "genotypes.htm#16_11_14at14_42_45_comp_sum_het"))
		insDoc(aux2, gLnk("R", "Theta(H)", "genotypes.htm#16_11_14at14_42_45_comp_sum_thetaH"))
		insDoc(aux2, gLnk("R", "No. of alleles", "genotypes.htm#16_11_14at14_42_45_comp_sum_numAll"))
		aux2 = insFld(aux1, gFld("Genetic structure (samp=pop)", "genotypes.htm#16_11_14at14_42_45_pop_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "genotypes.htm#16_11_14at14_42_45_pop_pairw_diff"))
		insDoc(aux2, gLnk("R", "Exact tests", "genotypes.htm#16_11_14at14_42_45_pop_exct_tests"))
