
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///C:/Users/Martin/Downloads/WinArl35_extracted/WinArl35/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (arp2_genotypes.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 13/11/14 at 14:10:28", "arp2_genotypes.htm#13_11_14at14_10_28"))
	insDoc(aux1, gLnk("R", "Settings", "arp2_genotypes.htm#13_11_14at14_10_28_run_information"))
		aux2 = insFld(aux1, gFld("Samples", ""))
		insDoc(aux2, gLnk("R", "Beach_1", "arp2_genotypes.htm#13_11_14at14_10_28_group0"))
		insDoc(aux2, gLnk("R", "Beach_2", "arp2_genotypes.htm#13_11_14at14_10_28_group1"))
		aux2 = insFld(aux1, gFld("Within-samples summary", ""))
		insDoc(aux2, gLnk("R", "Basic indices", "arp2_genotypes.htm#13_11_14at14_10_28_comp_sum_Basic"))
		insDoc(aux2, gLnk("R", "Heterozygosity", "arp2_genotypes.htm#13_11_14at14_10_28_comp_sum_het"))
		insDoc(aux2, gLnk("R", "Theta(H)", "arp2_genotypes.htm#13_11_14at14_10_28_comp_sum_thetaH"))
		insDoc(aux2, gLnk("R", "No. of alleles", "arp2_genotypes.htm#13_11_14at14_10_28_comp_sum_numAll"))
