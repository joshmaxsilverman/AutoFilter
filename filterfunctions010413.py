import csv
import string
import qMS
import sets
import numpy
import matplotlib.pyplot as plt

datapath = "/Users/joshsilverman/Dropbox/Research/AutoFilter/gogat_isocsv/gogat40_iso.csv"
data = list( csv.reader( open(datapath, 'rU') ) )
header = data[0]

# 1 - Find the indices for the quantities of interest using list comprehensions

ampu_index = [index for index, item in enumerate(header) if item == "AMP_U"][0]
ampl_index = [index for index, item in enumerate(header) if item == "AMP_L"][0]
isoz_charge_index = [index for index, item in enumerate(header) if item == "isoz_charge"][0]
isopep_index = [index for index, item in enumerate(header) if item == "isopep"][0]
protein_index = [index for index, item in enumerate(header) if item == "protein"][0]

# 2 - Declare the peptide list and the protein set

peptide_list = []
protein_set  = set()
peptide_set = set()

# 3 - Loop over data set, collect amplitudes, charge state, peptide sequence, protein id into protein_set

for line in data[1:]:
	
	ampu = float(line[ampu_index])
	ampl = float(line[ampl_index])
	isoz_charge = int(line[isoz_charge_index])
	isopep = line[isopep_index]
	protein = line[protein_index]
	
	identifier = [isopep, isoz_charge, protein, ampu/(ampu/ampl)]

	protein_set.add(protein)
	peptide_set.add(isopep)
	peptide_list.append(identifier)
	
# -- Scoring Functions 

# SF1 - Loop over proteins, calculate MAD for the peptides in each group

peptide_MAD_dict = {}

# print list(peptide_set)


# for peptide in list(peptide_set):
# 	temp_peptides = filter(lambda item: item[0] == peptide, peptide_list)
# 	
# 	
# 	temp_fraclabs = map(lambda item: item[3], temp_peptides)
# # 	print qMS.MAD(temp_fraclabs)
# 	peptide_MAD_dict.setdefault(peptide, qMS.MAD(temp_fraclabs) )	

for protein in list(protein_set):
	temp_peptides = filter(lambda item: item[2] == protein, peptide_list)
	temp_fraclabs = map(lambda item: item[3], temp_peptides)
	temp_protein = temp_peptides[0][2]

	peptide_MAD_dict.setdefault(temp_protein , qMS.MAD(temp_fraclabs) )	
	peptide_MAD_dict.setdefault(temp_protein , [qMS.MAD(temp_fraclabs),  len(temp_fraclabs)])



plt.hist(peptide_MAD_dict.values(), bins = 160)
plt.xlim([0, 0.4])
# plt.ylim([0, 500])
plt.show()

# print peptide_MAD_dict.values()