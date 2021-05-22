import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

proteome = pd.read_csv('prot (2).csv')

# print(len(proteome))

# analyzed_seq = ProteinAnalysis(str("NGEPRVGSLGRAFYSAPIQIWDNTTGTVASFATSFT"))
# print(analyzed_seq.molecular_weight())
# print(analyzed_seq.gravy())
# print(analyzed_seq.count_amino_acids())
# print(analyzed_seq.secondary_structure_fraction())

my_list = []
my_list_index = []

# 0
# 4
# 5
# 1466
# 3184

non_float = []
proteome["molecular_weight"] = 0

for i in range(len(proteome)):
  if not type(proteome["Protein"][i]) == str:
    non_float.append(i)

for i in range(len(proteome)):
  if i not in non_float:
    analyzed_seq = ProteinAnalysis(proteome["Protein"][i][:-1])
    my_list.append(analyzed_seq.molecular_weight())
  else:
    my_list.append(0)

for i in my_list:
  proteome["molecular_weight"] = my_list

print(proteome["molecular_weight"])
