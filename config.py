import os

version = "2.4"

module_dir = os.path.dirname(os.path.realpath(__file__))

cmd_bowtie = "bowtie"
cmd_R = os.path.join(module_dir, "command.R")
cmd_fourdigit = os.path.join(module_dir, "command_fourdigit.R")

bowtiebuild_class_I = os.path.join(module_dir, "references", "ClassIWithoutNQex2-3.plus75")
bowtiebuild_class_I_nonclass = os.path.join(module_dir, "references", "Class1NonClassical")
bowtiebuild_class_II = os.path.join(module_dir, "references", "HLA2.ex2.plus75")

fasta_class_I = os.path.join(module_dir, "references", "ClassIWithoutNQex2-3.plus75.fasta")
fasta_class_I_nonclass = os.path.join(module_dir, "references", "Class1NonClassical.fasta")
fasta_class_II = os.path.join(module_dir, "references", "HLA2.ex2.plus75.fasta")

dbmhc_I = os.path.join(module_dir, "HLAI.dbmhc")
dbmhc_II = os.path.join(module_dir, "HLAII.dbmhc")

class_I_list = ["A", "B", "C"]
class_I_nonclass_list = ["E", "F", "G", "H", "J", "K", "L", "P", "V"]
class_II_list = ["DQA1", "DQB1", "DRB1", "DRA", "DPA1", "DPB1"]


# HLA Class I list and according sizes
len_dict_I = {
    "A": 694,
    "B": 694,
    "C": 694
}

# HLA Class I non-class list and according sizes
len_dict_I_nonclass = {
    "E": 694,
    "F": 1089,
    "G": 694,
    "H": 1097,
    "J": 1095,
    "K": 1093,
    "L": 1096,
    "P": 732,
    "V": 538
}

# HLA Class II list and according sizes
len_dict_II = {
    "DQA1": 400,
    "DQB1": 421,
    "DRB1": 421,
    "DRA": 405,
    "DPA1": 396,
    "DPB1": 414
}
