from dataclasses import dataclass
import argparse
from tqdm import tqdm


def argument():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gene_file", required=True, type=str)
    parser.add_argument("-m", "--membership_file", required=True, type=str)

    args = parser.parse_args()
    file = args.gene_file
    membership_file = args.membership_file

    return file, membership_file


@dataclass
class Gene:
    gene_name_A: str
    gene_name_B: str
    exon_gap: int
    chr: int
    size: int
    pre_gene: str
    next_gene: str
    strand: str


def get_gene_list(membership_file):
    membership_gene_list_A = []
    membership_gene_list_B = []

    for line in tqdm(membership_file):

        if line.startswith("#"):
            continue
        line = line.strip("\n").split("\t")

        gene_A = line[4]
        gene_B = line[10]

        membership_gene_list_A.append(gene_A)
        membership_gene_list_B.append(gene_B)

    return membership_gene_list_A, membership_gene_list_B


def adjacent_ortholog_check(adjacent_gene_A, adjacent_gene_B, membership_gene_list_A, membership_gene_list_B):
    is_ortho = False  # False means NOT ORTHO
    A_is_in = False  # False means NOT IN
    B_is_in = False  # False means NOT IN
    overall_not_in = False  # False means IN

    if adjacent_gene_A in membership_gene_list_A:
        A_is_in = True  # True means IN
    if adjacent_gene_B in membership_gene_list_B:
        B_is_in = True  # True means IN
    if adjacent_gene_A not in membership_gene_list_A or adjacent_gene_B not in membership_gene_list_B:
        overall_not_in = True  # True means NOT IN

    if A_is_in and B_is_in:
        i = membership_gene_list_A.index(adjacent_gene_A)
        if membership_gene_list_B[i] == adjacent_gene_B:
            is_ortho = True  # True means ORTHO

    return is_ortho, A_is_in, B_is_in, overall_not_in


def parse(file, membership_file):
    file_handler = open(file, "r")
    line_count = 0
    membership_file_handler = open(membership_file, "r")
    global_dict = {}
    list_A, list_B = get_gene_list(membership_file_handler)

    for line in tqdm(file_handler):
        line_count += 1
        line = line.strip("\n").split("\t")
        chr_A = str(line[2])
        chr_B = str(line[3])
        size_A = int(line[4])
        size_B = int(line[5])
        exon_count_A = int(line[6])
        exon_count_B = int(line[7])
        exon_len_A = int(line[8])
        exon_len_B = int(line[9])
        pre_gene_A = str(line[10])
        pre_gene_B = str(line[11])
        nex_gene_A = str(line[12])
        nex_gene_B = str(line[13])
        strand_A = str(line[14])
        strand_B = str(line[15])

        if chr_A == chr_B:
            key = "chr"
            global_dict[key] = global_dict.get(key, 0) + 1

        if abs(size_A - size_B)/(size_A + size_B) <= 0.1:
            key = "gene_size"
            global_dict[key] = global_dict.get(key, 0) + 1

        if exon_count_A == exon_count_B:
            key = "exon_count"
            global_dict[key] = global_dict.get(key, 0) + 1

        if abs(exon_len_A - exon_len_B)/(exon_len_A + exon_len_B) <= 0.1:
            key = "exon_size"
            global_dict[key] = global_dict.get(key, 0) + 1

        if strand_A == strand_B:
            key = "strand"
            global_dict[key] = global_dict.get(key, 0) + 1

        is_ortho_pre, A_is_in_pre, B_is_in_pre, overall_not_in_pre = adjacent_ortholog_check(pre_gene_A, pre_gene_B, list_A, list_B)
        if is_ortho_pre:
            key = "pre_ortho"
            global_dict[key] = global_dict.get(key, 0) + 1
        if A_is_in_pre is False:  # if preA is not in membership file
            key = "pre_A_not_in_mem"
            global_dict[key] = global_dict.get(key, 0) + 1
        if B_is_in_pre is False:  # if preB is not in membership file
            key = "pre_B_not_in_mem"
            global_dict[key] = global_dict.get(key, 0) + 1
        if overall_not_in_pre:  # if preB is in membership file
            key = "pre_gene_not_in_mem"
            global_dict[key] = global_dict.get(key, 0) + 1

        is_ortho_nex, A_is_in_nex, B_is_in_nex, overall_not_in_nex = adjacent_ortholog_check(nex_gene_A, nex_gene_B, list_A, list_B)
        if is_ortho_nex:
            key = "nex_ortho"
            global_dict[key] = global_dict.get(key, 0) + 1
        if A_is_in_nex is False:
            key = "nex_A_not_in_mem"
            global_dict[key] = global_dict.get(key, 0) + 1
        if B_is_in_nex is False:
            key = "nex_B_not_in_mem"
            global_dict[key] = global_dict.get(key, 0) + 1
        if overall_not_in_nex:
            key = "nex_gene_not_in_mem"
            global_dict[key] = global_dict.get(key, 0) + 1

    file_handler.close()
    membership_file_handler.close()


    return global_dict, line_count


def add_not_ortho_gene(file, membership_file):
    global_dict, line_count = parse(file, membership_file)
    pre_not_ortho = line_count - global_dict["pre_gene_not_in_mem"] - global_dict["pre_ortho"]
    nex_not_ortho = line_count - global_dict["nex_gene_not_in_mem"] - global_dict["nex_ortho"]
    global_dict["pre_not_ortho"] = pre_not_ortho
    global_dict["nex_not_ortho"] = nex_not_ortho
    global_dict_percent = {}
    for evidence, agreement in global_dict.items():
        global_dict_percent[evidence] = int(agreement)/int(line_count) * 100
    print(global_dict_percent)
    return global_dict_percent



def main():
    file, membership_file = argument()
    add_not_ortho_gene(file, membership_file)


if __name__ == "__main__":
    main()
