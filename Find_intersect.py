import argparse


def argument():
    parser = argparse.ArgumentParser()
    parser.add_argument("-cgun_emac", "--cgun_vs_emac_indicator_file", required=True, type=str)
    parser.add_argument("-ceso_emac", "--ceso_vs_emac_indicator_file", required=True, type=str)
    parser.add_argument("-ceso_cgun", "--ceso_vs_cgun_indicator_file", required=True, type=str)
    parser.add_argument("-spp", "--target_species", required=True, choices=['cgun', 'emac', 'ceso'], help="can only be cgun, emac, or ceso")
    # parser.add_argument("-evi", "--evidence", type=str, help="Evidence (optional)")
    # parser.add_argument("-m", "--membership_file", required=True, type=str)

    args = parser.parse_args()
    cgun_emac = args.cgun_vs_emac_indicator_file
    ceso_emac = args.ceso_vs_emac_indicator_file
    ceso_cgun = args.ceso_vs_cgun_indicator_file
    spp = args.target_species
    # evi = args.evidence
    # membership_file = args.membership_file

    return cgun_emac, ceso_emac, ceso_cgun, spp,


# def get_gene_list(membership_file):
#     """
#     Get gene list from membership file in order to check adjacent genes
#     """
#     membership_gene_list_A = []
#     membership_gene_list_B = []
#
#     for line in tqdm(membership_file):
#
#         if line.startswith("#"):
#             continue
#         line = line.strip("\n").split("\t")
#
#         gene_A = line[4]
#         gene_B = line[10]
#
#         membership_gene_list_A.append(gene_A)
#         membership_gene_list_B.append(gene_B)
#
#     return membership_gene_list_A, membership_gene_list_B


# def adjacent_ortholog_check(adjacent_gene_A, adjacent_gene_B, membership_gene_list_A, membership_gene_list_B):
#     """
#     Check adjacent genes
#     """
#     is_ortho = False  # False means NOT ORTHO
#     A_is_in = False  # False means NOT IN
#     B_is_in = False  # False means NOT IN
#
#     if adjacent_gene_A in membership_gene_list_A:
#         A_is_in = True  # True means IN
#     if adjacent_gene_B in membership_gene_list_B:
#         B_is_in = True  # True means IN
#
#     if A_is_in and B_is_in:
#         i = membership_gene_list_A.index(adjacent_gene_A)
#         if membership_gene_list_B[i] == adjacent_gene_B:
#             is_ortho = True  # True means ORTHO
#
#     return is_ortho


def parse(file_path, index):
    gene_set = set()
    line_dict = {}

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            line_dict[parts[index]] = parts
            gene_set.add(parts[index])

    return gene_set, line_dict


def get_intersection(cgun_emac, ceso_emac, ceso_cgun, spp):
    if spp == "cgun":
        gene_set_a, line_dict_a = parse(cgun_emac, 0)
        gene_set_b, line_dict_b = parse(ceso_cgun, 1)
    elif spp == "emac":
        gene_set_a, line_dict_a = parse(cgun_emac, 1)
        gene_set_b, line_dict_b = parse(ceso_emac, 1)
    else:  # "ceso"
        gene_set_a, line_dict_a = parse(ceso_emac, 0)
        gene_set_b, line_dict_b = parse(ceso_cgun, 0)
    intersection_set = gene_set_a.intersection(gene_set_b)

    return intersection_set, line_dict_a, line_dict_b


def output(intersection_set, line_dict_a, line_dict_b, spp):
    result_dict = {}
    for gene in intersection_set:
        indicator_dict = test_indicator(gene, line_dict_a, line_dict_b, spp)
        result_dict[gene] = indicator_dict

    # print(result_dict)
    return result_dict


def test_indicator(gene, line_dict_a, line_dict_b, spp):
    indicator_dict = {}

    if spp == "cgun":
        other_spp = ["emac", "ceso"]
    elif spp == "emac":
        other_spp = ["cgun", "ceso"]
    else:
        other_spp = ["emac", "cgun"]

    chr_list = []
    size_list = []
    exon_count_list = []
    exon_len_list = []
    strand_list = []

    # chr
    if line_dict_a[gene][2] == line_dict_a[gene][3]:
        chr_list.append(other_spp[0])
    if line_digict_b[gene][2] == line_dict_b[gene][3]:
        chr_list.append(other_spp[1])

    # size
    if abs(int(line_dict_a[gene][4]) - int(line_dict_a[gene][5]))/(int(line_dict_a[gene][4]) + int(line_dict_a[gene][5])) <= 0.1:
        size_list.append(other_spp[0])
    if abs(int(line_dict_b[gene][4]) - int(line_dict_b[gene][5]))/(int(line_dict_b[gene][4]) + int(line_dict_b[gene][5])) <= 0.1:
        size_list.append(other_spp[1])

    # exon_count
    if line_dict_a[gene][6] == line_dict_a[gene][7]:
        exon_count_list.append(other_spp[0])
    if line_dict_b[gene][6] == line_dict_b[gene][7]:
        exon_count_list.append(other_spp[1])

    # exon_len
    if abs(int(line_dict_a[gene][8]) - int(line_dict_a[gene][9]))/(int(line_dict_a[gene][8]) + int(line_dict_a[gene][9])) <= 0.1:
        exon_len_list.append(other_spp[0])
    if abs(int(line_dict_b[gene][8]) - int(line_dict_b[gene][9]))/(int(line_dict_a[gene][8]) + int(line_dict_a[gene][9])) <= 0.1:
        exon_len_list.append(other_spp[1])

    # strand
    if line_dict_a[gene][14] == line_dict_a[gene][15]:
        strand_list.append(other_spp[0])
    if line_dict_b[gene][14] == line_dict_b[gene][15]:
        strand_list.append(other_spp[1])

    if chr_list:
        indicator_dict['chr'] = chr_list
    if size_list:
        indicator_dict['size'] = size_list
    if exon_count_list:
        indicator_dict['exon_count'] = exon_count_list
    if exon_len_list:
        indicator_dict['exon_len'] = exon_len_list
    if strand_list:
        indicator_dict['strand'] = strand_list
    # print(indicator_dict)

    return indicator_dict


def report_percentage(result_dict, intersection_set, line_dict_a, line_dict_b):
    intersection_size = len(intersection_set)
    intersection_size_a = len(intersection_set) / len(line_dict_a)
    intersection_size_b = len(intersection_set) / len(line_dict_b)

    a = 0
    b = 0
    c = 0

    for gene, indicator_dict in result_dict.items():
        indicator_count = len(indicator_dict)
        count = indicator_count
        for indicator, spp_list in indicator_dict.items():
            if len(spp_list) != 2:
                count -= 1

        if indicator_count == 5 and indicator_count == count:
            a += 1
            # print(gene)
            # print(indicator_dict)
        elif indicator_count == 0:
            c += 1
        else:
            b += 1
    A = a/intersection_size
    B = b/intersection_size
    C = c/intersection_size

    print(intersection_size)
    print(intersection_size_a)
    print(intersection_size_b)
    print("-----------------------")
    print(A)
    print(B)
    print(C)

def main():
    cgun_emac, ceso_emac, ceso_cgun, spp = argument()
    intersection_set, line_dict_a, line_dict_b = get_intersection(cgun_emac, ceso_emac, ceso_cgun, spp)
    result_dict = output(intersection_set, line_dict_a, line_dict_b, spp)
    for genes, indicator_dict in result_dict.items():
        print(f"{genes}: {indicator_dict}")
    report_percentage(result_dict, intersection_set, line_dict_a, line_dict_b)


if __name__ == "__main__":
    main()
