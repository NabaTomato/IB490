import gzip
import argparse
from dataclasses import dataclass
from tqdm import tqdm


def argument():
    parser = argparse.ArgumentParser()
    parser.add_argument("-M", "--membership_file", required=True, type=str, help="path to membership file")
    parser.add_argument("-GA", "--gtf_file_A", required=True, type=str)
    parser.add_argument("-GB", "--gtf_file_B", required=True, type=str)

    args = parser.parse_args()
    membership_file = args.membership_file
    gtf_file_A = args.gtf_file_A
    gtf_file_B = args.gtf_file_B

    return membership_file, gtf_file_A, gtf_file_B


@dataclass
class Gene:
    gene_name_A: str
    gene_name_B: str
    chr_A: str
    chr_B: str
    size_A: int
    size_B: int
    exon_count_A: int
    exon_count_B: int
    exon_len_A: int
    exon_len_B: int
    pre_gene_A: str
    pre_gene_B: str
    nex_gene_A: str
    nex_gene_B: str
    strand_A: str
    strand_B: str
    evalue: float
    gene_func_A: str
    gene_func_B: str


def membership_parse(membership_file):
    gene_info_list = []
    counts = 0
    for line in membership_file:
        if line.startswith("#"):
            continue
        line = line.strip("\n").split("\t")

        Type = line[1]
        gene_A = line[4]
        chr_A = line[6]
        size_A = int(line[9]) - int(line[8])
        gene_B = line[10]
        chr_B = line[12]
        size_B = int(line[15]) - int(line[14])

        if Type == "gene":
            gene_info = Gene(gene_A, gene_B, chr_A, chr_B, size_A, size_B, 0, 0, 0, 0, "none", "none",
                             "none", "none", "none", "none", 0, "none", "none")
            gene_info_list.append(gene_info)
            counts += 1
            if counts == 100:
                break

    return gene_info_list


def preprocess_gene_info_list(gene_info_list):
    gene_info_dict = {}
    for gene in gene_info_list:
        mylist = []
        for prefix in ('A', 'B'):
            gene_name = getattr(gene, f'gene_name_{prefix}')
            mylist.append(gene_name)
        gene_id_tuple = tuple(mylist)
        gene_info_dict[gene_id_tuple] = gene
    return gene_info_dict


def get_exon_size_revised(gtf_file, gene_info_dict, attr_prefix):
    gtf_gene_list = []
    gene_id_to_class = {}
    if attr_prefix == "A":
        idx = 0
    else:
        idx = 1

    for gene_id_tuple in gene_info_dict:
        gene_id = gene_id_tuple[idx]
        gene_class = gene_info_dict[gene_id_tuple]
        gene_id_to_class[gene_id] = gene_class

    for line in tqdm(gtf_file):
        if line.startswith("#"):
            continue
        line = line.strip("\n").split("\t")
        types, length, strand, description = line[2], int(line[4]) - int(line[3]), line[6], line[8]
        gtf_gene_id = description.split('gene_id "')[1].split('"')[0]

        gene_class = gene_id_to_class.get(gtf_gene_id)
        if gene_class and gtf_gene_id == getattr(gene_class, f'gene_name_{attr_prefix}'):
            if types == "exon":
                exon_count = getattr(gene_class, f'exon_count_{attr_prefix}') + 1
                setattr(gene_class, f'exon_count_{attr_prefix}', exon_count)
                exon_len = getattr(gene_class, f'exon_len_{attr_prefix}') + length
                setattr(gene_class, f'exon_len_{attr_prefix}', exon_len)
            if types == "gene":
                setattr(gene_class, f'gene_func_{attr_prefix}', description)
                setattr(gene_class, f'strand_{attr_prefix}', strand)

        if types == "gene":
            gtf_gene_list.append(description)

    return gtf_gene_list


def get_pre_nex_gene(gene_info_dict, gtf_gene_list, attr_prefix):
    """
    Get the previous gene and the next gene of A and B genes
    """
    gene_id_to_class = {}
    if attr_prefix == "A":
        idx = 0
    else:
        idx = 1

    for gene_id_tuple in gene_info_dict:
        gene_id = gene_id_tuple[idx]
        gene_class = gene_info_dict[gene_id_tuple]
        gene_id_to_class[gene_id] = gene_class

    for i, gene_description in enumerate(gtf_gene_list):
        gene_description = gene_description.split('gene_id "')
        gtf_gene_id = gene_description[1].split('"')[0]

        gene_class = gene_id_to_class.get(gtf_gene_id)
        if gene_class and gtf_gene_id == getattr(gene_class, f'gene_name_{attr_prefix}'):
            if i > 0:
                previous_gene = gtf_gene_list[i - 1]
                parts = previous_gene.split('gene_id "')
                pre_gene = parts[1].split('"')[0]
                setattr(gene_class, f'pre_gene_{attr_prefix}', pre_gene)
            else:
                setattr(gene_class, f'nex_gene_{attr_prefix}', None)

            if i < len(gtf_gene_list) - 1:
                next_gene = gtf_gene_list[i + 1]
                parts = next_gene.split('gene_id "')
                nex_gene = parts[1].split('"')[0]
                setattr(gene_class, f'nex_gene_{attr_prefix}', nex_gene)
            else:
                setattr(gene_class, f'pre_gene_{attr_prefix}', None)


def gtf(membership_file, gtf_file_A, gtf_file_B):
    """
    Update info
    """
    gene_info_list = membership_parse(membership_file)
    gene_info_dict = preprocess_gene_info_list(gene_info_list)
    gtf_gene_list_A = get_exon_size_revised(gtf_file_A, gene_info_dict, 'A')
    gtf_gene_list_B = get_exon_size_revised(gtf_file_B, gene_info_dict, 'B')
    get_pre_nex_gene(gene_info_dict, gtf_gene_list_A, 'A')
    get_pre_nex_gene(gene_info_dict, gtf_gene_list_B, 'B')
    # for gene in gene_info_list:
    #     print(gene.gene_name_A + "\t" + gene.gene_name_B + "\t" + gene.chr_A + "\t" + gene.chr_B + "\t" + str(gene.size_A) + "\t" + str(gene.size_B) +
    #           "\t" + str(gene.exon_count_A) + "\t" + str(gene.exon_count_B) + "\t" + str(gene.pre_gene_A) + "\t" + str(gene.pre_gene_B) + "\t" +
    #           str(gene.nex_gene_A) + "\t" + str(gene.nex_gene_B) + "\t" + str(gene.strand_A) + "\t" + str(gene.strand_A) + "\n")

    return gene_info_list


def output(membership_file, gtf_file_A, gtf_file_B):
    membership_file = open(membership_file, "r")
    gtf_file_A = gzip.open(gtf_file_A, "rt")
    gtf_file_B = gzip.open(gtf_file_B, "rt")
    # gtf(membership_file, gtf_file_A, gtf_file_B)
    with open("/home/hanyu/Python_Practice/cgun_vs_emac_first_100_gene.tsv", "w") as file:
        gene_info_list = gtf(membership_file, gtf_file_A, gtf_file_B)
        for gene in gene_info_list:
            file.write(gene.gene_name_A + "\t" + gene.gene_name_B + "\t" + gene.chr_A + "\t" + gene.chr_B + "\t" + str(gene.size_A) + "\t" + str(gene.size_B) +
                       "\t" + str(gene.exon_count_A) + "\t" + str(gene.exon_count_B) + "\t" + str(gene.exon_len_A) + "\t" + str(gene.exon_len_B) + "\t" +
                       str(gene.pre_gene_A) + "\t" + str(gene.pre_gene_B) + "\t" + str(gene.nex_gene_A) + "\t" + str(gene.nex_gene_B) + "\t" +
                       str(gene.strand_A) + "\t" + str(gene.strand_A) + "\n")
    membership_file.close()
    gtf_file_A.close()
    gtf_file_B.close()


def main():
    membership_file, gtf_file_A, gtf_file_B = argument()
    output(membership_file, gtf_file_A, gtf_file_B)


if __name__ == "__main__":
    main()
