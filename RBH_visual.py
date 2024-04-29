import cairo
import argparse
from dataclasses import dataclass


def argument():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--chromosome_file")
    parser.add_argument("-g", "--membership_file")
    parser.add_argument("-w", "--width", default=1, type=float)
    parser.add_argument("-l", "--length", default=5, type=float)
    parser.add_argument("-x", "--x", help="X coordinate of the left side of the square", default=10, type=float)
    parser.add_argument("-y", "--y", help="Y coordinate of the the top side of the square", default=10, type=float)
    parser.add_argument("-lw", "--line_width", help="a line width", default=0.5, type=float)
    parser.add_argument("-t", "--threshold", help="the threshold of chromosome length", default=10000000, type=int)
    myargs = parser.parse_args()
    chromosome_file = myargs.chromosome_file
    gene_file = myargs.membership_file
    width = myargs.width
    length = myargs.length
    x_start = myargs.x
    y_start = myargs.y
    line_width = myargs.line_width
    threshold = myargs.threshold

    return chromosome_file, gene_file, width, length, x_start, y_start, line_width, threshold


@dataclass
class Chromosome:
    org: str
    chr_id: str
    chr_length: int


@dataclass
class Gene:
    chr_id: str
    gene_id: str
    start_bp: int
    end_bp: int
    gene_center_x: float


def length_limit(length, threshold):
    """
    Threshold
    """
    return int(length) >= threshold


def chr_file_parse(chromosome_file, threshold):
    """
    parse chromosome file
    """
    chr_A_list = []
    chr_B_list = []
    filter_list = []
    total_bp_dict = dict()

    for line in chromosome_file:

        if line.startswith("#"):
            continue
        line = line.strip("\n").split("\t")

        org = line[0]  # Organism name.
        chr_name = line[1]  # Chromosome name.
        chr_length = int(line[2])  # Chromosome length.

        # Check if organism is targeted and length meets the threshold.
        if org == "ceso.ftc.fv8" and length_limit(chr_length, threshold):
            if org not in filter_list:
                total_bp_dict[org] = 0
                filter_list.append(org)
            if chr_name not in chr_A_list:
                chr = Chromosome(org, chr_name, chr_length)
                chr_A_list.append(chr)
                total_bp_dict[org] += chr_length
        if org == "cgun.ftc.fv8" and length_limit(chr_length, threshold):
            if org not in filter_list:
                total_bp_dict[org] = 0
                filter_list.append(org)
            if chr_name not in chr_B_list:
                chr = Chromosome(org, chr_name, chr_length)
                chr_B_list.append(chr)
                total_bp_dict[org] += chr_length


    return chr_A_list, chr_B_list, total_bp_dict


def gene_file_parse(gene_file, chr_A_list, chr_B_list):
    """
    parse gene file
    """
    geneA_list = []
    geneB_list = []
    for line in gene_file:
        if line.startswith("#"):
            continue
        line = line.strip("\n").split("\t")

        Type = line[1]
        gene_name_A = line[4]
        chr_name_A = line[6]
        gene_name_B = line[10]
        chr_name_B = line[12]
        start_bp_A = int(line[8])
        end_bp_A = int(line[9])
        start_bp_B = int(line[14])
        end_bp_B = int(line[15])

        flag_A = False
        flag_B = False

        for chromosomes in chr_A_list:
            if chr_name_A == chromosomes.chr_id and not flag_A and Type == "gene":
                geneA = Gene(chr_name_A, gene_name_A, start_bp_A, end_bp_A, 0)
                geneA_list.append(geneA)
                flag_A = True
        for chromosomes in chr_B_list:
            if chr_name_B == chromosomes.chr_id and not flag_B and Type == "gene":
                geneB = Gene(chr_name_B, gene_name_B, start_bp_B, end_bp_B, 0)
                geneB_list.append(geneB)
                flag_B = True

        if not flag_A:
            geneA_list.append(None)
        if not flag_B:
            geneB_list.append(None)

    return geneA_list, geneB_list


def sort_chr_list(geneA_list, geneB_list, chr_B_list):
    chr_dict = {}
    for i in range(len(geneA_list)):
        if geneA_list[i] is not None and geneB_list[i] is not None:
            if geneA_list[i].chr_id not in chr_dict:
                chr_dict[geneA_list[i].chr_id] = {}
            chr_dict[geneA_list[i].chr_id].setdefault(geneB_list[i].chr_id, 0)
            chr_dict[geneA_list[i].chr_id][geneB_list[i].chr_id] += 1

    sorted_chr = []
    for chr_A_id, chr_B_id_count_dict in chr_dict.items():
        current_max = 0
        matched_chr = None
        for chr in chr_B_id_count_dict.keys():
            if current_max <= chr_B_id_count_dict[chr]:
                matched_chr = chr
                current_max = chr_B_id_count_dict[chr]
        sorted_chr.append(matched_chr)

    sorted_chr_B_list = []
    for chr in sorted_chr:
        for chromosomes in chr_B_list:
            if chromosomes.chr_id == chr and chromosomes not in sorted_chr_B_list:
                sorted_chr_B_list.append(chromosomes)

    return sorted_chr_B_list


def get_color(num):
    color_list = [[0.8784, 0.0667, 0.3725],
                  [0.0588, 0.3216, 0.7294],
                  [0.3140, 0.7840, 0.4710],
                  [0.5569, 0.8667, 0.1333],
                  [0.3216, 0.8784, 0.4980],
                  [0.9333, 0.4353, 0.0902]]
    color = color_list[num % len(color_list)]

    return color


def draw_square(context, width, length, x_start, y_start, line_width, color):
    context.set_line_width(line_width)
    context.move_to(x_start, y_start)
    context.line_to(x_start + width, y_start)
    context.line_to(x_start + width, y_start + length)
    context.line_to(x_start, y_start + length)
    context.set_source_rgb(color[0], color[1], color[2])
    context.close_path()

    context.stroke()


def draw_text(context, text, x_start, y_start):
    context.set_source_rgb(0, 0, 0)
    context.move_to(x_start, y_start)
    context.show_text(text)

    context.stroke()


def calculate_gene_position(gene_start_bp, gene_end_bp, chr_start_x, chr_start_y, pixel_per_bp):
    gene_start_x = chr_start_x + (gene_start_bp * pixel_per_bp)
    gene_end_x = chr_start_x + (gene_end_bp * pixel_per_bp)
    gene_y = chr_start_y
    gene_center_x = (gene_start_x + gene_end_x) / 2
    return gene_center_x, gene_y


def draw_gene_position(context, gene_start_bp, gene_end_bp, chr_start_x, chr_start_y, pixel_per_bp, length, color):
    gene_center_x, gene_y = calculate_gene_position(gene_start_bp, gene_end_bp, chr_start_x, chr_start_y, pixel_per_bp)
    context.set_line_width(0.1)
    context.move_to(gene_center_x, gene_y)
    context.line_to(gene_center_x, gene_y + length)
    context.set_source_rgb(color[0], color[1], color[2])
    context.stroke()

    return gene_center_x


def gene_position(context, gene_list, chr_id, current_x, current_y, pixel_per_bp, length, color):
    for genes in gene_list:
        if genes is not None and genes.chr_id == chr_id:
            gene_start_x = genes.start_bp
            gene_end_x = genes.end_bp
            gene_center_x = draw_gene_position(context, gene_start_x, gene_end_x, current_x, current_y,
                                               pixel_per_bp, length, color)
            genes.gene_center_x = gene_center_x


def draw_connection(context, start_x, start_y, end_x, end_y):
    context.set_line_width(0.1)
    context.move_to(start_x, start_y)
    context.line_to(end_x, end_y)
    context.set_source_rgb(0, 0, 0)
    context.stroke()


def draw_multiple_squares(chromosome_file, gene_file, length, x_start, y_start, line_width, threshold):
    chr_file_handler = open(chromosome_file, "r")
    gene_file_handler = open(gene_file, "r")
    chr_A_list, chr_B_list, total_bp_dict = chr_file_parse(chr_file_handler, threshold)
    gene_A_list, gene_B_list = gene_file_parse(gene_file_handler, chr_A_list, chr_B_list)
    sorted_chr_B_list = sort_chr_list(gene_A_list, gene_B_list, chr_B_list)

    surface = cairo.SVGSurface("ceso_vs_cgun_21.svg", 1000, 1000)
    context = cairo.Context(surface)

    context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(10)

    current_x = x_start
    current_y = y_start
    spacing = line_width
    filter_list = []

    for i, chromosomes in enumerate(chr_A_list):
        if chromosomes.org not in filter_list:
            draw_text(context, chromosomes.org, current_x, current_y + 5)
            filter_list.append(chromosomes.org)
            current_x += 70
        if chromosomes.chr_id != "21":
            continue
        # total_bp = total_bp_dict[chromosomes.org]
        total_bp = 96460698
        pixel_per_bp = (1000 - 200) / total_bp
        pixel_per_chr = pixel_per_bp * chromosomes.chr_length - spacing
        width = pixel_per_chr

        color = get_color(i)

        draw_square(context, width, length, current_x, current_y, line_width, color)

        draw_text(context, chromosomes.chr_id, current_x + width/2, current_y - 2)

        gene_position(context, gene_A_list, chromosomes.chr_id, current_x, current_y, pixel_per_bp, length, color)

        current_x += width + spacing

    current_x = x_start
    current_y += 100

    for i, chromosomes in enumerate(sorted_chr_B_list):
        if chromosomes.org not in filter_list:
            draw_text(context, chromosomes.org, current_x, current_y + 5)
            filter_list.append(chromosomes.org)
            current_x += 70
        if chromosomes.chr_id != "21":
            continue
        # total_bp = total_bp_dict[chromosomes.org]
        total_bp = 96460698
        pixel_per_bp = (1000 - 200) / total_bp
        pixel_per_chr = pixel_per_bp * chromosomes.chr_length - spacing
        width = pixel_per_chr

        color = get_color(i)

        draw_square(context, width, length, current_x, current_y, line_width, color)

        draw_text(context, chromosomes.chr_id, current_x + width/2, current_y + 15)

        gene_position(context, gene_B_list, chromosomes.chr_id, current_x, current_y, pixel_per_bp, length, color)

        current_x += width + spacing

    for j in range(len(gene_A_list)):
        if gene_A_list[j] is not None and gene_B_list[j] is not None:
            draw_connection(context, gene_A_list[j].gene_center_x, y_start + length, gene_B_list[j].gene_center_x,
                            y_start + 100)

    surface.finish()
    surface.flush()

    gene_file_handler.close()
    chr_file_handler.close()


def main():
    chromosome_file, gene_file, width, length, x_start, y_start, line_width, threshold = argument()
    draw_multiple_squares(chromosome_file, gene_file, length, x_start, y_start, line_width, threshold)


if __name__ == "__main__":
    main()
