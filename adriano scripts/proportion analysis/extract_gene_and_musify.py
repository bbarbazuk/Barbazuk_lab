import sys
import os


#Takes a file output by the regulation finder script and line by line extracts the gene name. If musify parameter is 1 will change gene name from the acomys 
#version to the Mus version using the pasa output file. If musify parameter is 0 will just output the gene name from the regulation finder script.
#Usage: python extract_gene_and_musify.py <regulation_finder_output> <pasa_output> <musify>


def process_regulation_finder_output(regulation_finder_output):
    genes_list = []
    with open(regulation_finder_output, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            gene_name = "FC705_" + line.split('_')[1] if line.startswith("FC705") else line.split('_')[0]
            if gene_name not in genes_list:
                genes_list.append(gene_name)
    return genes_list
            

def musify_gene_names(genes_list, pasa_output):
    new_genes_list = []
    gene_map = {}
    with open(pasa_output, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            lines = line.split('\t')
            acomys_gene = "FC705_" + lines[3].split('_')[1] if lines[3].startswith("FC705") else lines[3].split('_')[0]
            mus_gene = lines[2].split('_')[1]
            gene_map[acomys_gene] = mus_gene
    for gene in genes_list:
        if gene in gene_map:
            new_genes_list.append(gene_map[gene])
    return new_genes_list

def write_output(genes_list, musify, regulation_finder_output):
    base_file_name = os.path.basename(regulation_finder_output)
    output_file = f"{os.path.splitext(base_file_name)[0]}_genes.txt" if musify == 0 else f"{os.path.splitext(base_file_name)[0]}_genes_musified.txt"
    with open(output_file, 'w') as file:
        for gene in genes_list:
            file.write(gene + '\n')
    print(f"Output written to {output_file}")


def main():
    if sys.argv[1] == "-h" or sys.argv[1] == "help":
        print("Usage: python extract_gene_and_musify.py <regulation_finder_output> <pasa_output> <musify>")
        print("Extracts the gene name from the regulation finder output and optionally changes it to the Mus version using the pasa output file.")
        return
    if len(sys.argv) != 4:
        print("Error: Incorrect number of arguments.")
        print("Usage: python extract_gene_and_musify.py <regulation_finder_output> <pasa_output> <musify>")
        return
    regulation_finder_output = sys.argv[1]
    pasa_output = sys.argv[2]
    musify = sys.argv[3]
    if musify not in ["0", "1"]:
        print("Error: Musify parameter must be 0 or 1.")
        return
    musify = int(musify)
    if not regulation_finder_output.endswith('.txt'):
        print("Error: Regulation finder output file must be a .txt file.")
        return
    genes_list = process_regulation_finder_output(regulation_finder_output)
    if musify == 1:
        genes_list = musify_gene_names(genes_list, pasa_output)
    write_output(genes_list, musify, regulation_finder_output)
    

if __name__ == "__main__":
    main()