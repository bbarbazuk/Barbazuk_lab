import sys
import csv


def get_genes(gene_list):
    #Gene list is a text file containing gene names separated by newlines
    with open(gene_list, 'r') as file:
        genes = file.read().splitlines()
        return genes

def process_output(genes, pasa_output):
    with open(pasa_output, 'r') as file:

        base_name = pasa_output.rsplit('.', 1)[0]
        output_file = f"{base_name}_filtered.txt"

        with open(output_file, 'w') as outfile:
            for line in file:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                gene_infos = fields[2].split(',')
                for gene_info in gene_infos:
                    gene_parts = gene_info.split('_')
                    if len(gene_parts) > 1:
                        gene_name = gene_parts[1]
                        if gene_name in genes:
                            outfile.write(f"{gene_info}\n")


#usage python sort_pasa_output_by_mousegene.py <gene_list.txt> <pasa_output.txt>
def main():
    if len(sys.argv) != 3:
        print("Error: Format is sort_pasa_output_by_mousegene.py gene_list.txt pasa_output.txt")
        sys.exit(1)
    
    gene_list = sys.argv[1]
    pasa_output = sys.argv[2]
    genes = get_genes(gene_list)
    process_output(genes, pasa_output)
    
if __name__ == '__main__':
    main()
