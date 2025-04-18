import sys
import csv


def get_genes(gene_list):
    #Gene list is a text file containing gene names separated by newlines
    with open(gene_list, 'r') as file:
        genes = file.read().splitlines()
        return genes

def process_output(genes, pasa_output, prefix, keep_prefix):
    with open(pasa_output, 'r') as file:
        debug_count = 0
        base_name = pasa_output.rsplit('.', 1)[0]
        output_file = f"{base_name}_filtered_{prefix}.txt"

        with open(output_file, 'w') as outfile:
            for line in file:
                if line.startswith('#'):
                    continue                
                fields = line.strip().split('\t')

                index = None
                if fields[2].startswith(prefix):
                    index = 2
                elif fields[3].startswith(prefix):
                    index = 3
                if index is None:
                    print("Error parsing gene name from PASA output")
                    sys.exit(1)

                gene_infos = fields[index].split(',')
                for gene_info in gene_infos:
                    gene_parts = gene_info.split('_')
                    if len(gene_parts) > 1:
                        gene_name = None
                        if keep_prefix:
                            gene_name = gene_parts[0] + "_" + gene_parts[1]
                        else:
                            gene_name = gene_parts[1]
                        if gene_name in genes:
                            outfile.write(f"{gene_info}\n")
                            debug_count += 1
        print(f"Isolated {debug_count} events from {pasa_output} that correspond to input genes")


#usage python sort_pasa_output_by_gene.py <gene_list.txt> <pasa_output.txt> <species_prefix>
def main():
    if len(sys.argv) == 2 and (sys.argv[1] == "-h" or sys.argv[1] == "--help"):
        print("given a list of newline-separated gene names and a pasa event file, this script will output the events in the pasa output file that correspond to the genes in the list of species corresponding to the prefix (ie, Mus or FC705). If keep_prefix is indicated, the prefix should be present in the input file and will be kept in the output file. Otherwise, the prefix will be removed. The output file will be named <pasa_output>_filtered.txt")
        print("Usage: python sort_pasa_output_by_gene.py <gene_list.txt> <pasa_output.txt> <species_prefix> <keep_prefix>")
        sys.exit(1)
    if len(sys.argv) != 5:
        print("Error: Format is sort_pasa_output_by_gene.py <gene_list.txt> <pasa_output.txt> <species_prefix> <keep_prefix (0 or 1)>")
        sys.exit(1)
    
    gene_list = sys.argv[1]
    pasa_output = sys.argv[2]
    prefix = sys.argv[3]
    keep_prefix = int(sys.argv[4])
    
    genes = get_genes(gene_list)
    process_output(genes, pasa_output, prefix, keep_prefix)
    
if __name__ == '__main__':
    main()
