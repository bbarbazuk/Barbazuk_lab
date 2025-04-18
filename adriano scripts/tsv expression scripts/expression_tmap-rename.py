import sys
import csv

#Tries to reconcile possible cufflink matches using a .tmap file (should get from pipeline) in expression data (should get from rsem)
#if one to one match, renames the gene/isoform in the .tsv file directly from cufflink name to reference name
#if one to many, renames the gene/isoform in the .tsv file directly from cufflink name to reference name + ?x_cufflinkname where x is priority code
#check cuffcompare documentation for priority codes, https://cole-trapnell-lab.github.io/cufflinks/cuffcompare/
def process_gene_tmap(pasa_cluster_map, tmap_file):
    with open(tmap_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')

        for row in reader:
            if row['class_code'] == '=':
                ref_gene_id = row['ref_gene_id']
                cuff_gene_id = row['cuff_gene_id']
                if cuff_gene_id and ref_gene_id:
                    pasa_cluster_map[cuff_gene_id] = ref_gene_id
            elif row['class_code'] == 'j':
                cuff_gene_id = row.get('cuff_gene_id')
                ref_gene_id = row.get('ref_gene_id')
                if cuff_gene_id and ref_gene_id:
                    pasa_cluster_map[cuff_gene_id] = ref_gene_id + '?j_' + cuff_gene_id
            elif row['class_code'] == 'o':
                cuff_gene_id = row.get('cuff_gene_id')
                ref_gene_id = row.get('ref_gene_id')
                if cuff_gene_id and ref_gene_id:
                    pasa_cluster_map[cuff_gene_id] = ref_gene_id + '?o_' + cuff_gene_id
            elif row['class_code'] == 'u':
                cuff_gene_id = row.get('cuff_gene_id')
                pasa_cluster_map[cuff_gene_id] = 'NEW?_' + cuff_gene_id

def process_isoform_tmap(pasa_cluster_map, tmap_file):
    with open(tmap_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')

        for row in reader:
            if row['class_code'] == '=':
                cuff_id = row.get('cuff_id')
                ref_id = row.get('ref_id')
                if cuff_id and ref_id:
                    pasa_cluster_map[cuff_id] = ref_id
            elif row['class_code'] == 'j':
                cuff_id = row.get('cuff_id')
                ref_id = row.get('ref_id')
                if cuff_id and ref_id:
                    pasa_cluster_map[cuff_id] = ref_id + '?_' + cuff_id
            elif row['class_code'] == 'o':
                cuff_id = row.get('cuff_id')
                ref_id = row.get('ref_id')
                if cuff_id and ref_id:
                    pasa_cluster_map[cuff_id] = ref_id + '?o_' + cuff_id
            elif row['class_code'] == 'u':
                cuff_id = row.get('cuff_id')
                pasa_cluster_map[cuff_id] = 'NEW?_' + cuff_id

def process_gene_tsv(pasa_cluster_map, tsv_file):
    debug_counter = 0
    with open(tsv_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')

        fieldnames = reader.fieldnames
        base_name = tsv_file.rsplit('.', 1)[0]
        output_file = f"{base_name}_gene_renamed.tsv"

        with open(output_file, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for row in reader:
                original_name = row['Name']

                if original_name in pasa_cluster_map:
                    row['Name'] = pasa_cluster_map[original_name]
                    debug_counter = debug_counter + 1
                writer.writerow(row)
    print(f"Renamed Gene Entries: {debug_counter}")

def process_isoform_tsv(pasa_cluster_map, tsv_file):
    debug_counter = 0
    with open(tsv_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t', quotechar='"')

        fieldnames = reader.fieldnames
        base_name = tsv_file.rsplit('.', 1)[0]
        output_file = f"{base_name}_isoform_renamed.tsv"

        with open(output_file, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for row in reader:
                original_name = row['Name']

                if original_name in pasa_cluster_map:
                    row['Name'] = pasa_cluster_map[original_name]
                    debug_counter = debug_counter + 1

                writer.writerow(row)
    print(f"Renamed Gene Entries: {debug_counter}")

def main():
    if len(sys.argv) == 2 and (sys.argv[1] == "-h" or sys.argv[1] == "help"):
        print("Usage: python tmap-rename.py <file.tmap> <file.tsv> <gene|isoform>\n")
        print("Tries to reconcile possible cufflink matches using a .tmap file (should get from pipeline) in expression data (should get from rsem)\n")
        print("if one to one match, renames the gene/isoform in the .tsv file directly from cufflink name to reference name\n")
        print("if one to many, renames the gene/isoform in the .tsv file directly from cufflink name to reference name + ?x_cufflinkname where x is priority code\n")
        print("check cuffcompare documentation for priority codes, https://cole-trapnell-lab.github.io/cufflinks/cuffcompare/")
        sys.exit(0)

    if len(sys.argv) != 4:
        print("Usage: python tmap-rename.py <file.tmap> <file.tsv> <gene|isoform>")
        sys.exit(1)

    tsv_file = None
    tmap_file = None
    level = None

    for arg in sys.argv[1:]:
        if arg.endswith('.tmap'):
            tmap_file = arg
        elif arg.endswith('.tsv'):
            tsv_file = arg
        elif arg.lower() in ['gene', 'isoform']:
            level = arg.lower()

    if not tmap_file:
        print("Error: A .tmap file must be provided.")
        sys.exit(1)
    if not tsv_file:
        print("Error: A .tsv file must be provided.")
        sys.exit(1)
    if not level:
        print("Error: The third argument must be 'gene' or 'isoform'.")
        sys.exit(1)


    pasa_cluster_map = {}
    if level == 'gene':
        process_gene_tmap(pasa_cluster_map, tmap_file)
        process_gene_tsv(pasa_cluster_map, tsv_file)
    elif level == 'isoform':
        process_isoform_tmap(pasa_cluster_map, tmap_file)
        process_isoform_tsv(pasa_cluster_map, tsv_file)

    return

if __name__ == '__main__':
    main()
