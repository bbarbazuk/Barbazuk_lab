import sys
import csv

#Filters for log2FoldChange values outside of the provided range
#eg. python filter_change.py file.tsv -1 1
#will filter out all values greater than -1 and less than 1 and write to file_log2filtered.tsv
def process_tsv(tsv_file, min, max):
    debug_count = 0
    with open(tsv_file, 'r') as infile:
        reader = csv.DictReader(infile, delimiter='\t')

        fieldnames = reader.fieldnames
        base_name = tsv_file.rsplit('.', 1)[0]
        output_file = f"{base_name}_log2filtered.tsv"

        with open(output_file, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

            for row in reader:
                change = float(row['log2FoldChange'])
                
                if change < min or change > max:
                    writer.writerow(row)
                    debug_count += 1
    
    print(f"Done filtering. New tsv contains {debug_count} entries.")
            
                    



def main():
    if len(sys.argv) == 2 and (sys.argv[1] == "-h" or sys.argv[1] == "help"):
        print("Usage: python filter_change.py <file.tsv> <min_to_filter> <max_to_filter>\n")
        print("Filters for log2FoldChange values outside of the provided range")
        print("eg. python filter_change.py file.tsv -1 1")
        print("will filter out all values greater than -1 and less than 1 and write to file_log2filtered.tsv")
        sys.exit(0)
    if len(sys.argv) != 4:
        print("Usage: python filter_change.py <file.tsv> <min_to_filter> <max_to_filter>")
        sys.exit(1)
    
    tsv_file = None
    min_change = None
    max_change = None
    
    if sys.argv[1].endswith('.tsv'):
        tsv_file = sys.argv[1]
    if isinstance(float(sys.argv[2]), float):
        min_change = float(sys.argv[2])
    if isinstance(float(sys.argv[3]), float):
        max_change = float(sys.argv[3])

    if not tsv_file:
        print("Error: A .tsv file must be provided.")
        sys.exit(1)

    process_tsv(tsv_file, min_change, max_change)

if __name__ == '__main__':
    main()