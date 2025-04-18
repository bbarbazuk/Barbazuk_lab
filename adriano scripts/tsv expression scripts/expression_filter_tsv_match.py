import sys
import csv

#Filters for 1-1 matches based on annotation
#Format from earlier script that tried to reconcile possible matches
def process_tsv(tsv_file):
    debug_count = 0
    with open(tsv_file, 'r') as infile:
        reader = csv.DictReader(infile, delimiter='\t')

        fieldnames = reader.fieldnames
        base_name = tsv_file.rsplit('.', 1)[0]
        output_file = f"{base_name}_matchfiltered.tsv"

        with open(output_file, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

            for row in reader:
                name = row['Name']
                
                if "align_id" not in name:
                    writer.writerow(row)
                    debug_count += 1
    
    print(f"Done filtering. New tsv contains {debug_count} entries.")
            
                    



def main():
    if len(sys.argv) != 2:
        print("Usage: python filter_match.py <file.tsv>")
        sys.exit(1)
    
    if sys.argv[1] == "-h" or sys.argv[1] == "help":
        print("Usage: python filter_match.py <file.tsv>")
        print("Filters for 1-1 matches based on annotation")
        print("Format from earlier script that tried to reconcile possible matches")
        sys.exit(0)
    
    tsv_file = None
    
    if sys.argv[1].endswith('.tsv'):
        tsv_file = sys.argv[1]

    if not tsv_file:
        print("Error: A .tsv file must be provided.")
        sys.exit(1)

    process_tsv(tsv_file)

if __name__ == '__main__':
    main()