import sys

#TSV format ex days = 2 samples = 2:
#Name, sample 1-1, sample 1-2, sample 2-1, sample 2-2
#Adds a third column with the average of the samples for each day and makes new file named "sampletsv_averaged.tsv"
#eg Name, sample 1-1, sample 1-2, average 1, sample 2-1, sample 2-2, average 2
def process_tsv(tsv_file, days, samples):
    with open(tsv_file, 'r') as infile:
        debug_counter = 0
        lines = infile.readlines()
        fieldnames = lines[0].strip().split('\t')
        base_name = tsv_file.rsplit('.', 1)[0]
        output_file = f"{base_name}_averaged.tsv"

        with open(output_file, 'w', newline='') as outfile:
            outfile.write(f"{fieldnames[0]}\t")
            for day in range(int(days)):
                for sample in range(samples):
                    index = 1 + day * samples + sample
                    outfile.write(f"{fieldnames[index]}\t")
                outfile.write(f"average {day + 1}\t")
            outfile.write("\n")

            for line in lines[1:]:
                row = line.strip().split('\t')
                outfile.write(f"{row[0]}\t")
                for day in range(int(days)):
                    temp = []
                    for sample in range(samples):
                        index = 1 + day * samples + sample
                        temp.append(float(row[index]))
                        outfile.write(f"{row[index]}\t")
                    average = sum(temp) / samples
                    outfile.write(f"{average}\t")
                outfile.write(f"\n")
                debug_counter += 1
    print(f"Done averaging. New tsv contains {debug_counter} entries.")

def main():

    if len(sys.argv) == 2 and (sys.argv[1] == "-h" or sys.argv[1] == "help"):
        print("Usage: python average_counts.py <file.tsv> <days> <samples/day>\n")
        print("TSV format ex days = 2 samples = 2: \n")
        print("Name, sample 1-1, sample 1-2, sample 2-1, sample 2-2\n")
        print("Adds a third column with the average of the samples for each day and makes new file named 'sampletsv_averaged.tsv'\n")
        print("eg Name, sample 1-1, sample 1-2, average 1, sample 2-1, sample 2-2, average 2")
        sys.exit(0)
    if len(sys.argv) != 4:
        print("Usage: python average_counts.py <file.tsv> <days> <samples/day>")
        sys.exit(1)
    
    tsv_file = None
    days = None
    samples = None
    
    if sys.argv[1].endswith('.tsv'):
        tsv_file = sys.argv[1]
    
    if sys.argv[2].isdigit():
        days = int(sys.argv[2])
    
    if sys.argv[3].isdigit():
        samples = int(sys.argv[3])

    if not tsv_file:
        print("Error: A .tsv file must be provided.")
        sys.exit(1)
    
    if not days or not samples:
        print("Error: Number of days and samples per day must be provided.")
        sys.exit(1)

    process_tsv(tsv_file, days, samples)

if __name__ == '__main__':
    main()