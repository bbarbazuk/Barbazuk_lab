import sys

#Takes in an expression change file in the format output by event_expression_proportions.py and outputs a file for upregulated events and another for downregulated events.
#for each set of days in the data (ie, for Days 0, 7, 14 will output six files for 0-7, 7-14, and 0-14).
#An event is considered differently regulated if the difference between the groups is greater than 30% of the first value. 
#For example, if the virst value is 0.5, the second must be greater than 0.5 + 0.15 = 0.65 to be considered upregulated and less than 0.35 to be considered downregulated.
#The output files will simply contain a list of events that are upregulated or downregulated.
#usage python find_proportion_regulation.py <expression_proportions.txt>


def find_regulation(expression_proportions, days):
    with open(expression_proportions, 'r') as file:
        lines = file.readlines()
        base_name = expression_proportions.rsplit('.', 1)[0]
        indices = [(i, j) for i in range(days) for j in range(i + 1, days)]
        for index in indices:
            start_day = index[0]
            end_day = index[1]
            up_count = 0
            down_count = 0
            upregulated_file = f"{base_name}_{start_day}-{end_day}_upregulated.txt"
            downregulated_file = f"{base_name}_{start_day}-{end_day}_downregulated.txt"
            with open(upregulated_file, 'w') as upfile, open(downregulated_file, 'w') as downfile:
                for line in lines[1:]:
                    if line.startswith('\t') and not line.startswith('\t\t'):
                        parts = line.strip().split()
                        event = parts[0]
                        values = [float(v.strip(')')) for v in parts[-days:]]
                        if all(v < 0.1 for v in values):
                            continue
                        start_value = values[index[0]]
                        end_value = values[index[1]]
                        threshold = 0.3 * start_value
                        if end_value > start_value + threshold:
                            upfile.write(event + '\n')
                            up_count += 1
                        elif end_value < start_value - threshold:
                            downfile.write(event + '\n')
                            down_count += 1
            print(f"Found {up_count} upregulated events and {down_count} downregulated events for indices {start_day} to {end_day}.")

def main():
    if len(sys.argv) == 2 and (sys.argv[1] == "-h" or sys.argv[1] == "--help"):
        print("Takes in event/transcript match list and count data (run through script that adds average columns) to output a file containing the proportion of each transcript that is associated with each event (using the average across samples).")
        print("for each set of days in the data (ie, for Days 0, 7, 14 will output six files for 0-7, 7-14, and 0-14).")
        print("An event is considered differently regulated if the difference between the groups is greater than 30 percent of the first value.")
        print("For example, if the virst value is 0.5, the second must be greater than 0.5 + 0.15 = 0.65 to be considered upregulated and less than 0.35 to be considered downregulated.")
        print("The output files will simply contain a list of events that are upregulated or downregulated.")
        print("Usage: python find_proportion_regulation.py <expression_proportions.txt> <days>")
        sys.exit(1)
    
    if len(sys.argv) != 3:
        print("Usage: python find_proportion_regulation.py <expression_proportions.txt> <days>")
        sys.exit(1)
    
    expression_proportions = sys.argv[1]
    days = int(sys.argv[2])
    find_regulation(expression_proportions, days)

if __name__ == '__main__':
    main()