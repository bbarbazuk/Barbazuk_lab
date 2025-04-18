import sys
from decimal import Decimal, getcontext

# Set the precision for Decimal operations
getcontext().prec = 10

class Gene:
    gene_name = ""
    transcripts = []
    def __init__(self, gene_name, transcripts):
        self.gene_name = gene_name
        self.transcripts = transcripts
    def get_transcripts(self):
        return self.transcripts

def process_data(event_transcript_list, event_transcript_outgroup, expression_data, days, samples_per_day, species_id):
    with open(event_transcript_list, 'r') as event_file:
        with open(expression_data, 'r') as expression_file:
            with open(event_transcript_outgroup, 'r') as event_outgroup_file:
                next(expression_file)  # Skip the header line
                data = {}
                data_complement = {}
                genes = {} 
                current_gene = None
                for line1, line2 in zip(event_file, event_outgroup_file):
                    fields_in = line1.strip().split(';')
                    fields_out = line2.strip().split(';')
                    event = fields_in[0]
                    gene = event.split(':')[0].rstrip('_')
                    transcripts_in = [transcript for transcript in fields_in[1:] if transcript]
                    transcripts_out = [transcript for transcript in fields_out[1:] if transcript]
                    if gene != current_gene:
                        current_gene = gene
                        data[current_gene] = {}
                        data_complement[current_gene] = {}
                        genes[gene] = Gene(gene, (transcripts_in + transcripts_out))
                    data[current_gene][event] = {transcript: [Decimal(0)] * days for transcript in transcripts_in if transcript}
                    data_complement[current_gene][event] = {transcript: [Decimal(0)] * days for transcript in transcripts_out if transcript}
                transcript_day_map = {}
                event_totals_map = {}
                event_totals_map_complement = {}
                for line in expression_file:
                    fields = line.strip().split('\t')
                    transcript = fields[0].strip('"')
                    gene = transcript.split('_')[0]
                    transcript_day_map[transcript] = [Decimal(fields[1 + samples_per_day + day * (samples_per_day + 1)]) for day in range(days)]
                for gene in data:
                    total_counts = [Decimal(0)] * days
                    for transcript in genes[gene].get_transcripts():
                        if transcript in transcript_day_map:
                            counts = transcript_day_map[transcript]
                            total_counts = [total + count for total, count in zip(total_counts, counts)]
                    for event in data[gene]:
                        event_totals_map[event] = [Decimal(0)] * days
                        event_totals_map_complement[event] = [Decimal(0)] * days
                        for transcript in data[gene][event]:
                            if transcript in transcript_day_map:
                                counts = transcript_day_map[transcript]
                                proportions = [count / total if total > 0 else Decimal(0) for count, total in zip(counts, total_counts)]

                                proportions = [proportion.quantize(Decimal('0.001')) for proportion in proportions]
                                data[gene][event][transcript] = proportions
                                event_totals_map[event] = [total + proportion for total, proportion in zip(event_totals_map[event], proportions)]
                    for event in data_complement[gene]:
                        event_totals_map_complement[event] = [Decimal(0)] * days
                        for transcript in data_complement[gene][event]:
                            if transcript in transcript_day_map:
                                counts = transcript_day_map[transcript]
                                proportions = [count / total if total > 0 else Decimal(0) for count, total in zip(counts, total_counts)]
                                proportions = [proportion.quantize(Decimal('0.001')) for proportion in proportions]
                                data_complement[gene][event][transcript] = proportions
                                event_totals_map_complement[event] = [total + proportion for total, proportion in zip(event_totals_map_complement[event], proportions)]
            with open(f"{species_id}_expression_proportions.txt", 'w') as outfile:
                for gene in data:
                    outfile.write(f"{gene}\n")
                    for event in data[gene]:
                        outfile.write(f"\t{event} (event_totals: {' '.join(map(str, event_totals_map[event]))})\n")
                        for transcript in data[gene][event]:
                            outfile.write(f"\t\t{transcript}: {' '.join(map(str, data[gene][event][transcript]))}\n")
            with open(f"{species_id}_expression_proportions_complement.txt", 'w') as outfile:
                for gene in data_complement:
                    outfile.write(f"{gene}\n")
                    for event in data_complement[gene]:
                        outfile.write(f"\t{event} (event_complement_totals: {' '.join(map(str, event_totals_map_complement[event]))})\n")
                        for transcript in data_complement[gene][event]:
                            outfile.write(f"\t\t{transcript}: {' '.join(map(str, data_complement[gene][event][transcript]))}\n")
    #print(f"Done processing data. Output in {species_id}_expression_proportions.txt")
    #print(f"Compare to proportions in {species_id}_expression_proportions_complement.txt and assure that the sum of the proportions for each event is 1.0.")

def main():
    if len(sys.argv) == 2 and (sys.argv[1] == "-h" or sys.argv[1] == "help"):
        print("Usage: python event_expression_proportions.py <event_transcript_ingroup> <event_transcript_outgroup> <expression_data> <days> <samples/day> <species_id>")
        print("Takes in event/transcript match list and count data (run through script that adds average columns) to output a file containing the proportion of each transcript that is associated with each event (using the average across samples).")
        print("Additionally separated by gene assuming standard pipeline output format.")
        print()
        print("Example grouping input:")
        print("Flna_event1;transcript1;transcript2;transcript3")
        print("Flna_event2;transcript1;transcript2")
        print()
        print("Count input should be in the format of the output from expression_average_counts.py (a tsv file with the average of each sample for each day) in format:")
        print("Transcript, sample 1-1, sample 1-2, average 1, sample 2-1, sample 2-2, average 2")
        print()
        print("Example output (species_id_expression_proportions.txt) assuming four days:")
        print("gene 1")
        print("      event 1:")
        print("            transcript 1: 0.5 0.7 1.0 0.3")
        print("            transcript 2: 0.5 0.3 0.0 0.7")
        print("      event 2:")
        print("            transcript 1: 0.25 0.25 0.5 0.75")
        print("            transcript 2: 0.75 0.75 0.5 0.25")
        print("gene 2")
        print("      event 1:")
        print("            transcript 1: 0.5 0.5 0.5 0.5")
        print("            transcript 2: 0.5 0.5 0.5 0.5")
        print("")
        print("Compare this to the complement output (species_id_expression_proportions_complement.txt) to assure that the sum of the proportions for each event is 1.0.")
        sys.exit(0)

    if len(sys.argv) != 7:
        print("Usage: python event_expression_proportions.py <event_transcript_ingroup> <event_transcript_outgroup> <expression_data> <days> <samples/day> <species_id>")
        sys.exit(1)
    
    event_transcript_in = sys.argv[1]
    event_transcript_out = sys.argv[2]
    expression_data = sys.argv[3]
    days = int(sys.argv[4])
    samples_per_day = int(sys.argv[5])
    species_id = sys.argv[6]
    
    if not event_transcript_in.endswith('.txt') or not event_transcript_out.endswith('.txt') or not expression_data.endswith('.tsv'):
        print("Usage: python event_expression_proportions.py <event_transcript_ingroup> <event_transcript_outgroup> <expression_data> <days> <samples/day>")
        sys.exit(1)
    
    process_data(event_transcript_in, event_transcript_out, expression_data, days, samples_per_day, species_id)

if __name__ == "__main__":
    main()
