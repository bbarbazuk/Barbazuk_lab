import sys
class exon:
    def __init__(self, chromosome, start, end, strand, format):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand
        self.format = format


def parse_pasa_output(exon_map, pasa_output):
    with open(pasa_output, 'r') as file:
        for line in file:
            parts = line.strip().split(':')
            if len(parts) > 1:
                gene_id = parts[0][:-1]
                exon_info = parts[1].split('_')
                chromosome = exon_info[1]
                num_1 = exon_info[-4]
                num_2 = exon_info[-3]
                if exon_info[-1] == '-':
                    start = int(num_2)
                    end = int(num_1)
                else:
                    start = int(num_1)
                    end = int(num_2)
                strand = exon_info[-1]
                if gene_id not in exon_map:
                    exon_map[gene_id] = []
                exon_map[gene_id].append(exon(chromosome, start, end, strand, line.strip()))


def parse_gff(exon_map, gff_file, species_id):
    with open(species_id + '_event_transcript_outgroups.txt', 'w') as outfile2:
        with open(species_id + '_event_transcript_ingroups.txt', 'w') as outfile:
            for gene_id, exons in exon_map.items():
                for exon in exons:
                    outfile.write(f"{exon.format};")
                    outfile2.write(f"{exon.format};")
                    exon_skip_flag = True
                    gene_flag = False
                    current_transcript = None
                    transcripts = []
                    transcripts_alt = []
                    with open(gff_file, 'r') as infile:
                        for line in infile:
                            if line.startswith('#'):
                                continue
                            fields = line.strip().split('\t')
                            metadata_fields = fields[8].split(';')
                            gene_name = ''
                            for field in metadata_fields:
                                if field.startswith(' gene_id'):
                                    gene_name = field.split('"')[1]
                                    break
                            if gene_name != gene_id and gene_flag:
                                gene_flag = False
                                break
                            elif gene_name != gene_id:
                                continue
                            elif gene_name == gene_id:
                                gene_flag = True
                            if fields[2] == 'transcript':
                                if current_transcript:
                                    if exon_skip_flag:
                                        transcripts.append(current_transcript)
                                    else:
                                        transcripts_alt.append(current_transcript)
                                    current_transcript = None
                                    exon_skip_flag = True
                                for field in metadata_fields:
                                    if 'transcript_id' in field:
                                        current_transcript = field.split('"')[1]
                                        break
                            if fields[2] == 'exon':
                                if int(fields[3]) == int(exon.start) and int(fields[4]) == int(exon.end) and fields[6] == exon.strand:
                                    exon_skip_flag = False
                    if current_transcript and exon_skip_flag:
                        transcripts.append(current_transcript)
                    elif current_transcript and not exon_skip_flag:
                        transcripts_alt.append(current_transcript)
                    if transcripts and transcripts_alt:
                        outfile.write(';'.join(transcripts) + ';\n')
                        outfile2.write(';'.join(transcripts_alt) + ';\n')
                    elif transcripts and not transcripts_alt:
                        outfile.write(';'.join(transcripts) + ';\n')
                        outfile2.write('\n')
                    elif not transcripts and transcripts_alt:
                        outfile.write('\n')
                        outfile2.write(';'.join(transcripts_alt) + ';\n')
                    else:
                        outfile.write('\n')
                        outfile2.write('\n')


#USAGE sort_events_transcripts.py <gtf_file> <pasa_output> <species_id>
def main():
    if len(sys.argv) == 2 and (sys.argv[1] == "-h" or sys.argv[1] == "help"):
        print("Usage: sort_events_transcripts.py <gtf_file> <pasa_output> <species_id>")
        print("Groups transcripts together by presence of exon skip event taking in a gff file, a pasa output file (exon skip), and a species id")
        print("Example output: event1;transcript1;transcript3;transcript4")
        print("Will output file species_event_transcript_ingroups.txt with the above format")
        print("Will also output species_event_transcript_outgroups.txt with the same format, but for the outgroup (transcripts that do not contain the exon skip event)")
        print("This version is tested with Acomys cahirinus and was modified to account for the unpublished gtf file")
        sys.exit(0)
    if len(sys.argv) != 4:
        print("Error: Format is sort_events_transcripts.py gtf_file.gtf event_clusters.txt Mus")
        sys.exit(1)
    
    pasa_file = sys.argv[2]
    gtf_file = sys.argv[1]
    species_id = sys.argv[3]
    exon_map = {}
    
    parse_pasa_output(exon_map, pasa_file)
    parse_gff(exon_map, gtf_file, species_id)

if __name__ == '__main__':
    main()
