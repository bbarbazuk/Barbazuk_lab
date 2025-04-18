import sys
import multiprocessing as mp

class exon:
    def __init__(self, gene, chromosome, start, end, strand, format):
        self.gene = gene
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand
        self.format = format

def parse_pasa_output(exons, pasa_output):
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
                exons.append(exon(gene_id, chromosome, start, end, strand, line.strip()))


def process_exon(exon_data):
    evexon, gff_file = exon_data
    gene_id = evexon.gene
    transcripts = []
    transcripts_alt = []
    
    with open(gff_file, 'r') as infile:
        gene_flag = False
        exon_skip_flag = True
        current_transcript = None
        for line in infile:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            metadata_fields = fields[8].split(';')
            
            if fields[2] == 'gene':
                if current_transcript:
                    break
                for field in metadata_fields:
                    if field.startswith(' gene_name'):
                        gene_name = field.split('"')[1]
                        break
                gene_flag = (gene_name == gene_id)
            
            if not gene_flag:
                continue
            
            if fields[2] == 'transcript':
                if current_transcript:
                    if exon_skip_flag:
                        transcripts.append(current_transcript)
                    else:
                        transcripts_alt.append(current_transcript)
                current_transcript = None
                exon_skip_flag = True
                
                for field in metadata_fields:
                    if field.startswith(' transcript_id'):
                        current_transcript = field.split('"')[1]
                        break
            
            if fields[2] == 'exon':
                    if int(fields[3]) == evexon.start and int(fields[4]) == evexon.end and fields[6] == evexon.strand:
                        exon_skip_flag = False
        if current_transcript:
            if exon_skip_flag:
                transcripts.append(current_transcript)
            else:
                transcripts_alt.append(current_transcript)
    return (evexon, transcripts, transcripts_alt)

def parse_gff(exons, gff_file, species_id):
    pool = mp.Pool(mp.cpu_count())
    results = pool.map(process_exon, [(evexon, gff_file) for evexon in exons])
    pool.close()
    pool.join()

    results.sort(key=lambda x: x[0].format)

    with open(species_id + '_event_transcript_outgroups.txt', 'w') as outfile2:
        with open(species_id + '_event_transcript_ingroups.txt', 'w') as outfile:
            for exon, transcripts, transcripts_alt in results:
                ingroup_line = f"{exon.format};" + ';'.join(transcripts) + ';\n' if transcripts else f"{exon.format};\n"
                outgroup_line = f"{exon.format};" + ';'.join(transcripts_alt) + ';\n' if transcripts_alt else f"{exon.format};\n"
                
                outfile.write(ingroup_line)
                outfile2.write(outgroup_line)

def main():
    if len(sys.argv) == 2 and (sys.argv[1] == "-h" or sys.argv[1] == "help"):
        print("Usage: sort_events_transcripts.py <gtf_file> <pasa_output> <species_id>")
        print("Groups transcripts together by presence of exon skip event taking in a gff file, a pasa output file (exon skip), and a species id")
        print("Example output: event1;transcript1;transcript3;transcript4")
        print("Will output file species_event_transcript_ingroups.txt with the above format")
        print("Will also output species_event_transcript_outgroups.txt with the same format, but for the outgroup (transcripts that do not contain the exon skip event)")
        print("Tested primarily with Mus musculus but should work with other species as long as the gff file is formatted correctly")
        sys.exit(0)
    if len(sys.argv) != 4:
        print("Error: Format is sort_events_transcripts.py gtf_file.gtf event_clusters.txt Mus")
        sys.exit(1)
    
    pasa_file = sys.argv[2]
    gtf_file = sys.argv[1]
    species_id = sys.argv[3]
    exons = []
    
    parse_pasa_output(exons, pasa_file)
    parse_gff(exons, gtf_file, species_id)

if __name__ == '__main__':
    main()