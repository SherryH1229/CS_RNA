import argparse, sys

parser = argparse.ArgumentParser(description = 'Match cufflinks output to corresponding RNA biotypes from Ensembl.')
parser.add_argument('fpkm_file', type = str, help = 'FPKM file ("isoforms.fpkm_tracking") from cufflinks')
parser.add_argument('gtffile', metavar = 'ensembl_gtf_file', help = 'Ensembl GTF Files')
#parser.add_argument('refToType', metavar = 'refseqID_to_type_file', help = 'A file with two columns: refseq ID and RNA biotype, can be generated via Biomart.')
args = parser.parse_args()

refIDToType = dict()
#geneIDToType = dict()

with open(args.gtffile, 'rU') as fIn:
    for line in fIn:
        if not line.startswith('#'):
            tokens = line.strip().split('\t')
            if len(tokens) < 9:
                print >> sys.stderr, line
            annotation = tokens[8]
            annoTokens = annotation.strip().split('; ')

            gene_id = '.'
            gene_biotype = '.'
            for annoToken in annoTokens:
                keyValuePair = annoToken.strip().split(' ')
                if len(keyValuePair) < 2:
                    keyValuePair = annoToken.strip().split('=')

                if keyValuePair[0] == 'gene_id':
                    gene_id = keyValuePair[1].strip('";')
                if keyValuePair[0] == 'gene_biotype':
                    gene_biotype = keyValuePair[1].strip('";')

            if gene_id not in refIDToType:
                refIDToType[gene_id] = gene_biotype

with open(args.fpkm_file, 'rU') as fIn:
    opening = True
    for line in fIn:
        if opening:
            print >> sys.stdout, line.strip() + "\tBiotype"
            opening = False
            continue
        tokens = line.strip().split()
        if tokens[0] in refIDToType:
            print >> sys.stdout, line.strip() + "\t" + refIDToType[tokens[0]]
        else:
            print >> sys.stdout, line.strip() + "\t."

