import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from chipseq_gene_annoatation_function_V5_5 import type1_downstream_geneannotation, type1_upstream_geneannotation
from chipseq_gene_annoatation_function_V5_5 import type2_geneannotation, type3_geneannotation, remove_type3_in_type2
from chipseq_gene_annoatation_function_V5_5 import two_type1_selection, summit_signal, add_description, extract_common_name
from chipseq_gene_annoatation_function_V5_5 import get_bdg_200bp_per_bins_promoter, get_bdg_200bp_per_bins_tail, get_bdg_gene_body ## working in V5.3 removing in V5.4 and V5.5
import argparse

start_time = time.time()
print('\n\n')
print('start time:', time.ctime())


##################################### arguments #####################################################

parser = argparse.ArgumentParser(description="This is a in-home script developed by Raiman (ruiwenchen@um.edu.mo) for annotating target genes based on a peak summit.\n\
    python package: pandas is required.\n\
    Before using this script, please ensure chipseq_gene_annoatation_function-VX.X.py is appear in the same folder.\n\
",\
formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument("-summit_file", help="sample_summit.bed output from macs2 callpeak or other tab delimited file which the first three columns are chr, start, end")

parser.add_argument("-gff", help = "gff_file")

parser.add_argument("-bdg_file", help = "The bdg file can be generated from a mapping file using MACS2.")

parser.add_argument("-chrom_size", help = "The chrom_size file which for generate bed file to get fasta, this file to ensure peak summit + 250 bp would not exceed the chromosome end.")

parser.add_argument("-o", help = "output prefix, default=./output", default='./output')

parser.add_argument("-type1_method", help = "use ATG or TSS to annotate type1 gene, defalut=ATG",\
 choices=['ATG', 'TSS'], default='ATG')

parser.add_argument("--d1", "-distance1", help = "max distance in type1 gene annotation, defalut=2000",  default=2000)

parser.add_argument("--d_r", "-distance_ratio", help = "When there are two genes within d1 bp at both ends of the peak summit,\
                     if d_far / d_close >= distance_ratio, only the closest one will be selected as the target gene,\
                     otherwise, both will be regarded as the target gene, defalut=1.5",  default=1.5)

parser.add_argument("--d2", "-distance2", help = "When there are two genes within d1 bp at both ends of the peak summit,\
                     and the distance between the TSS of the two genes is less than d2, both will be regarded as the target gene (ignoring the distance ratio),\
                     defalut=500",  default=500)

parser.add_argument("--h_t", "-high_quality_threshold", help = "summit signal greater than which will be considered as high quality summmit,\
                    defalut=20",  default=20)


parser.add_argument("-ncRNA_gene", help="Consider ncRNA_gene or not. Set to '1' to include ncRNA_gene,\
 or '0' to exclude ncRNA_gene. Default is '0' (unconsider).", choices=['0', '1'], default='0')

parser.add_argument("-mito_gene", help="Consider mitochondria or not. Set to '1' to include mitochondria,\
 or '0' to exclude ncRNA_gene. Default is '0' (unconsider). ONLY WORK FOR ASPERGILLUS NIDULANS.", choices=['0', '1'], default='0')

args = parser.parse_args()

chr_index = 0
summit_index = 2
d1 = int(args.d1)
d2 = int(args.d2)
d_r = float(args.d_r)
h_t = int(args.h_t)
type1_method = args.type1_method
ncRNA_gene = int(args.ncRNA_gene)
mito_gene = int(args.mito_gene)






#################################### read input files ###############################################
bdg_cols = ['chr', 'start', 'end', 'score']
bdg = pd.read_table(args.bdg_file, names=bdg_cols)

genome_gff_header = ['chr','source','Genome_attr','start','end','score','strand','phase','ID']
genome_gff = pd.read_csv(args.gff, sep='\t', comment='#', names=genome_gff_header)
CHIP_seq=pd.read_csv(args.summit_file, sep='\t', header = None)

chrom_size = pd.read_csv(args.chrom_size, sep="\t", names=['chr', 'length'])


CHIP_seq_output = CHIP_seq[[chr_index, summit_index]].copy()
if(mito_gene==0):
    CHIP_seq_output = CHIP_seq_output[CHIP_seq_output[chr_index]!='mito_A_nidulans_FGSC_A4']

genome_gff['Genome_attr'] = genome_gff['Genome_attr'].replace('protein_coding_gene','gene')
genome_gff['Genome_attr'] = genome_gff['Genome_attr'].replace('pseudogene','gene')
if(ncRNA_gene == 1):
    genome_gff['Genome_attr'] = genome_gff['Genome_attr'].replace('ncRNA_gene','gene')

################## sub_genome_gff for annotating type1, type2 and type3 genes  #####################
#type1
if(type1_method == 'ATG'):
    genome_gff_type1 = genome_gff[genome_gff['Genome_attr']=='CDS'].copy()
elif(type1_method == 'TSS'):
    genome_gff_type1 = genome_gff[genome_gff['Genome_attr']=='gene'].copy()
#type2
genome_gff_fiveUTR = genome_gff[genome_gff['Genome_attr']=='five_prime_UTR'].copy()
genome_gff_threeUTR = genome_gff[genome_gff['Genome_attr']=='three_prime_UTR'].copy()
#type3
genome_gff_gene = genome_gff[genome_gff['Genome_attr']=='gene'].copy()


###################################   main function   ################################################
log = open('%s_working_log.txt'%args.o, 'w+')
log.write("1. log for annotating target genes\n")
print(time.ctime())
print("calculating type3 target genes")
CHIP_seq_output['gene_body'] = CHIP_seq_output.apply(lambda X:
                                                      type3_geneannotation(genome_gff_gene, X[0], X[summit_index], log),
                                                     axis=1)

print(time.ctime())
print("calculating type2 target genes - 5'-UTR")
CHIP_seq_output['five_prime_UTR'] = CHIP_seq_output.apply(lambda X:
                                                      type2_geneannotation(genome_gff_fiveUTR, X[0], X[summit_index]),
                                                     axis=1)

print(time.ctime())
print("calculating type2 target genes - 3'-UTR")
CHIP_seq_output['three_prime_UTR'] = CHIP_seq_output.apply(lambda X:
                                                      type2_geneannotation(genome_gff_threeUTR, X[0], X[summit_index]),
                                                     axis=1)

if(type1_method == 'ATG'):
    print(time.ctime())
    print("type1_method: ATG")
    print("running function remove_type3_in_type2")
    CHIP_seq_output['gene_body'] = CHIP_seq_output.apply(lambda X:
                                                        remove_type3_in_type2(X['three_prime_UTR'],X['five_prime_UTR'], X['gene_body']),
                                                        axis=1)

print(time.ctime())
print("Annotating type1 downstream target genes")
result = CHIP_seq_output.apply(lambda X:
                                                      type1_downstream_geneannotation(genome_gff_type1,
                                                       X[chr_index], X[summit_index], d1, X['gene_body'], log), axis=1)
result = np.vstack([np.asarray(t) for t in result])
CHIP_seq_output[['type1_down', 'down_distance']] = result 

print(time.ctime())
print("Annotating type1 upstream target genes")
result = CHIP_seq_output.apply(lambda X:
                                                      type1_upstream_geneannotation(genome_gff_type1, X[chr_index], 
                                                      X[summit_index], d1, X['gene_body'], log), axis=1)
result = np.vstack([np.asarray(t) for t in result])
CHIP_seq_output[['type1_up', 'up_distance']] = result

#### change the distance to float type so that we can compare them.
CHIP_seq_output['down_distance'] = CHIP_seq_output['down_distance'].replace('None', float('nan'))
CHIP_seq_output['down_distance'] = CHIP_seq_output['down_distance'].astype(float)
CHIP_seq_output['up_distance'] = CHIP_seq_output['up_distance'].replace('None', float('nan'))
CHIP_seq_output['up_distance'] = CHIP_seq_output['up_distance'].astype(float)

#### If we want to analysis something, such as draw a scatter plot of up_distance and down_distance.
#### we can generate a tem file here.


print(time.ctime())
print("Two type1 selection")
result = CHIP_seq_output.apply(lambda X:
                                                      two_type1_selection(d_r, X['type1_down'], 
                                                      X['down_distance'], X['type1_up'], X['up_distance'], d2), axis=1)
result = np.vstack([np.asarray(t) for t in result])
CHIP_seq_output[['type1_down', 'down_distance', 'type1_up', 'up_distance']] = result



### create ID_description data frame
ID_description = pd.DataFrame()
ID = genome_gff_gene['ID'].tolist()
ID_list=[]
description_list=[]
for i in range(len(ID)):
    ANID = ID[i].split('ID=')[1].split(';')[0]
    description = ID[i]
    ID_list.append(ANID)
    description_list.append(description)
    
ID_description['ID'] = ID_list
ID_description['description'] = description_list

### add description
print("adding type3 description")
CHIP_seq_output['gene_body_description'] = CHIP_seq_output.apply(lambda X:
                                                      add_description(X['gene_body'], ID_description),
                                                     axis=1)
print("adding type1_down description")
CHIP_seq_output['type1_down_description'] = CHIP_seq_output.apply(lambda X:
                                                      add_description(X['type1_down'], ID_description),
                                                     axis=1)
print("adding type1_up description")
CHIP_seq_output['type1_up_description'] = CHIP_seq_output.apply(lambda X:
                                                      add_description(X['type1_up'], ID_description),
                                                     axis=1)

### rename chr start end summit
CHIP_seq_output = CHIP_seq_output.rename(columns={
    chr_index: 'chr',
    summit_index: 'summit'
})

### rearrange columns
### type2 not need, remove them.
CHIP_seq_output = CHIP_seq_output[['chr', 'summit', 'gene_body', 'gene_body_description',\
                                    'type1_down', 'down_distance', 'type1_down_description', 'type1_up', 'up_distance', 'type1_up_description']]


#### sorted
print(time.ctime())
print("adding signal and sorting")
Chrs = []
summits = []
binding_types = []
signals = []
geneIDs = []
signals2 = []
for i in range(len(CHIP_seq_output)):
    Chr = CHIP_seq_output.iloc[i,0]
    summit = CHIP_seq_output.iloc[i,1]
    signal = summit_signal(bdg, Chr, summit)
    signals2.append(signal)
    if(CHIP_seq_output.iloc[i,2]!=None):
        genes = CHIP_seq_output.iloc[i,2].split(", ")
        for gene in genes:
            Chrs.append(Chr)
            summits.append(summit)
            binding_types.append("gene body")
            signals.append(signal)
            geneIDs.append(gene)

    if(CHIP_seq_output.iloc[i,4]!=None):
        genes = CHIP_seq_output.iloc[i,4].split(", ")
        for gene in genes:
            Chrs.append(Chr)
            summits.append(summit)
            binding_types.append("promoter")
            signals.append(signal)
            geneIDs.append(gene)

    if(CHIP_seq_output.iloc[i,7]!=None):
        genes = CHIP_seq_output.iloc[i,7].split(", ")
        for gene in genes:
            Chrs.append(Chr)
            summits.append(summit)
            binding_types.append("promoter")
            signals.append(signal)
            geneIDs.append(gene)


CHIP_seq_output['signal'] = signals2
CHIP_seq_output = CHIP_seq_output[['chr', 'summit', 'signal', 'gene_body', 'gene_body_description',\
                                    'type1_down', 'down_distance', 'type1_down_description', 'type1_up', 'up_distance', 'type1_up_description']]

CHIP_seq_output_sorted = pd.DataFrame()
CHIP_seq_output_sorted['chr'] = Chrs
CHIP_seq_output_sorted['summit'] = summits
CHIP_seq_output_sorted['binding_type'] = binding_types
CHIP_seq_output_sorted['signal'] = signals
CHIP_seq_output_sorted['geneID'] = geneIDs

CHIP_seq_output_sorted = pd.merge(CHIP_seq_output_sorted, ID_description, left_on='geneID', right_on='ID', how='left')
CHIP_seq_output_sorted = CHIP_seq_output_sorted.sort_values('signal', ascending=False)
CHIP_seq_output_sorted = CHIP_seq_output_sorted.drop('ID', axis=1)
CHIP_seq_output_sorted['common_name'] = CHIP_seq_output_sorted['description'].apply(lambda x:extract_common_name(x))
CHIP_seq_output_sorted['description'] = CHIP_seq_output_sorted['description'].apply(lambda x:x.split("description=")[1].split(";")[0])
CHIP_seq_output_sorted = CHIP_seq_output_sorted[['chr', 'summit', 'binding_type', 'signal', 'geneID', 'common_name', 'description']]

### high quality peak target for GO
print(time.ctime())
print("selecting high quality peak")
hq_peaks = CHIP_seq_output_sorted[CHIP_seq_output_sorted['signal']>h_t]
hq_genes = hq_peaks['geneID'].tolist()
hq_genes = list(set(hq_genes))

### high quality summit flanking 500 for meme
def test_length(Chr, start, end, chrom_size):
    length = chrom_size[chrom_size['chr']==Chr]['length'].tolist()[0]
    if(start<1):
        return 1, 500
    elif(end > length):
        return length-500, length
    else:
        return start, end

hq_peak_meme = pd.DataFrame()
hq_peak_meme['chr'] = hq_peaks['chr']
hq_peak_meme['start'] = hq_peaks['summit'] - 249
hq_peak_meme['end'] = hq_peaks['summit'] + 250
hq_peak_meme = hq_peak_meme.drop_duplicates(subset=hq_peak_meme.columns.tolist(), keep='first')
result = hq_peak_meme.apply(lambda X : test_length(X['chr'], X['start'], X['end'], chrom_size), axis=1)
result = np.vstack([np.asarray(t) for t in result])
hq_peak_meme[['start', 'end']] = result






######################################## output result #########################################################
CHIP_seq_output.to_excel('%s_target_gene.xlsx'%args.o, index=False)
CHIP_seq_output_sorted.to_excel('%s_target_gene_sorted.xlsx'%args.o, index=False)
hq_peak_meme.to_csv("%s_hq_peak_summit_flanking_500bp_ht_%d.bed"%(args.o, h_t), index=None, header=None, sep='\t')

with open('%s_hq_genes_ht_%d.txt'%(args.o, h_t), 'w+') as file:
    for gene in hq_genes:
        file.write(str(gene) + '\n')

with open('%s_hq_genes_KEGG_ht_%d.txt'%(args.o, h_t), 'w+') as file:
    for gene in hq_genes:
        file.write(str(gene) +'.2' + '\n')

with open('%s_log.txt'%args.o, 'w+') as file:
    file.write("summit file:%s"%args.summit_file)
    file.write("\n")
    file.write("gff:%s"%args.gff)
    file.write("\n")
    file.write("distance1:%s, distance_ratio:%s, distance2:%s"%(d1, d_r, d2))
    file.write("\n")
    file.write("type1_method:%s"%type1_method)
    file.write("\n")
    file.write("ncRNA_gene:%s"%ncRNA_gene)
    file.write("\n")
    file.write("mito_gene:%s"%mito_gene)


log.close()
######################################## end ###################################################################
print('\n\n')
elapsed_time = time.time() - start_time
print('end time:', time.ctime())
print('running time = %.2fs'%elapsed_time)
print('end')

