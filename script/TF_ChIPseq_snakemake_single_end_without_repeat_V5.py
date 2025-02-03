# # config.yaml

# ######################### input files and reference files ###############################################
# FASTQ: /home/raiman/AflR_20241114/fastq/
# fqsuffix: fq
# gff: /home/raiman/genome_info/Anidulans/FungiDB-65_AnidulansFGSCA4.gff
# fasta: /home/raiman/genome_info/Anidulans/FungiDB-65_AnidulansFGSCA4_Genome.fasta
# gfolder: /home/raiman/genome_info/Anidulans/
# gaf: /home/raiman/reference/Anidulans/FungiDB-65_AnidulansFGSCA4_Curated_GO.gaf
# obo: /home/raiman/reference/go-basic.obo
# KEGG_relation: /home/raiman/reference/Anidulans/ani_KEGG_info.csv
# ############################ mapping parameters  #####################################################
# ncores : 4
# HT2index: /home/raiman/genome_info/Anidulans/Anidulans_hisat2_index
# trim5: 8
# ############################## call peaks parameters  #################################################
# genomesize: 3e7
# extsize: 200
# ############################## selectpeak parameters  #################################################
# h_t: 20 #high quality peak threshold, Peaks with summit bdg values greater than "h_t" will be classified as high quality.
# ############################### peak annotation parameters #############################################
# type1_method: ATG
# d1: 2000
# d2: 500
# d_r: 1.5
# ################################# meme parameters  ################################################
# db: /home/raiman/reference/motif_databases/JASPAR/JASPAR2022_CORE_fungi_non-redundant_v2.meme
# output: /home/raiman/AflR_20241114/output


configfile: "./config.yaml"

outputdir = config["output"]
FASTQdir = config["FASTQ"]
#samples = config["samples"]
gff_file = config["gff"]


samples= glob_wildcards("/home/raiman/AflR_20241114/fastq/{fname}.fq")[0]
samples = [sample for sample in samples]





rule all:
    input:
        ##QC
        expand(os.path.join(outputdir, "FastQC", "".join(["{sample}",  "_fastqc.zip"])), sample=samples),
        os.path.join(outputdir, "MultiQC", "multiqc_report.html"),
        
        ## reference use in bedGraphToBigWig 
        config["fasta"] + '.fai',
        os.path.join(config["gfolder"], "chrom.sizes"),
        
        #macs2 result
        expand(os.path.join(outputdir,"{sample}", "macs2", "".join(["{sample}", "_normalized.bdg"])), sample=samples),


        ## inhome script result
        expand(os.path.join(outputdir, "{sample}", "GO","".join(["{sample}", "_GO_result.xlsx"])), sample=samples),
        
        ## meme required .fasta file
        expand(os.path.join(outputdir, "{sample}", "meme","".join(["{sample}", "_hq_peak_summit_flanking_500bp_ht_", str(config["h_t"]), ".fasta"])), sample=samples),

        # write report
        expand(os.path.join(outputdir, "{sample}", "{sample}.html"), sample=samples)





## ------------------------------------------------------------------------------------ ##
## Quality control
## ------------------------------------------------------------------------------------ ##
## FastQC, original reads
rule fastqc:
    input:
        os.path.join(FASTQdir, "".join(["{sample}", ".", str(config["fqsuffix"])])),
    output:
        os.path.join(outputdir, "FastQC", "".join(["{sample}",  "_fastqc.zip"])),
    params:
        FastQC = lambda wildcards, output: os.path.dirname(output[0])  
    singularity:
        "docker://biocontainers/fastqc:v0.11.9_cv8"
    threads:
        config["ncores"]
    shell:
        "fastqc -o {params.FastQC} -t {threads} {input}"


## MultiQC
rule multiqc:
    input:
        expand(os.path.join(outputdir, "FastQC", "".join(["{sample}", "_fastqc.zip"])), sample=samples),
    output:
        os.path.join(outputdir, "MultiQC", "multiqc_report.html")
    params:
        inputdirs = os.path.join(outputdir, "FastQC"),
        MultiQCdir = lambda wildcards, output: os.path.dirname(output[0])  ## dirname of first output
    log:
        os.path.join(outputdir, "logs", "multiqc.log")
    singularity:
        "docker://ewels/multiqc:1.7"
    shell:
        "multiqc {params.inputdirs} -f -o {params.MultiQCdir} > {log} 2>&1"


rule samtools_faidx:
    input:
        fasta = config["fasta"]
    output:
        fastafai = config["fasta"] + '.fai'
    shell:
        "samtools faidx {input.fasta}"

rule cutfastafai:
    input:
        config["fasta"] + '.fai'
    output:
        os.path.join(config["gfolder"], "chrom.sizes")
    shell:
        "cut -f1,2 {input} > {output}"


rule mapping:
    input:
        a1 = os.path.join(FASTQdir, "".join(["{sample}", ".", str(config["fqsuffix"])])),
    output:
        os.path.join(outputdir, "{sample}", "hisat2", "".join(["{sample}", "_sorted.bam"]))
    params:
        index = config["HT2index"],
        trim5 = config["trim5"]
    singularity:
        "docker://zlskidmore/hisat2:2.1.0"
    shell:
        """
        hisat2 --trim5 {params.trim5} -x {params.index} -U {input.a1} | samtools view -@ 4 -bS - | samtools sort -@ 4  -O bam -o {output}
        samtools index {output}
        """
         



rule pileup_bdg:
    input:
        os.path.join(outputdir, "{sample}", "hisat2", "".join(["{sample}", "_sorted.bam"]))
    output:
        os.path.join(outputdir,"{sample}", "macs2", "".join(["{sample}", "_pileup.bdg"])) 
    shell:
        "macs2 pileup -f BAM --extsize 200 -i {input} -o {output}"

rule samtools_alignmentstat:
    input:
        os.path.join(outputdir, "{sample}", "hisat2", "".join(["{sample}", "_sorted.bam"]))
    output:
        os.path.join(outputdir, "{sample}", "hisat2", "".join(["{sample}", "_alignmentstat"]))
    shell:
        "samtools flagstat {input} > {output}"

rule normalized_bdg:
    input:
        pileup = os.path.join(outputdir,"{sample}", "macs2", "".join(["{sample}", "_pileup.bdg"])),
        alignment = os.path.join(outputdir, "{sample}", "hisat2", "".join(["{sample}", "_alignmentstat"]))
    output:
        tem1 = os.path.join(outputdir,"{sample}", "macs2", "".join(["tem_","{sample}", "_normalized.bdg"])) 
    params:
        tem2 = os.path.join(outputdir,"{sample}", "macs2", "".join(["tem_tem_","{sample}", "_normalized.bdg"])) 
    shell:
        """
        mappedReads=$(awk 'NR==8 {{print $1}}' {input.alignment})
        scale=$(bc -l <<< "scale=3; 1000000 / $mappedReads")
        echo " Normalizing {input.pileup} with factor $scale"
        macs2 bdgopt -i {input.pileup} -m multiply -p $scale -o {params.tem2}
        sed -n '2,${{p}}' {params.tem2} > {output.tem1}
        rm {params.tem2}
        """

rule bedGraphToBigWig:
    input:
        tem_bdg = os.path.join(outputdir,"{sample}", "macs2", "".join(["tem_","{sample}", "_normalized.bdg"])) ,
        chromsize =  os.path.join(config["gfolder"], "chrom.sizes")
    output:
        bdg = os.path.join(outputdir,"{sample}", "macs2", "".join(["{sample}", "_normalized.bdg"])),
        bw = os.path.join(outputdir,"{sample}", "macs2", "".join(["{sample}", "_normalized.bw"]))
    shell:
        """
        bedSort {input.tem_bdg} {output.bdg}
        bedGraphToBigWig {output.bdg} {input.chromsize} {output.bw} 
        """

rule call_peak:
    input:
        s = os.path.join(outputdir, "{sample}", "hisat2", "".join(["{sample}", "_sorted.bam",])),
    output:
        os.path.join(outputdir, "{sample}", "macs2", "withoutCtrl_narrow", "".join(["{sample}", "_peaks.narrowPeak"])),
        os.path.join(outputdir, "{sample}", "macs2", "withoutCtrl_narrow", "".join(["{sample}", "_peaks.xls"])),
        os.path.join(outputdir, "{sample}", "macs2", "withoutCtrl_narrow", "".join(["{sample}", "_summits.bed"]))
    params:
        name = "{sample}",
        outdir = os.path.join(outputdir, "{sample}", "macs2", "withoutCtrl_narrow"),
        genomesize = config["genomesize"],
        extsize = config["extsize"]
    shell:
        "macs2 callpeak -t {input.s} --name {params.name} --outdir {params.outdir} --nomodel -g {params.genomesize} --extsize {params.extsize} --call-summits"


rule peak_annotation:
    input:
        summit = os.path.join(outputdir, "{sample}", "macs2", "withoutCtrl_narrow", "".join(["{sample}", "_summits.bed"])),
        gff = gff_file,
        bdg = os.path.join(outputdir,"{sample}", "macs2", "".join(["{sample}", "_normalized.bdg"])),
    output:
        os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_target_gene.xlsx"])),
        os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_hq_genes_ht_", str(config["h_t"]), ".txt"])),
        os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_hq_genes_KEGG_ht_", str(config["h_t"]), ".txt"])),
        os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_hq_peak_summit_flanking_500bp_ht_", str(config["h_t"]), ".bed"]))
    params:
        o = os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}"])),
        type1_method = config["type1_method"],
        distance1 = config["d1"],
        distance2 = config["d2"],
        chromsize = os.path.join(config["gfolder"], "chrom.sizes"),
        d_r = config["d_r"],
        h_t = config["h_t"]
    shell:
        "python chipseq_gene_annotation_V5.5.py -summit_file {input.summit} -gff {input.gff} -o {params.o} --d1 {params.distance1} --d2 {params.distance2} --d_r {params.d_r} -type1_method {params.type1_method} -bdg_file {input.bdg} --h_t {params.h_t} -chrom_size {params.chromsize}"


rule GO_analysis:
    input:
        os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_hq_genes_ht_", str(config["h_t"]), ".txt"]))
    output:
        os.path.join(outputdir, "{sample}", "GO","".join(["{sample}", "_GO_result.xlsx"]))
    params:
        o = os.path.join(outputdir, "{sample}", "GO","".join(["{sample}", "_GO_result"])),
        gaf = config["gaf"],
        obo = config["obo"]
    shell:
        "python GO_analysis_V1.2.py -target_gene_file {input} -gaf {params.gaf} -obo {params.obo} -o {params.o}"

rule meme:
    input:
        bed = os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_hq_peak_summit_flanking_500bp_ht_", str(config["h_t"]), ".bed"])),
        db = config["db"]
    output:
        fasta = os.path.join(outputdir, "{sample}", "meme","".join(["{sample}", "_hq_peak_summit_flanking_500bp_ht_", str(config["h_t"]), ".fasta"])),
        meme = directory(os.path.join(os.path.join(outputdir, "{sample}", "meme", "meme_result")))
    params:
        fasta = config["fasta"]
    shell:
        """
        bedtools getfasta -fi {params.fasta} -bed {input.bed} -fo {output.fasta}
        meme-chip -oc {output.meme} -time 240 -ccut 100 -dna -order 2 -minw 6 -maxw 15 -db {input.db} -meme-mod zoops -meme-nmotifs 6 -meme-searchsize 100000 -meme-p 32 -streme-pvt 0.05 -streme-align center -streme-totallength 4000000 -centrimo-score 5.0 -centrimo-ethresh 10.0 {output.fasta}
        """

rule report:
    input:
        directory(os.path.join(os.path.join(outputdir, "{sample}", "meme", "meme_result"))),
        os.path.join(outputdir, "{sample}", "GO","".join(["{sample}", "_GO_result.xlsx"]))
    output:
        os.path.join(outputdir, "{sample}", "{sample}.html"),
    params:
        output_m_path = directory(os.path.join(outputdir, "{sample}"))
    shell:
        "python TF_write_html_report_V2.py -sample_name {wildcards.sample} -output_m_path {params.output_m_path}"
