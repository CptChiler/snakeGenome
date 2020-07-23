import sys
import os

shell.executable('/bin/bash')
t = config['threads']

IDS, = glob_wildcards('reads/{sample}_R1.fastq.bz2')

rule all:
    input:
        expand('results/{sample}/{sample}_proteome.gbk', sample=IDS),
        expand('results/{sample}/{sample}_proteome.faa', sample=IDS),
        expand('results/{sample}/{sample}_genome.fasta', sample=IDS),
        expand('results/{sample}/{sample}_genome.gfa', sample=IDS),
        expand('results/{sample}/{sample}_assembly_QC.pdf', sample=IDS),
        expand('results/{sample}/assembly/proteome_{sample}', sample=IDS),
        expand('results/{sample}/{sample}_Plasmid_abricate.csv', sample=IDS),
        expand('results/{sample}/{sample}_AMR_abricate.csv', sample=IDS),
        expand('results/{sample}/{sample}_QC.html', sample=IDS),
        expand('results/{sample}/{sample}_QC.json', sample=IDS),
        expand('results/{sample}/{sample}_species.html', sample=IDS),
        expand('reads/QC_reads/{sample}_QC_R2.fastq.gz', sample=IDS),
        expand('reads/QC_reads/{sample}_QC_R1.fastq.gz', sample=IDS),
        expand('reads/{sample}_R2.fastq.gz', sample=IDS),
        expand('reads/{sample}_R1.fastq.gz', sample=IDS),

rule bz_2_gz_R1:
    input: 'reads/{sample}_R1.fastq.bz2',
	output:
	    temp('reads/{sample}_R1.fastq.gz'),
	threads: t
	shell: "lbunzip2 {input} -dc | pigz -4 -c -p {threads} > {output}"

rule bz_2_gz_R2:
    input: 'reads/{sample}_R2.fastq.bz2',
	output:
	    temp('reads/{sample}_R2.fastq.gz'),
	threads: t
	shell: "lbunzip2 {input} -dc | pigz -4 -c -p {threads} > {output}"

rule fastp:
    input:
        s1 = 'reads/{sample}_R1.fastq.gz',
        s2 = 'reads/{sample}_R2.fastq.gz',
    output:
        s1 = 'reads/QC_reads/{sample}_QC_R1.fastq.gz',
        s2 = 'reads/QC_reads/{sample}_QC_R2.fastq.gz',
        s3 = 'results/{sample}/{sample}_QC.html',
        s4 = 'results/{sample}/{sample}_QC.json',
    threads: t
    shell:
        '''
        fastp -i {input.s1} -I {input.s2} -o {output.s1} -O {output.s2} -h {output.s3} -j {output.s4} --detect_adapter_for_pe  --thread {threads}
        '''

rule kraken2:
    input:
        s1 = 'reads/QC_reads/{sample}_QC_R1.fastq.gz',
        s2 = 'reads/QC_reads/{sample}_QC_R2.fastq.gz',
    output:
        temp('scripts/kraken2/{sample}_kraken2.txt'),
    threads: t
    params:
         db = config['db1']
    shell:
        '''
	    kraken2 --quick --threads {threads} --db {params.db} --gzip-compressed --paired {input.s1} {input.s2} --output {output}
	    '''

rule kraken_plot:
    input: 'scripts/kraken2/{sample}_kraken2.txt',
    output: 'results/{sample}/{sample}_species.html',
    shell: 'ktImportTaxonomy -q 2 -t 3 {input} -o {output}'

rule unicycler:
    input:
        s1 = 'reads/QC_reads/{sample}_QC_R1.fastq.gz',
        s2 = 'reads/QC_reads/{sample}_QC_R2.fastq.gz',
    output:
        directory('results/{sample}/assembly/{sample}_genome_unicycler'),
    threads: t
    params:
        md = config['md']
    shell:
        '''
        unicycler -1 {input.s1} -2 {input.s2} -o {output} --mode {params.md} --no_correct --threads {threads} --keep 0 --verbosity 1
        '''

rule abricate_AMR:
    input:
        'results/{sample}/assembly/{sample}_genome_unicycler',
    output:
        'results/{sample}/{sample}_AMR_abricate.csv',
    threads: t
    shell:
        '''
        abricate --db card --csv --threads {threads} {input}/assembly.fasta > {output}
        '''

rule abricate_Plasmid:
    input:
        'results/{sample}/assembly/{sample}_genome_unicycler',
    output:
        'results/{sample}/{sample}_Plasmid_abricate.csv',
    threads: t
    shell:
        '''
        abricate --db plasmidfinder --csv --threads {threads} {input}/assembly.fasta > {output}
        '''

rule prokka:
    input:
        'results/{sample}/assembly/{sample}_genome_unicycler',
    output:
        directory('results/{sample}/assembly/proteome_{sample}'),
    threads: t
    shell:
        '''
        prokka {input}/assembly.fasta -prefix proteome --addgenes --compliant --outdir {output} --cpu {threads}
        '''

rule ideel_1:
    input:
        'results/{sample}/assembly/{sample}_genome_unicycler',
    output:
        temp('scripts/proteins/{sample}.faa'),
    shell:
        '''
        prodigal -a {output} -q -i {input}/assembly.fasta
        '''

rule ideel_2:
    input:
        'scripts/proteins/{sample}.faa',
    output:
        temp('scripts/lengths/{sample}.data'),
    threads: t
    params:
        db = config['db2'],
        of = '6 qlen slen',
    shell:
        '''
        diamond blastp --threads {threads} --max-target-seqs 1 --db {params.db} --query {input} --outfmt {params.of} --out {output}
        ''' 

rule ideel_3:
    input: 'scripts/lengths/{sample}.data',
    output: 'results/{sample}/{sample}_assembly_QC.pdf',
    script: 'scripts/hist.R'
 
rule unicycler_rename:
    input:
        'results/{sample}/assembly/{sample}_genome_unicycler',
    output:
        s1 = 'results/{sample}/{sample}_genome.fasta',
        s2 = 'results/{sample}/{sample}_genome.gfa',
    shell:
        '''
        cp '{input}/assembly.fasta' '{output.s1}'
        cp '{input}/assembly.gfa' '{output.s2}'
        '''

rule prokka_rename:
    input:
        'results/{sample}/assembly/proteome_{sample}',
    output:
        s1 = 'results/{sample}/{sample}_proteome.faa',
        s2 = 'results/{sample}/{sample}_proteome.gbk',
    shell:
        '''
        cp '{input}/proteome.faa' '{output.s1}'
        cp '{input}/proteome.gbk' '{output.s2}'
        '''
