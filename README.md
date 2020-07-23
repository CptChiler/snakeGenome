## Short read Illumina genome pipe

- QC of raw reads
- Short read assmembly
- Genome Annotation
- Species identification on reads(NCBI_taxon)
- AMR detection on assembly with CARD
- Plasmid detection on assembly plasmidfinder

### Dependencies:

- snakemake
- unicycler
- fastp
- kraken2
- abricate
- prokka

### Setup databases

```bash
abricate --check
abricate --setupdb
kraken2-build --download-library bacteria --threads 100 --db NCBI_bac_tax
```

```bash
prokka --setupdb
prokka --listdb  
```
in `config.json` you can change path to your database

```bash
{
	"md": "normal",
	"db1": "kraken2_db",
	"db2": "diamond_uniprot_db",
	"threads": 100
}
```

### check pipe

Clone the repo.

Put your Illumina data in folder `reads/` inside of repo and name them like this.

>Illumina reads: `TP1234_R1.fastq.bz2` and `TP1234_R2.fastq.bz2` 

in `config.json` you can change the mode of the assemblie with the 3 possible modes (see more info _Unicycler_).

> conservative   
normal  -> default
bold

```bash
snakemake --configfile config.json --forceall --dag | dot -Tpdf > dag.pdf
```

### run pipe

```bash
cd ...Illumina_genome_pipe/
snakemake --configfile config.json --cores 100
```
Enjoy the results in `/results`


### Troubleshooting

If prokka gets error blastASN reinstall it:
```bash
conda install prokka --force --yes
```

If unicycler failes at spades downgrade it:
```bash
conda install spades=3.11.0 --force --yes
```
