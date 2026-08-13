# Lutiebacter paper analyses
Final workflow and code for manuscript

# Kansas precipitation gradient map
Total annual precipitation for each Kansas county was collected from NOAA National Centers for Environmental Information (NCEI)
- [Total annual precipitation data](<Kansas precipitation gradient map/PrecipTotal-2010-2020.csv>)
- [Code related to Figure 1a Kansas annual precipitation map](<Kansas precipitation gradient map/Kansas_map_precip.R>)

# Quantitative osmotic stress phenotyping of select Luteibacter strains.
- [Supplemental Table S3. CellTiterBlue assay data](<Strain osmotic stress phenotype analysis/Supp.Table.S3.CellTiterBlue.assay.measurements.that.passed.qualifty.filtering.xlsx>)
- [Number of technical reps per biological reps](<Strain osmotic stress phenotype analysis/Supp.Table.S2. Number of technical reps per biological rep after filtering.xlsx>)
- [Code for the analyses related to Figures 1b and 2a-d](<Strain osmotic stress phenotype analysis/growth.assay.analyses.R>)

# Whole genome assembly and annotation
- Raw sequences and assembled whole genoms are available on NCBI SRA/GenBank under [BioProject PRJNA1300453](<https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1300453>)

- Bash scripts (created by Brian Sanderson):
    - [fastp.quality.filtering.sh](<Genome assembly and annotation/01_fastp.sh>)
    - [SPAdes.de.novo.assembly.sh](<Genome assembly and annotation/02_shovill.sh>)
    - [Annotation.with.Bakta.sh](<Genome assembly and annotation/03_bakta.sh>)

# Genome Structure (96 genomes)
### Create working enviornment
```sh
cd /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS
mkdir -p structure/{assemblies,annotations,mash,ani,phylogeny,metadata,results,scripts,logs}
```
- /assemblies contains .fasta files for genomes
- /annotations contains .gff files for genomes
- /metadata contains [genome/strain metadata file](genome_metadata_with_files.tsv)
- /mash contains a list of the names of the assemblies files (assemblies_list.txt)
- /phylogeny contains a list of the names of the annotation files (gff_list.txt)

## MASH (strains vs. strains)
run [01_run_mash.sh](<Genome structure/01_run_mash.sh>)
```sh
sbatch 01_run_mash.sh
```
#### convert MASH distance file to a square matrix
run [02_mash_to_matrix.py](<Genome structure/02_mash_to_matrix.py>)
```sh
python 02_mash_to_matrix.py
```
#### download files and visualize in R on local device
- mash_dist_matrix.tsv
- genome_metadata_with_files.tsv
- [MASH_dist.plot.R](<Genome structure/MASH_dist.plot.R>)
    - This produced Figure 3b
    - [Supp.Table.S6.MinHash distance matrix.xlsx](<Genome structure/Supp.Table.S6.MinHash distance matrix.xlsx>)

## FastANI (strains vs. strains)
run [04_run_fastani_within.sh](<Genome structure/04_run_fastani_within.sh>)
summarize results: [05_summarize_ani.R](<Genome structure/05_summarize_ani.R>)
- [Supp.Table.S7. ANI matrix.xls](<Genome structure/Supp.Table.S7. ANI matrix.xlsx>)

## Comparison to reference Luteibacter species genomes
### Downloaded Luteibacter genomes from NCBI
NCBI RefSeq
- Luteibacter pinisoli MAH-14; GCF_006385595.1; from rhizosphere of Korean Pine Tree
- L. rhizovicinus DSM 16549; GCF_001887595.1; from rhizosphere of barley, Denmark
- L. mycovicinus DSM 112764; GCF_000745235.1; fungal endophyte
- L. jiangsuensis W1I16; GCF_050434905.1; wheat rhizosphere, Washington State
- L. anthropi SM7.4; GCF_023699965.1; Artificial bog, Chapel Hill, USA
- L. aegosomaticola 335; GCF_023078475.1; from longhorn beetle, South Korea
- L. yeojuensis DSM 17673; GCF_011742875.1; greenhouse soil, South Korea
- L. flocculans EIF3; GCF_023612255.1;from eutrophic pond water, Germany

```sh
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/006/385/595/GCF_006385595.1_ASM638559v1/GCF_006385595.1_ASM638559v1_genomic.fna.gz
# repeat to get all the genomes
gunzip *.gz
```
### Run mash with the reference genomes
```sh
module load mash/2.3
cd /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS
mkdir -p mash_ref
# sketch your genomes + references together
mash sketch assemblies/*.fasta ../reference_genomes/*.fna -o mash_ref/luteibacter_all -p 8
# distances
mash dist mash_ref/luteibacter_all.msh mash_ref/luteibacter_all.msh > mash_ref/dist_to_flocculans.tsv
```
Closest matches:
- SVR strains = 0.09 - 0.13 mash dist to L. yeojuensis DSM 17673; GCF_011742875.1
- HAY strain = 0.11 mash dist to L. yeojuensis DSM 17673; GCF_011742875.1
- TLI & KNZ strains = ~0.08 mash dist to L. jiangsuensis W1I16; GCF_050434905.1

### Run FastANI comparing ref genomes to strains
```sh
module load fastani/1.34
fastANI --ql assemblies_list.txt --rl ref_list.txt -o mash_ref/ani_vs_refs.tsv -t 8
```
Closest matches:
- Across all strains: ANI ranged from 81.6588 to 90.6446 to the ref genomes
- SVR strains = 86.9248 - 87.5447 ANI to L. yeojuensis DSM 17673; GCF_011742875.1
- HAY strain = 85.7533 ANI to L. yeojuensis DSM 17673; GCF_011742875.1
- TLI & KNZ strains = 90.4 - 90.6 ANI to L. jiangsuensis W1I16; GCF_050434905.1

This suggests that all of our strains are novel species not represented in NCBI, but are all in the Luteibacter genus

## Roary Pangenome
more information on roary here: https://sanger-pathogens.github.io/Roary/

- Roary requires .gff files that have the annotation and then the fasta sequence at the end. This is what the output from Prokka looks like. Since we annotated using bactka, we need to convert the .gff file (symlinks) to this compatable format
submit slurm job:
    - [06a_make_roary_gffs.sh](Pangenome/06a_make_roary_gffs.sh)
#### Now run Roary (split paralogs)
By default, Roary uses synteny to split identical or highly similar sequences into separate clusters

[06_run_roary.sh](Pangenome/06_run_roary.sh)
Output summary:
```sh
Core genes        496   (in 99-100% of strains)
Soft core         19    (in 95-99% of strains)
Shell             6675  (in 15-95% of strains)
Cloud             5248  (in <15% of strains) rare, accessory, strain specific
Total            12438 
```
#### Repeated Roary with no split paralogs option (gets gene copy number information)
- "-s" and save new files in structure/phylogeny/roary_nosplit
- this clusters paralogs into the same cluster and keeps all of the locus tags, which can be converted to counts

### Create Core gene Phylogeny using IQ-Tree
Submit slurm job: [07_run_iqtree.sh](Pangenome/07_run_iqtree.sh)
- The Best-fit model by BIC: GTR+F+G4
- 7 identical sequences in the core alignment were ignored (clones); 89 strains?
- 96 sequences; 494,616 alignment columns; 62,813 parsimony-informative sites
- chi2 test failed because the genomes are compositionally heterogeneous, deep structure, as we already know.
```sh
NOTE: KNZ13-10G is identical to KNZ12-12B but kept for subsequent analysis
NOTE: KNZ2-8B is identical to KNZ2-1B but kept for subsequent analysis
NOTE: SVR3-9D is identical to SVR3-8D but kept for subsequent analysis
NOTE: TLI6-5F is identical to TLI6-11F but kept for subsequent analysis
NOTE: TLI7-5F is identical to TLI6-1G but kept for subsequent analysis
NOTE: TLI8-10B is identical to TLI6-3F but kept for subsequent analysis
Checking for duplicate sequences: done in 0.00171804 secs using 770% CPU
Identifying sites to remove: done in 0.299021 secs using 489.3% CPU
NOTE: 7 identical sequences (see below) will be ignored for subsequent analysis
NOTE: TLI6-8H (identical to TLI6-11F) is ignored but added at the end
NOTE: TLI6-IE (identical to TLI6-11F) is ignored but added at the end
NOTE: TLI7-7B (identical to TLI6-11F) is ignored but added at the end
NOTE: TLI7-7G (identical to TLI6-11F) is ignored but added at the end
NOTE: TLI8-6C (identical to TLI6-11F) is ignored but added at the end
NOTE: TLI9-11B (identical to TLI6-11F) is ignored but added at the end
NOTE: TLI9-3A (identical to TLI6-1G) is ignored but added at the end
```
Output files
- luteibacter_core.contree (bootstrap consensus tree; plot this)
    - I used iToL to plot this: [luteibacter_core.contree](Pangenome/luteibacter_core.contree)
    - The colapsed version is Figure 3a

- luteibacter_core.treefile (ML tree with branch lengths)

- Convert gene_presence_absence file into a binary format
    - Run [make_gene_presence_absence_binary.py](Pangenome/make_gene_presence_absence_binary.py)

#### Convert gene_presence_absence.csv to a gene count table
- use the no split paralogs Roary output
Run [08_make_gene_count_matrix_nosplit.py](Pangenome/08_make_gene_count_matrix_nosplit.py)
```sh
module load miniconda3
python 08_make_gene_count_matrix_nosplit.py
```
 summary of output counts
 ```sh
 Value counts across all strain-gene cells:
0    827149 <- absent
1    355982 <- single copy
2       533 <- duplicated
6         5 (hypothetical protein group_753)
7         3 (hypothetical proteins group_21 & group_753)
8         6 (hypothetical proteins group_21 & group_753)
9         2 (hypothetical protein group_21)
```
Then blast uncharacterized groups sequence on InterProScan (https://www.ebi.ac.uk/interpro/) & UniProt (https://www.uniprot.org/tool-dashboard)

### Interactively plotting in R on local device
[copy_number_var.R](Pangenome/copy_number_var.R)
- Plot single copy genes in venn diagrams comparing site and lineages
- Analyzed and plot copy number variance
- Compare genome size

# Functional annotations for Roary gene clusters (KO, KEGG, & Pfams)
### KOfamscan
- Run:
  - KOfamScan → KO IDs + KEGG hierarchy
- using the split paralogs Roary output (single copy)
- Input: /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary/pan_genome_reference.fa

Make new directory for functional annotations
```sh
cd /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS
mkdir -p structure/functional_annotation/{kofamscan,pfam,scripts,logs}
```
- The UCR HPC kofamscan module is using a ko_list and profile that have not been updated since July 2023.
  - Therefore, I downloading the most up to date ko_list from https://www.genome.jp/ftp/db/kofam/ (last updated March 25, 2026)
  - and I downloading the most up to date profiles from https://www.genome.jp/ftp/db/kofam/profiles.tar.gz (last updated March 26, 2026)
```sh
cd /bigdata/roperlab/nginn001/databases
wget -O ko_list.03.25.26.gz https://www.genome.jp/ftp/db/kofam/ko_list.gz
gunzip -f ko_list.03.25.26.gz

wget -O ko_profiles.03.26.26.tar.gz https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
tar -xzf ko_profiles.03.26.26.tar.gz #Extract
mv profiles ko_profiles.03.26.26 # rename extracted profiles directory
```
#### Create file with Roary groups, but .faa (amino acid seqs) format
Run [00_make_roary_cluster_representatives.py](<Functional annotations/00_make_roary_cluster_representatives.py>)
```sh
conda create -n roary_kofam python=3.9 biopython -y #instal biopython
conda activate roary_kofam

python 00_make_roary_cluster_representatives.py \
  --clustered_proteins /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary/clustered_proteins \
  --faa_dir /bigdata/roperlab/nginn001/KU_Luteibacter/faa.protein_files_bakta/ \
  --out_faa /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/kofamscan/roary_cluster_representatives.faa \
  --out_map /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/kofamscan/roary_cluster_representatives.tsv \
  --out_missing /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/kofamscan/roary_cluster_representatives.missing.tsv \
  --pick longest
```
Output
```sh
Loaded protein sequences: 357222
Clusters parsed: 12438
Representative sequences written: 12435
Clusters with no members found in .faa files: 3
Total missing cluster-member IDs: 23
```

Run script [01_run_kofamscan_roary_reps.slurm](<Functional annotations/01_run_kofamscan_roary_reps.slurm>) using the most up to date databases and the newly created roary_cluster_representatives.faa.

Output:
```sh
Input 12,435 representative protein clusters
370,855 total hit lines
10,872 unique representative clusters with at least one KO hit
87% of the input protein clusters got assigned a KO
6282 * hits: scores exceed the KO threshold (meaning they have high confidence)
```
### Parse KO hits and map to pfams and KEGG -> add gene count information
#### Download KEGG mapping files
 - Need to download the latest lists to match the newest KOs
 ```sh
 mkdir -p /bigdata/roperlab/nginn001/databases/kegg_maps_2026_04_15
 cd /bigdata/roperlab/nginn001/databases/kegg_maps_2026_04_15
# download
 curl -L https://rest.kegg.jp/list/ko > list_ko.2026_04_15.tsv
 curl -L https://rest.kegg.jp/list/pathway > list_pathway.2026_04_15.tsv
 curl -L https://rest.kegg.jp/list/module > list_module.2026_04_15.tsv
 curl -L https://rest.kegg.jp/link/pathway/ko > link_pathway_ko.2026_04_15.tsv
 curl -L https://rest.kegg.jp/link/module/ko > link_module_ko.2026_04_15.tsv
 curl -L https://rest.kegg.jp/get/br:ko00001 > brite_ko00001.2026_04_15.txt

 mkdir -p /bigdata/roperlab/nginn001/databases/pfam
 cd /bigdata/roperlab/nginn001/databases/pfam
 wget -O Pfam-A.hmm.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
 ```
#### Run hmmscan (Pfam) and parse
```sh
cd /bigdata/roperlab/nginn001/databases/pfam
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm #indexes the HMM files for scanning
```
Run [run_pfam_hmmscan_roary_reps.slurm](<Functional annotations/run_pfam_hmmscan_roary_reps.slurm>)

Then run: [04_parse_pfam_domtblout_besthit.py](<Functional annotations/04_parse_pfam_domtblout_besthit.py>)

```sh
python 04_parse_pfam_domtblout_besthit.py \
  --domtblout /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/pfam/pfam.domtblout \
  --out_tsv /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/pfam/pfam_besthit.tsv
```
#### Parse the BRITE file
run [parse_brite_ko00001.py](<Functional annotations/parse_brite_ko00001.py>)
```sh
module load miniconda3
python parse_brite_ko00001.py \
  --brite_txt /bigdata/roperlab/nginn001/databases/kegg_maps_2026_04_15/brite_ko00001.2026_04_15.txt \
  --out_tsv /bigdata/roperlab/nginn001/databases/kegg_maps_2026_04_15/brite_ko00001.parsed.tsv
```
Output file head should look like this:
```sh
KO	 brite_A_id	 brite_A_name	    brite_B_id	   brite_B_name	          brite_C_id	brite_C_name	brite_ko_entry
K00844	09100	Metabolism	09101	Carbohydrate metabolism	ko00010	Glycolysis / Gluconeogenesis	HK; hexokinase [EC:2.7.1.1]
K12407	09100	Metabolism	09101	Carbohydrate metabolism	ko00010	Glycolysis / Gluconeogenesis	GCK; glucokinase [EC:2.7.1.2]
K00845	09100	Metabolism	09101	Carbohydrate metabolism	ko00010	Glycolysis / Gluconeogenesis	glk; glucokinase [EC:2.7.1.2]
K25026	09100	Metabolism	09101	Carbohydrate metabolism	ko00010	Glycolysis / Gluconeogenesis	glk; 
```
#### Parse kofamscan results and map to tables
Run [02_parse_kofamscan_besthit.py](<Functional annotations/02_parse_kofamscan_besthit.py>)
This will:
- keeps the top KO hit per representative cluster... whether or not it had * (met the threshold)
- adds a column saying whether that top hit was a trusted hit (had a *)
- uses the pathway, module, BRITE, and Pfam lookup tables downloaded above

To run
```sh
python 02_parse_kofamscan_besthit.py \
>   --kofam_results /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/kofamscan/kofamscan_results.tsv \
>   --cluster_map /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/kofamscan/roary_cluster_representatives.tsv \
>   --list_ko /bigdata/roperlab/nginn001/databases/kegg_maps_2026_04_15/list_ko.2026_04_15.tsv \
>   --link_module_ko /bigdata/roperlab/nginn001/databases/kegg_maps_2026_04_15/link_module_ko.2026_04_15.tsv \
>   --list_module /bigdata/roperlab/nginn001/databases/kegg_maps_2026_04_15/list_module.2026_04_15.tsv \
>   --link_pathway_ko /bigdata/roperlab/nginn001/databases/kegg_maps_2026_04_15/link_pathway_ko.2026_04_15.tsv \
>   --list_pathway /bigdata/roperlab/nginn001/databases/kegg_maps_2026_04_15/list_pathway.2026_04_15.tsv \
>   --out_besthit /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/kofamscan/kofamscan_besthit_annotations.tsvp /bigdata/roperlab/nginn001/databases/kegg_maps_2026_04_15/kegg_ko_to_pfam.tsv
```
#### Put it all TOGETHER
Create file annotation folder
```sh
mkdir -p /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/final_annotations
```
run [05_join_all_annotations_and_build_matrices.py](<Functional annotations/05_join_all_annotations_and_build_matrices.py>)
```sh
python 05_join_all_annotations_and_build_matrices.py \
  --roary_rtab /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary/gene_presence_absence.Rtab \
  --besthit_annotations /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/kofamscan/kofamscan_besthit_annotations.tsv \
  --brite_parsed /bigdata/roperlab/nginn001/databases/kegg_maps_2026_04_15/brite_ko00001.parsed.tsv \
  --pfam_besthit /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/pfam/pfam_besthit.tsv \
  --out_cluster_join /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/final_annotations/roary_clusters_with_all_annotations.tsv \
  --out_annotation_long /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/final_annotations/cluster_annotation_long.tsv \
  --out_ko_matrix /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/final_annotations/ko_presence_absence.Rtab \
  --out_module_matrix /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/final_annotations/module_presence_absence.Rtab \
  --out_pathway_matrix /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/final_annotations/pathway_presence_absence.Rtab \
  --out_brite_B_matrix /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/final_annotations/brite_B_presence_absence.Rtab \
  --out_brite_C_matrix /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/final_annotations/brite_C_presence_absence.Rtab \
  --out_pfam_matrix /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/final_annotations/pfam_presence_absence.Rtab
  ```
Output summary:
```sh
Roary clusters loaded: 12438
Best-hit annotation rows loaded: 10872
BRITE KO rows loaded: 28189
Pfam rows loaded: 11272
Joined cluster rows written: 12438
KO features written: 3319
Module features written: 348
Pathway features written: 844
BRITE B features written: 46
BRITE C features written: 263
Pfam features written: 2290
```
### Interactively plotting in R on local device





# Genome-wide association studies (GWAS)
## SNP calling

## GWAS using Pyseer
#




# Phage related analyses
