# Lutiebacter paper analyses
## *Not peer-reviewed yet*
### Workflow and code for draft manuscript to be shared with co-authors and reviewers.

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

# Roary Pangenome
more information on roary here: https://sanger-pathogens.github.io/Roary/

- Roary requires .gff files that have the annotation and then the fasta sequence at the end. This is what the output from Prokka looks like. Since we annotated using bactka, we need to convert the .gff file (symlinks) to this compatable format
submit slurm job:
    - [06a_make_roary_gffs.sh](Pangenome/06a_make_roary_gffs.sh)
### Now run Roary (split paralogs)
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
### Repeated Roary with no split paralogs option (gets gene copy number information)
- "-s" and save new files in structure/phylogeny/roary_nosplit
- this clusters paralogs into the same cluster and keeps all of the locus tags, which can be converted to counts

## Create Core gene Phylogeny using IQ-Tree
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

### Convert gene_presence_absence.csv to a gene count table
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
### Run hmmscan (Pfam) and parse
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
### Parse the BRITE file
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
### Parse kofamscan results and map to tables
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
### Put it all TOGETHER
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
## Annotate Roary no split paralog output clusters.
Input: /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary_nosplit/pan_genome_reference.fa
#### Create file with Roary groups, but .faa (amino acid seqs) format
Run [00_nosplit.py](<Functional annotations/00_nosplit.py>)
```sh
conda create -n roary_kofam python=3.9 biopython -y #instal biopython
conda activate roary_kofam

python 00_nosplit.py \
  --clustered_proteins /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary_nosplit/clustered_proteins \
  --faa_dir /bigdata/roperlab/nginn001/KU_Luteibacter/faa.protein_files_bakta/ \
  --out_faa /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/kofamscan/no_split/roary_no_split_cluster_representatives.faa \
  --out_map /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/kofamscan/no_split/roary_no_split_cluster_representatives.tsv \
  --out_missing /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/functional_annotation/kofamscan/no_split/roary_no_split_cluster_representatives.missing.tsv \
  --pick longest
```
Output
```sh
Loaded protein sequences: 357222
Clusters parsed: 12330
Representative sequences written: 12327
Clusters with no members found in .faa files: 3
Total missing cluster-member IDs: 23
```
Run script [01_nosplit](<Functional annotations/01_nosplit>) using the most up to date databases and the newly created roary_cluster_representatives.faa.
### Interactively plotting in R on local device
Files produced above and used in this analysis
- [roary_clusters_with_all_annotations.tsv](<Functional annotations/roary_clusters_with_all_annotations.tsv>) (Supplemental Table S12 in manuscript)
- [genome_metadata_with_files.tsv](<genome_metadata_with_files.tsv>)

Analyses and plotting code:
[functional.enrichment.R](<Functional annotations/functional.enrichment.R>)

Outfiles and figures
- produces Figure 4a-c and Supplemental Figure S5
- [sig.ko.sorted.in.cats.filt.5_1_26.unit.names.xlsx](<Functional annotations/sig.ko.sorted.in.cats.filt.5_1_26.unit.names.xlsx>)
- [Supp.Table.S13.BRITE.B.KW.results.xlsx](<Functional annotations/Supp.Table.S13.BRITE.B.KW.results.xlsx>)
- [Supp.Table.S14.KO.KW.results.xls](<Functional annotations/Supp.Table.S14.KO.KW.results.xlsx>)

# Genome-wide association studies (GWAS)
questions:
- Which genes/SNPs are associated with osmotic tolerance?
- Which genomic features differ across species/site groups?
- Which signals persist within species, not just across species?
- GWAS by site or region is too confounded by species
### Input files
- /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary/gene_presence_absence.Rtab
- near_clone_drop.txt
  - I created this file locally and uploaded to HPC
  - Removed strains that have 99.9-100% similarity, MASH dist of 0, and have very similar phenotype
  - SNP VCF file
### Strain derplication
- For the dereplicated strain GWAS, this is a phenotype-aware, maximal-retention filtering pass, so it is intentionally conservative. It removes only the most obvious redundancies with the following rules:
  - strains with >99.99% ANI were evaluated
  - For strain clusters >99.99% ANI, but had very different phenotypes, I kept strains with differening phenotypes
  - When possible I prioitized keeping strains that have phenotype data for both NaCl and PEG assays.
This resulted in the exclusion of 20 strains
exclusion strains:
```sh
KNZ2-5D
KNZ2-6D
KNZ12-10F
KNZ12-1A
KNZ12-7H
KNZ1-1H
TLI6-11F
TLI6-IE
TLI6-1G
TLI6-5H
TLI6-8H
TLI6-4A
TLI7-6A
TLI8-4H
TLI8-9D
TLI8-10B
TLI9-1C
TLI9-11B
TLI9-11C
SVR3-5F
```
## SNP calling (Create SNP VCF file)
- use snippy/4.6.0 module on the HPC
- From my genomes I selected KNZ12-1B and SVR3-8D as the reference genomes because:
  - I wanted to use specie 1 as a reference
  - This strain had the lowest number of contigs and had phenotyped for both NaCl and PEG within their sites
```sh
cd /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS
mkdir -p gwas/{inputs,scripts,logs,snippy/all,snippy/sp1,vcf}
```
run [06_make_snippy_sample_lists.py](GWAS/snippy/06_make_snippy_sample_lists.py)
```sh
module load miniconda3
conda create -n gwas_py python=3.9 pandas -y
conda activate gwas_py
python /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/scripts/06_make_snippy_sample_lists.py
```
output
```sh
Wrote:
/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/inputs/snippy_all_samples.tsv
/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/inputs/snippy_sp1_samples.tsv
All genomes: 96
sp1 genomes: 56
```
Now run: 
-  [7_run_snippy_all.sh](GWAS/snippy/07_run_snippy_all.sh) (KNZ reference)
-  [08_run_snippy_sp1.sh](GWAS/snippy/08_run_snippy_sp1.sh) (KNZ reference)
-  [09_run_snippy_all_SVRref.sh](GWAS/snippy/09_run_snippy_all_SVRref.sh) (SVR reference)

This produced snippy results foleder per strain for each run. Note: I used FASTA files, so the SNPs are not annotated. We can map to significant SNPs back to the .gff3 files after the GWAS.

- ** I decided to make the SVR reference just for SVR strains. So, I deleted the non-SVR folders produced by "09_run..."

Now run snippy core to combine all of the per strain SNPs into a single dataset to get a core (multi-genomes) SNP alignment. Note: it will only keep SNP positions where all genomes have a valid base call (core):
- [10_run_snippy_core_all.sh](GWAS/snippy/10_run_snippy_core_all.sh)
- [11_run_snippy_core_sp1.sh](GWAS/snippy/11_run_snippy_core_sp1.sh)
- [12_run_snippy_core_all_SVRref.sh](GWAS/snippy/12_run_snippy_core_all_SVRref.sh) <- this will be SVR strains ONLY

## GWAS using Pyseer
### GWAS plan (12 separate GWASs):
Genome subsets:
1. all genomes 
2. all genomes, minus the clones in the "clones_to_exclude.txt" 
3. sp1 genomes 
4. sp1 genomes, minus the clones in the "clones_to_exclude.txt" 
5. SVR site genomes
6. SVR site genomes, minus the clones in the "clones_to_exclude.txt"
Phenotypes:
1. log10_relative_growth_nacl
2. log10_relative_growth_peg
Type (both can be completed in one run/job):
1. Gene presence/absence
2. SNPs

### Input files:
```sh
# MASH distance (square): controls for genome structure
/bigdata/KU_Luteibacter/microGWAS/structure/mash/mash_dist_matrix.tsv
# Metadata
/bigdata/KU_Luteibacter/microGWAS/structure/metadata/genome_metadata_with_files.tsv
# Clone list
/bigdata/KU_Luteibacter/microGWAS/gwas/inputs/clones_to_exclude.txt
# Gene matrix
/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary/gene_presence_absence.Rtab
# Master VCFs from snippy-core
/bigdata/KU_Luteibacter/microGWAS/gwas/vcf/all_genomes_knz12_1b/core.vcf.gz
/bigdata/KU_Luteibacter/microGWAS/gwas/vcf/sp1_knz12_1b/core.vcf.gz
/bigdata/KU_Luteibacter/microGWAS/gwas/vcf/SVR_genomes_svr3_8d/core.vcf.gz
```
### Pipeline details
Each of the 12 runs will create a run directory with:
- phenotypes.tsv
- strains.txt
- genes.Rtab
- mash.tsv
- variants.vcf.gz
- variants.vcf.gz.csi

Then, a Slurm job runs:
pyseer --pres genes.Rtab ... > genes_pyseer.tsv
pyseer --vcf variants.vcf.gz ... > snps_pyseer.tsv

### Make run folders and subset input files
run [01_prepare_pyseer_runs.py](GWAS/pyseer/01_prepare_pyseer_runs.py)
```sh
module load bcftools
source /opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate gwas_py
python 01_prepare_pyseer_runs.py
```
output
```sh
Prepared /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/runs/01_all_nacl with 94 strains
Prepared /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/runs/01_all_peg with 71 strains
Prepared /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/runs/02_all_noclone_nacl with 75 strains
Prepared /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/runs/02_all_noclone_peg with 60 strains
Prepared /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/runs/03_sp1_nacl with 54 strains
Prepared /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/runs/03_sp1_peg with 38 strains
Prepared /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/runs/04_sp1_noclone_nacl with 36 strains
Prepared /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/runs/04_sp1_noclone_peg with 28 strains
Prepared /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/runs/05_SVR_nacl with 38 strains
Prepared /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/runs/05_SVR_peg with 31 strains
Prepared /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/runs/06_SVR_noclone_nacl with 37 strains
Prepared /bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/runs/06_SVR_noclone_peg with 30 strains
```
### Install pyseer using interactive job on HPC
```sh
srun -p short -t 02:00:00 --mem=32G --cpus-per-task=4 --pty bash
module load miniconda3
conda create -n pyseer_env -c conda-forge -c bioconda -c defaults pyseer -y
# once it completes
exit # leave the compute node interactive job
```
test that it worked
```sh
conda activate pyseer_env
which pyseer
pyseer --version
# pyseer 1.4.1
```
## Pyseer runs
Run [03_submit_all_pyseer.sh](GWAS/pyseer/03_submit_all_pyseer.sh), which will pass each of the 12 genomes subsets to [02_run_pyseer_subset.sh](GWAS/pyseer/02_run_pyseer_subset.sh). This will submit 12 jobs running 24 GWAS.
```sh
bash 03_submit_all_pyseer.sh
```
I tested the 02 scripts before launching all of the runs with the 03 script
```sh
sbatch 02_run_pyseer_subset.sh \
/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/gwas/runs/01_all_nacl
```
### Map SNP positions back to reference genome .gff3 file to get gene names
Currently the names are locus positions, e.g. contig_1_1144_C_G
We need to extract the locus tag / gene ID / gene name / product from the GFF file
Two references
- KNZ12-1B.gff3 for all and sp1 runs
- SVR3-8D.gff3 for SVR subset
Then this will map to our full final annotations to get KO, KEGG, and pfam annotations for the SNPs

RUN [05_submit_all_snp_annotation.sh](GWAS/pyseer/05_submit_all_snp_annotation.sh) which will submit a slurm job for each folder using [04_annotate_pyseer_snps.slurm.sh](GWAS/pyseer/04_annotate_pyseer_snps.slurm.sh), which passes to use this python script[04_annotate_pyseer_snps.py](GWAS/pyseer/04_annotate_pyseer_snps.py). 

```sh
bash 05_submit_snp_annotations.sh 
```
## GWAS hit analyses
- moved all of the genes_pyseer.tsv and snps_pyseer.tsv files to local device and work in Rstudio.
### Narrow down candidate genes
None of the FDR corrected p-values were significant, likely because my # of strains is pretty low. So, I focused on gene clusters where ltr-p-value was <0.10. Then if a gene cluster was a hit for 4 or more GWAS runs within either gene or SNPs focused runs I consider it a candidate. I also considered gene clusters that met that p-value cut off and was significant in the sp1 no clones run a candidate.
- Candates fall into the following catagories
  - osmotic
  - membrane_transport
  - motility
  - metabolism
  - phage
- Top candidates based on annotations
  - otsB (trehalose)
  - efflux pump
  - glycosyltransferase
  - trxA
  - motility gene

beta > 0 → allele increases phenotype
beta < 0 → allele decreases phenotype

### Sorting and plotting code
[pyseer.GWAS.results.R](GWAS/pyseer/pyseer.GWAS.results.R)

Outputs:
- [Supp.Table.S15.GWAS.candidate.stats.xlsx](GWAS/pyseer/Supp.Table.S15.GWAS.candidate.stats.xlsx)
- Figure 5
- Table 1

# Phage related analyses
Prophage regions were annotated from Luteibacter whole genome nucleotide sequences using PHASTER (PHAge Search Tool Enhanced Release) webserver
Citation: Arndt, D., Grant, J., Marcu, A., Sajed, T., Pon, A., Liang, Y., Wishart, D.S. (2016) PHASTER: a better, faster version of the PHAST phage search tool. Nucleic Acids Res., 44(Web Server issue): W16-W21
- Summary of results: [Supp.Table.S16.Prophage.information.xlsx](Phage/Supp.Table.S16.Prophage.information.xlsx)

## Creating prophage tree file
1. cluster phage regions: [cluster.phage.sh](Phage/cluster.phage.sh) # identified 8 clusters
2. Alignment phage sequences with MAFFT, trim with trimAI, and build tree with IQtree: [full_phylogeny.sh](Phage/full_phylogeny.sh)

## R analysis and plotting
Analysis completed in R using: [lute.phage.analysis.May2026.R](Phage/lute.phage.analysis.May2026.R)

### Input files
- [full.phage.region.treefile](Phage/full.phage.region.treefile)
- [luteibacter_core.treefile](Phage/luteibacter_core.treefile) (this is an output from Roary)
- [phaster_phage_prediction_results.xlsx](Phage/phaster_phage_prediction_results.xlsx)

### outputs
- Figure 6
- Supplemental Figure S6

