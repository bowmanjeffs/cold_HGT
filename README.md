cold_HGT
========

code for psychrophile HGT project

Scripts are run in the order given below.  Note that additional small scripts setting up directory structures and performing other simple tasks have not been included.

#### aquisition time analysis ####

## select_genomes.py: finds closest mesophilic relatives to psychrophiles

## get_select_genomes.py: pulls the genomes selected by previous script

make_genbank_cds_select.py: annotates fna files, using hmmscan and pfam-A to identify likely cds

## setup_gcamel.py: establishes directory structure and copies files for gcamel analysis

select_genomes_get_gc.py: identifies gc anomalies in the annotated genomes

gcamel_select_genomes.r: this is a modified version of Eric Collins gcamel program (see github repository)

gc_by_pos.py: a dependency of gcamel_select_genomes.r.  Does the same thing as gc_by_pos.pl, but in Python, for those of us who are deficient in Perl...

## combine_csv_files.py

# we have provided to *gz files with the output of the amerlioration dating estimate.  These files are the input for the following r script.

select_genomes_gc_amel_process.r

## select_genomes_GC_anomaly.r: evaluates GC anomalies in genomes

#### genome divergence ####

select_genomes_divergence.py: calculates compositional vectors and 16S distance matrix for genomes

parse_blast_xml.py: dependency of select_genomes_divergence.py

select_genomes_compositional_vector.r: evaluates distance between compositional vectors and compares with 16S distance matrix


