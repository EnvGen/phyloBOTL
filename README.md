# Requirements and installation

1. miniconda - Follow instruction on https://docs.conda.io/projects/miniconda/en/latest/

   make sure you have the conda version 23.7.3

         conda activate base
         conda install conda=23.7.3 -y

3. clone the repository:

         git clone https://github.com/lfdelzam/phyloBOTL

## Installation option 1 - step by step, in case you already have some programs installed or databases downloaded
Make sure you are using the same version as indicated hereafter:

Conda R_env:

         cd phyloBOTL
         conda env create -f conda_env/R_env.yaml -y
         conda activate R_env
         Rscript --vanilla support/Install_Rpackages.R <numberofcpus>

Snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html):

         conda create -n phylobotl_env -c bioconda snakemake=7.25.0 Python=3.11.4 -y
         conda activate phylobotl_env
         conda config --set channel_priority strict

Genomad and genomad_db (https://portal.nersc.gov/genomad/index.html)

          conda create -n genomad_env -c conda-forge -c bioconda genomad=1.6.1 -y
          conda activate genomad_env
          mkdir -p <directory_DB_path>/genomad_env/data
          genomad download-database <directory_DB_path>/genomad_env/data
          conda deactivate

EGGNOG: emapper.py and EGGNOG database (https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#user-content-v2112)

          conda create -n eggnog_mapper_env -c bioconda -c conda-forge eggnog-mapper=2.0.1 -y
          conda activate eggnog_mapper_env
          mkdir -p <directory_DB_path>/eggnog_mapper_env/data
          download_eggnog_data.py --data_dir <directory_DB_path>/eggnog_mapper_env/data -y
          conda deactivate

## Installation option 2 - all in one go, in case you want to install/download all the required programs and databases from scratch.

      cd phyloBOTL
      bash install_key_envs.sh <directory_DB_path> <cpus>

<directory_DB_path> a directory where you want to download the databases.
<cpus> number of cpus to download R packages when creating the R_env

## Optional

if you want to use gtdb or kSNP4 to build the phylogenetic tree

A. gtdbtk - source: https://ecogenomics.github.io/GTDBTk/installing/index.html#installing
          - conda: https://ecogenomics.github.io/GTDBTk/installing/bioconda.html

    conda create -n gtdb_env -c conda-forge -c bioconda gtdbtk=2.1.1 -y
    conda activate gtdb_env
    download-db.sh
    # Set the environment variable to the directory containing the GTDB-Tk reference data
    conda env config vars set GTDBTK_DATA_PATH="/path/to/unarchived/gtdbtk/data"


B. kSNP4  - source: https://sourceforge.net/projects/ksnp/files/latest/download
          - manual: https://sourceforge.net/projects/ksnp/files/kSNP4.1%20User%20Guide.pdf/download


# Usage

Modify config file:

      nano support/config_phylobotl.json

      "workdir":"/abs/path/to/phyloBOTL",
      "threads": 20,   - CPUs
      "input": {
                "File_list":"GENOME_LIST", - csv file <sample name>,<abs/path/to/contig.fasta[.gz]>,<group label>,<special_group_name>
                "group_1_label": "C", - your selection
                "group_1_name": "Clinical", - your selection
                "group_2_label": "E", - your selection
                "group_2_name": "Environmental", - your selection
                "Special_group_name": "None" - your selection or "None". For instance "Baltic Sea"
              },
      "tree_file": "Tree_rooted_boot.treefile", - tree file name. If you already have a tree, set the path to it here and the pipeline won't build the tree
      "tree_using": "iqtree",
      "iqtree": {
                "iqtree_params":"-m GTR+I+G -B 1000 -bnni",
                "iqtree_rooted": "TRUE" - or "FALSE"
                },
      "FastTree_params":"-gtr -nt -gamma",
      "GTDB": {
              "gtdb_params": "--bacteria",
              "taxa_filter":"p__Proteobacteria",
              "outgroup_taxon":"p__Firmicutes"
              },
      "KSNP": {
              "kSNP4_param_phylo_method":"ML", - Maximum Likelihood
              "kSNP4_param_genome_fraction": "core",
              "Path_to_kSNP4pkg": "/usr/local/kSNP4pkg"
              },
      "Prokka_params": "--rawproduct --quiet",
      "Proteome_dir": "proteome", - directory containing the protein files <sample name>.faa.
      "orthology":
                {
                "orthofinder_parameters": "-f proteome -a 4 -S diamond -og", - double check that name in -f <Proteom_dir> corresponds to "Proteome_dir" above
                "ortholog_count_table": "Orthogroups.GeneCount.tsv",
                "ortholog_table": "Orthogroups.tsv",
                "path_to_orthologs_sequences": "Orthologues/Results_dir/Orthogroup_Sequences"
                },if you provide a full path to existing "ortholog_count_table","ortholog_table" and "path_to_orthologs_sequences" files, the pipeline won't run orthofinder  
      "up_to_Eggnog_all_orthologs": False, - Set to True if you want to perform only the EggNOG annotation of orthologs groups
      "Eggnog_all_orthologs_selection": "longest", - Options : "longest", "random". criterium for selection one sequences from the groups of orthologs group, either the longest sequence or one randomly selected.  
      "path_to_eggnog_db":"directory_DB_path/eggnog_mapper_env/data", - abs path to eggnog database directory.        
      "PHYLOGLM": {
                  "phyloglm_Bootnumber":0,- phyloglm parameter, number of independent bootstrap replicates, 0 means no bootstrap.
                  "p_adj_value_cutoff":0.05,
                  "orthologue_ratio_in_genome_dataset":0.95,- In this case if the orthologue is present/absent in 95% of the genomes, it won't be considered in the analysis.
                  "phyloglm_btol_number":10, - phyloglm parameter, (logistic regression only) bound on the linear predictor to bound the searching space.
                  "phyloglm_outfiles_prefix":"Vv"
                  },
      "synteny_visualization": {
                            "present_in_at_least_n_genomes": 10, - parameter used when creating the Graph, representing the weigth cutoff to be included in graph
                            "locus_size_window": 40000, - bps, window used when looking for enriched orthologs locilised in the region
                            "max_n_genomes_per_fig": 22 - 22 representative genomes will be included in the .pdf figure
                          },
      "tax_specific_protein_annotation": {     - if you want to annotated Enriched and depleted orthologs using a specific database. For instance UniRef90 Vibrionacea
                                          "status": True,
                                          "path_to_database": "/abs/path/to/Your_Specific_database.fasta"
                                                },
      "genomad":{ "include": False, - Set to True in you want to predict plasmid and Virus on genomes using genomad
                  "path_to_genomad_db": "directory_DB_path/genomad_env/data/genomad_db", - abs path to genomad database directory
                  "Genomad_score_cut_off": 0.8, - Minimum score to classify a contig as plasmid (or Virus)
                  "params_genomad": "--cleanup --conservative --enable-score-calibration", - genomad parameters. The flag --enable-score-calibration will work if there are more than 1000 sequences in the fasta file (in a genome file)
                },
      "output_dir": "Results_with_IQTREE_rooted" - Output directory name




and save CTRL+x, y

Run:
    Option 1: Using a laptop or one node on a HPC.
      conda activate phylobotl_env
      snakemake -s phylobotl.smk --use-conda --conda-frontend conda --cores <numberofcpus>

    Option 2:

    If you want to run the pipeline using several nodes in a HPC:

		a. Download cookiecutter:

		conda create -n cookiecutter_env -c conda-forge cookiecutter -y
		conda activate  cookiecutter_env

		b. Create the profile directory and answer the questions:

		profile_dir="/abs/path/to/phyloBOTL/.config/snakemake"
		mkdir -p "$profile_dir"
		template="gh:Snakemake-Profiles/slurm"
		cookiecutter --output-dir "$profile_dir" "$template"
		conda deactivate

		c. Run the pipeline:

    conda activate phylobotl_env
    snakemake --profile .config/snakemake/phylobotl_slurm -s phylobotl.smk --use-conda --conda-frontend conda


# Output files

├── Genomes
│   ├── <isolate1>.fa

Annotations
│   ├── GBK_files
│   ├── GFF_files
│   ├── Orthologues
│   └── Rep_seq_Orth_groups

├── Genes
│   ├── <isolate1>.ffn

├── Orthologues
│   └── Results_dir
├── Pangenome_graph
│   ├── pangenome.h5
│   ├── MSA

├── phyloglm_input
│   ├── Annotations.txt
│   └── Orthogroups.tsv
├── Prokka_out

├── proteome
│   ├── <isolate1>.faa

├── GENOMAD
│   ├── Genomad_output_<isolate1>

├── Trees
│   ├── IQ_TREE.bionj


In <output_dir>:
├── Annotations
│   ├── Depleted_eggNOG
│   ├── Depleted_GBK_files
│   ├── Depleted_KEGG
│   ├── Depleted_KEGG_core
│   ├── Depleted_Loci
│   ├── Enriched_eggNOG
│   ├── Enriched_GBK_files
│   ├── Enriched_KEGG
│   ├── Enriched_KEGG_core
│   ├── Loci
│   ├── Specific_db
│   ├── Vv_Candidates_depleted_orthologues.tsv
│   ├── Vv_Candidates_enriched_orthologues.tsv
│   ├── Vv_core_depleted_orthologues.tsv
│   ├── Vv_core_enriched_orthologues.tsv
│   ├── Vv_depleted_core_orthologues_annotation.tsv
│   ├── Vv_Depleted_orthologues_annotation.tsv
│   ├── Vv_Enriched_core_orthologues_annotation.tsv
│   └── Vv_Enriched_orthologues_annotation.tsv
├── co_localization_figures
├── co_localization_figures_Depleted
├── Depleted_Orthologues_DNA_sequences
│   ├──<Ortholog1>.fna
├── Enriched_Orthologues_DNA_sequences
│   ├──<Ortholog1>.fna
├── Unsupervise
├── Unsupervise_core
├── Visualization
│   ├── Vv_Tree_with_core_orthologues_fan_branch.length.pdf
│   ├── Vv_Tree_with_depleted_core_orthologues_fan_branch.length.pdf
│   ├── Vv_Tree_with_enriched_core_orthologues_fan_branch.length.pdf
│   └── Vv_Tree_with_Selected_orthologues_fan_branch.length.pdf
├── Vv_Depleted_core.txt
├── Vv_Depleted.txt
├── Vv_Enriched_core.txt
└── Vv_Enriched.txt
