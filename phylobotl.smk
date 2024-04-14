import os

#Usage: snakemake -s phylobotl.smk --use-conda --conda-frontend conda --cores 20
#Optional: --use-envmodules --use-singularity

configfile: "./support/config_phylobotl.json"
workdir: config["workdir"]

def read_File_list(myfile):
  strs=[]
  fls=[]
  with open(myfile, "r") as fin:
    for line in fin:
      line=line.rstrip()
      line=line.split(",")
      strs.append(line[0])
      fls.append(line[1])

  return strs, fls


strains, files =  read_File_list(config["input"]["File_list"])

def files_rule_all():
    file_list=[]
    if config["up_to_Eggnog_all_orthologs"]:
      file_list=["Annotations/Rep_seq_Orth_groups/All_EggNOG/selected.emapper.annotations"]
      return file_list

    if config["genomad"]["include"]:
      file_list=[expand("GENOMAD/Genomad_output_{s}/{s}_summary.log",s = strains),"GENOMAD/Summary_genomad/Genomes_with_plasmid.tsv",
      expand("Genomes/{s}.fa", s = strains)] 

    commun_files=[config["output_dir"]+"/Enriched_Orthologues_DNA_sequences/report.txt",config["output_dir"]+"/Depleted_Orthologues_DNA_sequences/report.txt",
       config["output_dir"]+"/Annotations/Enriched_GBK_files/Genomes_without_selected_orthologs.txt",
       config["output_dir"]+"/Annotations/Depleted_GBK_files/Genomes_without_selected_orthologs.txt",
       config["output_dir"]+"/Annotations/Enriched_KEGG/Kegg_annotations.tsv",config["output_dir"]+"/Annotations/Depleted_KEGG/Kegg_annotations.tsv",
       "Pangenome_graph/pangenome.h5",
       config["output_dir"]+"/co_localization_figures/Leiden_clusters_using_max_loci_size_"+str(config["synteny_visualization"]["locus_size_window"])+"_orthologues.ft",
       config["output_dir"]+"/co_localization_figures_Depleted/Leiden_clusters_using_max_loci_size_"+str(config["synteny_visualization"]["locus_size_window"])+"_orthologues.ft",
       config["output_dir"]+"/Annotations/Enriched_KEGG_core/Kegg_annotations.tsv",
       config["output_dir"]+"/Annotations/Depleted_KEGG_core/Kegg_annotations.tsv"]

    if config["tax_specific_protein_annotation"]["status"]:
       file_list.append(commun_files)
       file_list.append([config["output_dir"]+"/Annotations/Specific_db/Enriched_orthologs_most_frequent_annotation.tsv",
       config["output_dir"]+"/Annotations/Specific_db/Depleted/Depleted_orthologs_best_annotation.tsv",
       config["output_dir"]+"/Annotations/Loci/Specific_db/cluster_most_frequent_annotations.tsv"])  
    else:
        file_list.append(commun_files)
    
    return file_list


rule all:
    input: files_rule_all()

rule Genomes:
    input: config["input"]["File_list"]
    output: expand("Genomes/{s}.fa", s = strains)
    threads: 1
    resources:
        runtime = lambda wildcards, attempt: attempt*60 * 4, mem_mb=6400	
    shell:  """
                cat {input} | while read line
                        do
                            if [[ "$line" != \#* ]]; then
                                s=$(echo $line | cut -d ',' -f1)
                                f=$(echo $line | cut -d ',' -f2)
                                if [[ ! -s Genomes/$s.fa ]]; then
                                  if [[ $f =~ \.gz$ ]]; then
                                    gunzip -cd $f > Genomes/$s.fa
                                  else
                                    cp $f Genomes/$s.fa
                                  fi
                                fi
                            fi
                        done
            """


if config["tree_using"] == "gtdb":
	rule tree:
  		input: expand("Genomes/{s}.fa", s = strains)
  		output: config["tree_file"]
  		params: txf=config["GTDB"]["taxa_filter"], og=config["GTDB"]["outgroup_taxon"], gt=config["GTDB"]["gtdb_params"]
  		threads: config["threads"]
		resources:
        		runtime = lambda wildcards, attempt: attempt*60 *48, mem_mb=6400*config["threads"]
  		conda: "gtdb_env"
  		shell: """
             		gtdbtk de_novo_wf --genome_dir Genomes --taxa_filter {params.txf} --outgroup_taxon {params.og} {params.gt} --out_dir GTDBTK_out --extension .fa --cpus {threads}
			        cp GTDBTK_out/gtdbtk.*decorated.tree {output}
         	   """

elif config["tree_using"] == "ani":

  rule fastANI:
      input: expand("Genomes/{s}.fa", s = strains)
      output: "fastANI_output/results"
      threads: config["threads"]
      resources:
                runtime = lambda wildcards, attempt: attempt*60 *48, mem_mb=6400*config["threads"]
      container: "docker://staphb/fastani"
      conda: "conda_envs/fastani.yaml"
      shell: """
			     ls Genomes/*.fa > list_temp
                 fastANI --ql list_temp --rl list_temp -t {threads} -o {output}
			     rm list_temp
             """

elif config["tree_using"] == "ksnp":

  rule kSNP_tree:
    input: expand("Genomes/{s}.fa", s = strains)
    output: config["tree_file"]
    params: m=config["KSNP"]["kSNP4_param_phylo_method"], f=config["KSNP"]["kSNP4_param_genome_fraction"], k=config["KSNP"]["Path_to_kSNP4pkg"]
    threads: config["threads"]
    resources:
              runtime = lambda wildcards, attempt: attempt*60 *48, mem_mb=6400*config["threads"]
    container: "docker://staphb/ksnp4"
    shell: """
               export PATH="{params.k}:$PATH"
               mkdir -p kSNP4_Results
               bash support/run_kSNP4.sh Genomes kSNP4_inFile kSNP4_Results {threads} -{params.m} -{params.f}
			   cp kSNP4_Results/tree.{params.f}_SNPs.{params.m}.tre {output}
           """


rule prokka:
  input: "Genomes/{s}.fa"
  output: "Prokka_out/{s}.txt"
  params: f="Genomes/{s}.fa", dout="Prokka_out", prx="{s}" ,p=config["Prokka_params"]
  threads: config["threads"]
  resources:
          runtime = lambda wildcards, attempt: attempt*60 *10, mem_mb=6400*config["threads"]
  container: "docker://staphb/prokka"
  conda: "conda_envs/prokka.yaml"
  shell:  """
            prokka --outdir {params.dout} --cpus {threads} --prefix {params.prx} {params.f} --force {params.p} 
          """

rule Proteome_dir:
  input: "Prokka_out/{s}.txt"
  output: config["Proteome_dir"]+"/{s}.faa"
  params: d=config["Proteome_dir"], s="{s}"
  threads: 1
  resources:
          runtime = lambda wildcards, attempt: attempt*60 *0.2, mem_mb=6400	
  shell: """
            mv Prokka_out/{params.s}.faa {params.d}/.
         """

rule Gene_dir:
  input: "Prokka_out/{s}.txt"
  output: "Genes/{s}.ffn"
  params: s="{s}"
  threads: 1
  resources:
           runtime = lambda wildcards, attempt: attempt*60 *0.2, mem_mb=6400
  shell: """
            mv Prokka_out/{params.s}.ffn Genes/.
         """

rule annot_dir:
  input: "Prokka_out/{s}.txt"
  output: "Annotations/GBK_files/{s}.gbk", "Annotations/GFF_files/{s}.gff"
  params: s="{s}"
  threads: 1
  resources:
          runtime = lambda wildcards, attempt: attempt*60 *0.5, mem_mb=6400
  shell: """
            sed '/^##FASTA/q' Prokka_out/{params.s}.gff > Annotations/GFF_files/{params.s}.gff
            egrep -v '^(ACCESSION|VERSION)' Prokka_out/{params.s}.gbk > Annotations/GBK_files/{params.s}.gbk
         """

rule pangenome_graph:
     input: expand("Annotations/GBK_files/{s}.gbk", s= strains)
     output: "Pangenome_graph/pangenome.h5"
     threads: config["threads"]
     conda: "conda_envs/ppanggolin.yaml"
     resources:
           runtime = lambda wildcards, attempt: attempt*60 *20, mem_mb=6400*config["threads"]
     shell:  """
               ls {input} | while read f; do name=$(basename $f); n=$( echo $name | sed s/'\.gbk'//); echo "$n\t$f"; done > annot_list
               ppanggolin all --anno annot_list -c {threads} -o Pangenome_graph -f
               rm annot_list
             """

if config["tree_using"] == "FastTree" or config["tree_using"] == "iqtree":
  rule ppanggolin_msa:
    input: "Pangenome_graph/pangenome.h5"
    output: "Pangenome_graph/MSA/core_genome_alignment.aln"
    threads: config["threads"]
    conda: "conda_envs/ppanggolin.yaml"
    resources:
           runtime = lambda wildcards, attempt: attempt*60 *20, mem_mb=6400*config["threads"]
    shell:  "ppanggolin msa -p Pangenome_graph/pangenome.h5 --partition core --source dna -o Pangenome_graph/MSA --phylo -c {threads} -f "

if config["tree_using"] == "FastTree":
  rule FastTree:
    input: "Pangenome_graph/MSA/core_genome_alignment.aln"
    output: config["tree_file"]
    threads: 1
    resources:
            runtime = lambda wildcards, attempt: attempt*60 *30, mem_mb=6400
    params: config["FastTree_params"]
    container: "docker://staphb/fasttree"
    conda: "conda_envs/fastTree2.yaml"
    envmodules:"bioinfo-tools/FastTree/2.1.10"
    shell: "FastTree {params} {input} > {output} "

if config["tree_using"] == "iqtree":
  rule IQTREE:
    input: "Pangenome_graph/MSA/core_genome_alignment.aln"
    output: config["tree_file"]
    threads: config["threads"]
    params: i=config["iqtree"]["iqtree_params"], r=config["iqtree"]["iqtree_rooted"], o="Trees"
    container: "docker://staphb/iqtree2"
    conda: "conda_envs/iqtree.yaml"
    resources:
             runtime = lambda wildcards, attempt: attempt*60 *60, mem_mb=6400*config["threads"]
    envmodules:"bioinfo-tools/iqtree"
    shell:  """
              mkdir -p {params.o}
              iqtree -s {input} -T {threads}  --seqtype DNA --prefix IQ_TREE {params.i}
              if [[ {params.r} == "TRUE" ]]; then
                iqtree2 -m 12.12 -s {input}  --seqtype DNA -T {threads} -te IQ_TREE.treefile --prefix IQTREE_rooted
                cp IQTREE_rooted.treefile {output}
                mv IQTREE_rooted.* {params.o}/.
		        mv IQ_TREE.* {params.o}/.
              else
                cp IQ_TREE.treefile {output}
                mv IQ_TREE.* {params.o}/.
              fi
            """



rule orthofinder:
  input: expand(config["Proteome_dir"]+"/{s}.faa", s = strains)
  output: f=config["orthology"]["ortholog_count_table"],g=config["orthology"]["ortholog_table"]
  threads: config["threads"]
  params: p=config["orthology"]["orthofinder_parameters"], f="Orthologues/Results_dir/Orthogroups/Orthogroups.GeneCount.tsv",
          g="Orthologues/Results_dir/Orthogroups/Orthogroups.tsv"
  container: "docker://davidemms/orthofinder"
  conda: "conda_envs/orthofinder.yaml"
  resources: slurm_extra="-C fat", runtime = lambda wildcards, attempt: attempt*60 *72, mem_mb=6400*config["threads"]
  shell:  """
              orthofinder -t {threads} -o Orthologues -n dir {params.p}
              cp {params.f} {output.f}
              cp {params.g} {output.g}
          """

# for bacteria total Ram memory/(#genomes*0.02) = a [ 128GBram/(428*0.02) = 14,95 --> a max = 14

if config["up_to_Eggnog_all_orthologs"]:
    rule selecte_seq_in_orth_group:
        input:i=config["orthology"]["ortholog_table"]
        output: "Annotations/Rep_seq_Orth_groups/report.txt"
        threads: 1 
        resources:
                runtime = lambda wildcards, attempt: attempt*60 *5, mem_mb=6400
        params:i=config["orthology"]["path_to_orthologs_sequences"],
               o="Annotations/Rep_seq_Orth_groups",
               c=config["Eggnog_all_orthologs_selection"]
        shell:  """
                   python support/Rep_seq_Orth_groups.py -i {params.i} -c {params.c} -o {params.o}
                """

    rule add_OG_in_header:
        input: "Annotations/Rep_seq_Orth_groups/report.txt"
        output: "Annotations/Rep_seq_Orth_groups/selected_seq_OGs.fasta"
        threads: 4
        resources:
              runtime = lambda wildcards, attempt: attempt*60 *5, mem_mb=6400*4
        params: i=config["workdir"]+"/Annotations/Rep_seq_Orth_groups/selection/"
        shell: 	"""
		              ls {params.i} | while read f
			            do
				             ort=$(echo $(basename $f) | cut -d"." -f1)
				             cat {params.i}$f | sed s/">"/">$ort."/ >> {output}
			            done
		        """

    rule eggNOG_all:
        input:  i="Annotations/Rep_seq_Orth_groups/selected_seq_OGs.fasta"
        output: "Annotations/Rep_seq_Orth_groups/All_EggNOG/selected.emapper.annotations"
        threads: config["threads"]
        params: db=config["path_to_eggnog_db"], o="Annotations/Rep_seq_Orth_groups/All_EggNOG"
        conda: "eggnog_mapper_env"
        resources:
            runtime = lambda wildcards, attempt: attempt*60 *72, mem_mb=6400*config["threads"]
        shell:  """
                    mkdir -p {params.o}

                    fileanot={params.o}/selected.emapper.annotations

                    if [[ ! -s $fileanot ]]; then
                          filehit={params.o}/selected.emapper.hits
                          if [[ -s $filehit ]]; then
                              echo "INFO: $filehit is present, annotation will be resumed"

                              emapper.py -i {input} -o selected --output_dir {params.o} --cpu {threads} --data_dir {params.db} -m diamond --resume
                          else
                              emapper.py -i {input} -o selected --output_dir {params.o} --cpu {threads} --data_dir {params.db} -m diamond
                          fi
                    else
                        echo "INFO: $fileanot already done"
                    fi

                """


rule converting_orthfiles:
  input: o=config["orthology"]["ortholog_table"], p=expand(config["Proteome_dir"]+"/{s}.faa", s = strains)
  output: o="phyloglm_input/Orthogroups.tsv", p="phyloglm_input/Annotations.txt"
  threads: 1
  resources:
           runtime = lambda wildcards, attempt: attempt*60 *5, mem_mb=6400
  shell: "bash support/converting.sh {input.o} {output.o} {output.p}"

if config["tree_using"] == "ani":
  rule comparative_with_ANI:
    input:  g=config["orthology"]["ortholog_count_table"], m="phyloglm_input/Annotations.txt",
            s="phyloglm_input/Orthogroups.tsv",n="fastANI_output/results", i=config["input"]["File_list"], k=expand("Annotations/GFF_files/{s}.gff", s= strains)
    output: c=config["output_dir"]+"/Annotations/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Candidates_enriched_orthologues.tsv",
            e=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Enriched.txt",
            e2=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Enriched_core.txt",
            dp=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Depleted.txt",
            dp2=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Depleted_core.txt",
            tr=config["output_dir"]+"/tree/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Cleaned_tree.txt"
    threads: config["threads"]
    params: b=config["PHYLOGLM"]["phyloglm_Bootnumber"], q=config["PHYLOGLM"]["p_adj_value_cutoff"],
            r=config["PHYLOGLM"]["orthologue_ratio_in_genome_dataset"],
            e=config["input"]["Special_group_name"],
            a=config["PHYLOGLM"]["phyloglm_btol_number"], l=config["PHYLOGLM"]["phyloglm_outfiles_prefix"],
            o=os.path.join(config["workdir"],config["output_dir"]), k="Annotations/GFF_files",
            g1=config["input"]["group_1_name"], l1=config["input"]["group_1_label"],
            g2=config["input"]["group_2_name"], l2=config["input"]["group_2_label"]
    container: "docker://lfdelzam/phylobotl_r_image"
    conda: "R_env"
    resources:
             runtime = lambda wildcards, attempt: attempt*60 *48, mem_mb=6400*config["threads"]
    shell:  """
              Rscript support/phylobotl.R -q {params.q} -l {params.l} -g {input.g} -N {input.n} \
              -b {params.b} -a {params.a} -o {params.o} -r {params.r} -s {input.s} -m {input.m} \
              -i {input.i} -k {params.k} -y {params.g1} -w {params.l1} \
              -z {params.g2} -v {params.l2} -e {params.e}
            """

else:
  rule comparative:
    input:  g=config["orthology"]["ortholog_count_table"], m="phyloglm_input/Annotations.txt",
            s="phyloglm_input/Orthogroups.tsv", t=config["tree_file"], i=config["input"]["File_list"], k=expand("Annotations/GFF_files/{s}.gff", s= strains)
    output: c=config["output_dir"]+"/Annotations/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Candidates_enriched_orthologues.tsv",
            e=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Enriched.txt",
            e2=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Enriched_core.txt",
            dp=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Depleted.txt",
            dp2=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Depleted_core.txt",
            tr=config["output_dir"]+"/tree/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Cleaned_tree.txt"
    threads: config["threads"]
    resources:
                 runtime = lambda wildcards, attempt: attempt*60 *48, mem_mb=6400*config["threads"]
    params: b=config["PHYLOGLM"]["phyloglm_Bootnumber"], q=config["PHYLOGLM"]["p_adj_value_cutoff"],
            r=config["PHYLOGLM"]["orthologue_ratio_in_genome_dataset"],
            e=config["input"]["Special_group_name"],
            a=config["PHYLOGLM"]["phyloglm_btol_number"], l=config["PHYLOGLM"]["phyloglm_outfiles_prefix"],
            o=os.path.join(config["workdir"],config["output_dir"]), k="Annotations/GFF_files",
            g1=config["input"]["group_1_name"], l1=config["input"]["group_1_label"],
            g2=config["input"]["group_2_name"], l2=config["input"]["group_2_label"]
    container: "docker://lfdelzam/phylobotl_r_image"
    conda: "R_env"
    shell:  """
              Rscript support/phylobotl.R -q {params.q} -l {params.l} -g {input.g} -t {input.t} -b {params.b} -a {params.a} \
              -o {params.o} -r {params.r} -s {input.s} -m {input.m} -i {input.i} -k {params.k} \
              -y {params.g1} -w {params.l1} -z {params.g2} -v {params.l2} -e {params.e}
            """

rule orthologues_gff:
    input: f=expand("Annotations/GFF_files/{s}.gff", s= strains), t=config["orthology"]["ortholog_table"]
    output: "Annotations/Orthologues/orthologues_gff.tsv"
    params: "Annotations/GFF_files"
    threads:1
    resources:
              runtime = lambda wildcards, attempt: attempt*60 *5, mem_mb=6400
    shell: "python support/parse_gff.py -i {params} -l {input.t} -o Annotations/Orthologues"

## Enriched section
rule vis_co:
    input: a="Annotations/Orthologues/orthologues_gff.tsv",
           e=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Enriched.txt",
           t=config["output_dir"]+"/tree/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Cleaned_tree.txt"
    output: config["output_dir"]+"/co_localization_figures/Leiden_clusters_using_max_loci_size_"+str(config["synteny_visualization"]["locus_size_window"])+"_orthologues.ft"
    params: m=config["synteny_visualization"]["present_in_at_least_n_genomes"],
            s=config["synteny_visualization"]["locus_size_window"],
            f=config["synteny_visualization"]["max_n_genomes_per_fig"],
            o=config["output_dir"]+"/co_localization_figures",
            l=" Enriched",i=config["input"]["File_list"]
    threads: config["threads"]
    resources:
             runtime = lambda wildcards, attempt: attempt*60 *24, mem_mb=6400*config["threads"]
    container: "docker://lfdelzam/phylobotl_r_image"
    conda: "R_env"
    shell:  """
                Rscript support/synteny_visual.R -a {input.a} -t {input.t} -e {input.e} -o {params.o} \
                -l {params.l} -m {params.m} -s {params.s} -f {params.f} -i {params.i}
            """

rule DNA_selected_orthologues:
  input: i=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Enriched.txt",
         g="phyloglm_input/Orthogroups.tsv", d=expand("Genes/{s}.ffn", s = strains)
  output: config["output_dir"]+"/Enriched_Orthologues_DNA_sequences/report.txt"
  params: o=config["output_dir"]+"/Enriched_Orthologues_DNA_sequences"
  threads: 1
  resources:
          runtime = lambda wildcards, attempt: attempt*60 *8, mem_mb=6400
  shell: "python support/DNA_selected_orthologues.py -i {input.i} -d Genes -g {input.g} -o {params.o}"


rule eggNOG_enriched:
    input:  i=config["orthology"]["ortholog_table"],
            e=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Enriched.txt"
    output: config["output_dir"]+"/Annotations/Enriched_KEGG/Kegg_annotations.tsv"
    threads: config["threads"]
    params: i=config["orthology"]["path_to_orthologs_sequences"], p="support/KEGG_Pathways_list.tsv", m="support/KEGG_Modules_list.tsv",
            o1=config["output_dir"]+"/Annotations/Enriched_eggNOG",
            o2=config["output_dir"]+"/Annotations/Enriched_KEGG",
            db=config["path_to_eggnog_db"]
    conda: "eggnog_mapper_env"
    resources:
               runtime = lambda wildcards, attempt: attempt*60 *40, mem_mb=6400*config["threads"]
    shell:  """
                mkdir -p {params.o1}
                grep -v '#' {input.e} | while read ort; do
                fileanot={params.o1}/$ort.emapper.annotations

                if [[ ! -s $fileanot ]]; then
                      filehit={params.o1}/$ort.emapper.hits
                      if [[ -s $filehit ]]; then
                          echo "INFO: $filehit is present, annotation will be resumed"
                          emapper.py -i {params.i}/$ort.fa -o $ort --output_dir {params.o1} --cpu {threads} --data_dir {params.db} -m diamond --resume
                      else
                          emapper.py -i {params.i}/$ort.fa -o $ort --output_dir {params.o1} --cpu {threads} --data_dir {params.db} -m diamond
                      fi
                else
                    echo "INFO: $fileanot already done"
                fi
                done
                python support/parse_eggNOGs.py -i {params.o1} -p {params.p} -m {params.m} -o {params.o2}

            """

rule eggNOG_loci:
    input:  e=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Enriched_core.txt",
            o=config["output_dir"]+"/Annotations/Enriched_KEGG/Kegg_annotations.tsv",
            l=config["output_dir"]+"/co_localization_figures/Leiden_clusters_using_max_loci_size_"+str(config["synteny_visualization"]["locus_size_window"])+"_orthologues.ft"
    output: config["output_dir"]+"/Annotations/Enriched_KEGG_core/Kegg_annotations.tsv"
    threads: config["threads"]
    params: p="support/KEGG_Pathways_list.tsv", #i2=config["output_dir"]+"/Annotations/eggNOG",
            m="support/KEGG_Modules_list.tsv", o2=config["output_dir"]+"/Annotations/Enriched_KEGG_core",
            c=config["output_dir"]+"/Annotations/Enriched_eggNOG_core",
            d=config["output_dir"]+"/Annotations/Loci",
            i=config["orthology"]["path_to_orthologs_sequences"],
            e=config["output_dir"]+"/Annotations/Enriched_eggNOG",
            db=config["path_to_eggnog_db"]
    conda: "eggnog_mapper_env"
    resources:
              runtime = lambda wildcards, attempt: attempt*60 *40, mem_mb=6400*config["threads"]
    shell:  """
                toprt="no"
                grep -v '#' {input.l} | while read line; do
                        if [ $(echo $line | grep -c '>') -eq 1  ]; then
                            if [ $toprt == "yes" ]; then
                                mkdir -p {params.d}/KEGG/$cl
                                python support/parse_eggNOGs.py -i $dirout -p {params.p} -m {params.m} -o {params.d}/KEGG/$cl
                            fi

                            cl=$(echo $line | sed s/'>'//)
                            dirout={params.d}/eggNOG/$cl
                            mkdir -p $dirout
                            toprt="yes"

                        else
                            fext={params.e}/$line.emapper.annotations
                            if [[ -s $fext ]]; then
                                cp $fext $dirout/.
                            fi

                            fileanot=$dirout/$line.emapper.annotations

                            if [[ ! -s $fileanot ]]; then
                                  filehit=$dirout/$line.emapper.hits
                                  if [[ -s $filehit ]]; then
                                      echo "INFO: $filehit is present, annotation will be resumed"
                                      emapper.py -i {params.i}/$ort.fa -o $ort --output_dir $dirout --cpu {threads} --data_dir {params.db} -m diamond --resume
                                  else
                                      emapper.py -i {params.i}/$ort.fa -o $ort --output_dir $dirout --cpu {threads} --data_dir {params.db} -m diamond
                                  fi
                            else
                                echo "INFO: $fileanot already done"
                            fi

                        fi
                    done

                        cl=$(grep -c '>' {input.l})
                        dirout={params.d}/eggNOG/$cl
                        mkdir -p {params.d}/KEGG/$cl
                        python support/parse_eggNOGs.py -i $dirout -p {params.p} -m {params.m} -o {params.d}/KEGG/$cl


                mkdir -p {params.c}
                grep -v '#' {input.e} | while read ort; do cp {params.e}/$ort.emapper.annotations {params.c}/ ; done
                python support/parse_eggNOGs.py -i {params.c} -p {params.p} -m {params.m} -o {params.o2}
                rm -r {params.c}
             """


rule enriched_gbk:
  input:  g=expand("Annotations/GBK_files/{s}.gbk", s= strains),
          a="Annotations/Orthologues/orthologues_gff.tsv",
          s=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Enriched.txt"
  output: config["output_dir"]+"/Annotations/Enriched_GBK_files/Genomes_without_selected_orthologs.txt"
  threads: 1
  resources:
            runtime = lambda wildcards, attempt: attempt*60 *5, mem_mb=6400
  params: i="Annotations/GBK_files",o=config["output_dir"]+"/Annotations/Enriched_GBK_files"
  shell:  "python support/parse_gbk.py -i {params.i} -e {input.a} -o {params.o} -s {input.s}"


rule vis_co_depleted:
    input: a="Annotations/Orthologues/orthologues_gff.tsv",
           e=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Depleted.txt",
           t=config["output_dir"]+"/tree/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Cleaned_tree.txt"
    output: config["output_dir"]+"/co_localization_figures_Depleted/Leiden_clusters_using_max_loci_size_"+str(config["synteny_visualization"]["locus_size_window"])+"_orthologues.ft"
    params: m=config["synteny_visualization"]["present_in_at_least_n_genomes"],
            s=config["synteny_visualization"]["locus_size_window"],
            f=config["synteny_visualization"]["max_n_genomes_per_fig"],
            o=config["output_dir"]+"/co_localization_figures_Depleted",
            l=" Depleted",i=config["input"]["File_list"]
    threads: config["threads"]
    resources:
              runtime = lambda wildcards, attempt: attempt*60 *48, mem_mb=6400*config["threads"]
    container: "docker://lfdelzam/phylobotl_r_image"
    conda: "R_env"
    shell:  """
                Rscript support/synteny_visual.R -a {input.a} -t {input.t} -e {input.e} -o {params.o} \
                -l {params.l} -m {params.m} -s {params.s} -f {params.f} -i {params.i}
            """

rule DNA_selected_orthologues_depleted:
  input: i=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Depleted.txt",
         g="phyloglm_input/Orthogroups.tsv", d=expand("Genes/{s}.ffn", s = strains)
  output: config["output_dir"]+"/Depleted_Orthologues_DNA_sequences/report.txt"
  params: o=config["output_dir"]+"/Depleted_Orthologues_DNA_sequences"
  threads: 1
  resources:
            runtime = lambda wildcards, attempt: attempt*60 *5, mem_mb=6400
  shell: "python support/DNA_selected_orthologues.py -i {input.i} -d Genes -g {input.g} -o {params.o}"


rule eggNOG_Depleted:
    input:  i=config["orthology"]["ortholog_table"],
            e=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Depleted.txt"
    output: config["output_dir"]+"/Annotations/Depleted_KEGG/Kegg_annotations.tsv"
    threads: config["threads"]
    params: i=config["orthology"]["path_to_orthologs_sequences"], p="support/KEGG_Pathways_list.tsv", m="support/KEGG_Modules_list.tsv",
            o1=config["output_dir"]+"/Annotations/Depleted_eggNOG",
            o2=config["output_dir"]+"/Annotations/Depleted_KEGG",
            db=config["path_to_eggnog_db"]
    resources:
              runtime = lambda wildcards, attempt: attempt*60 *48, mem_mb=6400*config["threads"]
    conda: "eggnog_mapper_env"
    shell:  """
                mkdir -p {params.o1}
                grep -v '#' {input.e} | while read ort; do
                fileanot={params.o1}/$ort.emapper.annotations

                if [[ ! -s $fileanot ]]; then
                      filehit={params.o1}/$ort.emapper.hits
                      if [[ -s $filehit ]]; then
                          echo "INFO: $filehit is present, annotation will be resumed"
                          emapper.py -i {params.i}/$ort.fa -o $ort --output_dir {params.o1} --cpu {threads} --data_dir {params.db} -m diamond --resume
                      else
                          emapper.py -i {params.i}/$ort.fa -o $ort --output_dir {params.o1} --cpu {threads} --data_dir {params.db} -m diamond
                      fi
                else
                    echo "INFO: $fileanot already done"
                fi
                done
                python support/parse_eggNOGs.py -i {params.o1} -p {params.p} -m {params.m} -o {params.o2}

            """

rule eggNOG_loci_Depleted:
    input:  e=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Depleted_core.txt",
            o=config["output_dir"]+"/Annotations/Depleted_KEGG/Kegg_annotations.tsv",
            l=config["output_dir"]+"/co_localization_figures_Depleted/Leiden_clusters_using_max_loci_size_"+str(config["synteny_visualization"]["locus_size_window"])+"_orthologues.ft"
    output: config["output_dir"]+"/Annotations/Depleted_KEGG_core/Kegg_annotations.tsv"
    threads: config["threads"]
    params: p="support/KEGG_Pathways_list.tsv", #i2=config["output_dir"]+"/Annotations/eggNOG",
            m="support/KEGG_Modules_list.tsv", o2=config["output_dir"]+"/Annotations/Depleted_KEGG_core",
            c=config["output_dir"]+"/Annotations/Depleted_eggNOG_core",
            d=config["output_dir"]+"/Annotations/Depleted_Loci",
            i=config["orthology"]["path_to_orthologs_sequences"],
            e=config["output_dir"]+"/Annotations/Depleted_eggNOG",
            db=config["path_to_eggnog_db"]
    conda: "eggnog_mapper_env"
    resources:
            runtime = lambda wildcards, attempt: attempt*60 *48, mem_mb=6400*config["threads"]
    shell:  """
                toprt="no"
                grep -v '#' {input.l} | while read line; do
                        if [ $(echo $line | grep -c '>') -eq 1  ]; then
                            if [ $toprt == "yes" ]; then
                                mkdir -p {params.d}/KEGG/$cl
                                python support/parse_eggNOGs.py -i $dirout -p {params.p} -m {params.m} -o {params.d}/KEGG/$cl
                            fi

                            cl=$(echo $line | sed s/'>'//)
                            dirout={params.d}/eggNOG/$cl
                            mkdir -p $dirout
                            toprt="yes"

                        else
                            fext={params.e}/$line.emapper.annotations
                            if [[ -s $fext ]]; then
                                cp $fext $dirout/.
                            fi

                            fileanot=$dirout/$line.emapper.annotations

                            if [[ ! -s $fileanot ]]; then
                                  filehit=$dirout/$line.emapper.hits
                                  if [[ -s $filehit ]]; then
                                      echo "INFO: $filehit is present, annotation will be resumed"
                                      emapper.py -i {params.i}/$ort.fa -o $ort --output_dir $dirout --cpu {threads} --data_dir {params.db} -m diamond --resume
                                  else
                                      emapper.py -i {params.i}/$ort.fa -o $ort --output_dir $dirout --cpu {threads} --data_dir {params.db} -m diamond
                                  fi
                            else
                                echo "INFO: $fileanot already done"
                            fi

                        fi
                    done

                        cl=$(grep -c '>' {input.l})
                        dirout={params.d}/eggNOG/$cl
                        mkdir -p {params.d}/KEGG/$cl
                        python support/parse_eggNOGs.py -i $dirout -p {params.p} -m {params.m} -o {params.d}/KEGG/$cl


                mkdir -p {params.c}
                grep -v '#' {input.e} | while read ort; do cp {params.e}/$ort.emapper.annotations {params.c}/ ; done
                python support/parse_eggNOGs.py -i {params.c} -p {params.p} -m {params.m} -o {params.o2}
                rm -r {params.c}
             """


rule Depleted_gbk:
  input:  g=expand("Annotations/GBK_files/{s}.gbk", s= strains),
          a="Annotations/Orthologues/orthologues_gff.tsv",
          s=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Depleted.txt"
  output: config["output_dir"]+"/Annotations/Depleted_GBK_files/Genomes_without_selected_orthologs.txt"
  threads: 1
  resources:
            runtime = lambda wildcards, attempt: attempt*60 *5, mem_mb=6400
  params: i="Annotations/GBK_files",o=config["output_dir"]+"/Annotations/Depleted_GBK_files"
  shell:  "python support/parse_gbk.py -i {params.i} -e {input.a} -o {params.o} -s {input.s}"

if config["genomad"]["include"]:
    rule genomad:
        input: "Genomes/{s}.fa"
        output: "GENOMAD/Genomad_output_{s}/{s}_summary.log"
        threads: config["threads"]
        params: db=config["genomad"]["path_to_genomad_db"], i="Genomes/{s}.fa", o="GENOMAD/Genomad_output_{s}", p=config["genomad"]["params_genomad"]
        conda: "genomad_env"
        resources:
                runtime = lambda wildcards, attempt: attempt*60 *48, mem_mb=6400*config["threads"]
        shell: "genomad end-to-end {params.p} -t {threads} {params.i} {params.o} {params.db}"

    rule parse_genomad:
        input: 	i=expand("GENOMAD/Genomad_output_{s}/{s}_summary.log",s = strains),
            f="Annotations/Orthologues/orthologues_gff.tsv",
            e=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Enriched.txt"
        output: "GENOMAD/Summary_genomad/Genomes_with_plasmid.tsv"
        threads: 1
        resources:
                    runtime = lambda wildcards, attempt: attempt*60 *5, mem_mb=6400
        params: i="GENOMAD", o="GENOMAD/Summary_genomad", c=config["genomad"]["Genomad_score_cut_off"]
        shell: "python support/parse_genomad.py -i {params.i} -f {input.f} -e {input.e} -o {params.o} -c {params.c}"

if config["tax_specific_protein_annotation"]["status"]:
    rule specif_annot:
        input: d=config["tax_specific_protein_annotation"]["path_to_database"],
               e=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Enriched.txt"
        output: config["output_dir"]+"/Annotations/Specific_db/Enriched_orthologs_best_annotation.tsv",
                config["output_dir"]+"/Annotations/Specific_db/Enriched_orthologs_most_frequent_annotation.tsv"
        params: p="--max-accept 1 --start-sens 4 --sens-steps 3 -s 7 -v 1 --min-seq-id 0.9",
                s=config["orthology"]["path_to_orthologs_sequences"],
		o=config["output_dir"]+"/Annotations/Specific_db",
                d=config["output_dir"]+"/Annotations/Specific_db/Hits_tables"
        threads: config["threads"]
        resources:
                 runtime = lambda wildcards, attempt: attempt*60 *48, mem_mb=6400*config["threads"]
        conda: "conda_envs/mmseqs.yaml"
        shell:  """
		    mkdir -p {params.d}
                    cat {input.e} | while read ort
                    do
			if [[ $ort =~ ^[^#] ]]; then
                          filein={params.s}/$ort.fa
                          fileout={params.d}/$ort.hits.txt
                          mmseqs easy-search $filein {input.d} $fileout tmpo1_dir --threads {threads} {params.p}
                          rm -r tmpo1_dir
			fi
                    done
                    python support/parse_mmseqs.py -i {params.d} -r {input.d} -o {params.o} -p Enriched
                """

    rule specif_annot_clusters:
        input:  i=config["output_dir"]+"/Annotations/Specific_db/Enriched_orthologs_most_frequent_annotation.tsv",
                c=config["output_dir"]+"/co_localization_figures/Leiden_clusters_using_max_loci_size_"+str(config["synteny_visualization"]["locus_size_window"])+"_orthologues.ft"
        output: config["output_dir"]+"/Annotations/Loci/Specific_db/cluster_most_frequent_annotations.tsv"
        params: o=config["output_dir"]+"/Annotations/Loci/Specific_db/"
        threads: 1
        resources:
                  runtime = lambda wildcards, attempt: attempt*60 *5, mem_mb=6400
        shell:  """
                    python support/parse_Uniprot.py -i {input.i} -c {input.c} -o {params.o}
                """

    rule specif_annot_dep:
        input: d=config["tax_specific_protein_annotation"]["path_to_database"],
               e=config["output_dir"]+"/"+config["PHYLOGLM"]["phyloglm_outfiles_prefix"]+"_Depleted.txt"
        output: config["output_dir"]+"/Annotations/Specific_db/Depleted/Depleted_orthologs_best_annotation.tsv",
                config["output_dir"]+"/Annotations/Specific_db/Depleted/Depleted_orthologs_most_frequent_annotation.tsv"
        params: p="--max-accept 1 --start-sens 4 --sens-steps 3 -s 7 -v 1 --min-seq-id 0.9",
                s=config["orthology"]["path_to_orthologs_sequences"],
                o=config["output_dir"]+"/Annotations/Specific_db/Depleted",
                d=config["output_dir"]+"/Annotations/Specific_db/Depleted/Hits_tables"
        threads: config["threads"]
        resources:
                  runtime = lambda wildcards, attempt: attempt*60 *48, mem_mb=6400*config["threads"]
        conda: "conda_envs/mmseqs.yaml"
        shell:  """
                    mkdir -p {params.d}
		    cat {input.e} | while read ort
                    do
                      	if [[ $ort =~ ^[^#] ]]; then
                          filein={params.s}/$ort.fa
                          fileout={params.d}/$ort.hits.txt
                          mmseqs easy-search $filein {input.d} $fileout tmpo2_dir --threads {threads} {params.p}
                          rm -r tmpo2_dir
                        fi
                    done
                    python support/parse_mmseqs.py -i {params.d} -r {input.d} -o {params.o} -p Depleted
                """

