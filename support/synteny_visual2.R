rm(list = ls())
options(warn=-1)

list_of_packages <- c("ape","tidyverse","gggenomes", "igraph","tidytree","gridExtra","readr","ggtree", "ggnewscale", "randomcoloR",
                      "future.apply","foreach","doParallel","dplyr", "furrr", "future")

suppressPackageStartupMessages(library(argparser))
##something is not working with the fugires (some of the are empties)
# arguments
p <- arg_parser("synteny_visual.R")
p <- add_argument(p, "-a", help="Path to orthologues_gff.tsv", default="orthologues_gff.tsv")
p <- add_argument(p, "-m", help="Minimum Weigth in graph to be included", default=4)
p <- add_argument(p, "-t", help="Path to cleaned tree, .txt", default="Borrelo_407_with_alpha_phyloglm_and_logistf/tree/Vv_Cleaned_tree.txt") #must bee the cleaned one
p <- add_argument(p, "-s", help="loci size to look for co-location", default=40000) 
p <- add_argument(p, "-f", help="max number of genomes per figure", default=22)
p <- add_argument(p, "-o", help="output directory", default="Borrelo_407_with_alpha_phyloglm_and_logistf/synteny")
p <- add_argument(p, "-l", help="Label of Selected Orthologues in figures", default=" Enriched")
p <- add_argument(p, "-e", help="Path to enriched/depleted orthologues list", default="Borrelo_407_with_alpha_phyloglm_and_logistf/Vv_Enriched.txt")
p <- add_argument(p, "-i", help="genomes metadata", default="GENOME_LIST")
p <- add_argument(p, "-c", help="Number of CPUs", default=8, type="integer")
argv <- parse_args(p)

#libraries
for (s in list_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE))}

#functions

visualizing_ortholgues_in_cluster <- function(dfP, clster, prefixT, tree) {

  SIZE <- 3
  palette <- distinctColorPalette(length(unique(dfP$Cluster))-1)
  
  p1<-ggtree(tree, layout="fan", open.angle = 0, ladderize=TRUE, branch.length="branch.length") +
    geom_tiplab(size=1.4, align=T, linesize=.05, offset=0.003)+
    geom_text2(aes(subset = !isTip, label=label), size=SIZE,angle = 0,vjust = "inward", hjust = "inward",check_overlap = TRUE)+
    geom_nodepoint(color="black", alpha=1/3, size=0.5 )+ geom_rootedge(rootedge = 0.005) 
  
  NA_vector=rep(NA,(length(p1[["data"]][["label"]]))-length(dfP$Cluster))
  new_labels=c(gsub("Absent", "", dfP$Cluster), NA_vector)
  Pa <- p1 +geom_tiplab2(size=1.5, align=T,geom="text",linesize=0.0, offset=0.0021,aes(label=new_labels))
  
  P2<-gheatmap(Pa, dfP[,2, drop=FALSE], offset=0.001, width=0.02, colnames=F) +
    scale_fill_manual(values = c(palette,"white"),
                      name=paste0("Cluster ",clster, " subgroups")) +
    theme(legend.position="bottom", legend.box = "vertical",hjust = 0.5) +
    guides(fill = guide_legend(ncol = 22, title.position="top", title.hjust = 0.5,override.aes = list(color = "black")))

  
  P3 <- P2 + new_scale_fill() 
  
  P4<-gheatmap(P3, dfP[,1, drop=FALSE] , offset=0, width=0.02, colnames=F) +scale_fill_manual(
    values = c("black","blue"), name="Isolate")+
    theme(legend.position="bottom", legend.box = "vertical",hjust = 0.5) +
    guides(fill = guide_legend(ncol = 2, title.position="top", title.hjust = 0.5,override.aes = list(color = "black")))

  
  ggsave(paste(prefixT,paste("cluster",clster, "tree_fan_ladderized_branchlength.pdf", sep="_"), sep="/"), P4, width = 42, height = 60, units = "cm", device = "pdf")
}

region_in_genome <- function(geno, tDF, O1) {

  # Filter tDF for the given genome once, instead of repeatedly
  geno_tDF <- tDF %>% filter(bin_id == geno) #%>% arrange(seq_id, start) #new
  
  # Get the unique contigs for the genome
  contig_where_loci <- unique(geno_tDF$seq_id)
  
  # Filter O1 once for the given genome
  Orts_IN_genome <- O1 %>% filter(bin_id == geno)
  
  # Initialize an empty list to collect results 
  result_list <- vector("list", length(contig_where_loci))
  
  for (i in seq_along(contig_where_loci)) {
    c <- contig_where_loci[i]
    
    # Filter for the current contig once
    contig_tDF <- geno_tDF %>% filter(seq_id == c)
    
    # Calculate the min and max positions in one go
    start_p <- min(contig_tDF$start)
    end_p <- max(contig_tDF$end)
    
    # Filter Orts_IN_genome for the current contig and relevant start/end positions
    Orts_in_contig <- Orts_IN_genome %>% 
      filter(seq_id == c, start >= start_p, end <= end_p) #%>% arrange(seq_id, start) #new
    
    # Store the result in the list
    result_list[[i]] <- Orts_in_contig
  }
  
  # Combine the list into a single data frame at the end
  outdf <- bind_rows(result_list)
  
  return(outdf)
}

plot_synteny<-function(c1,loci1,O1, texto, pout) {

  if (length(loci1) >0) {
    
    cat("Orthologues in cluster", c1, ":", loci1, "\n", file = paste(argv$o,paste("Leiden_clusters_using_max_loci_size",max_loci,"orthologues.txt" ,sep="_"), sep="/"),
        append = T)
    if (pout) {
      cat(">",c1,"\n", file = paste(argv$o,paste("Leiden_clusters_using_max_loci_size",max_loci,"orthologues.ft",sep="_"), sep="/"),
          append = T, sep="")
      for (k in loci1) {
        cat(k,"\n", file = paste(argv$o,paste("Leiden_clusters_using_max_loci_size",max_loci,"orthologues.ft",sep="_"), sep="/"),
            append = T)
      }
    }
    
    # Filter for loci1 in feat_id and select necessary columns
    O2d <- O1 %>% filter(feat_id %in% loci1) #%>% select(bin_id, feat_id, strand, seq_id)
    
    # Extract unique bin_id as a vector
    list_of_genomesd <- O2d %>% distinct(bin_id) %>% pull(bin_id)
    
    # Initialize references and similars list
    references <- vector("character", length(list_of_genomesd))  # preallocate memory
    similars <- vector("list", length(list_of_genomesd))
    names(similars) <- list_of_genomesd
    
    # Pre-compute regions for all genomes to avoid redundant calls

    # Set up parallel processing using all available cores
    plan(multisession, workers = argv$c)  # Or 'multicore' for Linux/macOS
    
    # Parallelize the lapply using future_lapply
    genome_regions <- future_lapply(list_of_genomesd, function(genome) {
      region_in_genome(genome, O2d, O1)
    })
    plan(sequential)
    # Name the list elements for easier reference
    names(genome_regions) <- list_of_genomesd
#select ref the one with less contigs
    k <- 1
    while (length(list_of_genomesd) > 0) {
      # Use precomputed region for the first genome in the list
      aver1d <- genome_regions[[list_of_genomesd[1]]]
      order1 <- as.vector(aver1d$feat_id)
      strand1 <- ifelse(aver1d$strand == "+", 1, -1)
      rf=as.vector(list_of_genomesd[1])
      c_rf=n_distinct(as.factor(aver1d$seq_id))
      
      cluster1 <- c(rf,rep("NA", length(list_of_genomesd)-1))
      # Only enter the loop if more than 1 genome is left
      if (length(list_of_genomesd) > 1) {

        # Vectorized checking for similarities across genomes
        for (i in 2:length(list_of_genomesd)) {#i=2
          current_genome <- list_of_genomesd[i]
          if (!current_genome %in% cluster1 && !is.na(current_genome)) {
            gr2 <- genome_regions[[as.character(current_genome)]]
            order2 <- as.vector(gr2$feat_id)
            strand2 <- ifelse(gr2$strand == "+", 1, -1)
            
            # Check for equivalence or reverse equivalence in one condition
            if ((setequal(order1, order2) && setequal(strand1, strand2)) ||
                (setequal(rev(order1), order2) && setequal(rev(strand1 * (-1)), strand2))) {
              cluster1[i] <- as.vector(current_genome)
              p_c_rf=n_distinct(as.factor(gr2$seq_id))
              if (p_c_rf < c_rf) {rf=as.vector(current_genome)
              c_rf = p_c_rf} #reference, the one with least contigs
            }
          }
          
        }
       
      }

      # Clean up the cluster and update references
      cluster1 <- cluster1[cluster1 != "NA"]
      references[k] <- rf
      similars[[rf]] <- cluster1
      list_of_genomesd <- setdiff(list_of_genomesd, cluster1)
      k <- k + 1

    }

    # Filter non-empty references and similars
    references <- references[references != ""]
    similars <- similars[lapply(similars, length) > 0]
    
    sink(paste(argv$o,paste("group_of_genomes",texto,"in_figure.txt", sep="_"),sep="/"))
    print(similars)
    sink()
    
    # Pre-compute phyloorder and list_of_bins
    phyloorder <- names(similars[order(sapply(similars, length), decreasing = TRUE)])
    list_of_bins <- vector("character", length(phyloorder))
    
    # Initialize lists to store results, to be converted to data frames later
    fdfd_list <- vector("list", length(phyloorder))
    emale_nagF_list <- vector("list", length(phyloorder))
    
    # Loop over phyloorder
    for (j in seq_along(phyloorder)) {
      g <- phyloorder[j]
      ngeno <- length(similars[[g]])
      
      # Process tempdf
      tempdf <- region_in_genome(g, O2d, O1)
      tempdf$bin_id <- paste(tempdf$bin_id, " (", ngeno, ")", sep = "")
      fdfd_list[[j]] <- tempdf  # Store in list
      
      # Process temnag
      temnag <- emale_nag %>% filter(bin_id == g)
      temnag$bin_id <- paste(temnag$bin_id, " (", ngeno, ")", sep = "")
      emale_nagF_list[[j]] <- temnag  # Store in list
      
      # Store bin_id
      list_of_bins[j] <- paste(g, " (", ngeno, ")", sep = "")
    }
    
    # Combine lists into data frames
    fdfd <- do.call(rbind, fdfd_list)
    emale_nagF <- do.call(rbind, emale_nagF_list)
    
    plot_gggenomes(fdfd, emale_nagF, references, phyloorder, texto)
  }  
  return(similars)
}

plot_gggenomes <- function(fdfd, emale_nagF, references, phyloorder, texto) {
  
  fdfd <- fdfd %>% mutate(
      feat_id = factor(feat_id, levels = unique(feat_id)),
      strand = factor(strand, levels = unique(strand)),
      bin_id = factor(bin_id, levels = unique(bin_id))
    )

  emale_nagF <- emale_nagF %>% mutate(bin_id = factor(bin_id, levels = unique(bin_id)))
  
  if (length(references) > ste) {
    f_numb <- 1
    for (counter in seq(1, length(references), ste)) { #counter=1
      end_v <- min(counter + ste - 1, length(references))
      flocsub <- gggenomes(genes = fdfd[gsub(" (.*.)","", fdfd$bin_id) %in% phyloorder[counter:end_v],], spacing = 0.05, theme = "clean") %>% 
        add_feats(ngaros = emale_nagF) +
        geom_seq() + geom_bin_label(size = 3, stat = "identity") + 
        geom_gene(aes(fill = feat_id), position = "strand", stat = "identity", na.rm = TRUE, show.legend = FALSE) + 
        geom_feat_tag(aes(label = feat_id, color = type), key_glyph = "rect", angle = 20, nudge_y = 0.2, size = 2, check_overlap = FALSE) +
        scale_color_manual(values = c("darkblue", "grey", "red", "black", "purple", "darkgreen", "orange", "brown", "blue", "yellow"), name = "Orthologue label type") +
        theme(legend.position = "bottom", legend.key.size = unit(0.2, 'cm'))
      
      ggsave(paste(argv$o, paste(texto, f_numb, "by_ref_genomes", "figure.pdf", sep = "_"), sep = "/"), flocsub, width = 32, height = 60, units = "cm", device = "pdf")
      f_numb <- f_numb + 1
    }
  } else {
    floc <- gggenomes(genes = fdfd, spacing = 0.05, theme = "clean") %>% 
      add_feats(ngaros = emale_nagF) +
      geom_seq() + geom_bin_label(size = 3, stat = "identity") + 
      geom_gene(aes(fill = feat_id), position = "strand", stat = "identity", na.rm = TRUE, show.legend = FALSE) + 
      geom_feat_tag(aes(label = feat_id, color = type),key_glyph = "rect", angle = 20, nudge_y = 0.2, size = 2, check_overlap = FALSE) +
      scale_color_manual(values = c("darkblue", "grey", "red", "black", "purple", "darkgreen", "orange", "brown", "blue", "yellow"), name = "Orthologue label type") +
      theme(legend.position = "bottom", legend.key.size = unit(0.2, 'cm'))
    
    ggsave(paste(argv$o, paste(texto, "by_ref_genomes", "figure.pdf", sep = "_"), sep = "/"), floc, width = 32, height = 60, units = "cm", device = "pdf")
  }
}

find_co_locations <- function(genomesDF, max_loci = 40000, num_cores = parallel::detectCores() - 1) {
  cat("INFO: Finding synteny in genomes\n") 
  # Set up parallel backend
  plan(multisession, workers = num_cores)
  
  # Parallelize the genome processing
  co_locations <- future_lapply(unique(genomesDF$genome), function(g) {
    sdf_g <- genomesDF %>% filter(genome == g)
    co_loc_list <- list()
    bloq <- 1
    
    for (c in unique(sdf_g$Contig)) {
      
      sdf_c <- sdf_g %>% filter(Contig == c)

      if (length(sdf_c$Orthologue) > 1) {
        i <- 1
        while (i < length(sdf_c$Orthologue)) {
          o <- sdf_c$Orthologue[i]
          end1 <- max(sdf_c$end[sdf_c$Orthologue == o])
          sdf_bloq <- sdf_c %>% filter(end >= end1 & end <= max_loci + end1)
          
          if (length(sdf_bloq$Orthologue) > 1) {
            co_loc_list[[bloq]] <- unique(sdf_bloq$Orthologue)
            bloq <- bloq + 1
          }
          i <- i + 1
        } # close while
      } # close if
    } # close contig loop
    
    return(co_loc_list)
  })
  
  # Flatten the list of co-locations and remove empty elements
  co_locations <- do.call(c, co_locations)
  co_locations <- co_locations[lapply(co_locations, length) > 0]
  
  # Reset the parallel plan
  plan(sequential)
  
  return(co_locations)
}

create_network_table <- function(co_locations, num_cores = parallel::detectCores() - 1) {
  
  # Set up parallel backend
  plan(multisession, workers = num_cores)
  
  # Parallelize the combination creation for each co_location
  co_list <- future_lapply(co_locations, function(loc) {
    loc_length <- length(loc)
    
    # Generate combinations of orthologues if more than one orthologue is present
    if (loc_length > 1) {
      combs <- combn(loc, 2)
      data.frame(Ortholog1 = combs[1,], Ortholog2 = combs[2,], stringsAsFactors = FALSE)
    } else {
      NULL
    }
  })
  
  # Combine all data frames into one and remove NULL elements
  co_tables <- do.call(rbind, co_list)
  
  # Reset the parallel plan
  plan(sequential)
  
  return(co_tables)
}

create_graph_input <- function(co_tabl, m_threshold) {
  cat("INFO: Generating graph input\n")
  # Summarize counts for each unique pair of Ortholog1 and Ortholog2
  co_table_counts <- co_tabl %>%
    group_by(Ortholog1, Ortholog2) %>%
    summarize(count_coloc = n(), .groups = 'drop') %>%
    filter(count_coloc > m_threshold) %>%   # Filter by the threshold
    rename(weigth = count_coloc)            # Rename for consistency
  
  return(as.data.frame(co_table_counts))
}

#Reading info
list_enriched=read.csv(argv$e)
names(list_enriched)="Enriched"
list_enriched$Enriched=gsub(" ","",list_enriched$Enriched)

genomeDFWall <- read_tsv(argv$a, col_names = TRUE, show_col_types = FALSE, num_threads = argv$c) 
# Filter, rename, and drop columns in one step
genomeDFall <- genomeDFWall %>%
  filter(Orthologue %in% list_enriched$Enriched) %>%
  rename(genome = Genome, strand = Direction, start = Start, end = Stop, attributes = Description) %>%
  select(-Type)

dir.create(argv$o, showWarnings = FALSE)
#### Find co-location of enriched orthologues
max_loci=argv$s
if (!file.exists(paste(argv$o,paste(max_loci,"graph.rds", sep="_"), sep="/"))) {

  co_location <- find_co_locations(genomeDFall,num_cores = argv$c)
  cat("INFO: Creating network input table\n")
  co_table <- create_network_table(co_location,num_cores = argv$c)
  co_table_count <- create_graph_input(co_table,argv$m)

  cat("INFO: Generating graph and Leiden clustering\n")
  gra=graph_from_data_frame(co_table_count, directed = FALSE)
  saveRDS(gra, file = paste(argv$o,paste(max_loci,"graph.rds", sep="_"), sep="/"))
} else {
gra=readRDS(paste(argv$o,paste(max_loci,"graph.rds", sep="_"), sep="/"))
}
set.seed(1)
cls_leiden=cluster_leiden(
  gra,objective_function="modularity", # objective_function = c("CPM", "modularity"), #  weights = NULL,
  resolution_parameter = 1, beta = 0.01, n_iterations = 100, # vertex_weights = NULL, #  initial_membership = NULL
)


pdf(paste(argv$o,paste("Leiden_clusters_using_max_loci_size",max_loci, sep="_"), sep="/"))
plot(cls_leiden, gra,layout=layout_nicely, 
     vertex.label.cex = 0.6, 
     vertex.size=8,vertex.label.dist=0.5, edge.arrow.size=0.5, rescale=TRUE,
     main = paste(argv$l,"Orthologues"), sub= "Co-localization clusters")

dev.off()
####
for (c in 1:cls_leiden$nb_clusters) {
  if (length(cls_leiden[[c]]) > 9) {
          sg = induced_subgraph(gra, cls_leiden[[c]])
          sub=cluster_leiden(sg,objective_function="modularity")
          pdf(paste(argv$o,paste("Leiden_clusters_using_max_loci_size",max_loci,"subgraph_cluster",c, sep="_"), sep="/"))
          plot(sub,sg, cex=0.6)
          dev.off()
  }  
}
###
tree=read.tree(argv$t)
Y=as_tibble(tree)

emale_genes_all <- genomeDFWall %>% rename(seq_id=Contig, bin_id=Genome, feat_id=Orthologue, strand=Direction,
                                           start=Start, end=Stop, attributes=Description, type=Type)
emale_nag <- emale_genes_all %>%
  mutate(type = ifelse(feat_id %in% unique(genomeDFall$Orthologue), argv$l, type))

# Calculate the number of unique 'bin_id' values once
ngen <- n_distinct(emale_genes_all$bin_id)
# Store the value from argv
ste <- argv$f
##

# INFO: Reading genome/sample information
cat("INFO: Reading genome/sample information \n")
genome_table <- read.csv(argv$i, sep = ",", na.strings = "NA", header = FALSE, col.names = c("Strain", "File", "Type", "Special")) %>% mutate_all(~ gsub("\\s+", " ", .))

# Clinical Isolates names
PATHOGENIC_STRAINS <- genome_table$Strain[genome_table$Type == "C"]

# Vectorized approach to assign 1 or 0 for Clinical/Environmental isolates
y <- as.integer(tree$tip.label %in% PATHOGENIC_STRAINS)

# Check if all pathogenic strains have matching labels
if (sum(y) != length(PATHOGENIC_STRAINS)) {
  stop("Strain names and tree labels don't match, please check spelling in GENOME_LIST file")
}

# Map 1/0 to Clinical/Environmental types
tipo <- ifelse(y == 1, "Clinical", "Environmental")

# Create the data frame
DFP <- data.frame(Isolates = as.character(tipo), row.names = tree$tip.label)

cat("INFO: Co-location - Output per similar genomes\n")

# Store emale_genes_all$bin_id in a variable
bin_ids <- emale_genes_all$bin_id
# Use filter and mutate from dplyr for efficiency
O1 <- emale_genes_all %>%
  filter(bin_ids %in% Y$label[1:ngen]) %>%
  mutate(
    feat_id = factor(feat_id, levels = unique(feat_id)),
    strand = factor(strand, levels = unique(strand)),
    bin_id = factor(bin_id, levels = unique(bin_id))
  )

vis_outdir=paste(argv$o, "Tree_clusters", sep="/")
dir.create(vis_outdir)

dirgff=paste(argv$o, "GFF", sep="/")
dir.create(dirgff, showWarnings = FALSE)

GFF_DF <- genomeDFWall %>% rename(genome = Genome, strand = Direction, start = Start, end = Stop, attributes = Description)

for (c in 1:cls_leiden$nb_clusters) { #c=3
  cat("INFO: Cluster -",c,"\n")
  loci=cls_leiden[[c]]
  cluster_group <-plot_synteny(c,loci,O1,paste("grouped_cluster_leiden", c, sep="_"), pout=TRUE)
  DFP_c=DFP
  DFP_c$Cluster="Absent"
  for (z in 1:length(cluster_group)  ) {
    DFP_c$Cluster[ row.names(DFP_c) %in% cluster_group[[z]] ]=z
  }
  #Exporting loci gff
  Gff_cluster <- GFF_DF %>% filter(Orthologue %in% loci)
  fileout=paste(dirgff,paste0("cluster",c,".gff"),sep="/")
  write_tsv(Gff_cluster,fileout)
  
  cat("INFO: Visualising orthologues in Cluster -",c,"\n")
  visualizing_ortholgues_in_cluster(DFP_c, c, vis_outdir,tree) 
}

cat("INFO: co-location - all genomes\n")

dirall=paste(argv$o, "All_genomes", sep="/")
dir.create(dirall, showWarnings = FALSE)
#ploting all the genomes
for (c in 1:cls_leiden$nb_clusters) {
  loci=cls_leiden[[c]]
  if (length(loci) >1){ 
    texto=paste("cluster_leiden", c, sep="_")
    fnumb=1
    for (counter in seq(1,ngen,ste)) { #counter=1
      end_v=(counter+ste-1)
      if (end_v > ngen) {end_v = ngen}
      
      # Store emale_genes_all$bin_id and O1$feat_id in variables to avoid repeated lookups
      O1 <- emale_genes_all %>%
        filter(bin_id %in% Y$label[counter:end_v]) %>%
        mutate(feat_id = factor(feat_id, levels = unique(feat_id)))
      
      # Filter O1 for loci in feat_id
      O2 <- O1 %>% filter(feat_id %in% loci)
      
      if (nrow(O2)>0) {
        phyloorder=Y$label[counter:end_v]

        # Pre-allocate a list to store results
        fdfa_list <- vector("list", length(phyloorder))
        # Loop through phyloorder, store results in the list
        for (i in seq_along(phyloorder)) {
          fdfa_list[[i]] <- region_in_genome(phyloorder[i], O2, O1)
        }
        # Combine all list elements into a single data frame
        fdfa <- bind_rows(fdfa_list)
        # Convert columns to factors in a single mutate step
        fdfa <- fdfa %>%
          mutate(
            bin_id = factor(bin_id, levels = unique(bin_id)),
            feat_id = factor(feat_id, levels = unique(feat_id)),
            strand = factor(strand, levels = unique(strand))
          )

        floca=gggenomes(genes=fdfa,spacing = 0.05, theme ="clean") %>% 
          add_feats(ngaros=emale_nag)+ 
          geom_seq() +  geom_bin_label(size=3,stat="identity") + 
          geom_gene(aes(fill=feat_id),position="strand",stat="identity",  na.rm=TRUE, show.legend = F) + 
          geom_feat_tag(aes(label=feat_id, color=type), key_glyph = "rect",angle = 20, nudge_y=0.2, check_overlap = TRUE, size=2, na.rm=TRUE)+
          scale_color_manual(values=c("darkblue","grey", "red", "black", "purple"),
                             name="Orthologue label type")+ 
          theme(legend.position="bottom", legend.key.size = unit(1, 'cm')) 
        
        ggsave(paste(dirall,paste("genomes",texto,fnumb,"figure.pdf", sep="_"),sep="/"), floca, width = 32, height = 60, units = "cm", device = "pdf")
        fnumb=fnumb+1
      }      
      
    } # close  for (counter in seq(1,ngen,ste))
    
  }  
}


###########
#Extract gff sections of selected clusters
#genomeDFWall <- read_tsv(argv$a, col_names = TRUE, show_col_types = FALSE, num_threads = argv$c) 
# Filter, rename, and drop columns in one step
#3,4,13 on chr1 and clusters 9,5,12,11,10 

#dirgff=paste(argv$o, "GFF", sep="/")
#dir.create(dirgff, showWarnings = FALSE)

#for (c in 1:cls_leiden$nb_clusters) {

#loci=c(cls_leiden[[5]],cls_leiden[[9]],cls_leiden[[10]],cls_leiden[[11]],cls_leiden[[12]])  
#loci=cls_leiden[[c]]

#Gff_cluster <- genomeDFWall %>%
#  filter(Orthologue %in% loci) %>%
#  rename(genome = Genome, strand = Direction, start = Start, end = Stop, attributes = Description) #%>% select(-Type)
#fileout=paste(dirgff,paste0("cluster",c,".gff"),sep="/")
#fileout=paste(dirgff,paste0("clusters5-9-10-11-12.gff"),sep="/")
#write_tsv(Gff_cluster,fileout)
#}
