# Clear workspace and suppress startup messages
rm(list = ls())
options(warn=-1)
suppressPackageStartupMessages(library(argparser))

# arguments
p <- arg_parser("phylocomgen")
p <- add_argument(p, "-t", help="phylogenetic tree file", default="/Users/luisdelgado/Documents/Phyloglm_folder/Results_407_with_IQTREE_rooted/tree/Vv_Cleaned_tree.txt")
p <- add_argument(p, "-g", help="gene/orthologue count file", default="Orthogroups.GeneCount.tsv")
p <- add_argument(p, "-i", help="genomes metadata", default="GENOME_LIST")
p <- add_argument(p, "-e", help="Special group", default="Baltic Sea")
p <- add_argument(p, "-s", help="Orthogroups list file", default="phyloglm_input/Orthogroups.tsv")
p <- add_argument(p, "-m", help="list of protein annotation from prokka file", default="phyloglm_input/Annotations.txt")
p <- add_argument(p, "-l", help="outfiles prefix", default="Vv")
p <- add_argument(p, "-o", help="output directory", default="Borrelo_407_with_alpha_phyloglm_and_logistf")
p <- add_argument(p, "-r", help="FILTER: ratio of orthologue presence/absence in the genome dataset,
                  if higher than this value, the orthologues is not considered in the analysis", default=0.95)
p <- add_argument(p, "-b", help="phyloglm Bootnumber", default=0)
p <- add_argument(p, "-a", help="phyloglm btol number", default=10)
p <- add_argument(p, "-q", help="p.adj.value cutoff", default=0.05)
p <- add_argument(p, "-c", help="number of clusters", default=2)
p <- add_argument(p, "-N", help="FastANI output file", default="")
p <- add_argument(p, "-p", help="number of cpus", default=8, type="integer")
p <- add_argument(p, "-k", help="Path to gff files", default="Annotations/GFF_files")
p <- add_argument(p, "-w", help="Group 1 name", default="Clinical")
p <- add_argument(p, "-y", help="Group 1 label", default="C")
p <- add_argument(p, "-z", help="Group 2 name", default="Environmental")
p <- add_argument(p, "-d", help="Group 2 label", default="E")

argv <- parse_args(p)

##extra parameter - on development
core_fraction=0.97

#libraries
list_of_packages <- c("dplyr","ggtree","ggplot2", "phylolm", "factoextra", "cluster",
                      "FactoMineR","umap","tidyr","ape","phangorn","pheatmap",
                      "vegan","ggrepel","tidyverse", "data.table", "RColorBrewer",
                      "logistf", "ips", "gridExtra", "readr", "parallel", "foreach", "doParallel","utils")

# Load libraries
invisible(lapply(list_of_packages, function(pkg) suppressPackageStartupMessages(library(pkg, character.only = TRUE))))


#functions

clustersout<-function(clust, vec, list_pathog_str) {
  cat("\n ---------------- \n")
  cat("******** Cluster", clust, "\n")
  PAM_patho=names(vec$clustering[vec$clustering == clust])
  unknown_in_cluster=rep("Na", length(PAM_patho))
  patho_in_cluster=rep("Na", length(PAM_patho))
  cat("Total strains in cluster",clust, ":" ,length(PAM_patho), "\n")
  for (w in PAM_patho) {
    if (w %in% list_pathog_str) {
      patho_in_cluster=c(patho_in_cluster,w)
    } else {
      unknown_in_cluster=c(unknown_in_cluster,w)
    }
  }
  patho_in_cluster=patho_in_cluster[patho_in_cluster != "Na"]
  unknown_in_cluster=unknown_in_cluster[unknown_in_cluster != "Na"]
  
  cp=length(patho_in_cluster)
  Np=length(unknown_in_cluster)
  cat("Known Pathogenic:", cp, "->", cp*100/length(list_pathog_str), "% of the known pathogenic strains", "\n List:", patho_in_cluster,"\n")
  cat("Unknown:", Np, "\n List ",unknown_in_cluster,"\n")
}

print_out_annot <- function(vect_features, Orthogroups, annots) {
  # Convert data frames to data.table for faster lookups
  Orthogroups <- as.data.table(Orthogroups)
  annots <- as.data.table(annots)
  
  # Function to process a single feature
  process_feature <- function(s) {
    # Extract orthologues for the current feature
    lista <- unlist(Orthogroups[Orthogroups[[1]] == s,])
    lista <- lista[! lista %in% c("", s)]
    lista <- lista[!is.na(lista)]
    
    # Retrieve annotations
    Annotations <- sapply(lista, function(i) {
      if (length(annots[[2]][annots[[1]] == i]) != 0) {
        return(annots[[2]][annots[[1]] == i])
      } else {
        return("NA")
      }
    }, USE.NAMES = FALSE)
    
    # Filter out empty annotations
    Annotations <- Annotations[Annotations != "NA"]
    
    # Return data if annotations exist
    if (length(Annotations) != 0) {
      return(data.frame(Orthologue = s, as.data.frame(table(Annotations))))
    } else {
      return(NULL)
    }
  }
  
  # Parallel processing using multiple cores
  results <- mclapply(vect_features, process_feature, mc.cores = argv$p)
  
  # Combine all results into one data frame
  anot_table <- rbindlist(results, fill = TRUE)
  
  return(as.data.frame(anot_table))
}

unsupervised_learning<-function(texto, features, mDF ) {
  dirU=paste(argv$o, texto, sep="/")
  dir.create(dirU, showWarnings = FALSE)
  prefixU=paste(dirU,argv$l, sep="/")
  test=mDF[,features, drop=FALSE]
  test[test>1]=1
  
  dist_mat <- dist(test, method = "binary")
  try(hclust_avg <- hclust(dist_mat, method = 'average'), silent = TRUE)
  if (exists("hclust_avg")) {
    
    pdf(paste(prefixU,"hclust.pdf", sep = "_"))
    plot(hclust_avg, hang = -1, cex = 0.3)
    dev.off()
    
  }
  # Cutting tree by no. of clusters
  try(Members_in_clusters <- cutree(hclust_avg, k = 2 ), silent = TRUE)
  if (exists("Members_in_clusters")) {
    sink(paste(prefixU,"hclust_table.txt", sep = "_"))
    print(table(Members_in_clusters))
    cat("\nMembership: \n")
    print(Members_in_clusters)
    sink()
  }
  
  try(pheatmap(test,
               clustering_distance_rows = "binary",
               clustering_distance_cols = "binary",
               cutree_rows = 2,
               cutree_cols = 2,fontsize = 1.8, fontsize_row=2,filename= paste(prefixU,"heatmap_orthologs_binary.png", sep = "_"), main="Orthologues clustering using binary distances"
  ), silent = TRUE)
  
  ### Beta-diversity (Bray-Curtis):
  
  jaccard_dist = as.matrix(vegdist(t(test), method = "jaccard")) #Orthologues "bray"
  jaccard_distst = as.matrix(vegdist(test, method = "jaccard")) #Strains
  
  try(pheatmap(
    jaccard_dist,
    clustering_distance_rows = as.dist(jaccard_dist),
    clustering_distance_cols = as.dist(jaccard_dist),
    cutree_rows = 2,
    cutree_cols = 2,fontsize = 1.8, filename= paste(prefixU,"heatmap_orthologs_jaccard.png", sep = "_"), main="Orthologs clustering using jaccard dissimilarity"
  ), silent = TRUE)
  try(pheatmap(
    jaccard_distst,
    clustering_distance_rows = as.dist(jaccard_distst),
    clustering_distance_cols = as.dist(jaccard_distst),
    cutree_rows = 2,
    cutree_cols = 2,fontsize = 1.8,filename= paste(prefixU,"heatmap_strains_jaccard.png", sep = "_"), main="Strains clustering using jaccard dissimilarity"
  ), silent = TRUE)
  
  
  try(pheatmap(test,
               clustering_distance_rows = as.dist(jaccard_distst),
               clustering_distance_cols = as.dist(jaccard_dist),
               cutree_rows = 2,
               cutree_cols = 2,fontsize = 1.8, fontsize_row=1.8,filename= paste(prefixU,"heatmap_orthologs_jaccard_strain_otholg.png", sep = "_"), main="Orthologues and Strain clustering using Jaccard distances"
  ), silent = TRUE)
  

  try(pheatmap(test,
               clustering_distance_rows = as.dist(1-cor(t(test), method = "spearman")),
               clustering_distance_cols = as.dist(1-cor(test, method = "spearman")),
               cutree_rows = 2,
               cutree_cols = 2,fontsize = 1.8, fontsize_row=1.8,filename= paste(prefixU,"heatmap_orthologs_spearman_strain_otholg.png", sep = "_"), main="Orthologues and Strain clustering using Spearman distances"
  ), silent = TRUE)
  

  # run the mds algorithm
  try( mds <- metaMDS(jaccard_distst), silent = TRUE)
  if (exists("mds")) {
    color_type = rep("white", ncol(jaccard_distst))
    color_type[colnames(jaccard_distst) %in% PATHOGENIC_STRAINS]="red"
    
    pch_type = rep(21,  ncol(jaccard_distst)) # Circle symbols
    pch_type[colnames(jaccard_distst) %in% PATHOGENIC_STRAINS] = 23 # Diamond symbols
    
    df_mds=data.frame( X1=mds$points[,1], X2=mds$points[,2])
    df_mds$Isolate="Env"
    df_mds$Isolate[row.names(mds$points) %in% PATHOGENIC_STRAINS]="Clin"
    #plot the results
    set.seed(42)

    gmds=ggplot(data=df_mds, aes(x = X1, y = X2, color = Isolate))+
      geom_point()+ labs(x = "MDS1", y = "MDS2", subtitle = "MDS using Jaccard distance")+
      geom_text_repel(label=row.names(df_mds), #mds$points
                      size = 2.5,
                      box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.3, "lines")
      )+
      theme_bw()
    ggsave(paste(prefixU,"MDS_sel_features.pdf", sep="_"), gmds, width = 42, height = 60, units = "cm", device = "pdf")
  }
  
  res.pca3 <- PCA(test,  graph = FALSE)
  A<-fviz_pca_ind(res.pca3, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE, ggtheme = theme_minimal())
  
  ggsave(paste(prefixU,"PCA.pdf", sep="_"), A, width = 32, height = 32, units = "cm", device = "pdf")
  ####
  #UMAP
  
  df.umapSEL = umap(test, n_components = 2, random_state = 15)  #df[,selected_features]
  layoutSEL <- df.umapSEL[["layout"]]
  finalSEL <- data.frame(layoutSEL)
  finalSEL$Isolate="Env"
  finalSEL$Isolate[row.names(layoutSEL) %in% PATHOGENIC_STRAINS]="Clin"

  UMFSEL=ggplot(data=finalSEL, aes(x = X1, y = X2, color = Isolate))+
    geom_point()+
    labs(x = "UMAP1", y = "UMAP2", subtitle = "UMAP plot")+
    geom_text_repel(label=row.names(finalSEL),max.overlaps = 30,
                    size = 2,
                    box.padding = 0.25,
                    point.padding = 0.25 )+theme_bw()
  
  ggsave(paste(prefixU,"UMAP_sel_features.pdf", sep="_"), UMFSEL, width = 32, height = 32, units = "cm", device = "pdf")
  
  ####
  # Control variable colors using their contributions
  B<-fviz_pca_var(res.pca3, col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
  ggsave(paste(prefixU,"PCA_variables.pdf", sep="_"), B, width = 32, height = 32, units = "cm", device = "pdf")
  
  cl.res <- clara(test,argv$c)
  C<-fviz_cluster(cl.res, stand=FALSE, repel = TRUE, show.clust.cent=FALSE, ggtheme = theme_bw(), main="")
  ggsave(paste(prefixU,"Clara_clusters.pdf", sep="_"), C, width = 32, height = 32, units = "cm", device = "pdf")
  
  
  sink(paste(prefixU,"Clusters.txt", sep="_"))
  
  for (u in 1:argv$c) { clustersout(u, cl.res, PATHOGENIC_STRAINS ) }
  sink()
}

visualizing_ortholgues<-function(feature, mdf, parte,tree) {

  SELDF=mdf[,feature, drop=FALSE]
  SELDF$Isolate=0
  SELDF$Isolate[row.names(SELDF) %in% PATHOGENIC_STRAINS]=1
  SELDF=SELDF[order(SELDF$Isolate),]
  subSELDF=SELDF[row.names(SELDF) %in% tree$tip.label,]
  
  for (t in names(subSELDF)) {
    if (t == "Isolate") {
      subSELDF[[t]][subSELDF[[t]] > 0] = paste(argv$w,"isolate", sep=" ")
      subSELDF[[t]][subSELDF[[t]] == 0] = paste(argv$z,"isolate", sep=" ")
    } else {
      subSELDF[[t]][subSELDF[[t]] > 0] = "Orthologue present"
      subSELDF[[t]][subSELDF[[t]] == 0] =  "Orthologue absent"
    }
  }
  if (argv$e != "") {
    Special_group=genome_table[,1][genome_table[,4] == argv$e]
    subSELDF[[argv$e]]="."
    subSELDF[[argv$e]][row.names(subSELDF) %in% Special_group]=argv$e
  }
  
  N2<-ggtree(tree, layout="fan", open.angle = 12, ladderize=TRUE, branch.length="branch.length", size=0.5) +
    geom_tiplab(size=1.3, align=T, linesize=.05, offset =0.02*length(feature)/20 )+ #0.125
    geom_text2(aes(subset = !isTip, label=label), ,na.rm = TRUE, size=1.0,angle = 0,vjust = "inward", hjust = "inward",check_overlap = TRUE)+
    geom_nodepoint(color="black", alpha=1/3, size=0.5,na.rm = TRUE )+ geom_rootedge(rootedge = 0.005)


  TOR2<-gheatmap(N2, subSELDF, offset=0, width=0.02*length(feature),legend_title = "", colnames=T, colnames_position="top",
                 colnames_angle=80,hjust=0, font.size=0.8) +scale_fill_manual(
                   values = c("white", "blue", "red", "purple", "gray", "black", "darkgreen","brown","orange","yellow"), name="")

  ggsave(paste(prefixV,"Tree_with",parte,"orthologues_fan_branch.length.pdf", sep="_"), TOR2, width = 32, height = 32, units = "cm", device = "pdf")
    
  #Importance
  
  try(fitSEL <- phyloglm(Isolate~.,phy=tree,data=SELDF, boot=argv$b, method = meth, btol=argv$a), silent = TRUE)
  
  if (exists("fitSEL")) {
    sink(paste(prefix,parte,"feature_importance.txt", sep="_"))
    print(summary(fitSEL))
    sink()
  }
}

############# set input/output files
dir.create(argv$o, showWarnings = FALSE)
dirtree=paste(argv$o, "tree", sep="/")
dir.create(dirtree, showWarnings = FALSE)

prefix=paste(argv$o,argv$l, sep="/")
prefixT=paste(dirtree,argv$l, sep="/")

cat("INFO: Reading data\n")
cat("      Reading gene count table\n")

gene_count <- read_tsv(argv$g, col_names =TRUE, show_col_types = FALSE,num_threads = argv$p)

# Transpose and clean the gene count dataframe
df_A <- gene_count %>% 
  t() %>% 
  as.data.frame() %>%
  setNames(.[1, ]) %>%
  slice(-1, -n()) %>%  # Removing first and last rows (header and total count row)
  mutate_all(as.numeric)

###### read tree
cat("      Reading phylogenetic tree\n")
if (argv$N != "") {
  cat("      FASTANI output will be used to create the phylogenetic tree\n")
  table <- read_tsv(argv$N, col_names = F,num_threads = argv$p)
  table<-table[!(table$V1==table$V2),]
  table$V3 <- 1 - (table$V3/100)
  table$V4 <- NULL
  table$V5 <- NULL
  
  matrix <- spread(data = table, key = V1, value = V3 )
  row.names(matrix) <- matrix$V2
  my_tree <- as.phylo(upgma(dist(matrix)))
  
  my_tree$tip.label=gsub("Genomes/","",my_tree$tip.label)
  my_tree$tip.label=gsub(".fa","",my_tree$tip.label)
  my_tree$tip.label=gsub("'","",my_tree$tip.label)
  
  write.tree(my_tree, paste(prefixT,"Cleaned_tree.txt", sep="_"))
  tree=my_tree
} else {
  
  tree = ape::read.tree(argv$t)
  ix1 = c(grep("^RS_GCF_", tree$tip.label), grep("^GB_GCA_", tree$tip.label), grep("^GCA_", tree$tip.label))
  tree= drop.tip(tree, ix1)
  
  tree$tip.label=gsub("'","",tree$tip.label)
  
  write.tree(tree, paste(prefixT,"Cleaned_tree.txt", sep="_"))
  
}
# read strain/genome information
cat("      Reading genome/sample information \n")

genome_table <- read.csv(argv$i, sep = ",", na.strings = "NA", header = FALSE, col.names = c("Strain", "File", "Type", "Special")) %>%
  mutate_all(~ gsub("\\s+", " ", .))  # Efficient string replacement in all columns

# Clinical Isolates names
PATHOGENIC_STRAINS <- genome_table$Strain[genome_table$Type == argv$y]

#remove orthologues present/absent in more than user defined percentage (default 95%) of the total strains in the dataset
if (argv$r == "") {
  NC <- length(PATHOGENIC_STRAINS)
  NE <- length(genome_table$Strain[genome_table$Type == argv$d])
  NT <- NC + NE
  CUT_OFF_ORT=0.5*(1+ NC/NT)
} else { CUT_OFF_ORT=argv$r }

# Efficient removal of orthologues based on the cutoff
present_ratios <- colSums(df_A != 0) / nrow(df_A)
toremove <- which(present_ratios > CUT_OFF_ORT | (1 - present_ratios) > CUT_OFF_ORT)

# Remove columns that meet the cutoff condition
df_A <- df_A[, -toremove, drop = FALSE]

# Visualizing tree
y <- ifelse(tree$tip.label %in% PATHOGENIC_STRAINS, 1, 0)
# Check if strain names match tree labels
if (sum(y) != length(PATHOGENIC_STRAINS)) {
  stop("Strain names and tree labels don't match, please check spelling in GENOME_LIST file")
}

# Create a dataframe for visualization
dfP <- data.frame(Isolates = ifelse(y == 1, argv$w, argv$z), row.names = tree$tip.label)

# Mark special group if defined
if (argv$e != "") {
  Special_group <- genome_table$Strain[genome_table$Special == argv$e]
  dfP[[argv$e]] <- " ISOLATES"
  dfP[[argv$e]][row.names(dfP) %in% Special_group] <- argv$e
}

###
sink(paste(prefixT,"labels_tree.tsv", sep="_"))
print(data.frame(Pathogenicity=y,Strain=tree$tip.label), quote = FALSE, row.names = FALSE)
sink()

############# Phylogenetic Generalized Linear Model

df=df_A[match(tree$tip.label,row.names(df_A)),]

if (all(row.names(df)==tree$tip.label)) cat("INFO: Comparative genomics starts\n") else stop("Strain names in Gene count table and tree labels don't match, please check spelling in GENOME_LIST file")

m=mean(tipHeights(tree))
aplhaH=exp(4)*0.95/m #4 is default value in phyloglm IGL and is suggested by IG paper. 

# Setup output directory
dirR <- paste(argv$o, "R_objects_adjusted", sep = "/")
if (!dir.exists(dirR)) dir.create(dirR, showWarnings = FALSE)

# Check if results already exist
res_files <- c("RES_R", "COEFs_R", "ALPHAs_R")
existing_files <- file.exists(paste(dirR, res_files, sep = "/"))
if (!any(existing_files)) {
  meth <- "logistic_IG10"
  set.seed(123)
  # Preallocate lists for results
  RES <- vector("list", ncol(df))
  COEFs <- vector("list", ncol(df))
  ALPHAs <- vector("list", ncol(df))
  STAT_List <- vector("list", ncol(df))
  names(RES) <- names(df)
  names(COEFs) <- names(df)
  names(ALPHAs) <- names(df)
  names(STAT_List) <- names(df)
  
  # Set up progress bar
  cat("      It may take some minutes before the bar progress is updated\n")
  pb <- txtProgressBar(min = 0, max = ncol(df), style = 3, width = 50, char = "-")
  
  # Parallelize using mclapply
  results <- mclapply(1:ncol(df), function(i) {
    set.seed(123)
    skip_pv <- FALSE
    dat <- data.frame(row.names = tree$tip.label, trait = y, predictor = df[[i]])
    dat$predictor[dat$predictor > 1] <- 1
    
    # Try phyloglm fitting
    fit <- tryCatch({
      phyloglm(trait ~ predictor, phy = tree, data = dat, boot = argv$b, method = meth, btol = argv$a)
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      sta <- summary(fit)
      
      if (sta$alpha >= aplhaH) {
        # Firth's penalized likelihood approach
        fitC0 <- tryCatch({
          logistf(trait ~ predictor, data = dat)
        }, error = function(e) NULL)
        
        if (!is.null(fitC0)) {
          sta_firth <- summary(fitC0)
          pvalue <- sta_firth$prob[2]
          return(list(pvalue = pvalue, coef = sta_firth$coefficients[2], alpha = paste0(sta$alpha, "*"), stat = sta_firth, skip_pv = TRUE))
        }
      }
      
      # Default case if Firth's penalized likelihood is not used
      if (!skip_pv) {
        pvalue <- if (argv$b == 0) sta$coefficients[2, 4] else sta$coefficients[2, 6]
        if (pvalue <= 1) {
          return(list(pvalue = pvalue, coef = sta$coefficients[2, 1], alpha = sta$alpha, stat = sta, skip_pv = FALSE))
        }
      }
    }
    
    return(NULL)  # In case of no fit
  }, mc.cores = argv$p)  # Adjust core number if needed
  
  # Process the results
  for (i in 1:ncol(df)) {
    res <- results[[i]]
    if (!is.null(res)) {
      RES[[names(df)[i]]] <- res$pvalue
      COEFs[[names(df)[i]]] <- res$coef
      ALPHAs[[names(df)[i]]] <- res$alpha
      STAT_List[[names(df)[i]]] <- res$stat
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Filter out empty results
  RES <- RES[lapply(RES, length) > 0]
  COEFs <- COEFs[lapply(COEFs, length) > 0]
  ALPHAs <- ALPHAs[lapply(ALPHAs, length) > 0]
  STAT_List <- STAT_List[lapply(STAT_List, length) > 0]
  
  # Save results
  saveRDS(RES, file = paste(dirR, "RES_R", sep = "/"))
  saveRDS(COEFs, file = paste(dirR, "COEFs_R", sep = "/"))
  saveRDS(ALPHAs, file = paste(dirR, "ALPHAs_R", sep = "/"))
  saveRDS(STAT_List, file = paste(dirR, "STAT_List_R", sep = "/"))
  
} else {
  RES=readRDS(file=paste(dirR,"RES_R", sep="/"))
  COEFs=readRDS(file=paste(dirR,"COEFs_R", sep="/"))
  ALPHAs=readRDS(file=paste(dirR,"ALPHAs_R", sep="/"))
  meth="logistic_IG10"
}

res2=p.adjust(RES, method = "fdr")
res3=res2[res2<argv$q]
#### Order by pvalue
res3=res3[order(res3)]

gene_count_clinical = subset(df, row.names(df) %in% PATHOGENIC_STRAINS)
gene_count_others = subset(df, !row.names(df) %in% PATHOGENIC_STRAINS)

Freqs_orth=list()
for (o in names(res3)) {
  
  w1=gene_count_clinical[[o]]
  w1[w1>0]<-1
  w2=gene_count_others[[o]]
  w2[w2>0]<-1
  
  Freqs_orth[[o]][["Clin_Present"]]=length(w1[w1>0])
  Freqs_orth[[o]][["Clin_Abs"]]=length(w1[w1==0])
  Freqs_orth[[o]][["Env_Present"]]=length(w2[w2>0])
  Freqs_orth[[o]][["Env_Abs"]]=length(w2[w2==0])
}

dirA=paste(argv$o, "Annotations", sep="/") 
dir.create(dirA, showWarnings = FALSE)
prefixA=paste(dirA,argv$l, sep="/")

cat("Ortolog","Present in Clincial","Absent in Clincial",paste0("Present in", argv$z),paste0("Absent in", argv$z),"PadjValue","Coeffcient", "alpha\n", sep="\t",
    file=paste(prefixA,"Candidates_enriched_orthologues.tsv", sep="_"))

cat("Ortolog","Present in Clincial","Absent in Clincial",paste0("Present in", argv$z),paste0("Absent in", argv$z),"PadjValue","Coeffcient", "alpha\n", sep="\t",
    file=paste(prefixA,"Candidates_depleted_orthologues.tsv", sep="_"))

abs_markets=rep("0", length(names(res3)))
present_markets=rep("0", length(names(res3)))

k=1
for (o in names(res3)) {
  
  if (COEFs[[o]] < 0 ) {
    abs_markets[k]=o
    cat(o,Freqs_orth[[o]][["Clin_Present"]],
        Freqs_orth[[o]][["Clin_Abs"]],Freqs_orth[[o]][["Env_Present"]],
        Freqs_orth[[o]][["Env_Abs"]],res3[[o]], COEFs[[o]], ALPHAs[[o]], "\n", sep="\t",
        file=paste(prefixA,"Candidates_depleted_orthologues.tsv", sep="_"), append = TRUE)
  } else {
    present_markets[k]=o
    cat(o,Freqs_orth[[o]][["Clin_Present"]],
        Freqs_orth[[o]][["Clin_Abs"]],Freqs_orth[[o]][["Env_Present"]],
        Freqs_orth[[o]][["Env_Abs"]],res3[[o]], COEFs[[o]], ALPHAs[[o]], "\n", sep="\t",
        file=paste(prefixA,"Candidates_enriched_orthologues.tsv", sep="_"), append = TRUE)
  }
  k=k+1
}
abs_markets=abs_markets[abs_markets != "0" ]
present_markets=present_markets[present_markets != "0" ]

cat("Ortolog","Present in Clincial","Absent in Clincial",paste0("Present in", argv$z),paste0("Absent in", argv$z),"PadjValue","Coeffcient", "alpha\n", sep="\t",
    file=paste(prefixA,"core_enriched_orthologues.tsv", sep="_"))
cat("Ortolog","Present in Clincial","Absent in Clincial",paste0("Present in", argv$z),paste0("Absent in", argv$z),"PadjValue","Coeffcient", "alpha\n", sep="\t",
    file=paste(prefixA,"core_depleted_orthologues.tsv", sep="_"))

enriched_core=rep("0", length(present_markets))
depleted_core=rep("0", length(abs_markets))

j=1
l=1
for (o in names(res3)) {
  
  if (o %in% present_markets & Freqs_orth[[o]][["Clin_Present"]]/(Freqs_orth[[o]][["Clin_Abs"]]+Freqs_orth[[o]][["Clin_Present"]]) > core_fraction) {
    cat(o,Freqs_orth[[o]][["Clin_Present"]],
        Freqs_orth[[o]][["Clin_Abs"]],Freqs_orth[[o]][["Env_Present"]],
        Freqs_orth[[o]][["Env_Abs"]],res3[[o]], COEFs[[o]], ALPHAs[[o]], "\n", sep="\t",
        file=paste(prefixA,"core_enriched_orthologues.tsv", sep="_"), append = TRUE)
    enriched_core[j]=o
    j=j+1
  }
  if (o %in% abs_markets & Freqs_orth[[o]][["Clin_Abs"]]/(Freqs_orth[[o]][["Clin_Abs"]]+Freqs_orth[[o]][["Clin_Present"]]) > core_fraction) {
    cat(o,Freqs_orth[[o]][["Clin_Present"]],
        Freqs_orth[[o]][["Clin_Abs"]],Freqs_orth[[o]][["Env_Present"]],
        Freqs_orth[[o]][["Env_Abs"]],res3[[o]], COEFs[[o]], ALPHAs[[o]], "\n", sep="\t",
        file=paste(prefixA,"core_depleted_orthologues.tsv", sep="_"), append = TRUE)
    depleted_core[l]=o
    l=l+1
  }
  
}
enriched_core=enriched_core[enriched_core != "0" ]
depleted_core=depleted_core[depleted_core != "0" ]

selected_features=c()
selected_core_features=c()
if (length(present_markets) > 0 ) {
  cat("\n      Total Enriched orthologs in", argv$w, "Isolates",length(present_markets), "\n")
  selected_features=c(selected_features,present_markets)
  
  cat("#List of Enriched Orthologues\n",file=paste(prefix,"Enriched.txt", sep="_"))
  for (E in present_markets) {cat(E,"\n", file=paste(prefix,"Enriched.txt", sep="_"), append=TRUE)}
  
} else { cat("No enriched Ortholog found") }

if (length(enriched_core) > 0 ) {
  cat("      Total Enriched core orthologs in", argv$w, "Isolates",length(enriched_core), "\n")
  selected_core_features=c(selected_core_features,enriched_core)
  
  cat("#List of Enriched core Orthologues\n",file=paste(prefix,"Enriched_core.txt", sep="_"))
  for (Ec in enriched_core) {cat(Ec,"\n", file=paste(prefix,"Enriched_core.txt", sep="_"), append=TRUE)}
  
} else { cat("No enriched core Ortholog found\n") }

if (length(abs_markets) > 0 ) {
  cat("      Total depleted orthologs in", argv$w,"Isolates",length(abs_markets),"\n")
  selected_features=c(selected_features,abs_markets)
  
  cat("#List of depleted Orthologues\n",file=paste(prefix,"Depleted.txt", sep="_"))
  for (Dp in abs_markets) {cat(Dp,"\n", file=paste(prefix,"Depleted.txt", sep="_"), append=TRUE)}
  
} else {
  cat("No depleted Ortholog found\n")
}

if (length(depleted_core) > 0 ) {
  cat("      Total depleted core orthologs in", argv$w,"Isolates",length(depleted_core),"\n")
  selected_core_features=c(selected_core_features,depleted_core)
  
  cat("#List of Depleted core Orthologues\n",file=paste(prefix,"Depleted_core.txt", sep="_"))
  for (Dpc in depleted_core) {cat(Dpc,"\n", file=paste(prefix,"Depleted_core.txt", sep="_"), append=TRUE)}
  
} else {
  cat("No depleted core Ortholog found\n")
}

if (length(selected_features) == 0) stop("No enriched/depleted Ortholog found\n")

dirV=paste(argv$o, "Visualization", sep="/")
dir.create(dirV, showWarnings = FALSE)
prefixV=paste(dirV,argv$l, sep="/")

if (length(selected_features) >1) {
  cat("INFO: Unsupervise learning steps \n")
  invisible(capture.output(unsupervised_learning("Unsupervise", selected_features, df)))
  cat("INFO: Visualization steps - Selected Orthologues\n")
  visualizing_ortholgues(selected_features, df, "Selected",tree)
}

if (length(selected_core_features) >1) {
  cat("INFO: Unsupervise learning steps using core enriched/depleted orthologs \n")
  invisible(capture.output(unsupervised_learning("Unsupervise_core", selected_core_features, df)))
  cat("INFO: Visualization steps - core Orthologues\n")
  visualizing_ortholgues(selected_core_features, df, "core",tree)
  cat("INFO: Visualization steps - core enriched Orthologues\n")
  visualizing_ortholgues(enriched_core, df, "enriched_core",tree)

}
####
#annotations of selected orthologs
cat("INFO: Annotation steps\n")

Orthogroups_g <- read_tsv(argv$s, col_names =FALSE, show_col_types = FALSE, num_threads = argv$p)

annots_g <- read_delim(argv$m, delim = ";", col_names=FALSE, show_col_types = FALSE, num_threads = argv$p)

if (length(present_markets) >0) {
    cat("INFO: Annotation Enriched orthologues\n")
    pres_anot_table <- print_out_annot(present_markets, Orthogroups_g, annots_g)
    fwrite(pres_anot_table, file=paste(prefixA, "Enriched_orthologues_annotation.tsv", sep="_"), quote=FALSE, sep='\t', row.names = FALSE)
}
if (length(abs_markets) >0) {
  cat("INFO: Annotation Depleted orthologues\n")
  abs_anot_table <- print_out_annot(abs_markets, Orthogroups_g, annots_g)
  fwrite(abs_anot_table, file=paste(prefixA, "Depleted_orthologues_annotation.tsv", sep="_"), quote=FALSE, sep='\t', row.names = FALSE)
  
}
cat("Done!\n")
