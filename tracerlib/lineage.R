library(ggplot2)
library(alakazam)
library(igraph)
library(dplyr)


#Parse arguments
myArgs <- commandArgs(trailingOnly = TRUE)
outdir <- myArgs[1]
dnapars_exec <- myArgs[2]

### Read in Change-O db from file ###
#Filenames must be joint with output directory path

H_dbfile <- paste(outdir, "/igblast_H_db-modified_germ-pass.tab", sep="")
K_dbfile <- paste(outdir, "/igblast_K_db-modified_germ-pass.tab", sep="")
L_dbfile <- paste(outdir, "/igblast_L_db-modified_germ-pass.tab", sep="")

H_db <- readChangeoDb(H_dbfile, select = NULL, drop = NULL, seq_upper = TRUE)
K_db <- readChangeoDb(K_dbfile, select = NULL, drop = NULL, seq_upper = TRUE)
L_db <- readChangeoDb(L_dbfile, select = NULL, drop = NULL, seq_upper = TRUE)

v_cutoff <- 4

# Set background colours and isotypes- MUST MAKE DIFFERENT FOR MOUSE
colors <- c('#e6f7ff', '#e5ffcc',
          '#f7dede', '#efbfbf', '#ffffcc',
          '#f1e6ff', '#e2ccff', '#d4b3ff',
                '#c599ff', '#b3ffff', 'grey', "cornsilk3")
isotypes <- c("IgD", "IgM", "IgA1", "IgA2", "IgE", "IgG1", "IgG2",
                    "IgG3", "IgG4", "IgD/IgM", "Multiple isotypes", "Unknown isotype")

#HEAVY CHAIN

# Make Change-O clone objects
clones <- H_db %>%
    group_by(CLONE) %>%
    do(CHANGEO=makeChangeoClone(., text_fields=c("CELL", "ISOTYPE"), 
                                add_count=TRUE))


graphs <- lapply(clones$CHANGEO, buildPhylipLineage, 
                 dnapars_exec=dnapars_exec, rm_temp=TRUE)


# Note, clones with only a single sequence will not be processed.
# A warning will be generated and NULL will be returned by buildPhylipLineage
# These entries may be removed for clarity
graphs[sapply(graphs, is.null)] <- NULL

# The set of tree may then be subset by node count for further 
# analysis, if desired.
graphs <- graphs[sapply(graphs, vcount) >= v_cutoff]

# Plot graphs
for (graph in graphs){
    CLONE <- graph$clone
    plot(graph)
    output <- paste(outdir, "/lineage_trees/H", sep="")

    pdf(paste0(output, CLONE, ".pdf"))
     
    V(graph)$color <- "grey"
    V(graph)$color[V(graph)$name == "Germline"] <- "black"
    V(graph)$color[V(graph)$ISOTYPE == "IGHD"] <- colors[1]
    V(graph)$color[V(graph)$ISOTYPE == "IGHM"] <- colors[2]

    V(graph)$color[V(graph)$ISOTYPE == "IGHA1"] <- colors[3]
    V(graph)$color[V(graph)$ISOTYPE == "IGHA2"] <- colors[4]
    V(graph)$color[V(graph)$ISOTYPE == "IGHE"] <- colors[5]
    V(graph)$color[V(graph)$ISOTYPE == "IGHG1"] <- colors[6]
    V(graph)$color[V(graph)$ISOTYPE == "IGHG2"] <- colors[7]
    V(graph)$color[V(graph)$ISOTYPE == "IGHG3"] <- colors[8]
    V(graph)$color[V(graph)$ISOTYPE == "IGHG4"] <- colors[9]
    V(graph)$color[V(graph)$ISOTYPE == "IGHDM"] <- colors[10]
    V(graph)$color[V(graph)$ISOTYPE == "None"] <- colors[12]

    V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
    V(graph)$label <- V(graph)$CELL
    V(graph)$label.cex <- 1.0
    #If many cells in a node, somehow format so cell names fit within node?


    # Set node sizes
    V(graph)$size <- 60
    V(graph)$size <- 30 + (V(graph)$COLLAPSE_COUNT) + (V(graph)$COLLAPSE_COUNT-1)*5
    V(graph)$size[V(graph)$name == "Germline"] <- 20
    V(graph)$size[grepl("Inferred", V(graph)$name)] <- 15 

    # Remove large default margins
    par(mar=c(0, 0, 0, 0) + 0.1)

    # Plot graph
    plot(graph, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
         vertex.label.color="black")
    # Add legend
    legend_names <- c("Germline", "Inferred", isotypes)
    legend_fill <- c("black", "white", colors)
      
    # Add isotype legend
    legend("topright", legend_names, 
           fill=legend_fill, cex=0.8)
       
    
    dev.off()

    }

# KAPPA CHAIN

# Make Change-O clone objects
clones <- K_db %>%
    group_by(CLONE) %>%
    do(CHANGEO=makeChangeoClone(., text_fields=c("CELL"), 
                                add_count=TRUE))

graphs <- lapply(clones$CHANGEO, buildPhylipLineage, 
                 dnapars_exec=dnapars_exec, rm_temp=TRUE)
                 
graphs[sapply(graphs, is.null)] <- NULL
graphs <- graphs[sapply(graphs, vcount) >= v_cutoff]
    
    
# Plot graphs

for (graph in graphs){
    CLONE <- graph$clone
    plot(graph)
    output <- paste(outdir, "/lineage_trees/K", sep="")
    pdf(paste0(output, CLONE, ".pdf"))
     
    V(graph)$color <- "grey"
    V(graph)$color[V(graph)$name == "Germline"] <- "black"
    V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
    V(graph)$label <- V(graph)$CELL
    V(graph)$label.cex <- 1.0

    V(graph)$size <- 60
    V(graph)$size <- 30 + (V(graph)$COLLAPSE_COUNT) + (V(graph)$COLLAPSE_COUNT-1)*5
    V(graph)$size[V(graph)$name == "Germline"] <- 20
    V(graph)$size[grepl("Inferred", V(graph)$name)] <- 15 

    par(mar=c(0, 0, 0, 0) + 0.1)

    # Plot graph
    plot(graph, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
         vertex.label.color="black")
    # Add legend
    legend_names <- c("Germline", "Inferred", "Sample")
    legend_fill <- c("black", "white", "grey")
      
    # Add isotype legend
    legend("topright", legend_names, 
           fill=legend_fill, cex=0.8)
       
    
    dev.off()

    }
    
# LAMBDA CHAIN

# Make Change-O clone objects
clones <- L_db %>%
    group_by(CLONE) %>%
    do(CHANGEO=makeChangeoClone(., text_fields=c("CELL"), 
                                add_count=TRUE))

graphs <- lapply(clones$CHANGEO, buildPhylipLineage, 
                 dnapars_exec=dnapars_exec, rm_temp=TRUE)
                 
graphs[sapply(graphs, is.null)] <- NULL
graphs <- graphs[sapply(graphs, vcount) >= v_cutoff]
    
    
# Plot graphs

for (graph in graphs){
    CLONE <- graph$clone
    plot(graph)
    output <- paste(outdir, "/lineage_trees/L", sep="")
    pdf(paste0(output, CLONE, ".pdf"))
     
    V(graph)$color <- "grey"
    V(graph)$color[V(graph)$name == "Germline"] <- "black"
    V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
    V(graph)$label <- V(graph)$CELL
    V(graph)$label.cex <- 1.0

    V(graph)$size <- 60
    V(graph)$size <- 30 + (V(graph)$COLLAPSE_COUNT) + (V(graph)$COLLAPSE_COUNT-1)*5
    V(graph)$size[V(graph)$name == "Germline"] <- 20
    V(graph)$size[grepl("Inferred", V(graph)$name)] <- 15 

    par(mar=c(0, 0, 0, 0) + 0.1)

    # Plot graph
    plot(graph, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
         vertex.label.color="black")
    # Add legend
    legend_names <- c("Germline", "Inferred", "Sample")
    legend_fill <- c("black", "white", "grey")
      
    # Add isotype legend
    legend("topright", legend_names, 
           fill=legend_fill, cex=0.8)
       
    
    dev.off()

    }

