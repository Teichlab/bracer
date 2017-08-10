library(ggplot2)
library(alakazam)
library(igraph)
library(dplyr)


#Parse arguments
myArgs <- commandArgs(trailingOnly = TRUE)
outdir <- myArgs[1]
dnapars_exec <- myArgs[2]
species <- myArgs[3]


### Read in Change-O db from file if file exists ###
#Filenames must be joint with output directory path

cat_dbfile <- paste(outdir, "/concatenated_lineage_input.tab", sep="")
H_dbfile <- paste(outdir, "/igblast_H_db-modified_germ-pass.tab", sep="")
K_dbfile <- paste(outdir, "/igblast_K_db-modified_germ-pass.tab", sep="")
L_dbfile <- paste(outdir, "/igblast_L_db-modified_germ-pass.tab", sep="")

cat <- FALSE
H <- FALSE
K <- FALSE
L <- FALSE

if (file.exists(cat_dbfile)){
    cat_db <- readChangeoDb(cat_dbfile, select = NULL, drop = NULL, seq_upper = TRUE)
    cat <- TRUE
} 
    
if (cat!=TRUE){
    if (file.exists(H_dbfile)){
        H_db <- readChangeoDb(H_dbfile, select = NULL, drop = NULL, seq_upper = TRUE)
        H <- TRUE
    } 
    if (file.exists(K_dbfile)){
        K_db <- readChangeoDb(K_dbfile, select = NULL, drop = NULL, seq_upper = TRUE)
        K <- TRUE
    } 
    if (file.exists(L_dbfile)){
        L_db <- readChangeoDb(L_dbfile, select = NULL, drop = NULL, seq_upper = TRUE)
        L <- TRUE
    } 
}

v_cutoff <- 4

# Set background colours and isotypes
if (species=="Hsap"){
    colors <- c('#e6f7ff', '#e5ffcc', '#f7dede', '#efbfbf', '#ffffcc', '#f1e6ff', 
            '#e2ccff', '#d4b3ff', '#c599ff', '#b3ffff', 'grey', "cornsilk3")
    isotypes <- c("IgD", "IgM", "IgA1", "IgA2", "IgE", "IgG1", "IgG2", "IgG3", 
            "IgG4", "IgD and IgM", "Multiple isotypes", "Unknown isotype")
} else if (species=="Mmus"){
    colors <- c('#e6f7ff', '#e5ffcc', '#f7dede', '#f1e6ff', '#e2ccff', 
            '#d4b3ff', '#c599ff', '#a866ff', '#b3ffff', 'grey', "cornsilk3")
    isotypes <- c("IgD", "IgM", "IgA", "IgE", "IgG1", "IgG2A", "IgG2B",
                "IgG2C", "IgG3", "IgD and IgM", "Multiple isotypes", "Unknown isotype")    
} else {
    # Set main isotypes for all other species
    colors <- c('#e6f7ff', '#e5ffcc', '#f7dede', '#ffffcc',
                '#e2ccff', '#b3ffff', 'grey', "cornsilk3")
    isotypes <- c("IgD", "IgM", "IgA", "IgE", "IgG", "IgD and IgM", 
                "Multiple isotypes", "Unknown isotype")
}

# CONCATENATED HEAVY AND LIGHT CHAIN / ONLY HEAVY CHAIN

# Make Change-O clone objects

if (cat==TRUE){
    db <- cat_db
} else if (H==TRUE){
    db <- H_db
}


if (H==TRUE | cat==TRUE){
    clones <- db %>%
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
        if (cat==TRUE){
            output <- paste(outdir, "/lineage_trees/", sep="")
        } else {
            output <- paste(outdir, "/lineage_trees/H", sep="")
        }

        pdf(paste0(output, CLONE, ".pdf"))
     
        V(graph)$color <- "grey"
        V(graph)$color[V(graph)$name == "Germline"] <- "black"
        V(graph)$color[V(graph)$ISOTYPE == "IGHD"] <- '#e6f7ff'
        V(graph)$color[V(graph)$ISOTYPE == "IGHM"] <- '#e5ffcc'
        V(graph)$color[V(graph)$ISOTYPE == "IGHA1"] <- '#f7dede'
        V(graph)$color[V(graph)$ISOTYPE == "IGHA"] <- '#f7dede'
        V(graph)$color[V(graph)$ISOTYPE == "IGHA2"] <- '#efbfbf'
        if (species!="Hsap"){
            if (species!="Mmus"){
                V(graph)$color[V(graph)$ISOTYPE == "IGHA2"] <- '#f7dede'
            }
        }   
        V(graph)$color[V(graph)$ISOTYPE == "IGHE"] <- '#ffffcc'
        V(graph)$color[V(graph)$ISOTYPE == "IGHG1"] <- '#f1e6ff'
        V(graph)$color[V(graph)$ISOTYPE == "IGHG2"] <- '#e2ccff'
        if (species=="Hsap") {
            V(graph)$color[V(graph)$ISOTYPE == "IGHG3"] <- '#d4b3ff'
            V(graph)$color[V(graph)$ISOTYPE == "IGHG4"] <- '#c599ff'
        } else if (species=="Mmus") {
            V(graph)$color[V(graph)$ISOTYPE == "IGHG2A"] <- '#e2ccff'
            V(graph)$color[V(graph)$ISOTYPE == "IGHG2B"] <- '#d4b3ff'
            V(graph)$color[V(graph)$ISOTYPE == "IGHG2C"] <- '#c599ff'
            V(graph)$color[V(graph)$ISOTYPE == "IGHG3"] <- '#a866ff'
        } else {
            V(graph)$color[V(graph)$ISOTYPE == "IGHG1"] <- '#e2ccff'
            V(graph)$color[V(graph)$ISOTYPE == "IGHG2A"] <- '#e2ccff'
            V(graph)$color[V(graph)$ISOTYPE == "IGHG2B"] <- '#e2ccff'
            V(graph)$color[V(graph)$ISOTYPE == "IGHG2C"] <- '#e2ccff'
            V(graph)$color[V(graph)$ISOTYPE == "IGHG3"] <- '#e2ccff'
            V(graph)$color[V(graph)$ISOTYPE == "IGHG4"] <- '#e2ccff'              
        }
        V(graph)$color[V(graph)$ISOTYPE == "IGHDM"] <- '#b3ffff'
        V(graph)$color[V(graph)$ISOTYPE == "None"] <- "cornsilk3"

        V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
        V(graph)$label <- V(graph)$CELL
        V(graph)$label.cex <- 1.0

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

}

# KAPPA CHAIN

# Make Change-O clone objects
if (K==TRUE){
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
    }
    
# LAMBDA CHAIN

# Make Change-O clone objects
if (L==TRUE){
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
    }

