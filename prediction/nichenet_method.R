NicheNet_method <- function(dataFile,metaFile,defined_receiver_cell = "false",genesetFile="empty",groupFile = "empty",referencePath,
                            gene_names = "symbol",organism = "human",q_cutoff=0.05,databaseFile = "nichenet_LR_database",outPath){
  suppressMessages({
    library(tidyverse)
    library(Seurat)
    library(jsonlite)
    #########
    library(dplyr)
    library(ROCR)
    library(caTools)
    library(data.table)
  })
  
  stopMessage <- NULL
  
  UserUploadData <- tail(unlist(strsplit(dataFile,"/")),n=1)
  webPath <- gsub(paste0("/",UserUploadData),"",dataFile)
  ID_new <- paste(tail(unlist(strsplit(webPath,"/")),n=1),collapse = "/")
  ## updata Job_ID.txt, append a new job and state is run
  jobFile <- paste0(referencePath,"/Job_ID.txt")
  if(file.info(jobFile)$size==0){
    Task_ID_new <- cbind(ID_new,"run")
    write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
  }else{
    Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
    Task_ID_new <- rbind(Task_ID,c(ID_new,"run"))
    write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
  }
  
  tryCatch({
    
    if(file.info(dataFile)$size==0|file.info(metaFile)$size==0){
      stopMessage <- "wrong data."
    }else{
      data <- read.table(dataFile,sep='\t',header = T,as.is = T,fill = T,strip.white = T,check.names=F)
      if((class(data[,1]) == "character") | (class(data[,1]) == "factor")){
        rownames(data) <- as.matrix(data[,1])
        data <- data[,-1]
      }
      meta <- read.table(metaFile,sep='\t',header = T,as.is = T,fill = T,strip.white = T,check.names=F)
      
      if(organism == "human"){
        ligand_target_matrix <- readRDS(paste(referencePath,"ligand_target_matrix.rds",sep='/'))
        lr_network <- readRDS(paste(referencePath,"lr_network.rds",sep='/'))
        weighted_networks_lr <- readRDS(paste(referencePath,"weighted_networks_lr.rds",sep='/'))
      }else{
        ligand_target_matrix <- readRDS(paste(referencePath,"mouse_ligand_target_matrix.rds",sep='/'))
        lr_network <- readRDS(paste(referencePath,"mouse_lr_network.rds",sep='/'))
        weighted_networks_lr <- readRDS(paste(referencePath,"mouse_weighted_networks_lr.rds",sep='/'))
      }
      l <- load(paste(referencePath,"/",databaseFile,".RData",sep=''))
      database <- eval(parse(text = l))
      #database <- read.table(paste(referencePath,"/",databaseFile,".txt",sep=''),sep='\t',header = T,as.is=T,check.names = F)
      ligands0 <- database %>% pull(ligand) %>% unique()
      ligands <- intersect(colnames(ligand_target_matrix),ligands0)
      
      if(length(ligands)==0){
        stopMessage <- "no ligands."
      }else{
        receptors <- database %>% pull(receptor) %>% unique()

        if(defined_receiver_cell == "false"&genesetFile == "empty"){
          
          if(file.info(groupFile)$size==0){
            stopMessage <- "wrong groupFile."
          }else{
            group <- read.table(groupFile,sep='\t',header=T,as.is = T,fill=T, strip.white = TRUE)
            
            if(sum("cell" %in% colnames(meta))==0|sum("cell_type" %in% colnames(meta))==0|sum("cell" %in% colnames(group))==0|sum("compare_group" %in% colnames(group))==0){
              stopMessage <- "error setting column names in metaFile."
            }else{
              condition <- unique(group$compare_group)
              if(length(condition) < 2){
                stopMessage <- "the number of sample groups is less than 2."
              }else{
                meta <- merge(meta,group,by="cell")
                inter_ID <- intersect(meta$cell,colnames(data))
                if(length(inter_ID)==0){
                  stopMessage <- "no matched sample ID."
                }else{
                  meta <- meta[meta$cell %in% inter_ID,]
                  data <- data[,meta$cell]
                  data <- ensmbl_entrez_to_symbol(data,gene_names,organism,ID_transfer_Path = referencePath)
                  ##running
                  rownames(meta) <- meta$cell
                  ser.obj <- CreateSeuratObject(counts = data, meta.data = meta, project = "scRNA_seq")
                  Idents(ser.obj) <- "cell_type"
                  cell <- unique(ser.obj@meta.data$cell_type)
                  ########running method
                  result <- lapply(cell,function(x){
                    cat(x)
                    cat("\n")
                    receiver_cell <- x
                    seurat_obj_receiver <- subset(ser.obj, idents = receiver_cell)
                    seurat_obj_receiver <- SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["compare_group"]])
                    condition <- unique(group$compare_group)
                    
                    DE_table_receiver <- FindMarkers(object = seurat_obj_receiver, ident.1 = condition[1], ident.2 = NULL) %>% rownames_to_column("gene")
                    colnames(DE_table_receiver) <- c("gene","p_val","avg_logFC","pct.1","pct.2","p_val_adj")
                    
                    target_gene = DE_table_receiver %>% filter(p_val_adj <= as.numeric(q_cutoff) & abs(avg_logFC) >= 0.25) %>% pull(gene)
                    target_gene <- target_gene %>% .[. %in% rownames(ligand_target_matrix)]
                    if (length(target_gene) == 0) {
                      return(NULL)
                    }
                    sender_res <- calculate_ligand_activities(all_cell = cell,receiver_cell,ser.obj,receptors,ligands,ligand_target_matrix,database,target_gene)
                  })
                  result <- do.call(rbind,result)
                }
              }
            }
          }
        }else if(defined_receiver_cell == "true"&groupFile == "empty"){
          
          if(sum("cell" %in% colnames(meta))==0|sum("cell_type" %in% colnames(meta))==0|sum("class" %in% colnames(meta))==0){
            stopMessage <- "error setting column names in metaFile."
          }else{
            inter_ID <- intersect(meta$cell,colnames(data))
            if(length(inter_ID)==0){
              stopMessage <- "no matched sample ID."
            }else{
              meta <- meta[meta$cell %in% inter_ID,]
              data <- data[,meta$cell]
              data <- ensmbl_entrez_to_symbol(data,gene_names,organism,ID_transfer_Path = referencePath)
              
              if(file.info(genesetFile)$size==0){
                stopMessage <- "wrong genesetFile."
              }else{
                target_gene <- read.table(genesetFile,sep='\t',header=F,as.is = T,fill=T, strip.white = TRUE)
                target_gene <- as.matrix(target_gene)
                ##running
                rownames(meta) <- meta$cell
                ser.obj <- CreateSeuratObject(counts = data, meta.data = meta, project = "scRNA_seq")
                Idents(ser.obj) <- "cell_type"
                
                cell <- unique(ser.obj@meta.data$cell_type)
                receiver_cell <- unique(as.matrix(meta[meta$class=='receiver cell',"cell_type"]))
                result <- calculate_ligand_activities(all_cell = cell,receiver_cell,ser.obj,receptors,ligands,ligand_target_matrix,database,target_gene)
              }
            }
          }
        }else{
          stopMessage <- "wrong choices."
          write.table(paste0("{",'"error_attention" :','"',stopMessage,'"}'),paste(outPath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
          
          Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
          Task_ID_new <- Task_ID
          Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "dead"
          write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
        }
        
        if(is.null(result)){
          stopMessage <- "no results."
        }else if(nrow(result)==0){
          stopMessage <- "no results."
        }else{
          gene_pairs <- paste(result[,"ligand"],result[,"receptor"],sep='_')
          com_type <- database[match(gene_pairs,database$pair),"classification",drop=F]
          res_output <- cbind(result[,c(1:4,7),drop=F],com_type)
          ####============================================output table
          colnames(res_output)[6] <- "Function"
          res_output[,5] <- round(as.numeric(res_output[,5]),4)
          
          rownames(res_output) <- NULL
          res_output_df <- as.data.frame(res_output)
          pp <- paste(res_output_df$ligand,res_output_df$receptor,sep='_')
          res_output_df$curated <- database[match(pp,database$pair),,drop=F]$Evidence

          output <- list(data = res_output_df)
          jsoncars <- toJSON(output, pretty=TRUE)
          cat(jsoncars, file = paste(outPath,'NicheNet_method_r_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
          ####==========================================barplot
          res_output2 <- cbind(gene_pairs,res_output)
          cell_count <- aggregate(res_output2[,"gene_pairs"],list(res_output2[,"cell_from"],res_output2[,"cell_to"]),length)
          colnames(cell_count) <- c("cell_from","cell_to","number")

          xAxis0 <- paste(cell_count[,"cell_from"],cell_count[,"cell_to"],sep='_to_')
          xAxis <- paste(xAxis0,collapse = "','")
          xAxis_s <- paste("'",xAxis,"'",sep='')
          yAxis <- paste(cell_count[,"number"],collapse=',')

          write.table(t(c(xAxis_s,yAxis)),paste(outPath,"NicheNet_method_r_re_barplot.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
          ####============================================plot riverplot
          ligand_pearson <- unique(res_output[,c("ligand","pearson")])
          ligand_pearson_sort = ligand_pearson %>% arrange(-pearson)
          
          if(nrow(ligand_pearson_sort)>=50){
            ligand_plot <- ligand_pearson_sort[1:50,]
          }else{
            ligand_plot <- ligand_pearson_sort
          }
          res_output_plot <- res_output[res_output[,"ligand"] %in% as.matrix(ligand_plot[,1]),]
          
          river_jsoncars <- plot_riverplot(plot_data = res_output_plot)
          cat(river_jsoncars, file = paste(outPath,'NicheNet_method_r_re_riverplot.json',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
          #####================================================plot circles
          cell_count_sort <- cell_count[order(as.numeric(cell_count[,"number"]),decreasing = T),]
          #network
          res_output_plot_c <- res_output[res_output[,"cell_from"] %in% as.matrix(cell_count_sort[1,1])&res_output[,"cell_to"] %in% as.matrix(cell_count_sort[1,2]),,drop=F]
          if(nrow(res_output_plot_c)>=40){
            res_output_plot_c2 <- res_output_plot_c[1:40,]
          }else{
            res_output_plot_c2 <-res_output_plot_c
          }
          #node
          node <- union(res_output_plot_c2[,"ligand"],res_output_plot_c2[,"receptor"])
          res_output_plot_c2 <- as.matrix(res_output_plot_c2)
          degree <- table(c(res_output_plot_c2[,"ligand"],res_output_plot_c2[,"receptor"]))
          degree2 <- as.matrix(degree)
          degree3 <- degree2[node,]
          
          ID <- 0:(length(node)-1)
          symbolSize <- rep(30,times = length(node))
          ligand <- unique(res_output_plot_c2[,"ligand"])
          receptor <- unique(res_output_plot_c2[,"receptor"])
          cell_from <- unname(res_output_plot_c2[1,"cell_from"])
          cell_to <- unname(res_output_plot_c2[1,"cell_to"])
          inter <- intersect(ligand,receptor)
          cell_from_label <- paste(cell_from,"(ligand)",sep='')
          cell_to_label <- paste(cell_to,"(receptor)",sep='')
          
          if(length(inter)!=0){
            ligand_only <- setdiff(ligand,inter)
            receptor_only <- setdiff(receptor,inter)
            
            node_attribute <- rep(cell_from_label,times=length(node))
            node_attribute[node %in% receptor_only] <- cell_to_label
            node_attribute[node %in% inter] <- "both"
          }else{
            ligand_only <- ligand
            receptor_only <- receptor
            
            node_attribute <- rep(cell_from_label,times=length(node))
            node_attribute[node %in% receptor_only] <- cell_to_label
          }
          node2 <- cbind(ID,node,symbolSize,degree3,node_attribute)
          colnames(node2) <- c("id","name","symbolSize","value","category")
          
          category <- t(cbind(c(0,cell_from_label),c(1,cell_to_label),c(2,"both")))
          for(i in 1:nrow(category)){
            node2[node2[,"category"]==category[i,2],"category"] <- category[i,1]
          }
          category2 <- category[,2,drop=F]
          colnames(category2) <- "name"
          links_final <- cbind(node2[match(res_output_plot_c2[,"ligand"],node2[,"name"]),"id"],
                               node2[match(res_output_plot_c2[,"receptor"],node2[,"name"]),"id"])
          colnames(links_final) <- c("source","target")
          #output
          rownames(node2) <- NULL
          node2 <- as.data.frame(node2)
          node2[,"symbolSize"] <- as.numeric(as.matrix(node2[,"symbolSize"]))
          node2[,"value"] <- as.numeric(as.matrix(node2[,"value"]))
          node2[,"category"] <- as.numeric(as.matrix(node2[,"category"]))
          
          rownames(links_final) <- NULL
          links_final <- as.data.frame(links_final)
          rownames(category2) <- NULL
          category2 <- as.data.frame(category2)
          
          output <- list(nodes = node2,links = links_final,categories = category2)
          jsoncars <- toJSON(output, pretty=TRUE)
          cat(jsoncars, file = paste(outPath,"NicheNet_method_r_re_circle.json",sep='/'), fill = FALSE, labels = NULL, append = FALSE)
        }
      }
    }
    if(is.null(stopMessage)){
      download_table="NicheNet_method_r_re.txt"
      Res_barplot_file="NicheNet_method_r_re_barplot.txt"
      Res_riverplot_file="NicheNet_method_r_re_riverplot.json"
      Res_circle_file="NicheNet_method_r_re_circle.json"
      resut_merge=paste0("{",
                         '"download_table" :','"',download_table,'",',
                         '"Res_barplot_file" :','"',Res_barplot_file,'",',
                         '"Res_riverplot_file" :','"',Res_riverplot_file,'",',
                         '"Res_circle_file" :','"',Res_circle_file,'",',
                         '"error_attention" :','"no',
                         '"}')
      write.table(resut_merge,paste(outPath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
      
      Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
      Task_ID_new <- Task_ID
      Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "NicheNet_method+success"
      write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
      
      result_message <- paste0("{",'"error_attention" :','"success"',"}")
      return(result_message)
    }else{
      write.table(paste0("{",'"error_attention" :','"',stopMessage,'"}'),paste(outPath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
      
      Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
      Task_ID_new <- Task_ID
      Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "dead"
      write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
      
      result_message <- paste0("{",'"error_attention" :','"',stopMessage,'"}')
      return(result_message)
    }
  },error = function(e){
    stopMessage <- "task quit."
    write.table(paste0("{",'"error_attention" :','"',stopMessage,'"}'),paste(outPath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
    
    Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
    Task_ID_new <- Task_ID
    Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "dead"
    write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
    
    result_message <- paste0("{",'"error_attention" :','"',stopMessage,'"}')
    return(result_message)
  })
}
###############
ensmbl_entrez_to_symbol <- function(data,gene_names,organism,ID_transfer_Path=NULL){
  gene <- rownames(data)
  if(organism == "human"){
    if(gene_names == "symbol"){
      if(grepl("ENSG",gene[1])){
        gene_names = "ensembl"
      }else if(!grepl(paste("[",paste(toupper(letters),collapse ='|'),"]",sep=''),gene[1])){
        gene_names = "entrez"
      }
    }else if(gene_names == "ensembl"){
      if(!grepl(paste("[",paste(toupper(letters),collapse ='|'),"]",sep=''),gene[1])){
        gene_names = "entrez"
      }else if(grepl(paste("[",paste(toupper(letters),collapse ='|'),"]",sep=''),gene[1])& (!grepl("ENSG",gene[1]))){
        gene_names = "symbol"
      }
    }else if(gene_names == "entrez"){
      if(grepl("ENSG",gene[1])){
        gene_names = "ensembl"
      }else if(grepl(paste("[",paste(toupper(letters),collapse ='|'),"]",sep=''),gene[1])& (!grepl("ENSG",gene[1]))){
        gene_names = "symbol"
      }
    }
  }else if(organism == "mouse"){
    if(gene_names == "symbol"){
      if(grepl("ENSMUSG",gene[1])){
        gene_names = "ensembl"
      }else if(!grepl(paste("[",paste(letters,collapse ='|'),"]",sep=''),gene[1],ignore.case=TRUE)){
        gene_names = "entrez"
      }
    }else if(gene_names == "ensembl"){
      if(!grepl(paste("[",paste(letters,collapse ='|'),"]",sep=''),gene[1],ignore.case=TRUE)){
        gene_names = "entrez"
      }else if(grepl(paste("[",paste(letters,collapse ='|'),"]",sep=''),gene[1])){
        gene_names = "symbol"
      }
    }else if(gene_names == "entrez"){
      if(grepl("ENSMUSG",gene[1])){
        gene_names = "ensembl"
      }else if(grepl(paste("[",paste(letters,collapse ='|'),"]",sep=''),gene[1])){
        gene_names = "symbol"
      }
    }
  }
  if(organism == "human"){
    #ID_transfer_data <- read.table(paste(ID_transfer_Path,"ID_transfer_file.txt",sep='/'),sep='\t',header = T,as.is = T,fill=T,strip.white = T,quote = "",check.names = F)
    l2 <- load(paste(ID_transfer_Path,"ID_transfer_file.RData",sep='/'))
    ID_transfer_data <- eval(parse(text = l2))
  }else{
    #ID_transfer_data <- read.table(paste(ID_transfer_Path,"mouse_ID_transfer_file.txt",sep='/'),sep='\t',header = T,as.is = T,fill=T,strip.white = T,quote = "",check.names = F)
    l2 <- load(paste(ID_transfer_Path,"mouse_ID_transfer_file.RData",sep='/'))
    ID_transfer_data <- eval(parse(text = l2))
  }
  if(gene_names == "ensembl"){
    if(grepl("\\.",gene[1])){
      gene <- gsub("\\..*","",gene)
    }
    res <- unique(ID_transfer_data[ID_transfer_data[,2] %in% gene,c(2,3)])
    res2 <- res[match(rownames(data),res[,1]),]
    res2 <- na.omit(res2)
    frenq <- table(res2[,2])
    
    if(max(frenq)>1){
      exp_gene <- data[res2[,1],,drop=F]
      exp_gene2 <- cbind(res2,exp_gene)
      
      uni_exp <- exp_gene2[match(names(frenq[frenq==1]),exp_gene2[,2]),-c(1,2),drop=F]
      rownames(uni_exp) <- names(frenq[frenq==1])
      
      more_exp <- exp_gene2[exp_gene2[,2] %in% names(frenq[frenq>1]),,drop=F]
      more_exp2 <- apply(more_exp[,-c(1,2)],2,function(x) {
        tapply(x, factor(more_exp[,2]), function(x) mean(as.numeric(x))) 
      })
      data2 <- rbind(uni_exp,more_exp2)
    }else{
      data2 <- data[rownames(data) %in% res2[,1],,drop=F]
      rownames(data2) <- res2[match(rownames(data2),res2[,1]),2]
    }
  }else if(gene_names == "entrez"){
    res <- unique(ID_transfer_data[ID_transfer_data[,1] %in% gene,c(1,3)])
    res2 <- res[match(rownames(data),res[,1]),]
    res2 <- na.omit(res2)
    res3 <- res2[match(unique(res2[,2]),res2[,2]),]
    data2 <- data[res3[,1],]
    rownames(data2) <- res3[,2]
  }else if(gene_names == "symbol"){
    data2 <- data
  }
  return(data2)
}
###############
calculate_ligand_activities <- function(all_cell,receiver_cell,ser.obj,receptors,ligands,ligand_target_matrix,database,target_gene){
  sender_cell <- all_cell
  expressed_genes_receiver <- get_expressed_genes(receiver_cell, ser.obj, pct = 0.1)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)	
  if (length(expressed_receptors) == 0) {
    return(NULL)
  }
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  sender_res <- lapply(sender_cell,function(x){
    
    expressed_genes_sender <- get_expressed_genes(x, ser.obj, pct = 0.1)
    expressed_ligands <- intersect(ligands, expressed_genes_sender)
    if (length(expressed_ligands) == 0) {
      return(NULL)
    }
    potential_ligands = database %>% filter(ligand %in% expressed_ligands & 
                                              receptor %in% expressed_receptors) %>% pull(ligand) %>% unique()
    if (length(potential_ligands) == 0) {	
      return(NULL)
    }
    ligand_activities = predict_ligand_activities(geneset = target_gene, 
                                                  background_expressed_genes = background_expressed_genes,
                                                  ligand_target_matrix = ligand_target_matrix,potential_ligands = potential_ligands)
    best_upstream_ligands = ligand_activities %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
    lr_network_top = database %>% filter(ligand %in% best_upstream_ligands & receptor %in% expressed_receptors) %>% distinct(ligand,receptor)
    if(nrow(lr_network_top)==0){
      return(NULL)
    }
    LR <- cbind(lr_network_top[,1],cell_from = x,lr_network_top[,2],cell_to = receiver_cell)
    colnames(LR) <- c("ligand","cell_from","receptor","cell_to")
    ligand_activities <- as.data.frame(ligand_activities)
    LR2 <- merge(LR,ligand_activities,by.x="ligand",by.y= "test_ligand")
    return(LR2)
  })
  sender_res <- do.call(rbind,sender_res)
  return(sender_res)
}
###############
plot_riverplot <- function(plot_data){
  
  plot_data <- as.matrix(plot_data) 
  plot_data[,"cell_from"] <- paste(plot_data[,"cell_from"],"(from)",sep='')
  plot_data[,"cell_to"] <- paste(plot_data[,"cell_to"],"(to)",sep='')
  plot_data[,"ligand"] <- paste(plot_data[,"ligand"],"(ligand)",sep='')
  plot_data[,"receptor"] <- paste(plot_data[,"receptor"],"(receptor)",sep='')
  node <- matrix(unique(c(plot_data[,"cell_from"], plot_data[,"ligand"],plot_data[,"Function"],
                          plot_data[,"receptor"],plot_data[,"cell_to"])),ncol=1)
  colnames(node) <- "name"
  node <- as.data.frame(node)
  index <- rep(1,times = nrow(plot_data))
  plot_data <- cbind(index,plot_data)
  data_cell_from <- aggregate(plot_data[,"index"],list(plot_data[,c("cell_from")],plot_data[,c("ligand")]),length)
  colnames(data_cell_from) <- c("source","target","value")
  data_ligand <- aggregate(plot_data[,"index"],list(plot_data[,c("ligand")],plot_data[,c("Function")]),length)
  colnames(data_ligand) <- c("source","target","value")
  data_receptor <- aggregate(plot_data[,"index"],list(plot_data[,c("Function")],plot_data[,c("receptor")]),length)
  colnames(data_receptor) <- c("source","target","value")
  data_to <- aggregate(plot_data[,"index"],list(plot_data[,c("receptor")],plot_data[,c("cell_to")]),length)
  colnames(data_to) <- c("source","target","value")
  link <- Reduce(rbind,list(data_cell_from,data_ligand,data_receptor,data_to))
  output <- list(nodes=node,links=link)
  jsoncars <- toJSON(output, pretty=TRUE)
}
####==================================nichechetr=========================================####
get_expressed_genes = function(ident, seurat_obj, pct = 0.1, assay_oi = NULL){
  requireNamespace("Seurat")
  requireNamespace("dplyr")
  
  # input check
  
  
  if (!"RNA" %in% names(seurat_obj@assays)) {
    if ("Spatial" %in% names(seurat_obj@assays)) {
      if (class(seurat_obj@assays$Spatial@data) != "matrix" &
          class(seurat_obj@assays$Spatial@data) != "dgCMatrix") {
        warning("Spatial Seurat object should contain a matrix of normalized expression data. Check 'seurat_obj@assays$Spatial@data' for default or 'seurat_obj@assays$SCT@data' for when the single-cell transform pipeline was applied")
      }
      if (sum(dim(seurat_obj@assays$Spatial@data)) == 0) {
        stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$Spatial@data'")
      }
    }
  }
  else {
    if (class(seurat_obj@assays$RNA@data) != "matrix" &
        class(seurat_obj@assays$RNA@data) != "dgCMatrix") {
      warning("Seurat object should contain a matrix of normalized expression data. Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$integrated@data' for integrated data or seurat_obj@assays$SCT@data for when the single-cell transform pipeline was applied")
    }
    if ("integrated" %in% names(seurat_obj@assays)) {
      if (sum(dim(seurat_obj@assays$RNA@data)) == 0 & sum(dim(seurat_obj@assays$integrated@data)) ==
          0)
        stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$integrated@data' for integrated data")
    }
    else if ("SCT" %in% names(seurat_obj@assays)) {
      if (sum(dim(seurat_obj@assays$RNA@data)) == 0 & sum(dim(seurat_obj@assays$SCT@data)) ==
          0) {
        stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$SCT@data' for data corrected via SCT")
      }
    }
    else {
      if (sum(dim(seurat_obj@assays$RNA@data)) == 0) {
        stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data'")
      }
    }
  }
  if (sum(ident %in% unique(Idents(seurat_obj))) != length(ident)) {
    stop("One or more provided cell clusters is not part of the 'Idents' of your Seurat object")
  }
  
  if(!is.null(assay_oi)){
    if(! assay_oi %in% Seurat::Assays(seurat_obj)){
      stop("assay_oi should be an assay of your Seurat object")
    }
  }
  
  # Get cell identities of cluster of interest
  
  
  cells_oi = Idents(seurat_obj) %>% .[Idents(seurat_obj) %in%
                                        ident] %>% names()
  
  # Get exprs matrix: from assay oi or from most advanced assay if assay oi not specifcied
  
  if(!is.null(assay_oi)){
    cells_oi_in_matrix = intersect(colnames(seurat_obj[[assay_oi]]@data), cells_oi)
    exprs_mat = seurat_obj[[assay_oi]]@data %>% .[, cells_oi_in_matrix]
  } else {
    if ("integrated" %in% names(seurat_obj@assays)) {
      warning("Seurat object is result from the Seurat integration workflow. The expressed genes are now defined based on the integrated slot. You can change this via the assay_oi parameter of the get_expressed_genes() functions. Recommended assays: RNA or SCT")
      cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$integrated@data),
                                     cells_oi)
      if (length(cells_oi_in_matrix) != length(cells_oi))
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$integrated@data). Please check that the expression matrix contains cells in columns and genes in rows.")
      exprs_mat = seurat_obj@assays$integrated@data %>% .[,
                                                          cells_oi_in_matrix]
    }
    else if ("SCT" %in% names(seurat_obj@assays) & !"Spatial" %in%
             names(seurat_obj@assays)) {
      warning("Seurat object is result from the Seurat single-cell transform workflow. The expressed genes are defined based on the SCT slot. You can change this via the assay_oi parameter of the get_expressed_genes() functions. Recommended assays: RNA or SCT")
      cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$SCT@data),
                                     cells_oi)
      if (length(cells_oi_in_matrix) != length(cells_oi))
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$SCT@data). Please check that the expression matrix contains cells in columns and genes in rows.")
      exprs_mat = seurat_obj@assays$SCT@data %>% .[, cells_oi_in_matrix]
    }
    else if ("Spatial" %in% names(seurat_obj@assays) &
             !"SCT" %in% names(seurat_obj@assays)) {
      warning("Seurat object is result from the Seurat spatial object. The expressed genes are defined based on the Spatial slot. If the spatial data is spot-based (mixture of cells) and not single-cell resolution, we recommend against directly using nichenetr on spot-based data (because you want to look at cell-cell interactions, and not at spot-spot interactions! ;-) )")
      cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$Spatial@data),
                                     cells_oi)
      if (length(cells_oi_in_matrix) != length(cells_oi))
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$Spatial@data). Please check that the expression matrix contains cells in columns and genes in rows.")
      exprs_mat = seurat_obj@assays$Spatial@data %>% .[, cells_oi_in_matrix]
    }
    else if ("Spatial" %in% names(seurat_obj@assays) &
             "SCT" %in% names(seurat_obj@assays)) {
      warning("Seurat object is result from the Seurat spatial object, followed by the SCT workflow. If the spatial data is spot-based (mixture of cells) and not single-cell resolution, we recommend against directly using nichenetr on spot-based data (because you want to look at cell-cell interactions, and not at spot-spot interactions! The expressed genes are defined based on the SCT slot, but this can be changed via the assay_oi parameter.")
      cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$SCT@data),
                                     cells_oi)
      if (length(cells_oi_in_matrix) != length(cells_oi))
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$Spatial@data). Please check that the expression matrix contains cells in columns and genes in rows.")
      exprs_mat = seurat_obj@assays$SCT@data %>% .[, cells_oi_in_matrix]
    }
    else {
      if (sum(cells_oi %in% colnames(seurat_obj@assays$RNA@data)) ==
          0)
        stop("None of the cells are in colnames of 'seurat_obj@assays$RNA@data'. The expression matrix should contain cells in columns and genes in rows.")
      cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$RNA@data),
                                     cells_oi)
      if (length(cells_oi_in_matrix) != length(cells_oi))
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$RNA@data). Please check that the expression matrix contains cells in columns and genes in rows.")
      exprs_mat = seurat_obj@assays$RNA@data %>% .[, cells_oi_in_matrix]
    }
    
  }
  
  # use defined cells and exprs matrix to get expressed genes
  
  n_cells_oi_in_matrix = length(cells_oi_in_matrix)
  if (n_cells_oi_in_matrix < 5000) {
    genes = exprs_mat %>% apply(1, function(x) {
      sum(x > 0)/n_cells_oi_in_matrix
    }) %>% .[. >= pct] %>% names()
  }
  else {
    splits = split(1:nrow(exprs_mat), ceiling(seq_along(1:nrow(exprs_mat))/100))
    genes = splits %>% lapply(function(genes_indices, exprs,
                                       pct, n_cells_oi_in_matrix) {
      begin_i = genes_indices[1]
      end_i = genes_indices[length(genes_indices)]
      exprs = exprs[begin_i:end_i, ]
      genes = exprs %>% apply(1, function(x) {
        sum(x > 0)/n_cells_oi_in_matrix
      }) %>% .[. >= pct] %>% names()
    }, exprs_mat, pct, n_cells_oi_in_matrix) %>% unlist() %>%
      unname()
  }
  return(genes)
}
##############################
predict_ligand_activities = function(geneset,background_expressed_genes,ligand_target_matrix, potential_ligands, single = TRUE,...){
  setting = list(geneset) %>%
    lapply(convert_gene_list_settings_evaluation, name = "gene set", ligands_oi = potential_ligands, background = background_expressed_genes)
  if (single == TRUE){
    settings_ligand_prediction = setting %>%
      convert_settings_ligand_prediction(all_ligands = potential_ligands, validation = FALSE, single = TRUE)
    ligand_importances = settings_ligand_prediction %>% lapply(get_single_ligand_importances,ligand_target_matrix = ligand_target_matrix, known = FALSE) %>% bind_rows()
    
  } else {
    settings_ligand_prediction = setting %>%
      convert_settings_ligand_prediction(all_ligands = potential_ligands, validation = FALSE, single = FALSE)
    ligand_importances = settings_ligand_prediction %>% lapply(get_multi_ligand_importances,ligand_target_matrix = ligand_target_matrix, known = FALSE, ...) %>% bind_rows()
    
  }
  return(ligand_importances %>% select(test_ligand,auroc,aupr,pearson))
}
###########################
convert_gene_list_settings_evaluation = function(gene_list, name, ligands_oi, background) {
  # input check
  if(!is.character(gene_list))
    stop("gene_list should be character vector")
  if(!is.character(name) | length(name) > 1)
    stop("name should be character vector of length 1")
  if(!is.character(ligands_oi))
    stop("ligands_oi should be character vector")
  if(!is.character(background))
    stop("background should be character vector")
  
  requireNamespace("dplyr")
  
  background = background[(background %in% gene_list) == FALSE]
  
  background_logical = rep(FALSE,times = length(background))
  names(background_logical) = background
  gene_list_logical = rep(TRUE,times = length(gene_list))
  names(gene_list_logical) = gene_list
  response = c(background_logical,gene_list_logical)
  
  return(list(name = name, from = ligands_oi, response = response))
}
#################################
convert_settings_ligand_prediction = function(settings,all_ligands,validation = TRUE, single = TRUE){
  
  # input check
  if(!is.list(settings))
    stop("settings should be a list")
  if(!is.character(all_ligands))
    stop("all_ligands should be a character vector")
  if(!is.logical(validation) | length(validation) != 1)
    stop("validation should be TRUE or FALSE")
  if(!is.logical(single) | length(single) != 1)
    stop("single should be TRUE or FALSE")
  
  requireNamespace("dplyr")
  
  new_settings = list()
  if (validation == TRUE && single == TRUE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      for (k in 1:length(all_ligands)){
        test_ligand = all_ligands[[k]]
        new_settings[[length(new_settings) + 1]] = list(make_new_setting_ligand_prediction_single_validation(setting,test_ligand))
      }
    }
  } else if (validation == TRUE && single == FALSE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      new_settings[[length(new_settings) + 1]] = list(make_new_setting_ligand_prediction_multi_validation(setting,all_ligands))
    }
  } else if (validation == FALSE && single == TRUE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      for (k in 1:length(all_ligands)){
        test_ligand = all_ligands[[k]]
        new_settings[[length(new_settings) + 1]] = list(make_new_setting_ligand_prediction_single_application(setting,test_ligand))
      }
    }
  } else if (validation == FALSE && single == FALSE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      new_settings[[length(new_settings) + 1]] = list(make_new_setting_ligand_prediction_multi_application(setting,all_ligands))
    }
  }
  return(new_settings %>% unlist(recursive = FALSE))
}
###################################
make_new_setting_ligand_prediction_single_application = function(setting,test_ligand){
  new_setting = list()
  new_setting$name = setting$name
  new_setting$from = test_ligand
  new_setting$response = setting$response
  return(new_setting)
}
################################
get_single_ligand_importances = function(setting,ligand_target_matrix, ligands_position = "cols", known = TRUE){
  
  if(!is.logical(known) | length(known) > 1)
    stop("known should be a logical vector: TRUE or FALSE")
  
  requireNamespace("dplyr")
  
  metrics = evaluate_target_prediction(setting, ligand_target_matrix, ligands_position)
  
  colnames(metrics)[2] <- "test_ligand"
  #metrics = metrics %>% rename(ligand = test_ligand)
  if (known == TRUE){
    true_ligand = setting$ligand
    metrics_meta = metrics %>% select(setting,test_ligand) %>% bind_cols(tibble(ligand = true_ligand))
    metrics = inner_join(metrics_meta, metrics, by = c("setting","test_ligand"))
  }
  return(metrics)
}
###################################
evaluate_target_prediction = function(setting,ligand_target_matrix, ligands_position = "cols"){
  ## still make evaluation multiple ligands possible
  # input check
  if (!is.list(setting))
    stop("setting must be a list")
  if(!is.character(setting$from) | !is.character(setting$name))
    stop("setting$from and setting$name should be character vectors")
  if(!is.logical(setting$response) | is.null(names(setting$response)))
    stop("setting$response should be named logical vector containing class labels of the response that needs to be predicted ")
  if(!is.matrix(ligand_target_matrix))
    stop("ligand_target_matrix should be a matrix")
  if(!is.double(ligand_target_matrix) & !is.logical(ligand_target_matrix))
    stop("ligand_target matrix should be of type double if it contains numeric probabilities as predictions; or of type logical when it contains categorical target predictions (TRUE or FALSE)")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")
  
  requireNamespace("dplyr")
  
  if (length(setting$from) == 1){
    ligand_oi = setting$from
  } else {
    ligand_oi = paste0(setting$from,collapse = "-")
  }
  if (ligands_position == "cols"){
    if((ligand_oi %in% colnames(ligand_target_matrix)) == FALSE)
      stop("ligand should be in ligand_target_matrix")
    prediction_vector = ligand_target_matrix[,ligand_oi]
    names(prediction_vector) = rownames(ligand_target_matrix)
  } else if (ligands_position == "rows") {
    if((ligand_oi %in% rownames(ligand_target_matrix)) == FALSE)
      stop("ligand should be in ligand_target_matrix")
    prediction_vector = ligand_target_matrix[ligand_oi,]
    names(prediction_vector) = colnames(ligand_target_matrix)
  }
  response_vector = setting$response
  
  if(sd(prediction_vector) == 0)
    warning("all target gene probability score predictions have same value")
  if(sd(response_vector) == 0)
    stop("all genes have same response")
  performance = evaluate_target_prediction_strict(response_vector,prediction_vector,is.double(prediction_vector))
  output = bind_cols(tibble(setting = setting$name, ligand = ligand_oi), performance)
  
  return(output)
}
##############################
evaluate_target_prediction_strict = function(response,prediction,continuous = TRUE, prediction_response_df = FALSE){
  response_df = tibble(gene = names(response), response = response)
  prediction_df = tibble(gene = names(prediction), prediction = prediction)
  combined = inner_join(response_df,prediction_df, by = "gene")
  if (nrow(combined) == 0)
    stop("Gene names in response don't accord to gene names in ligand-target matrix (did you consider differences human-mouse namings?)")
  prediction_vector = combined$prediction
  names(prediction_vector) = combined$gene
  response_vector = combined$response
  names(response_vector) = combined$gene
  if (continuous == TRUE){
    performance = classification_evaluation_continuous_pred(prediction_vector,response_vector)
    
  } else{
    performance = classification_evaluation_categorical_pred(prediction_vector,response_vector)
  }
  if (prediction_response_df == TRUE){
    output = list(performance = performance, prediction_response_df = combined)
    return(output)
  } else {
    return(performance)
  }
  
}
#######################
classification_evaluation_continuous_pred = function(prediction,response, iregulon = TRUE){
  
  if ((sd(response) == 0 & sd(prediction) == 0) | is.null(prediction) | is.null(response)){ # problems can occur otherwise in these leave-one-in models
    return(dplyr::tibble(auroc = NA,
                         aupr = NA,
                         aupr_corrected = NA,
                         sensitivity_roc = NA,
                         specificity_roc = NA,
                         mean_rank_GST_log_pval = NA,
                         auc_iregulon = NA,
                         auc_iregulon_corrected = NA,
                         pearson = NA,
                         spearman = NA))
  }
  prediction_ROCR = ROCR::prediction(prediction, response)
  performance1 = ROCR::performance(prediction_ROCR, measure="tpr", x.measure="fpr")
  
  performance2 = ROCR::performance(prediction_ROCR, measure="prec", x.measure="rec")
  performance = tibble(tpr = performance1@y.values[[1]], fpr=performance1@x.values[[1]], precision=performance2@y.values[[1]], recall=performance2@x.values[[1]])
  
  performance = performance %>% replace_na(list(recall=0, precision=1))
  
  aupr = caTools::trapz(performance$recall, performance$precision)
  pos_class = sum(response)
  total = length(response)
  aupr_random = pos_class/total
  
  
  metrics = get_split_auroc(prediction, response)
  
  sensitivity = metrics$auroc_sensitivity
  specificity = metrics$auroc_specificity
  auroc = metrics$auroc
  
  cor_p = cor(prediction, response)
  cor_s = cor(prediction, response, method = "s")
  
  cor_p_pval = suppressWarnings(cor.test(as.numeric(prediction), as.numeric(response))) %>% .$p.value
  cor_s_pval = suppressWarnings(cor.test(as.numeric(prediction), as.numeric(response), method =  "s")) %>% .$p.value
  
  mean_rank_GST = limma::wilcoxGST(response, prediction)
  #### now start calculating the AUC-iRegulon
  tbl_perf = tibble(auroc = auroc,
                    aupr = aupr,
                    aupr_corrected = aupr - aupr_random,
                    sensitivity_roc = sensitivity,
                    specificity_roc = specificity,
                    mean_rank_GST_log_pval = -log(mean_rank_GST),
                    pearson_log_pval = -log10(cor_p_pval),
                    spearman_log_pval = -log10(cor_s_pval),
                    pearson = cor_p,
                    spearman = cor_s)
  if (iregulon == TRUE){
    output_iregulon = calculate_auc_iregulon(prediction,response)
    tbl_perf = tibble(auroc = auroc,
                      aupr = aupr,
                      aupr_corrected = aupr - aupr_random,
                      sensitivity_roc = sensitivity,
                      specificity_roc = specificity,
                      mean_rank_GST_log_pval = -log(mean_rank_GST),
                      auc_iregulon = output_iregulon$auc_iregulon,
                      auc_iregulon_corrected = output_iregulon$auc_iregulon_corrected,
                      pearson_log_pval = -log10(cor_p_pval),
                      spearman_log_pval = -log10(cor_s_pval),
                      pearson = cor_p,
                      spearman = cor_s)
  }
  return(tbl_perf)
}
########################
get_split_auroc = function(observed, known) {
  prediction = ROCR::prediction(observed, known)
  performance = ROCR::performance(prediction, measure="spec", x.measure="rec")
  metrics = tibble(tpr = performance@x.values[[1]], spec = performance@y.values[[1]])
  meetingpoint = which(-(1-metrics$spec) + 1 <= metrics$tpr)[[1]] # < or <= ?
  specs = c((1-metrics$spec)[seq_len(meetingpoint)-1],1-metrics$tpr[meetingpoint], 1)
  recs = c(metrics$tpr[seq_len(meetingpoint)], 0)
  auroc_specificity = caTools::trapz(specs, recs)
  auroc = caTools::trapz(1-metrics$spec, metrics$tpr)
  auroc_sensitivity = auroc - auroc_specificity
  tibble(auroc=auroc, auroc_sensitivity=auroc_sensitivity, auroc_specificity=auroc_specificity)
}
#################
calculate_auc_iregulon = function(prior,response){
  
  genes_prior = names(prior)
  dim(prior) = c(1,length(prior))
  colnames(prior) = genes_prior
  rownames(prior) = "ligand"
  
  prior_rank = apply(prior,1,rank_desc)
  rankings = tibble(prior=prior_rank[,1], rn = rownames(prior_rank))
  
  fake_rankings = rankings %>% mutate(rn = sample(rn))
  rankings = data.table::data.table(rankings)
  fake_rankings = data.table::data.table(fake_rankings)
  
  aucMaxRank = 0.03*nrow(rankings)
  
  # calculate enrichment over the expression settings
  
  geneSet = response[response == TRUE] %>% names()
  geneSet = unique(geneSet)
  nGenes = length(geneSet)
  geneSet = geneSet[which(geneSet %in% rankings$rn)]
  
  missing = nGenes-length(geneSet)
  
  gSetRanks = subset(rankings, rn %in% geneSet)[,-"rn", with=FALSE] # gene names are no longer needed
  
  aucThreshold = round(aucMaxRank)
  maxAUC = aucThreshold * nrow(gSetRanks)
  
  auc_iregulon = sapply(gSetRanks, .calcAUC, aucThreshold, maxAUC)
  
  gSetRanks_fake = subset(fake_rankings, rn %in% geneSet)[,-"rn", with=FALSE] # gene names are no longer needed
  auc_iregulon_fake = sapply(gSetRanks_fake, .calcAUC, aucThreshold, maxAUC)
  
  auc_iregulon_corrected = auc_iregulon - auc_iregulon_fake
  return(list(auc_iregulon = auc_iregulon, auc_iregulon_corrected = auc_iregulon_corrected))
}
#################
rank_desc = function(x){rank(desc(x), ties.method = "max")}
#################
.calcAUC = function(oneRanking, aucThreshold, maxAUC)
{
  x = unlist(oneRanking)
  x = sort(x[x<aucThreshold])
  y = 1:length(x)
  a = diff(c(x, aucThreshold)) * y
  return(sum(a)/maxAUC)
}