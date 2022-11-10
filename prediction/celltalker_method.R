#find significant LR pairs with celltalker
#need compare_group to identify differentially expressed genes between groups
# test.use  wilcox, bimod, roc, t,negbinom, poisson,LR,DESeq2
#organism = "human"
#databaseFile can chose from (celltalkDB_LR_database,icellnet_LR_database,iTALK_LR_database,nichenet_LR_database,ramilowski_pairs,singlecellsignalR_LR_database)
#organism = "mouse"
#databaseFile can chose from (celltalkDB_mouse_LR_database,RNAMagnet_mouse_LR_database)
# each replicate sample group should contain same cell types
# groupFile: need group and replicate file
celltalker_method <-function(dataFile,metaFile,group_replicateFile,referencePath,gene_names = "symbol",organism = "human",cells.reqd = 0,
	freq.pos.reqd = 0,freq.group.in.cluster = 0,databaseFile = "ramilowski_LR_database",outPath){
	suppressMessages({
		library(purrr)
		library(tibble)
		library(Seurat)
		library(tidyverse)
	  library(jsonlite)
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
    
    l <- load(paste(referencePath,"/",databaseFile,".RData",sep=''))
    database <- eval(parse(text = l))
    #ligand,receptor	
    if(file.info(dataFile)$size==0|file.info(metaFile)$size==0|file.info(group_replicateFile)$size==0){
      stopMessage <- "wrong data."
    }else{
      ##format error
      data <- read.table(dataFile,sep='\t',header=T,as.is = T,fill=T, strip.white = TRUE,check.names = F)
      if((class(data[,1]) == "character") | (class(data[,1]) == "factor")){
        rownames(data) <- as.matrix(data[,1])
        data <- data[,-1]
      }
      meta <- read.table(metaFile,sep='\t',header=T,as.is = T,fill=T, strip.white = TRUE,check.names = F)
      group <- read.table(group_replicateFile,sep='\t',header=T,as.is = T,fill=T, strip.white = TRUE,check.names = F)
      
      if(sum("cell" %in% colnames(meta))==0|sum("cell_type" %in% colnames(meta))==0|sum("cell" %in% colnames(group))==0|sum("compare_group" %in% colnames(group))==0|sum("replicate" %in% colnames(group))==0){
        stopMessage <- "error setting column names in metaFile and groupFile."
      }else{
        meta <- merge(meta,group,by="cell")
        inter_ID <- intersect(meta$cell,colnames(data))		
        ####
        if(length(inter_ID)==0){
          stopMessage <- "no matched sample ID."
        }else{
          meta <- meta[meta$cell %in% inter_ID,]
          data <- data[,meta$cell]
          data <- ensmbl_entrez_to_symbol(data,gene_names,organism,ID_transfer_Path = referencePath)
          ##running
          rownames(meta) <- as.matrix(meta$cell)
          ser.obj <- CreateSeuratObject(counts = data, meta.data = meta, project = "scRNA_seq")
          Idents(ser.obj) <- "compare_group"
          
          ligs <- as.character(unique(database$ligand))
          recs <- as.character(unique(database$receptor))
          ligs.present <- rownames(ser.obj)[rownames(ser.obj) %in% ligs]
          recs.present <- rownames(ser.obj)[rownames(ser.obj) %in% recs]
          genes.to.use <- base::union(ligs.present,recs.present)
          
          if(length(genes.to.use)==0){
            stopMessage <- "no matched ligands or receptors in expression data."
          }else{
            markers <- FindAllMarkers(ser.obj,assay="RNA",features=genes.to.use,only.pos=TRUE)
            if(nrow(markers)==0){
              stopMessage <- "no differential expressed ligands or receptors."
            }else{
              ligs.recs.use <- unique(markers$gene)
              interactions.forward1 <- database[as.character(database$ligand) %in% ligs.recs.use,]
              interactions.forward2 <- database[as.character(database$receptor) %in% ligs.recs.use,]
              interact.for <- rbind(interactions.forward1,interactions.forward2)
              
              if(nrow(interact.for)==0){
                stopMessage <- "no differential expressed ligand-receptor interactions."
              }else{
                expr.mat <- GetAssayData(ser.obj,slot="counts")
                defined.clusters <- ser.obj@meta.data$cell_type
                defined.groups <- ser.obj@meta.data$compare_group
                defined.replicates <- ser.obj@meta.data$replicate
                
                reshaped.matrices <- reshape_matrices(count.matrix=expr.mat,clusters=defined.clusters,groups=defined.groups,replicates = defined.replicates,
                                                      ligands.and.receptors=interact.for)
                
                consistent.lig.recs <- create_lig_rec_tib(exp.tib=reshaped.matrices,
                                                          clusters=defined.clusters,groups=defined.groups,
                                                          replicates=defined.replicates,cells.reqd= as.numeric(cells.reqd),
                                                          freq.pos.reqd=as.numeric(freq.pos.reqd),ligands.and.receptors=interact.for)

                put.int <- my_putative_interactions(ligand.receptor.tibble=consistent.lig.recs,
                                                    clusters=defined.clusters,groups=defined.groups,
                                                    freq.group.in.cluster=as.numeric(freq.group.in.cluster),ligands.and.receptors=interact.for,database=database)
                
                group <- unique(meta$compare_group)
                if(length(group)==2){
                  unique.ints <- unique_interactions(put.int,group1 = group[1],group2 = group[2],interact.for)
                  LR <- celltalker_LRS_to_df(unique.ints = unique.ints,put.int = put.int,group1 = group[1],group2 = group[2])
                }else{
                  combination<-combn(unique(meta$compare_group),2)
                  LR <-NULL
                  for(i in 1:ncol(combination)){
                    unique.ints <- unique_interactions(put.int,group1 = combination[1,i],group2 = combination[2,i],interact.for)
                    LR<- rbind(LR,celltalker_LRS_to_df(unique.ints = unique.ints,put.int = put.int,group1 = combination[1,i],group2 = combination[2,i]))
                  }
                  LR <- unique(LR)
                }
                
                if(is.null(LR)){
                  stopMessage <- "no results."
                }else{
                  pair <- paste(LR$ligand,LR$receptor,sep='_')
                  family <- database[match(pair,database$pair),"classification"]
                  LR <- cbind(LR,family)
                  #####============================================plot barplot
                  gene_pairs <- paste(LR[,"ligand"],LR[,"receptor"],sep='_')
                  LR2 <- cbind(gene_pairs,LR)
                  cell_count <- aggregate(LR2[,"gene_pairs"],list(LR2[,"cell_from"],LR2[,"cell_to"],LR2[,"sample_type"]),length)
                  colnames(cell_count) <- c("cell_from","cell_to","sample_type","number")
                  
                  cell_count_pairs <- paste(cell_count[,"cell_from"],cell_count[,"cell_to"],sep='_to_')
                  yAxis0 <-unique(cell_count_pairs)
                  
                  yAxis <- paste(yAxis0,collapse = "','")
                  yAxis_s <- paste("'",yAxis,"'",sep='')
                  
                  cell_count_bar <- cbind(cell_count_pairs,cell_count)
                  sample_type <- unique(cell_count_bar[,"sample_type"])
                  
                  legend <- paste(sample_type,collapse = "','")
                  legend2 <- paste("'",legend,"'",sep = "")
                  
                  bar_res <- lapply(sample_type,function(x){
                    bar_count <- cell_count_bar[cell_count_bar[,"sample_type"]==x,c(1,4,5)]
                    bar_count <- as.matrix(bar_count)
                    
                    name <- paste("name: '",x,"'",sep='')
                    type <- "type: 'bar'"
                    stack <- "stack: 'total'"
                    label <- "label: {show: true}"
                    emphasis <- "emphasis: {focus: 'series'}"
                    
                    data <- matrix(0,ncol=length(yAxis0),nrow=1)
                    colnames(data) <- yAxis0
                    data[,bar_count[,1]] <- bar_count[,3]
                    data2 <- paste(data,collapse = ",")
                    data_s <- paste("data: [",data2,"]",sep='')
                    
                    ss <- paste(name,type,stack,label,emphasis,data_s,sep = ",")
                    ss2 <-paste("{",ss,"}",sep = "")
                  })
                  
                  bar_res2 <- unlist(bar_res)
                  series <- paste(bar_res2,collapse=',')
                  write.table(t(c(legend2,yAxis_s,series)),paste(outPath,"celltalker_method_r_re_barplot.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
                  #####=====================================plot riverplot
                  cell_pairs <- paste(LR2[,"cell_from"],LR2[,"cell_to"],sep='_')
                  LR3 <- cbind(cell_pairs,LR2)
                  cell_count2 <- aggregate(LR3[,"cell_pairs"],list(LR3[,"gene_pairs"]),length)
                  colnames(cell_count2) <- c("gene_pairs","number")
                  #
                  cell_count2_sort <- cell_count2[order(cell_count2[,"number"],decreasing = T),,drop=F]
                  if(nrow(cell_count2_sort)>=50){
                    plot_pairs <- cell_count2_sort[1:50,]
                  }else{
                    plot_pairs <- cell_count2_sort
                  }
                  LR_plot <- LR3[LR3[,"gene_pairs"] %in% plot_pairs[,1],,drop=F] 
                  
                  river_jsoncars <- plot_riverplot(plot_data = LR_plot)
                  cat(river_jsoncars, file = paste(outPath,'celltalker_method_r_re_riverplot.json',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
                  #######======================================output table
                  res_output <- merge(LR3,cell_count2,by = "gene_pairs")
                  
                  res_output <- res_output[,c("ligand","cell_from","receptor","cell_to","sample_type","family","number")]
                  colnames(res_output)[c(6,7)] <- c("Function","score")
                  res_output_df <- as.data.frame(res_output)
                  pp <- paste(res_output_df$ligand,res_output_df$receptor,sep='_')
                  res_output_df$curated <- database[match(pp,database$pair),,drop=F]$Evidence
                  
                  output <- list(data = res_output_df)
                  jsoncars <- toJSON(output, pretty=TRUE)
                  cat(jsoncars, file = paste(outPath,'celltalker_method_r_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
                  #######======================================plot circles
                  cell_count_no_sample <- aggregate(LR3[,"gene_pairs"],list(LR3[,"cell_pairs"]),function(x) length(unique(x)))
                  cell_count_no_sample_sort <- cell_count_no_sample[order(cell_count_no_sample[,2],decreasing = T),]
                  #network
                  LR3_plot <- LR3[LR3[,"cell_pairs"] %in% as.matrix(cell_count_no_sample_sort[1,1]),,drop=F]
                  LR3_plot <- as.matrix(LR3_plot)
                  if(nrow(LR3_plot)>=40){
                    LR3_plot2 <- LR3_plot[1:40,]
                  }else{
                    LR3_plot2 <- LR3_plot
                  }

                  cell_from <- LR3_plot2[1,"cell_from"]
                  cell_to <- LR3_plot2[1,"cell_to"]
                  #node attribute
                  node <- union(LR3_plot2[,"ligand"],LR3_plot2[,"receptor"])
                  
                  ligand <- unique(LR3_plot2[,"ligand"])
                  receptor <- unique(LR3_plot2[,"receptor"])
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
                  
                  ID <- 0:(length(node)-1)
                  degree <- table(c(LR3_plot2[,"ligand"],LR3_plot2[,"receptor"]))
                  degree2 <- as.matrix(degree)
                  degree3 <- degree2[node,]
                  
                  symbolSize <- rep(50,times=length(node))
                  node2 <- cbind(ID,node,symbolSize,degree3,node_attribute)
                  colnames(node2) <- c("id","name","symbolSize","value","category")
                  
                  category <- cbind(0:(length(unique(node2[,"category"]))-1),unique(node2[,"category"]))
                  for(i in 1:nrow(category)){
                    node2[node2[,"category"]==category[i,2],"category"] <- category[i,1]
                  }
                  category2 <- category[,2,drop=F]
                  colnames(category2) <- "name"
                  
                  links_final <- cbind(node2[match(LR3_plot2[,"ligand"],node2[,"name"]),"id"],
                                       node2[match(LR3_plot2[,"receptor"],node2[,"name"]),"id"])
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
                  cat(jsoncars, file = paste(outPath,"celltalker_method_r_re_circle.json",sep='/'), fill = FALSE, labels = NULL, append = FALSE)    
                }
              }
            }
          }
        }
      }
    }
    if(is.null(stopMessage)){
      download_table="celltalker_method_r_re.txt"
      Res_barplot_file="celltalker_method_r_re_barplot.txt"
      Res_riverplot_file="celltalker_method_r_re_riverplot.json"
      Res_circle_file="celltalker_method_r_re_circle.json"
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
      Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "celltalker_method+success"
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
#######################function
celltalker_LRS_to_df <- function(unique.ints,put.int,group1,group2){
  
  group1.to.plot <- pull(unique.ints[1,2])[[1]]
  
  if(length(group1.to.plot)!=0){
  	for.circos.group1 <-pull(put.int[put.int$group %in% group1,2])[[1]][group1.to.plot]
  	LR.group1 <- names(for.circos.group1)
  	output <- lapply(LR.group1,function(x) LRS_to_df(for.circos.group1[[x]],LR = x))
  	output <- unique(do.call(rbind,output))
  	output$sample_type <- group1
  }else{
  	output <- NULL
  }
  group2.to.plot <- pull(unique.ints[2,2])[[1]]
  
  if(length(group2.to.plot)!=0){
  	for.circos.group2 <-  pull(put.int[put.int$group %in% group2,2])[[1]][group2.to.plot]
  	LR.group2 <- names(for.circos.group2)
  	output2 <- lapply(LR.group2,function(x) LRS_to_df(for.circos.group2[[x]],LR = x))
  	output2 <- unique(do.call(rbind,output2))
  	output2$sample_type <- group2
  	
  }else{
  	output2 <- NULL
  }
  commom.to.plot <- pull(unique.ints[3,2])[[1]]
  
  if(length(commom.to.plot)!=0){
  	for.circos.common1 <- pull(put.int[put.int$group %in% group1,2])[[1]][commom.to.plot]
  	LR.commom1 <- names(for.circos.common1)
 	  output3_1 <- lapply(LR.commom1,function(x) LRS_to_df(LR_list = for.circos.common1[[x]],LR = x))
 	  output3_1 <- do.call(rbind,output3_1)
 	  
 	  for.circos.common2 <- pull(put.int[put.int$group %in% group2,2])[[1]][commom.to.plot]
 	  LR.commom2 <- names(for.circos.common2)
 	  output3_2 <- lapply(LR.commom2,function(x) LRS_to_df(for.circos.common2[[x]],LR = x))
  	output3_2 <- do.call(rbind,output3_2)
  	output3 <- unique(rbind(output3_1,output3_2))
  	output3$sample_type <- paste(group1,group2,sep='_')
  }else{
  	output3 <- NULL
  }
  res_cat <- rbind(output,output2,output3)
  return(res_cat)
}

LRS_to_df <- function(LR_list,LR){
  ligand <- unlist(strsplit(LR,'_'))[1]
  receptor <- unlist(strsplit(LR,'_'))[2]
  for(i in 1:length(LR_list$receptor.cells)){
    LR_df <- data.frame(ligand = as.character(ligand),
                        receptor = as.character(receptor), 
                        cell_from = as.character(LR_list$ligand.cells), 
                        cell_to = as.character(LR_list$receptor.cells[i]))
  }
  return(LR_df)
}
##############################
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
      data2 <- data[as.character(res2[,1]),]
      rownames(data2) <- res2[,2]

  }else if(gene_names == "symbol"){
      data2 <- data
  }
  return(data2)
}
##overwrite reshape_matrices
reshape_matrices <- function(count.matrix,clusters,groups,replicates,ligands.and.receptors) {
	#Filter by ligands and receptors
	union.lig.rec <- base::union(ligands.and.receptors$ligand,ligands.and.receptors$receptor)
	mat.fil <- count.matrix[rownames(count.matrix) %in% union.lig.rec,]
	#Create combined metadata
	comb.meta <- data.frame("replicate.id"=replicates,"group.id"=groups,"cluster.id"=clusters)
	#Add metadata to count matrix
	mat.tib <- data.frame(as.matrix(mat.fil))
	mat.tib <- t(mat.tib)
	mat.tib <- cbind(mat.tib,comb.meta)
	mat.tib <- rownames_to_column(mat.tib,var="cell.names")
	mat.tib <- as_tibble(mat.tib)
	#Create nested data.frames by splitting and combining
	sp.rep.id <- mat.tib %>% split(.$replicate.id)
	sp.clust.id <- sp.rep.id %>% map(~.x %>% split(.$cluster.id))
	sp.clust.id.nest <- enframe(sp.clust.id,name="sample",value="expr.matrices")
	ref.table <- table(groups,replicates)
	res <- data.frame("replicate.id"=sp.clust.id.nest$sample,"group"=NA)
	if(nrow(ref.table)==1){
	  fun2 <- function(x) {rownames(ref.table[which.max(ref.table[,colnames(ref.table) %in% sp.clust.id.nest$sample[x]]),,drop=F])}
	  for (i in 1:ncol(ref.table)) {
	    res[i,2] <- fun2(i)
	  }
	}else{
	  fun <- function(x) { names(which.max(ref.table[,colnames(ref.table) %in% sp.clust.id.nest$sample[x]])) }
	  for (i in 1:ncol(ref.table)) {
	    res[i,2] <- fun(i)
	  }
	}
	groups.list <- res %>% split(.$group)
	listing <- vector('list',length=length(groups.list))
	for (i in 1:length(groups.list)) {
		listing[[i]] <- sp.clust.id.nest[sp.clust.id.nest$sample %in% groups.list[[i]]$replicate.id,]
	}
	names(listing) <- names(groups.list)
	groups.samples.nested <- enframe(listing,name="group",value="samples")
}
##overwrite create_lig_rec_tib
create_lig_rec_tib <- function(exp.tib,clusters,groups,replicates,cells.reqd,freq.pos.reqd,ligands.and.receptors) {
	clusters <- as.factor(clusters)
	replicate.tab <- vector("list",length=length(levels(as.factor(groups))))
	names(replicate.tab) <- levels(as.factor(groups))
	
	for (a in 1:length(levels(as.factor(groups)))) {
		pid.layer <- unnest(exp.tib[a,2],cols="samples")
		pid.cell.num <- vector("list",length=nrow(pid.layer))
		lig.rec.res <- vector("list",length=length(levels(clusters)))
		names(lig.rec.res) <- levels(clusters)

		for (z in 1:length(levels(clusters))) {
			for (i in 1:nrow(pid.layer)) {
				pid.cell.num[[i]] <- Reduce(rbind,unnest(pid.layer[i,2],cols="expr.matrices") %>% transmute(n.rows=map(expr.matrices,nrow)) %>% pull(n.rows))
			}
			n.cells.cluster <- Reduce(rbind,lapply(pid.cell.num,function(x) x[z,]))
			
			if (all(n.cells.cluster<cells.reqd)) {
				lig.rec.res[[z]] <- list(ligands=NA,receptors=NA)
			} else if (sum(n.cells.cluster>cells.reqd)/nrow(pid.layer)<=freq.pos.reqd) {
				lig.rec.res[[z]] <- list(ligands=NA,receptors=NA)
			} else {
				cluster.pos <- pid.layer %>% transmute(sample,pos=map(expr.matrices,~.x[[z]])) %>% pull(pos)
				cluster.pos <- cluster.pos[n.cells.cluster>cells.reqd]
				cols.to.drop <- c("cell.names","replicate.id","group.id","cluster.id")
				cluster.pos <- lapply(cluster.pos,function(x) dplyr::select(x,-all_of(cols.to.drop)))
				genes.pos <- lapply(cluster.pos,function(x) apply(x,2,function(y) sum(y>0)>freq.pos.reqd))
				genes.pos.vec <- Reduce(rbind,genes.pos)

			if (is.null(dim(genes.pos.vec))) {
				if (1/nrow(pid.layer)>freq.pos.reqd) {
					genes.to.include <- names(genes.pos.vec)[genes.pos.vec]
					########
					ligands <- genes.to.include[genes.to.include %in% ligands.and.receptors$ligand]
					receptors <- genes.to.include[genes.to.include %in% ligands.and.receptors$receptor]
					lig.rec.res[[z]] <- list(ligands=ligands,receptors=receptors)
				} else {lig.rec.res[[z]] <- list(ligands=NA,receptors=NA)}
				} else {
				genes.to.include <- apply(genes.pos.vec,2,function(x) sum(x)/nrow(pid.layer)>freq.pos.reqd)
				genes.to.include <- names(genes.to.include)[genes.to.include]
				ligands <- genes.to.include[genes.to.include %in% ligands.and.receptors$ligand]
				receptors <- genes.to.include[genes.to.include %in% ligands.and.receptors$receptor]
				lig.rec.res[[z]] <- list(ligands=ligands,receptors=receptors)
				}
			}
		}
		replicate.tab[[a]] <- tibble(cluster.id=levels(clusters),ligands.and.receptors=lig.rec.res)
	}
	lig.rec.tab <- enframe(replicate.tab,name="group",value="lig.rec.exp")
}
###
my_putative_interactions <- function(ligand.receptor.tibble,clusters,groups,freq.group.in.cluster,ligands.and.receptors,database) {
	group.list <- vector("list",length=length(levels(as.factor(groups))))
	names(group.list) <- levels(as.factor(groups))
	clusters.to.include <- vector("list",length=length(group.list))
	names(clusters.to.include) <- names(group.list)
	clusters.per.group <- table(clusters,groups)/rowSums(table(clusters,groups))
	for (q in 1:length(group.list)) {
		interactions.list <- vector("list",length=nrow(ligands.and.receptors))
		names(interactions.list) <- ligands.and.receptors$pair
		sub.list <- list("ligand.cells"=NULL,"receptor.cells"=NULL)
		for (i in 1:length(interactions.list)) {
			interactions.list[[i]] <- sub.list
		}
		group.unnest <- unnest(ligand.receptor.tibble[q,2],cols = c(lig.rec.exp))
		clusters.to.use <- rownames(clusters.per.group)[clusters.per.group[,q]>freq.group.in.cluster]
		group.unnest <- group.unnest[group.unnest$cluster.id %in% clusters.to.use,]
		
		for (z in 1:nrow(group.unnest)) {
			cell.ligs <- pull(group.unnest[z,2])[[1]]$ligands
			rec.list <- lapply(pull(group.unnest),function(x) x[[2]])
			names(rec.list) <- group.unnest$cluster.id

			for (a in 1:length(rec.list)) {
				interactions <- NULL
				if (any(is.null(rec.list[[a]]),is.na(rec.list[[a]]))) {
				} else {
					expanded <- expand.grid(cell.ligs,rec.list[[a]])
					expanded$pair <- paste(expanded[,1],expanded[,2],sep="_")
					interactions <- expanded$pair[expanded$pair %in% database$pair]
				}
				if (length(interactions)==0) {
				} else {
					for (i in 1:length(interactions.list[interactions])) {
						interactions.list[interactions][[i]]$ligand.cells <- c(interactions.list[interactions][[i]]$ligand.cells,names(rec.list)[z])
						interactions.list[interactions][[i]]$ligand.cells <- unique(interactions.list[interactions][[i]]$ligand.cells)
						interactions.list[interactions][[i]]$receptor.cells <- c(interactions.list[interactions][[i]]$receptor.cells,names(rec.list)[a])
						interactions.list[interactions][[i]]$receptor.cells <- unique(interactions.list[interactions][[i]]$receptor.cells)
					}
				}
			}
		}
		group.list[[q]] <- interactions.list
	}
	group.tab <- enframe(group.list,name="group",value="lig_rec_list")
}
##overwrite unique_interactions
unique_interactions <- function(putative.interactions.tib,group1,group2,ligands.and.receptors) {
	group1.exp <- pull(putative.interactions.tib[putative.interactions.tib$group %in% group1,2])[[1]]
	group2.exp <- pull(putative.interactions.tib[putative.interactions.tib$group %in% group2,2])[[1]]
	non.null1 <- group1.exp[sapply(group1.exp,function(x) !any(is.null(x[[1]])))]
	non.null2 <- group2.exp[sapply(group2.exp,function(x) !any(is.null(x[[1]])))]
	unique1 <- setdiff(names(non.null1),names(non.null2))
	unique2 <- setdiff(names(non.null2),names(non.null1))
	common.interactions <- intersect(names(non.null1),names(non.null2))
	enframe(list("unique1v2"=unique1,"unique2v1"=unique2,"common"=common.interactions),name="comparison",value="ligands.and.receptors")
}
#######################plot riverplot
plot_riverplot <- function(plot_data){
  plot_data <- as.matrix(plot_data)
  plot_data[,"cell_from"] <- paste(plot_data[,"cell_from"],"(from)",sep='')
  plot_data[,"cell_to"] <- paste(plot_data[,"cell_to"],"(to)",sep='')
  plot_data[,"ligand"] <- paste(plot_data[,"ligand"],"(ligand)",sep='')
  plot_data[,"receptor"] <- paste(plot_data[,"receptor"],"(receptor)",sep='')
  cell_to_sample <- paste(plot_data[,"sample_type"],"(to)",sep='')
  plot_data <- cbind(plot_data,cell_to_sample)

  node <- matrix(unique(c(plot_data[,"sample_type"],plot_data[,"cell_from"], plot_data[,"ligand"],plot_data[,"family"],
                          plot_data[,"receptor"],plot_data[,"cell_to"],plot_data[,"cell_to_sample"])),ncol=1)
  colnames(node) <- "name"
  node <- as.data.frame(node)
  index <- rep(1,times = nrow(plot_data))
  plot_data <- cbind(index,plot_data)
  
  data_sample_type <- aggregate(plot_data[,"index"],list(plot_data[,c("sample_type")],plot_data[,c("cell_from")]),length)
  colnames(data_sample_type) <- c("source","target","value")
  
  data_cell_from <- aggregate(plot_data[,"index"],list(plot_data[,c("cell_from")],plot_data[,c("ligand")]),length)
  colnames(data_cell_from) <- c("source","target","value")
  
  data_ligand <- aggregate(plot_data[,"index"],list(plot_data[,c("ligand")],plot_data[,c("family")]),length)
  colnames(data_ligand) <- c("source","target","value")
  
  data_family <- aggregate(plot_data[,"index"],list(plot_data[,c("family")],plot_data[,c("receptor")]),length)
  colnames(data_family) <- c("source","target","value")
  
  data_receptor <- aggregate(plot_data[,"index"],list(plot_data[,c("receptor")],plot_data[,c("cell_to")]),length)
  colnames(data_receptor) <- c("source","target","value")
  
  data_to <- aggregate(plot_data[,"index"],list(plot_data[,c("cell_to")],plot_data[,c("cell_to_sample")]),length)
  colnames(data_to) <- c("source","target","value")
  
  link <- Reduce(rbind,list(data_sample_type,data_cell_from,data_ligand,data_family,data_receptor,data_to))
  
  output <- list(nodes=node,links=link)
  jsoncars <- toJSON(output, pretty=TRUE)
}