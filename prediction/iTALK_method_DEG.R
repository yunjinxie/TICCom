# find significant LR pairs with iTALK,do DEG between sample groups
# organism can chose human or mouse
# min_gene_expressed: Genes expressed in minimum number of cells,to find expressed_genes
# min_valid_cells: Minimum number of genes detected in the cell,to find valid_cells
# organism (human,mouse)
# databaseFile (iTALK_LR_database,icellnet_LR_database,RNAMagnet_mouse_LR_database)
# DEG_method
# Wilcox,DESeq2,edgeR
iTALK_method_DEG <- function(dataFile,metaFile,groupFile,referencePath,gene_names = "symbol",organism = "human",min_valid_cells = 0,min_gene_expressed = 0,
	DEG_method = "Wilcox",q_cutoff = 0.05,databaseFile= "iTALK_LR_database",outPath){
  
	suppressMessages({
		library(iTALK)
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
    #database <- read.table(paste(referencePath,"/",databaseFile,".txt",sep=''),sep='\t',header=T,as.is=T)
    #ligand,receptor,classification
    ##format error
    if(file.info(dataFile)$size==0|file.info(metaFile)$size==0|file.info(groupFile)$size==0){
      stopMessage <- "wrong data."
    }else{
      data <- read.table(dataFile,sep='\t',header = T,as.is = T,fill = T,strip.white = T,check.names=F)
      if((class(data[,1]) == "character") | (class(data[,1]) == "factor")){
        rownames(data) <- as.matrix(data[,1])
        data <- data[,-1]
      }
      meta <- read.table(metaFile,sep='\t',header = T,as.is = T,fill = T,strip.white = T,check.names=F)
      group <- read.table(groupFile,sep='\t',header = T,as.is = T,fill = T,strip.white = T,check.names=F)
      
      if(sum("cell" %in% colnames(meta))==0|sum("cell_type" %in% colnames(meta))==0|sum("cell" %in% colnames(group))==0|sum("compare_group" %in% colnames(group))==0){
        stopMessage <- "error setting column names in metaFile and groupFile."
      }else{
        meta <- merge(meta,group,by="cell")
        inter_ID <- intersect(meta$cell,colnames(data))
        
        if(length(inter_ID)==0){
          stopMessage <- "no matched sample ID."
        }else{
          meta <- meta[meta$cell %in% inter_ID,]
          data <- data[,meta$cell]
          
          #running the method 
          #running
          data <- ensmbl_entrez_to_symbol(data,gene_names,organism,ID_transfer_Path = referencePath)	
          data <- as.data.frame(t(data),stringsAsFactors=F)
          data$cell_type <- meta$cell_type
          data$compare_group <- meta$compare_group
          cell <- unique(data$cell_type)
          ###########each cell type at least has two groups
          n_count <- unlist(lapply(cell,function(x) length(unique(data[data$cell_type==x,]$compare_group))))
          
          if(sum(n_count<2)!=0){
            stopMessage <- "each cell type need have at least two sample groups."
          }else{
            ##########each group in each cell type need at least two samples
            index <- 1:nrow(data)
            index_group <- cbind(index,as.matrix(data$cell_type),as.matrix(data$compare_group))
            sam_count <- aggregate(index_group[,"index"],list(index_group[,2],index_group[,3]),function(x) length(x))
            
            len <- length(cell)*length(unique(data$compare_group))
            if(sum(sam_count[,3]>=2)!=len){
              stopMessage <- "each sample group in each cell type need have at least two samples."
            }else{
              #cat("DEG_method can chose from Wilcox,DESeq2,SCDE,monocle,edgeR,DESingle,MAST")

              DEG_gene <- do.call(
                rbind,lapply(cell,function(x) DEG(data %>% dplyr::filter(cell_type == x), method = DEG_method, 
                                                  min_gene_expressed = as.numeric(min_gene_expressed), 
                                                  min_valid_cells = as.numeric(min_valid_cells), 
                                                  q_cut = as.numeric(q_cutoff))))
              
              DEG_gene <- na.omit(DEG_gene)
              if(nrow(DEG_gene)==0){
                stopMessage <- "no differential expressed genes among sample groups."
              }else{
                comm_list<-unique(database$classification)
                res_cat <- do.call(rbind,lapply(comm_list,function(x) my_FindLR(DEG_gene,datatype='DEG',comm_type=x,database = database)))
                
                if(is.null(res_cat)){
                  stopMessage <- "no results."
                }else{
                  #####rank based on logFC,sort by mean(rank)
                  ligand_FC <- abs(res_cat[,3])
                  ligand_rank <- rank(ligand_FC,ties.method = "min")
                  
                  receptor_FC <- abs(res_cat[,6])
                  receptor_rank <- rank(receptor_FC,ties.method = "min")
                  
                  mean_rank <- rowMeans(cbind(ligand_rank,receptor_rank))
                  res_cat$score <- mean_rank
                  #######output table
                  res_output <- res_cat[,c("ligand","cell_from","cell_from_logFC",
                                           "receptor","cell_to","cell_to_logFC","comm_type","score")]
                  colnames(res_output)[c(3,6,7)] <- c("logFC1","logFC2","Function")
                  res_output[,3] <- round(as.numeric(res_output[,3]),4)
                  res_output[,6] <- round(as.numeric(res_output[,6]),4)
                  
                  res_output[res_output[,3]=="-Inf"|res_output[,3]=="Inf",3] <- "-"
                  res_output[res_output[,6]=="-Inf"|res_output[,6]=="Inf",6] <- "-"
                  
                  res_output_df <- as.data.frame(res_output)
                  pp <- paste(res_output_df$ligand,res_output_df$receptor,sep='_')
                  res_output_df$curated <- database[match(pp,database$pair),,drop=F]$Evidence
                  
                  output <- list(data = res_output_df)
                  jsoncars <- toJSON(output, pretty=TRUE)
                  cat(jsoncars, file = paste(outPath,'iTALK_method_DEG_r_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
                  
                  #########plot barplot
                  gene_pairs <- paste(res_cat[,"ligand"],res_cat[,"receptor"],sep='_')
                  res_cat2 <- cbind(gene_pairs,res_cat)
                  cell_count <- aggregate(res_cat2[,"gene_pairs"],list(res_cat2[,"cell_from"],res_cat2[,"cell_to"]),length)
                  colnames(cell_count) <- c("cell_from","cell_to","number")
                  
                  xAxis0 <- paste(cell_count[,"cell_from"],cell_count[,"cell_to"],sep='_to_')
                  xAxis <- paste(xAxis0,collapse = ",")
                  #xAxis_s <- paste("'",xAxis,"'",sep='')
                  
                  yAxis <- paste(cell_count[,"number"],collapse=',')
                  
                  #write.table(t(c(xAxis_s,yAxis)),paste(outPath,"iTALK_method_DEG_r_re_barplot.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
                  ########plot riverplot
                  res_cat_sort <- res_cat[order(res_cat[,"score"],decreasing = F),]
                  if(nrow(res_cat_sort)>=50){
                    res_cat_plot <- res_cat_sort[1:50,c("cell_from","ligand","comm_type","receptor","cell_to")]
                  }else{
                    res_cat_plot <- res_cat_sort[,c("cell_from","ligand","comm_type","receptor","cell_to"),drop=F]
                  }
                  river_jsoncars <- plot_riverplot(plot_data = res_cat_plot)
                  cat(river_jsoncars, file = paste(outPath,'iTALK_method_DEG_r_re_riverplot.json',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
                  
                  ######################################################plot circles
                  cell_count_sort <- cell_count[order(cell_count[,3],decreasing = T),]
                  
                  #network
                  res_cat_plot <- res_cat[res_cat[,"cell_from"]==cell_count_sort[1,1]&res_cat[,"cell_to"]==cell_count_sort[1,2],,drop=F]
                  if(nrow(res_cat_plot)>=40){
                    res_cat_plot2 <- res_cat_plot[1:40,]
                  }else{
                    res_cat_plot2 <- res_cat_plot
                  }
                  #node
                  node <- union(res_cat_plot2[,"ligand"],res_cat_plot2[,"receptor"])
                  degree <- table(c(res_cat_plot2[,"ligand"],res_cat_plot2[,"receptor"]))
                  degree2 <- as.matrix(degree)
                  degree3 <- degree2[node,]
                  
                  ID <- 0:(length(node)-1)
                  
                  symbolSize <- rep(50,times=length(node))
                  
                  ligand <- unique(res_cat_plot2[,"ligand"])
                  receptor <- unique(res_cat_plot2[,"receptor"])
                  cell_from <- unname(res_cat_plot2[1,"cell_from"])
                  cell_to <- unname(res_cat_plot2[1,"cell_to"])
                  
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
                  
                  links_final <- cbind(node2[match(res_cat_plot2[,"ligand"],node2[,"name"]),"id"],
                                       node2[match(res_cat_plot2[,"receptor"],node2[,"name"]),"id"])
                  colnames(links_final) <- c("source","target")
                  #######output
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
                  
                  cat(jsoncars, file = paste(outPath,"iTALK_method_DEG_r_re_circle.json",sep='/'), fill = FALSE, labels = NULL, append = FALSE)
                }
              }
            }
          }
        }
      }
    }
    if(is.null(stopMessage)){
      download_table="iTALK_method_DEG_r_re.txt"
      Res_barplot_x=xAxis
      Res_barplot_value=yAxis
      Res_riverplot_file="iTALK_method_DEG_r_re_riverplot.json"
      Res_circle_file="iTALK_method_DEG_r_re_circle.json"
      resut_merge=paste0("{",
                         '"download_table" :','"',download_table,'",',
                         '"Res_barplot_x" :','"',Res_barplot_x,'",',
                         '"Res_barplot_value" :','"',Res_barplot_value,'",',
                         '"Res_riverplot_file" :','"',Res_riverplot_file,'",',
                         '"Res_circle_file" :','"',Res_circle_file,'",',
                         '"error_attention" :','"no',
                         '"}')
      write.table(resut_merge,paste(outPath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
      
      Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
      Task_ID_new <- Task_ID
      Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "iTALK_method_DEG+success"
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
######################function
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
##change from FindLR of iTALK
my_FindLR<-function(data_1,data_2=NULL,datatype,comm_type,database){
  database<-database[database$classification==comm_type,]
  if(datatype=='mean count'){
    gene_list_1<-data_1
    if(is.null(data_2)){
      gene_list_2<-gene_list_1
    }else{
      gene_list_2<-data_2
    }
    ligand_ind<-which(database$ligand %in% gene_list_1$gene)
    receptor_ind<-which(database$receptor %in% gene_list_2$gene)
    ind<-intersect(ligand_ind,receptor_ind)
    FilterTable_1<-database[ind,c('ligand','receptor')] %>%
      left_join(gene_list_1[,c('gene','exprs','cell_type')],by=c('ligand'='gene')) %>%
      dplyr::rename(cell_from_mean_exprs=exprs,cell_from=cell_type) %>%
      left_join(gene_list_2[,c('gene','exprs','cell_type')],by=c('receptor'='gene')) %>%
      dplyr::rename(cell_to_mean_exprs=exprs,cell_to=cell_type)
    ligand_ind<-which(database$ligand %in% gene_list_2$gene)
    receptor_ind<-which(database$receptor %in% gene_list_1$gene)
    ind<-intersect(ligand_ind,receptor_ind)
    FilterTable_2<-database[ind,c('ligand','receptor')] %>%
      left_join(gene_list_2[,c('gene','exprs','cell_type')],by=c('ligand'='gene')) %>%
      dplyr::rename(cell_from_mean_exprs=exprs,cell_from=cell_type) %>%
      left_join(gene_list_1[,c('gene','exprs','cell_type')],by=c('receptor'='gene')) %>%
      dplyr::rename(cell_to_mean_exprs=exprs,cell_to=cell_type)
    FilterTable<-rbind(FilterTable_1,FilterTable_2)
  }else if(datatype=='DEG'){
    gene_list_1<-data_1
    if(is.null(data_2)){
      gene_list_2<-gene_list_1
    }else{
      gene_list_2<-data_2
    }
    ligand_ind<-which(database$ligand %in% gene_list_1$gene)
    receptor_ind<-which(database$receptor %in% gene_list_2$gene)
    ind<-intersect(ligand_ind,receptor_ind)
    
    if(length(ind)==0){
      return(NULL)
    }
    
    FilterTable_1<-database[ind,c('ligand','receptor')] %>%
      left_join(gene_list_1[,c('gene','logFC','q.value','cell_type')],by=c('ligand'='gene')) %>%
      dplyr::rename(cell_from_logFC=logFC,cell_from_q.value=q.value,cell_from=cell_type) %>%
      left_join(gene_list_2[,c('gene','logFC','q.value','cell_type')],by=c('receptor'='gene')) %>%
      dplyr::rename(cell_to_logFC=logFC,cell_to_q.value=q.value,cell_to=cell_type)
    
    ligand_ind<-which(database$ligand %in% gene_list_2$gene)
    receptor_ind<-which(database$receptor %in% gene_list_1$gene)
    ind<-intersect(ligand_ind,receptor_ind)
    
    if(length(ind)==0){
      return(NULL)
    }
    
    FilterTable_2<-database[ind,c('ligand','receptor')] %>%
      left_join(gene_list_2[,c('gene','logFC','q.value','cell_type')],by=c('ligand'='gene')) %>%
      dplyr::rename(cell_from_logFC=logFC,cell_from_q.value=q.value,cell_from=cell_type) %>%
      left_join(gene_list_1[,c('gene','logFC','q.value','cell_type')],by=c('receptor'='gene')) %>%
      dplyr::rename(cell_to_logFC=logFC,cell_to_q.value=q.value,cell_to=cell_type)
    FilterTable<-rbind(FilterTable_1,FilterTable_2)
  }else{
    stop('Error: invalid data type')
  }

  FilterTable<-FilterTable[!duplicated(FilterTable),]
  res<-as.data.frame(FilterTable) %>% dplyr::rename(ligand=ligand,receptor=receptor)
  if(datatype=='DEG'){
    res<-res[!(res$cell_from_logFC==0.0001 & res$cell_to_logFC==0.0001),]
  }
  
  if(nrow(res)==0){
    return(NULL)
  }
  res<-res %>% mutate(comm_type=comm_type)
  return(res)
}

#######################plot riverplot
plot_riverplot <- function(plot_data){
  
  plot_data$cell_from <- paste(plot_data$cell_from,"(from)",sep='')
  plot_data$cell_to <- paste(plot_data$cell_to,"(to)",sep='')
  plot_data$ligand <- paste(plot_data$ligand,"(ligand)",sep='')
  plot_data$receptor <- paste(plot_data$receptor,"(receptor)",sep='')

  node <- matrix(unique(c(plot_data$cell_from,plot_data$ligand,plot_data$comm_type,
                          plot_data$receptor,plot_data$cell_to)),ncol=1)
  colnames(node) <- "name"
  node <- as.data.frame(node)
  plot_data$index <- 1
  
  data_cell_from <- aggregate(plot_data[,"index"],list(plot_data[,c("cell_from")],plot_data[,c("ligand")]),length)
  colnames(data_cell_from) <- c("source","target","value")
  
  data_ligand <- aggregate(plot_data[,"index"],list(plot_data[,c("ligand")],plot_data[,c("comm_type")]),length)
  colnames(data_ligand) <- c("source","target","value")
  
  data_receptor <- aggregate(plot_data[,"index"],list(plot_data[,c("comm_type")],plot_data[,c("receptor")]),length)
  colnames(data_receptor) <- c("source","target","value")
  
  data_to <- aggregate(plot_data[,"index"],list(plot_data[,c("receptor")],plot_data[,c("cell_to")]),length)
  colnames(data_to) <- c("source","target","value")
  
  link <- Reduce(rbind,list(data_cell_from,data_ligand,data_receptor,data_to))
  
  output <- list(nodes=node,links=link)
  jsoncars <- toJSON(output, pretty=TRUE)
  
}

