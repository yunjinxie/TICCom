#find significant LR pairs with ICELLNET,do DEG between sample groups
#organism = "human",can not change
ICELLNET_method <- function(dataFile,metaFile,referencePath,gene_names = "symbol",organism = "human",top = 4,direction = "out",databaseFile = "icellnet_LR_database",
                            outPath){
    suppressMessages({
      library(psych)
      library(icellnet)
      library(tidyverse)
      library(plyr)
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
    if(file.info(dataFile)$size==0|file.info(metaFile)$size==0){
      stopMessage <- "wrong data."
    }else{
      data <- read.table(dataFile,sep='\t',header = T,as.is = T,fill = T,strip.white = T,check.names=F)
      if((class(data[,1]) == "character") | (class(data[,1]) == "factor")){
        rownames(data) <- as.matrix(data[,1])
        data <- data[,-1]
      }
      meta <- read.table(metaFile,sep='\t',header = T,as.is = T,fill = T,strip.white = T,check.names=F)
      
      if(sum("cell" %in% colnames(meta))==0|sum("cell_type" %in% colnames(meta))==0){
        stopMessage <- "error setting column names in metaFile."
      }else{
        inter_ID <- intersect(meta$cell,colnames(data))
        ####
        if(length(inter_ID)==0){
          stopMessage <- "no matched sample ID."
        }else{
          meta <- meta[meta$cell %in% inter_ID,]
          data <- data[,meta$cell]
          colnames(meta)[colnames(meta)=="cell"] <- "ID"
          colnames(meta)[colnames(meta)=="cell_type"] <- "Class"
          rownames(meta) <- meta$ID
          data <- ensmbl_entrez_to_symbol(data = data,gene_names = gene_names,organism = organism,ID_transfer_Path = referencePath)
          
          if(databaseFile=="icellnet_LR_database"){
            #database0 <- read.table(paste(referencePath,"/",databaseFile,".txt",sep=''),sep='\t',header = T,as.is=T)
            # l0 <- load(paste(referencePath,"/",databaseFile,".RData",sep=''))
            # database0 <- eval(parse(text = l0))
            l <- load(paste(referencePath,"/icellnet_LR_database4",".RData",sep=''))
            database <- eval(parse(text = l))
            #database <- read.table(paste(referencePath,"/","icellnet_LR_database4.txt",sep=''),sep='\t',header = T,as.is=T,check.names = F)
            int = intersect(as.matrix(rownames(data)), as.matrix(database[, 1:5]))
            n <- as.numeric(top)
            if(length(int)==0 | (length(int) >0 & length(int) < n)){
              stopMessage <- "'top' is too much large or there are no expressed ligands or receptors."
            }else{
              #running the method
              data <- data[int,,drop=F]
              all_cells <- unique(meta$Class)
              scores <- lapply(all_cells,function(x) all_cells_score(cc_cell=x,meta,data,n,database,direction))  
              names(scores) <- all_cells
              
              lr_score <- lapply(all_cells,function(x){
                sc <- scores[[x]][[2]]
                return(sc)
              })
              lr_score <- do.call(rbind,lr_score)
              lr_score <- lr_score[lr_score[,"score"]!=0,,drop=F]
              
              cell_score <- lapply(all_cells,function(x){
                sc <- scores[[x]][[1]]
                return(sc)
              })
              cell_score <- do.call(rbind,cell_score)
            }
          }else{
            #database <- read.table(paste(referencePath,"/",databaseFile,".txt",sep=''),sep='\t',header = T,as.is=T)
            l <- load(paste(referencePath,"/",databaseFile,".RData",sep=''))
            database <- eval(parse(text = l))
            database <- database[,c("ligand","receptor","pair","classification","Evidence")]
            #database0 <- database
            int = intersect(as.matrix(rownames(data)), as.matrix(database[, 1:2]))
            
            n <- as.numeric(top)
            if(length(int)==0 | (length(int) >0 & length(int) < n)){
              stopMessage <- "'top' is too much large or there are no expressed ligands or receptors."
            }else{
              data <- data[int,,drop=F]
              all_cells <- unique(meta$Class)
              scores <- lapply(all_cells,function(x) other_all_cells_score(cc_cell=x,meta,data,n,database,direction))  
              names(scores) <- all_cells
              
              lr_score <- lapply(all_cells,function(x){
                sc <- scores[[x]][[2]]
                return(sc)
              })
              lr_score <- do.call(rbind,lr_score)
              lr_score <- lr_score[lr_score[,"score"]!=0,,drop=F]
              
              cell_score <- lapply(all_cells,function(x){
                sc <- scores[[x]][[1]]
                return(sc)
              })
              cell_score <- do.call(rbind,cell_score)
            }
          }
          #######======================================output table
          if(!is.null(lr_score)){
            if(nrow(lr_score)==0){
              stopMessage <- "no results."
            }else{
              res_output <- lr_score[,c("ligand","cell_from",
                                        "receptor","cell_to","family","score")]
              colnames(res_output)[5] <-"Function"
              res_output[,"score"] <- round(as.numeric(as.matrix(res_output[,"score"])),4)

              rownames(res_output) <- NULL
              res_output_df <- as.data.frame(res_output)
              if(databaseFile=="icellnet_LR_database"){
                res_output_df$Curated <- "manually curated"
              }else{
                pp <- paste(res_output_df$ligand,res_output_df$receptor,sep='_')
                res_output_df$curated <- database[match(pp,database$pair),,drop=F]$Evidence
              }
              output <- list(data = res_output_df)
              jsoncars <- toJSON(output, pretty=TRUE)
              cat(jsoncars, file = paste(outPath,'ICELLNET_method_r_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
              #####==================================cell score heatmap
              cell_score <- as.matrix(cell_score)
              cell_score[,"score"] <- round(as.numeric(cell_score[,"score"]),4)
              
              ucf <- unique(cell_score[,"cell_from"])
              ucfy <- cbind(ucf,c(0:(length(ucf)-1)))
              
              uct <- unique(cell_score[,"cell_to"])
              uctx <- cbind(uct,c(0:(length(uct)-1)))
              
              cell_score2 <- cbind(cell_score,x=0,y=0)
              ####
              for(i in 1:nrow(ucfy)){
                cell_score2[which(cell_score2[,"cell_from"] %in% ucfy[i,1]),"y"] <- ucfy[i,2]
              }
              for(i in 1:nrow(uctx)){
                cell_score2[which(cell_score2[,"cell_to"] %in% uctx[i,1]),"x"] <- uctx[i,2]
              }
              ####[y,x,value]
              xyAxis0 <- paste("[",cell_score2[,"y"],",",cell_score2[,"x"],",",cell_score2[,"score"],"]",sep='')
              xyAxis1 <- paste(xyAxis0,collapse =',')
              xyAxis <- paste("[",xyAxis1,"]",sep = "")
              y0<- paste(ucf,collapse="','")
              y <- paste("['",y0,"']",sep='')
              x0<- paste(uct,collapse="','")
              x <- paste("['",x0,"']",sep='')
              
              value_min <- min(as.numeric(cell_score2[,"score"]))
              value_max <- max(as.numeric(cell_score2[,"score"]))
              write.table((t(c(x,y,xyAxis,value_min,value_max))),paste(outPath,"ICELLNET_method_r_re_heatmap.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
              #####=======================================plot barplot
              gene_pairs <- paste(lr_score[,"ligand"],lr_score[,"receptor"],sep='_')
              lr_score2 <- cbind(gene_pairs,lr_score)
              cell_count <- aggregate(lr_score2[,"gene_pairs"],list(lr_score2[,"cell_from"],lr_score2[,"cell_to"]),length)
              colnames(cell_count) <- c("cell_from","cell_to","number")
              
              xAxis0 <- paste(cell_count[,"cell_from"],cell_count[,"cell_to"],sep='_to_')
              xAxis <- paste(xAxis0,collapse = "','")
              xAxis_s <- paste("'",xAxis,"'",sep='')
              yAxis <- paste(cell_count[,"number"],collapse=',')
              write.table(t(c(xAxis_s,yAxis)),paste(outPath,"ICELLNET_method_r_re_barplot.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
              #####========================================plot riverplot
              lr_score_sort <- lr_score[order(as.numeric(lr_score[,"score"]),decreasing = T),]
              if(nrow(lr_score_sort)>=50){
                lr_score_sort <- lr_score_sort[1:50,c("cell_from","ligand","family","receptor","cell_to")]
              }else{
                lr_score_sort <- lr_score_sort[,c("cell_from","ligand","family","receptor","cell_to"),drop=F]
              }
              river_jsoncars <- plot_riverplot(plot_data = lr_score_sort)
              cat(river_jsoncars, file = paste(outPath,'ICELLNET_method_r_re_riverplot.json',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
              #####========================================plot circles
              cell_score_sort <- cell_score[order(as.numeric(cell_score[,"score"]),decreasing = T),]
              #network
              lr_score_plot <- lr_score[lr_score[,"cell_from"] %in% as.matrix(cell_score_sort[1,1])&lr_score[,"cell_to"] %in% as.matrix(cell_score_sort[1,2]),,drop=F]
              if(nrow(lr_score_plot)>=40){
                lr_score_plot2 <- lr_score_plot[1:40,]
              }else{
                lr_score_plot2 <- lr_score_plot
              }
              #node
              node <- union(lr_score_plot2[,"ligand"],lr_score_plot2[,"receptor"])
              lr_score_plot2 <- as.matrix(lr_score_plot2)
              
              degree <- table(c(lr_score_plot2[,"ligand"],lr_score_plot2[,"receptor"]))
              degree2 <- as.matrix(degree)
              degree3 <- degree2[node,]
              
              ID <- 0:(length(node)-1)
              symbolSize <- rep(30,times = length(node))
              
              ligand <- unique(lr_score_plot2[,"ligand"])
              receptor <- unique(lr_score_plot2[,"receptor"])
              cell_from <- unname(lr_score_plot2[1,"cell_from"])
              cell_to <- unname(lr_score_plot2[1,"cell_to"])
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
              links_final <- cbind(node2[match(lr_score_plot2[,"ligand"],node2[,"name"]),"id"],
                                   node2[match(lr_score_plot2[,"receptor"],node2[,"name"]),"id"])
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
              cat(jsoncars, file = paste(outPath,"ICELLNET_method_r_re_circle.json",sep='/'), fill = FALSE, labels = NULL, append = FALSE)
            }
          }else{
            stopMessage <- "no results."
          }
        }
      }
    }
    if(is.null(stopMessage)){
      download_table="ICELLNET_method_r_re.txt"
      Res_heatmap_file="ICELLNET_method_r_re_heatmap.txt"
      Res_barplot_file="ICELLNET_method_r_re_barplot.txt"
      Res_riverplot_file="ICELLNET_method_r_re_riverplot.json"
      Res_circle_file="ICELLNET_method_r_re_circle.json"
      resut_merge=paste0("{",
                         '"download_table" :','"',download_table,'",',
                         '"Res_heatmap_file" :','"',Res_heatmap_file,'",',
                         '"Res_barplot_file" :','"',Res_barplot_file,'",',
                         '"Res_riverplot_file" :','"',Res_riverplot_file,'",',
                         '"Res_circle_file" :','"',Res_circle_file,'",',
                         '"error_attention" :','"no',
                         '"}')
      write.table(resut_merge,paste(outPath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
      
      Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
      Task_ID_new <- Task_ID
      Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "ICELLNET_method+success"
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
#################
all_cells_score <- function(cc_cell,meta,data,n,database,direction){
  
  cc <- meta[meta$Class==cc_cell,]
  cc_data <- data[,rownames(cc),drop=F]
  cc_data= gene.scaling(data = cc_data, n= n, db = database)
  cc_data$Symbol=rownames(cc_data)
  
  pc <- meta[meta$Class!=cc_cell,]
  pc_data <- data[,rownames(pc)]
  pc_data= gene.scaling(data = pc_data, n= n, db = database)
  pc_data$Symbol=rownames(pc_data)
  
  my_selection <- unique(as.matrix(pc$Class))

  score_computation= my_icellnet.score(direction=direction, PC.data=pc_data, CC.data= cc_data,PC.target = pc, PC=my_selection, db = database)
  
  score <-score_computation[[1]]
  
  score <- score[score[,1]!=0&score[,1]!="NaN",,drop=F]
  if(nrow(score)==0){
    return(list(score2 = NULL,lr_score = NULL))
  }
  
  if(direction=="out"){
      score2 <- data.frame(cell_from=cc_cell,cell_to=rownames(score),score=as.numeric(score),stringsAsFactors=F)
  }else if(direction=="in"){
      score2 <- data.frame(cell_from=rownames(score),cell_to=cc_cell,score=as.numeric(score),stringsAsFactors=F)
  }
  
  lr <- score_computation[[2]]
  lr <- lr[,rownames(score),drop=F]
  is_na <- apply(lr,1,function(x){sum(is.na(x))==length(x)})
  lr2 <- lr[!is_na,,drop=F]
  
  if(is.null(lr2)){
    return(list(score2 = score2,lr_score = NULL))
  }
  LR <- as.matrix(rownames(lr2))
  
  pair_family <- name.lr.couple (db = database, type = "Family")
  LR_family <- pair_family[match(LR,pair_family[,1]),2]
  
  ligand <- apply(LR,1,function(x) unlist(strsplit(x," / "))[1])
  receptor <- apply(LR,1,function(x) unlist(strsplit(x," / "))[2])
  if(direction == 'out'){
      lr_score <- data.frame(ligand = rep(ligand,ncol(lr2)),receptor = rep(receptor,ncol(lr2)),cell_from = cc_cell,
                             cell_to = rep(colnames(lr2),each = nrow(lr2)),family = LR_family,score = as.numeric(as.matrix(lr2)),stringsAsFactors=F)
  }else if(direction == 'in'){
      lr_score <- data.frame(ligand = rep(ligand,ncol(lr2)),receptor = rep(receptor,ncol(lr2)),cell_from = rep(colnames(lr2),each = nrow(lr2)),
                             cell_to = cc_cell,family = LR_family,score = as.numeric(as.matrix(lr2)),stringsAsFactors=F)
  }
  return(list(score2 = score2,lr_score = lr_score))
}

#This function rescale each gene expression across condition to range from minimum 0 to maximum 10.
gene.scaling <- function(data, n, db) {
  dim2 = dim(data)[2]
  dim1 = dim(data)[1]
  int = intersect(as.matrix(rownames(data)), as.matrix(db[, 1:5]))
  data = data[which(rownames(data) %in% int),,drop=F]
  
  if(ncol(data)==1){
    for(i in 1:length(int)){
      max <- data[i,]
      
      if (max > 0) {
        data[i, ] = data[i, ] /max  * 10
        data[i, which(data[i, ] > 10)] <- 10
      }else (data [i, ] <- 0)
    }
    
  }else{
    for (i in 1:length(int)) {
      progressBar(i, max = length(int))
      sorted = sort(data[i, ], decreasing = TRUE)
      
      if(n > ncol(data)){
        max = sum(sorted[1, 1:ncol(data)])/ncol(data)
      }else{
        max = sum(sorted[1, 1:n])/n
      }
      
      if (max > 0) {
        data[i, ] = data[i, ] /max  * 10
        data[i, which(data[i, ] > 10)] <- 10
      }else (data [i, ] <- 0)
      
    }
  }
  
  
  data = (data[complete.cases(data),,drop=F])
  return(data)
}

#PC Vector selecting a list of peripheral cell for the cell-cell communication analysis
my_icellnet.score = function(direction = c("out", "in"),
                          CC.data = CC.data,
                          PC.data = PC.data,
                          PC = PC,
                          PC.target = PC.target,
                          db = db)
  {
  if(dim(CC.data)[2]==1| dim(PC.data)[2] ==1){note("Check that PC.data and/or CC.data contains rownames. Ignore this note if this is the case")}
  #Creation of output matrix
  score_global = matrix(nrow = length(PC), ncol = 1)
  rownames(score_global) = PC
  colnames(score_global) = direction

  lr_global = matrix(nrow = dim(db)[1] , ncol = length(PC))
  colnames(lr_global) = PC
  rownames(lr_global) = name.lr.couple (db, type = "Family")[, 1]

  #Compute the score
  if (direction == "out") {
    for (cell in PC) {
      # check all cell types defined in PC exists in PC.data and select PC.target$ID
      cell.IDs=PC.target$ID[grepl(cell, PC.target$Class)]
      if (length(cell.IDs)==0){
        warning (paste0( "Cell type ", cell, " is not found in the PC.target file"))
        score_global[cell, ] = "NaN"
      }else{
          rc.PC = my_receptor.average.RNAseq(db = db,data = as.data.frame(PC.data[, cell.IDs], row.names = rownames(PC.data)))
          lg.CC = my_ligand.average.RNAseq(db = db,data = CC.data)
          lr_global[, cell] = lg.CC * rc.PC
          score_global[cell, ] = sum(lg.CC * rc.PC, na.rm = T)
      }
    }

  } else if (direction == "in") {
    for (cell in PC) {
      # check all cell types defined in PC exists in PC.data and select PC.target$ID
      cell.IDs=PC.target$ID[grepl(cell, PC.target$Class)]
      if (length(cell.IDs)==0){
          warning (paste0( "Cell type ", cell, " is not found in the PC.target file"))
          score_global[cell, ] = "NaN"
      } else {
        lg.PC = my_ligand.average.RNAseq(db = db,data = as.data.frame(PC.data[, cell.IDs ], row.names = rownames(PC.data)))
        rc.CC = my_receptor.average.RNAseq(db = db,data = CC.data)
        lr_global[, cell] = lg.PC * rc.CC
        score_global[cell, ] = sum(lg.PC * rc.CC, na.rm = T)
      }
    }
  } else {
    stop('Error : Direction of the communication ("in" or "out") must be specified ')
      }

  main <- list(score_global, lr_global)
  return(main)
}
#################
my_ligand.average.RNAseq <-function(db = db,data = data) {
    SYMBOL = rownames(data)
    data = as.data.frame(data)
    x.lg = vector(length = dim(db)[1])
    # for RNAseq expression data
    if (is.null(data$Symbol)) {
      data$Symbol = SYMBOL
    }
    for (mol in seq(1, dim(db)[1])) {
      # In ligand-provider data
      if (db$`Ligand 1`[mol] %in% data$Symbol) {
        if (is.na(db$`Ligand 2`[mol])) {
          x.lg[mol] = mean(data.matrix(data[which(data$Symbol == db$`Ligand 1`[mol]), which(colnames(data) !=
                                                                                              "Symbol")]))
        } else{
          x.lg[mol] = psych::geometric.mean(c(mean(data.matrix(data[which(data$Symbol == db$`Ligand 1`[mol]), which(colnames(data) !=
                                                                                                               "Symbol")])),
                                       mean(data.matrix(data[which(data$Symbol ==
                                                                     db$`Ligand 2`[mol]), which(colnames(data) != "Symbol")]))))
        }
      } else{
        #warning(paste0('Warning: ', db$`Ligand 1`[mol], ' is not found in data matrix'))
        x.lg[mol] = NA
      }
    }
    return(x.lg)
  }
######################
my_receptor.average.RNAseq <-function(db = db,data = data) {
    SYMBOL = rownames(data)
    # check type of data
    x.rc = vector(length = dim(db)[1])
    data = as.data.frame(data)
    # for RNAseq expression data
    if (is.null(data$Symbol))
      data$Symbol = rownames(data)
    for (mol in seq(1, dim(db)[1])) {
      if (db$`Receptor 1`[mol] %in% data$Symbol) {
        # In receptor-provider data
        if (is.na(db$`Receptor 2`[mol]) & is.na(db$`Receptor 3`[mol])) {
          x.rc[mol] = mean(data.matrix(data[which(data$Symbol == db$`Receptor 1`[mol]), which(colnames(data) !=
                                                                                                "Symbol")]))
        } else if (is.na(db$`Receptor 3`[mol])) {
          x.rc[mol] = psych::geometric.mean(c(mean(data.matrix(data[which(data$Symbol == db$`Receptor 1`[mol]), which(colnames(data) !=
                                                                                                                 "Symbol")])),
                                       mean(data.matrix(data[which(data$Symbol ==
                                                                     db$`Receptor 2`[mol]), which(colnames(data) != "Symbol")]))))
        } else if (is.na(db$`Receptor 2`[mol])) {
          x.rc[mol] = psych::geometric.mean(c(mean(data.matrix(data[which(data$Symbol == db$`Receptor 1`[mol]), which(colnames(data) !=
                                                                                                                 "Symbol")])),
                                       mean(data.matrix(data[which(data$Symbol ==
                                                                     db$`Receptor 3`[mol]), which(colnames(data) != "Symbol")]))))
        } else{
          x.rc[mol] = psych::geometric.mean(c(
            mean(data.matrix(data[which(data$Symbol == db$`Receptor 1`[mol]), which(colnames(data) !=
                                                                                      "Symbol")])),
            mean(data.matrix(data[which(data$Symbol ==
                                          db$`Receptor 3`[mol]), which(colnames(data) != "Symbol")])),
            mean(data.matrix(data[which(data$Symbol ==
                                          db$`Receptor 2`[mol]), which(colnames(data) != "Symbol")]))
          ))
        }
      } else{
        #warning(paste0('Warning: ', db$`Ligand 1`[mol], ' is not found in data matrix'))
        x.rc[mol] = NA
      }
    }
    return(x.rc)
  }

###use other databases
other_all_cells_score <- function(cc_cell,meta,data,n,database,direction){

  cc <- meta[meta$Class==cc_cell,]
  cc_data <- data[,rownames(cc),drop=F]
  cc_data= other_gene_scaling(data = cc_data, n= n, db = database)
  cc_data$Symbol=rownames(cc_data)
  
  pc <- meta[meta$Class!=cc_cell,]
  pc_data <- data[,rownames(pc),drop=F]
  pc_data= other_gene_scaling(data = pc_data, n= n, db = database)
  pc_data$Symbol=rownames(pc_data)
  
  my_selection <- unique(as.matrix(pc$Class))
  score_computation= other_icellnet.score(direction=direction, PC.data=pc_data, CC.data= cc_data,PC.target = pc, PC=my_selection, db = database)
  score <-score_computation[[1]]
  
  score <- score[score[,1]!=0&score[,1]!="NaN",,drop=F]
  if(nrow(score)==0){
    return(list(score2 = NULL,lr_score = NULL))
  }
  
  if(direction=="out"){
    score2 <- data.frame(cell_from=cc_cell,cell_to=rownames(score),score=as.numeric(score),stringsAsFactors=F)
  }else if(direction=="in"){
    score2 <- data.frame(cell_from=rownames(score),cell_to=cc_cell,score=as.numeric(score),stringsAsFactors=F)
  }
  
  lr <- score_computation[[2]]
  lr <- lr[,rownames(score),drop=F]
  is_na <- apply(lr,1,function(x){sum(is.na(x))==length(x)})
  lr2 <- lr[!is_na,,drop=F]
  
  if(is.null(lr2)){
    return(list(score2 = NULL,lr_score = NULL))
  }
  
  LR <- as.matrix(rownames(lr2))
  LR_family <- database[match(LR,database[,"pair"]),"classification"]

  ligand <- apply(LR,1,function(x) unlist(strsplit(x,"_"))[1])
  receptor <- apply(LR,1,function(x) unlist(strsplit(x,"_"))[2])
  if(direction == 'out'){
      lr_score <- data.frame(ligand = rep(ligand,ncol(lr2)),receptor = rep(receptor,ncol(lr2)),cell_from = cc_cell,
                             cell_to = rep(colnames(lr2),each = nrow(lr2)),family = LR_family,score = as.numeric(as.matrix(lr2)),stringsAsFactors=F)
  }else if(direction == 'in'){
      lr_score <- data.frame(ligand = rep(ligand,ncol(lr2)),receptor = rep(receptor,ncol(lr2)),cell_from = rep(colnames(lr2),each = nrow(lr2)),
                             cell_to = cc_cell,family = LR_family,score = as.numeric(as.matrix(lr2)),stringsAsFactors=F)
  }
  
  return(list(score2 = score2,lr_score = lr_score))
}
#######################
other_gene_scaling <- function(data, n, db) {
  dim2 = dim(data)[2]
  dim1 = dim(data)[1]
  int = intersect(as.matrix(rownames(data)), as.matrix(db[, 1:2]))
  data = data[which(rownames(data) %in% int),,drop=F]
  
  if(ncol(data)==1){
    for(i in 1:length(int)){
      
      max <- data[i,]
      
      if(max > 0){
        data[i,] = data[i,] /max  * 10
        data[i, which(data[i, ] > 10)] <- 10
      }else (data [i,] <- 0)
      
    }
  }else{
    for (i in 1:length(int)) {
      progressBar(i, max = length(int))
      sorted = sort(data[i,], decreasing = TRUE)
      
      if(n > ncol(data)){
        max = sum(sorted[1, 1:ncol(data)])/ncol(data)
      }else{
        max = sum(sorted[1, 1:n])/n
      }
      
      if (max > 0) {
        data[i,] = data[i,] /max  * 10
        data[i, which(data[i, ] > 10)] <- 10
      }else (data [i,] <- 0)
      
    }
  }
  
  data = (data[complete.cases(data),,drop=F])
  return(data)
}
#######################
other_icellnet.score = function(direction = c("out", "in"),
                          CC.data = CC.data,
                          PC.data = PC.data,
                          PC = PC,
                          PC.target = PC.target,
                          db = db)
  {
  if(dim(CC.data)[2]==1| dim(PC.data)[2] ==1){note("Check that PC.data and/or CC.data contains rownames. Ignore this note if this is the case")}
  #Creation of output matrix
  score_global = matrix(nrow = length(PC), ncol = 1)
  rownames(score_global) = PC
  colnames(score_global) = direction

  lr_global = matrix(nrow = dim(db)[1] , ncol = length(PC))
  colnames(lr_global) = PC
  rownames(lr_global) = db$pair

  #Compute the score
  if (direction == "out") {
    for (cell in PC) {
      # check all cell types defined in PC exists in PC.data and select PC.target$ID
      cell.IDs=PC.target$ID[grepl(cell, PC.target$Class)]
      if (length(cell.IDs)==0){
        warning (paste0( "Cell type ", cell, " is not found in the PC.target file"))
        score_global[cell, ] = "NaN"
      }else{
          rc.PC = other_receptor_average(db = db,data = as.data.frame(PC.data[, cell.IDs], row.names = rownames(PC.data)))
          lg.CC = other_ligand_average(db = db,data = CC.data)
          lr_global[, cell] = lg.CC * rc.PC
          score_global[cell, ] = sum(lg.CC * rc.PC, na.rm = T)
      }
    }

  } else if (direction == "in") {
    for (cell in PC) {
      # check all cell types defined in PC exists in PC.data and select PC.target$ID
      cell.IDs=PC.target$ID[grepl(cell, PC.target$Class)]
      if (length(cell.IDs)==0){
        warning (paste0( "Cell type ", cell, " is not found in the PC.target file"))
        score_global[cell, ] = "NaN"
      } else {
        lg.PC = other_ligand_average(db = db,data = as.data.frame(PC.data[, cell.IDs ], row.names = rownames(PC.data)))
        rc.CC = other_receptor_average(db = db,data = CC.data)
        lr_global[, cell] = lg.PC * rc.CC
        score_global[cell, ] = sum(lg.PC * rc.CC, na.rm = T)
      }
    }
  } else {
    stop('Error : Direction of the communication ("in" or "out") must be specified ')
      }

  main <- list(score_global, lr_global)
  return(main)
}

###################
other_receptor_average <- function(db,data){
  
  x.rc = matrix(NA,nrow = dim(db)[1])
  data = as.data.frame(data)
  
  inter <- intersect(rownames(data),db$receptor)
  if(ncol(data[,which(colnames(data)!="Symbol"),drop=F])==1){
    inter_mean<- data[inter,which(colnames(data)!="Symbol")]
  }else{
    inter_mean<- rowMeans(data[inter,which(colnames(data)!="Symbol")])
  }
  
  for(i in 1:length(inter)){
    index <- which(db$receptor %in% inter[i])
    x.rc[index] <- inter_mean[i]
  }
  x.rc <- as.vector(x.rc)
  return(x.rc)
}
###################
other_ligand_average <- function(db,data){
  
  x.rc = matrix(NA,nrow = dim(db)[1])
  data = as.data.frame(data)
  
  inter <- intersect(rownames(data),db$ligand)
  
  if(ncol(data[,which(colnames(data)!="Symbol"),drop=F])==1){
    inter_mean<- data[inter,which(colnames(data)!="Symbol")]
  }else{
    inter_mean<- rowMeans(data[inter,which(colnames(data)!="Symbol")])
  }
  
  for(i in 1:length(inter)){
    index <- which(db$ligand %in% inter[i])
    x.rc[index] <- inter_mean[i]
  }
  x.rc <- as.vector(x.rc)
  return(x.rc)
}
#######################plot riverplot
plot_riverplot <- function(plot_data){
  
  plot_data <- as.matrix(plot_data) 
  plot_data[,"cell_from"] <- paste(plot_data[,"cell_from"],"(from)",sep='')
  plot_data[,"cell_to"] <- paste(plot_data[,"cell_to"],"(to)",sep='')
  plot_data[,"ligand"] <- paste(plot_data[,"ligand"],"(ligand)",sep='')
  plot_data[,"receptor"] <- paste(plot_data[,"receptor"],"(receptor)",sep='')
  
  node <- matrix(unique(c(plot_data[,"cell_from"], plot_data[,"ligand"],plot_data[,"family"],
                          plot_data[,"receptor"],plot_data[,"cell_to"])),ncol=1)
  colnames(node) <- "name"
  node <- as.data.frame(node)
  
  index <- rep(1,times = nrow(plot_data))
  plot_data <- cbind(index,plot_data)
  
  data_cell_from <- aggregate(plot_data[,"index"],list(plot_data[,c("cell_from")],plot_data[,c("ligand")]),length)
  colnames(data_cell_from) <- c("source","target","value")
  
  data_ligand <- aggregate(plot_data[,"index"],list(plot_data[,c("ligand")],plot_data[,c("family")]),length)
  colnames(data_ligand) <- c("source","target","value")
  
  data_receptor <- aggregate(plot_data[,"index"],list(plot_data[,c("family")],plot_data[,c("receptor")]),length)
  colnames(data_receptor) <- c("source","target","value")
  
  data_to <- aggregate(plot_data[,"index"],list(plot_data[,c("receptor")],plot_data[,c("cell_to")]),length)
  colnames(data_to) <- c("source","target","value")
  
  link <- Reduce(rbind,list(data_cell_from,data_ligand,data_receptor,data_to))
  
  output <- list(nodes=node,links=link)
  jsoncars <- toJSON(output, pretty=TRUE)
}

