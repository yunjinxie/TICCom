TCGA_regulators <- function(IDPath=NULL,cancer = NULL,immune_cell = NULL,referencePath,outpath){
  
  suppressMessages({library(jsonlite)})
  
  stopMessage <- NULL
  #ID
  ID_new <- tail(unlist(strsplit(IDPath,"/")),n=1)
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
    l <- load(paste(referencePath,"detail.RData",sep='/'))
    pairs <- eval(parse(text = l))
    pairs <- as.matrix(pairs)
    
    human <- pairs[pairs[,"Species"]=="human",,drop=F]
    human <- as.matrix(human)
    if(sum(cancer %in% "empty")==0){
      cancer <- unlist(strsplit(cancer,","))
      human2 <- human[human[,"Cancer"] %in% cancer,,drop=F]
    }else{
      human2 <- human
    }
    
    if(sum(immune_cell %in% "empty")==0){
      immune_cell <- unlist(strsplit(immune_cell,","))
      human3 <- human2[human2[,"Cell"] %in% immune_cell|human2[,"Cell1"] %in% immune_cell,,drop=F]
    }else{
      human3 <- human2
    }
    
    if(nrow(human3)==0){
      stopMessage <- "no interactions."
    }else{
      human3 <- unique(human3[,c("Gene1Name","Gene1Id","Gene2Name","Gene2Id"),drop=F])
      
      l2 <- load(paste(referencePath,"starBaseV3_hg19_miRNA_gene_programNum5.RData",sep='/'))
      starbase <- eval(parse(text = l2))
      starbase <- as.matrix(starbase)
      all_miRNA <- unique(starbase[,"miRNAname2"])
      ###TF
      l3 <- load(paste(referencePath,"Low_experimental_proved_TF_all_ensembl.RData",sep='/'))
      TF <- eval(parse(text = l3))
      TF <- as.matrix(TF)
      all_TF <- unique(TF[,"TF"])
      ################function()
      result <- lapply(1:nrow(human3),function(x) plot_one_pair(LR1 = human3[x,],starbase = starbase,all_miRNA = all_miRNA,TF = TF,all_TF = all_TF,outpath = outpath))
      result <- do.call(rbind,result)
      result <- result[as.numeric(result[,11])!=0,,drop=F]
      
      if(nrow(result)==0){
        stopMessage <- "no results."
      }else{
        colnames(result) <- c("Gene1Name","Gene1Id","MG(num)_1","TG(num)_1",
                              "Gene2Name","Gene2Id","MG(num)_2","TG(num)_2",
                              "MG(Total)","TG(Total)","Total num")
        ####============================================output table
        result_df <- as.data.frame(result)
        rownames(result_df) <- NULL
        output <- list(data = result_df)
        jsoncars <- toJSON(output, pretty=TRUE)
        cat(jsoncars, file = paste(outpath,'TCGA_regulators_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
        ####=============================================plot circle
        all_gene1 <- unique(human3[,"Gene1Id"])
        all_gene2 <- unique(human3[,"Gene2Id"])
        gene1_miRNA <- unique(starbase[starbase[,"geneID"] %in% all_gene1,c("miRNAname2","geneName"),drop=F])
        gene2_miRNA <- unique(starbase[starbase[,"geneID"] %in% all_gene2,c("miRNAname2","geneName"),drop=F])
        
        union_res <- unique(rbind(unname(human3[,c("Gene1Name","Gene2Name")]),
                                  unname(gene1_miRNA),unname(gene2_miRNA))) 
        
        gene1_TF <- unique(TF[TF[,"Ensembl"] %in% all_gene1,c("TF","Gene"),drop=F])
        gene2_TF <- unique(TF[TF[,"Ensembl"] %in% all_gene2,c("TF","Gene"),drop=F])
        
        TF_union_res <- unique(rbind(unname(gene1_TF),unname(gene2_TF))) 
        
        links <- rbind(union_res,TF_union_res)
        
        links_sort <- t(apply(links,1,function(x) sort(x)))
        links_unique <- links[!duplicated(links_sort),]
        node <- union(links_unique[,1],links_unique[,2])
        
        write.table(links_unique,paste(outpath,"TCGA_regulators_circle_download.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = c("node1","node2"))
        
        degree <- table(c(links_unique[,1],links_unique[,2]))
        degree2 <- as.matrix(degree)
        degree3 <- degree2[node,]
        
        node_attribute <- rep("Partner1",times=length(node))
        node_attribute[node %in% setdiff(all_miRNA,c(human3["Gene1Name"],human3["Gene2Name"]))] <- "miRNA"
        node_attribute[node %in% setdiff(all_TF,c(human3["Gene1Name"],human3["Gene2Name"]))] <- "TF"
        node_attribute[node %in% human3[,"Gene2Name"]] <- "Partner2"
        
        ID <- 0:(length(node)-1)
        
        symbolSize <- rep(30,times = length(node))
        symbolSize[node_attribute %in% c("Partner1","Partner2")] <- 50
        
        node2 <- cbind(ID,node,symbolSize,degree3,node_attribute)
        colnames(node2) <- c("id","name","symbolSize","value","category")
        
        category <- t(cbind(c(0,"Partner1"),c(1,"Partner2"),c(2,"miRNA"),c(3,"TF")))
        for(i in 1:nrow(category)){
          node2[node2[,"category"]==category[i,2],"category"] <- category[i,1]
        }
        
        category2 <- category[,2,drop=F]
        colnames(category2) <- "name"
        
        links_final <- cbind(node2[match(links_unique[,1],node2[,"name"]),"id"],
                             node2[match(links_unique[,2],node2[,"name"]),"id"])
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
        
        cat(jsoncars, file = paste(outpath,"TCGA_regulators_circle.json",sep='/'), fill = FALSE, labels = NULL, append = FALSE) 
      }
    }
    if(is.null(stopMessage)){
      download_table="TCGA_regulators_re.txt"
      Res_circle_file="TCGA_regulators_circle.json"
      resut_merge=paste0("{",
                         '"download_table" :','"',download_table,'",',
                         '"Res_circle_file" :','"',Res_circle_file,'",',
                         '"error_attention" :','"no',
                         '"}')
      write.table(resut_merge,paste(outpath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
      
      Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
      Task_ID_new <- Task_ID
      Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "TCGA_regulators+success"
      write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
      
      result_message <- paste0("{",'"error_attention" :','"success"',"}")
      return(result_message)
    }else{
      write.table(paste0("{",'"error_attention" :','"',stopMessage,'"}'),paste(outpath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
      
      Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
      Task_ID_new <- Task_ID
      Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "dead"
      write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
      
      result_message <- paste0("{",'"error_attention" :','"',stopMessage,'"}')
      return(result_message)
    }
  },error = function(e){
    stopMessage <- "task quit."
    write.table(paste0("{",'"error_attention" :','"',stopMessage,'"}'),paste(outpath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
    
    Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
    Task_ID_new <- Task_ID
    Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "dead"
    write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
    
    result_message <- paste0("{",'"error_attention" :','"',stopMessage,'"}')
    return(result_message)
  })
}
#########function
plot_one_pair <- function(LR1,starbase,all_miRNA,TF,all_TF,outpath){
  
  gene1_miRNA <- unique(starbase[starbase[,"geneID"] %in% LR1["Gene1Id"],c("miRNAname2","geneName"),drop=F])
  gene2_miRNA <- unique(starbase[starbase[,"geneID"] %in% LR1["Gene2Id"],c("miRNAname2","geneName"),drop=F])

  if(nrow(gene1_miRNA)==0){
    gene1_miRNA <- NULL
    gene1_num1 <- 0
  }else{
    gene1_miRNA <- gene1_miRNA[which(!gene1_miRNA[,"miRNAname2"] %in% LR1["Gene2Name"]),,drop=F]

    if(nrow(gene1_miRNA)==0){
      gene1_miRNA <- NULL
      gene1_num1 <- 0
    }else{
      gene1_num1 <- nrow(gene1_miRNA)
    }
  }
 
  if(nrow(gene2_miRNA)==0){
    gene2_miRNA <- NULL
    gene2_num1 <- 0
  }else{
    gene2_miRNA <- gene2_miRNA[which(!gene2_miRNA[,"miRNAname2"] %in% LR1["Gene1Name"]),,drop=F]
    
    if(nrow(gene2_miRNA)==0){
      gene2_miRNA <- NULL
      gene2_num1 <- 0
    }else{
      gene2_num1 <- nrow(gene2_miRNA)
    }

  }
 
  union_res <- unique(rbind(unname(LR1[c("Gene1Name","Gene2Name")]),
                            unname(gene1_miRNA),unname(gene2_miRNA))) 
  
  ################################################
  gene1_TF <- unique(TF[TF[,"Ensembl"] %in% LR1["Gene1Id"],c("TF","Gene"),drop=F])
  gene2_TF <- unique(TF[TF[,"Ensembl"] %in% LR1["Gene2Id"],c("TF","Gene"),drop=F])

  if(nrow(gene1_TF)==0){
    gene1_TF <- NULL
    gene1_TF_num1 <- 0
  }else{
    gene1_TF <- gene1_TF[which(!gene1_TF[,"TF"] %in% LR1["Gene2Name"]),,drop=F]
    
    if(nrow(gene1_TF)==0){
      gene1_TF <- NULL
      gene1_TF_num1 <- 0
    }else{
      gene1_TF_num1 <- nrow(gene1_TF)
    }
  }
  
  if(nrow(gene2_TF)==0){
    gene2_TF <- NULL
    gene2_TF_num1 <- 0
  }else{
    gene2_TF <- gene2_TF[which(!gene2_TF[,"TF"] %in% LR1["Gene1Name"]),,drop=F]
    
    if(nrow(gene2_TF)==0){
      gene2_TF <- NULL
      gene2_TF_num1 <- 0
    }else{
      gene2_TF_num1 <- nrow(gene2_TF)
    }
  }
  TF_union_res <- unique(rbind(unname(gene1_TF),unname(gene2_TF))) 
  
  links <- rbind(union_res,TF_union_res)
  
  if(nrow(links)==1){
    result <- cbind(LR1["Gene1Name"],LR1["Gene1Id"],gene1_num1,gene1_TF_num1,
                    LR1["Gene2Name"],LR1["Gene2Id"],gene2_num1,gene2_TF_num1,
                    0,0,0)
    
    return(result)
  }
  
  links_sort <- t(apply(links,1,function(x) sort(x)))
  links_unique <- links[!duplicated(links_sort),]
  node <- union(links_unique[,1],links_unique[,2])
#####################
  degree <- table(c(links_unique[,1],links_unique[,2]))
  degree2 <- as.matrix(degree)
  degree3 <- degree2[node,]

  node_attribute <- rep("Partner1",times=length(node))
  node_attribute[node %in% setdiff(all_miRNA,c(LR1["Gene1Name"],LR1["Gene2Name"]))] <- "miRNA"
  node_attribute[node %in% setdiff(all_TF,c(LR1["Gene1Name"],LR1["Gene2Name"]))] <- "TF"
  node_attribute[node %in% LR1["Gene2Name"]] <- "Partner2"

  ID <- 0:(length(node)-1)

  symbolSize <- rep(30,times = length(node))
  symbolSize[node_attribute %in% c("Partner1","Partner2")] <- 50

  node2 <- cbind(ID,node,symbolSize,degree3,node_attribute)
  colnames(node2) <- c("id","name","symbolSize","value","category")

  category <- t(cbind(c(0,"Partner1"),c(1,"Partner2"),c(2,"miRNA"),c(3,"TF")))
  for(i in 1:nrow(category)){
    node2[node2[,"category"]==category[i,2],"category"] <- category[i,1]
  }

  category2 <- category[,2,drop=F]
  colnames(category2) <- "name"

  links_final <- cbind(node2[match(links_unique[,1],node2[,"name"]),"id"],
                       node2[match(links_unique[,2],node2[,"name"]),"id"])
  colnames(links_final) <- c("source","target")
######################output
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

  file_name <- paste(LR1["Gene1Name"],"_",LR1["Gene2Name"],".json",sep='')
  cat(jsoncars, file = paste(outpath,file_name,sep='/'), fill = FALSE, labels = NULL, append = FALSE)
  
  if(is.null(TF_union_res)){
    result <- cbind(LR1["Gene1Name"],LR1["Gene1Id"],gene1_num1,gene1_TF_num1,
                    LR1["Gene2Name"],LR1["Gene2Id"],gene2_num1,gene2_TF_num1,
                    length(which(node_attribute=="miRNA")),0,length(which(node_attribute=="miRNA")))
  }else{
    result <- cbind(LR1["Gene1Name"],LR1["Gene1Id"],gene1_num1,gene1_TF_num1,
                    LR1["Gene2Name"],LR1["Gene2Id"],gene2_num1,gene2_TF_num1,
                    length(which(node_attribute=="miRNA")),length(which(node_attribute=="TF")),
                    length(which(node_attribute %in% c("miRNA","TF"))))
  }
  
  return(result)
}
