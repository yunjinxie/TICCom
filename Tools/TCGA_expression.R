TCGA_expression <- function(IDPath=NULL,expPath=NULL,gene_listPath=NULL,cancer=NULL,referencePath=NULL,outpath){
  
  stopMessage <- NULL
  
  ###
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
    if(gene_listPath!="empty"){
      gene_list <- read.table(gene_listPath,sep='\t',header = F,as.is = T,fill = T,strip.white = T)
      gene_list <- as.matrix(gene_list)
      ############ABL1-CSF1;LCAM-SDCBP
      gene_list2 <- Reduce(rbind,strsplit(gene_list,"_|:|&"))
      
      if(is.null(ncol(gene_list2))){
        if(length(gene_list2)==2){
          gene_list2 <- matrix(gene_list2,nrow=1,ncol=2)
          colnames(gene_list2) <- c("UserID1","UserID2")
          gene1 <- ID_transfer(gene_list = gene_list2[,1],ID_transfer_Path = referencePath)
          gene2 <- ID_transfer(gene_list = gene_list2[,2],ID_transfer_Path = referencePath)
          if(nrow(gene1)==0|nrow(gene2)==0){
            stopMessage <- "no matched gene symbols."
          }else{
            colnames(gene1) <- c("UserID1","Gene1Name","Gene1Id")
            gene_list3 <- merge(gene_list2,gene1,by='UserID1')
            colnames(gene2) <- c("UserID2","Gene2Name","Gene2Id")
            gene_list4 <- merge(gene_list3,gene2,by='UserID2')
            if(nrow(gene_list4)==0){
              stopMessage <- "no matched gene symbols."
            }else{
              pairs2 <- gene_list4[,c("Gene1Name","Gene1Id","Gene2Name","Gene2Id"),drop=F]
            }
          }
        }else{
          stopMessage <- "too few gene interactions input."
        }
      }else{
        colnames(gene_list2) <- c("UserID1","UserID2")
        gene1 <- ID_transfer(gene_list = gene_list2[,1],ID_transfer_Path = referencePath)
        gene2 <- ID_transfer(gene_list = gene_list2[,2],ID_transfer_Path = referencePath)
        if(nrow(gene1)==0|nrow(gene2)==0){
          stopMessage <- "no matched gene symbols."
        }else{
          colnames(gene1) <- c("UserID1","Gene1Name","Gene1Id")
          gene_list3 <- merge(gene_list2,gene1,by='UserID1')
          colnames(gene2) <- c("UserID2","Gene2Name","Gene2Id")
          gene_list4 <- merge(gene_list3,gene2,by='UserID2')
          if(nrow(gene_list4)==0){
            stopMessage <- "no matched gene symbols."
          }else{
            pairs2 <- gene_list4[,c("Gene1Name","Gene1Id","Gene2Name","Gene2Id"),drop=F]
          }
        }
      }
    }else{
      l <- load(paste(referencePath,"immune_cancer_genePairs.RData",sep='/'))
      pairs <- eval(parse(text = l))
      pairs2 <- pairs
    }
    
    if(is.null(stopMessage)){
      pairs3 <- as.matrix(pairs2)
      colnames(pairs3) <- NULL
      gene_name <- unique(rbind(pairs3[,c(1,2)],pairs3[,c(3,4)]))
      
      if(sum(cancer %in% "empty")==0){
        cancer <- unlist(strsplit(cancer,","))
        #run the method 
        TCGAFile <- paste(referencePath,"/",cancer,"_cancer_genes_expression_30_05_log2.RData",sep='')
        exp <- lapply(TCGAFile,function(x){
          expl <- load(x)
          exp_data <- eval(parse(text = expl))
          return(exp_data)
        })
        
        exp_result <- lapply(1:length(cancer),function(x) Expression_analysis(exp = exp[[x]],pairs = pairs2,gene_name = gene_name,cancer = cancer[x]))
        exp_result <- do.call(rbind,exp_result)
        
        if(is.null(exp_result)){
          stopMessage <- "no results."
        }else{
          exp_result[,"Exp_mean"] <- round(as.numeric(exp_result[,"Exp_mean"]),4)
          write.table(exp_result[,c(1:4),drop=F],paste(outpath,"TCGA_expression_heatmap_download.txt",sep='/'),sep='\t',quote = F,row.names = F)
          #pheatmap
          heatmap_res <- plot_heatmap_data(exp_result = exp_result,choose_cancer = TRUE) 
          all_data <- heatmap_res[[1]]
          plot_data <- heatmap_res[[2]]
          write.table(all_data,paste(outpath,"TCGA_expression_heatmap_all.txt",sep='/'),sep='\t',row.names = F,col.names = F,quote = F)
          write.table(plot_data,paste(outpath,"TCGA_expression_heatmap.txt",sep='/'),sep='\t',row.names = F,col.names = F,quote = F)
          #boxplot
          boxplot_res <- plot_box_data(all = exp_result,choose_cancer = TRUE)
          write.table(boxplot_res,paste(outpath,"TCGA_expression_boxplot.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
        }
      }else if(expPath!="empty"){
        exp0 <- read.table(expPath,sep='\t',header=T,as.is = T,fill=T,strip.white = T,check.names=F)
        if((class(exp0[,1]) == "character") | (class(exp0[,1]) == "factor")){
          rownames(exp0) <- as.matrix(exp0[,1])
          exp0 <- exp0[,-1]
        }
        gene <- rownames(exp0)
        gene2 <- ID_transfer(gene_list = gene,ID_transfer_Path = referencePath)
        if(nrow(gene2)==0){
          stopMessage <- "no matched gene ENSGs."
        }else{
          gene3 <- gene2[match(rownames(exp0),gene2[,1]),]
          gene3 <- na.omit(gene3)
          frenq <- table(gene3[,3])
          
          if(max(frenq)>1){
            exp_gene <- exp0[gene3[,1],,drop=F]
            exp_gene2 <- cbind(gene3,exp_gene)
            
            uni_exp <- exp_gene2[match(names(frenq[frenq==1]),exp_gene2[,3]),-c(1,2,3),drop=F]
            rownames(uni_exp) <- names(frenq[frenq==1])
            
            more_exp <- exp_gene2[exp_gene2[,3] %in% names(frenq[frenq>1]),,drop=F]
            more_exp2 <- apply(more_exp[,-c(1,2,3)],2,function(x) {
              tapply(x, factor(more_exp[,3]), function(x) mean(as.numeric(x))) 
            })
            exp <- rbind(uni_exp,more_exp2)
          }else{
            exp <- exp0[rownames(exp0) %in% gene3[,1],,drop=F]
            rownames(exp) <- gene3[match(rownames(exp),gene3[,1]),3]
          }
          exp_result <- Expression_analysis(exp = exp,pairs = pairs2,gene_name = gene_name)
          if(is.null(exp_result)){
            stopMessage <- "no results."
          }else{
            exp_result[,"Exp_mean"] <- round(as.numeric(exp_result[,"Exp_mean"]),4)
            write.table(exp_result[,c(1:4),drop=F],paste(outpath,"TCGA_expression_heatmap_download.txt",sep='/'),sep='\t',quote = F,row.names = F)
            #pheatmap
            heatmap_res <- plot_heatmap_data(exp_result = exp_result,choose_cancer = FALSE) 
            all_data <- heatmap_res[[1]]
            plot_data <- heatmap_res[[2]]
            write.table(all_data,paste(outpath,"TCGA_expression_heatmap_all.txt",sep='/'),sep='\t',row.names = F,col.names = F,quote = F)
            write.table(plot_data,paste(outpath,"TCGA_expression_heatmap.txt",sep='/'),sep='\t',row.names = F,col.names = F,quote = F)
            #boxplot
            boxplot_res <- plot_box_data(all = exp_result,choose_cancer = FALSE)
            write.table(boxplot_res,paste(outpath,"TCGA_expression_boxplot.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
          }
        }
      }
    }
    if(is.null(stopMessage)){
      heat_download_table="TCGA_expression_heatmap_download.txt"
      Res_heat_all_file="TCGA_expression_heatmap_all.txt"
      Res_heat_file="TCGA_expression_heatmap.txt"
      Res_boxplot_file="TCGA_expression_boxplot.txt"
      resut_merge=paste0("{",
                         '"heat_download_table" :','"',heat_download_table,'",',
                         '"Res_heat_all_file" :','"',Res_heat_all_file,'",',
                         '"Res_heat_file" :','"',Res_heat_file,'",',
                         '"Res_boxplot_file" :','"',Res_boxplot_file,'",',
                         '"error_attention" :','"no',
                         '"}')
      write.table(resut_merge,paste(outpath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
      
      Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
      Task_ID_new <- Task_ID
      Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "TCGA_expression+success"
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
#################
Expression_analysis <- function(exp,pairs,gene_name,cancer = NULL){
  gene1 <- rownames(exp)[rownames(exp) %in% pairs$Gene1Id]
  gene2 <- rownames(exp)[rownames(exp) %in% pairs$Gene2Id]
  
  index1 <- which(pairs$Gene1Id %in% gene1)
  index2 <- which(pairs$Gene2Id %in% gene2)
  index <- intersect(index1,index2)
  
  if(length(index)==0){
    return(NULL)  
  }
  pairs_gene <- unique(c(pairs[index,]$Gene1Id,pairs[index,]$Gene2Id))
  pairs_gene_name <- gene_name[match(pairs_gene,gene_name[,2]),1]
  
  gene1_exp <- matrix(as.numeric(as.matrix(exp[pairs_gene,,drop=F])),nrow=nrow(exp[pairs_gene,,drop=F]))
  pairs_gene_exp_mean <- rowMeans(gene1_exp)
  #random expression
  r <- sample(1:nrow(exp),length(pairs_gene))
  r_exp <- matrix(as.numeric(as.matrix(exp[r,,drop=F])),nrow=nrow(exp[r,,drop=F]))
  r_exp_mean <- rowMeans(r_exp)
  
  if(!is.null(cancer)==T){
    result <- cbind(cancer,pairs_gene,pairs_gene_name,pairs_gene_exp_mean)
    colnames(result) <- c("Cancer","GeneID","GeneName","Exp_mean")
  }else{
    result <- cbind("-",pairs_gene,pairs_gene_name,pairs_gene_exp_mean)
    colnames(result) <- c("Cancer","GeneID","GeneName","Exp_mean")
  }
  all <- cbind(result,r_exp_mean)
  return(all)
}
#######################################
ID_transfer <- function(gene_list,ID_transfer_Path=NULL){
  l2 <- load(paste(ID_transfer_Path,"ID_transfer_file.RData",sep='/'))
  ID_transfer_data <- eval(parse(text = l2))
  #ID_transfer_data <- read.table(paste(ID_transfer_Path,"ID_transfer_file.txt",sep='/'),sep='\t',header = T,as.is = T,fill=T,strip.white = T,quote = "",check.names = F)
  gene <- gene_list[1]
  
  if(grepl("ENSG",gene[1])){
    gene_names = "ensembl"
  }else if(!grepl(paste("[",paste(toupper(letters),collapse ='|'),"]",sep=''),gene[1])){
    gene_names = "entrez"
  }else{
    gene_names = "symbol"
  }
  
  if(gene_names == "ensembl"){
    res<- unique(ID_transfer_data[ID_transfer_data[,2] %in% gene_list,c(2,3,2)])
  }else if(gene_names == "entrez"){
    res <- unique(ID_transfer_data[ID_transfer_data[,1] %in% gene_list,c(1,3,2)])
  }else{
    res <- unique(ID_transfer_data[ID_transfer_data[,3] %in% gene_list,c(3,3,2)])
  }
  return(res)
}
###############
plot_heatmap_data <- function(exp_result,choose_cancer = TRUE){
  
  all_data <- plot_all_heatmap(plot_exp = exp_result,choose_cancer = choose_cancer)
  if(length(unique(as.matrix(exp_result[,"Cancer"])))==1){
      exp_result <- exp_result[order(as.numeric(exp_result[,"Exp_mean"]),decreasing = T),,drop=F]
      
      if(nrow(exp_result)>= 25){
        plot_exp <- exp_result[1:25,]
      }else{
        plot_exp <- exp_result
      }
  }else{
    gene_count <- aggregate(exp_result[,"Cancer"], list(gene = exp_result[,"GeneName"]), function(x) length(unique(x)))
    gene_count_order <- gene_count[order(gene_count[,2],decreasing = T),,drop=F]
    
    if(nrow(gene_count_order)>=25){
      plot_gene <- gene_count_order[1:25,]
    }else{
      plot_gene <- gene_count_order
    }
    plot_exp <- exp_result[exp_result[,"GeneName"] %in% plot_gene[,1],,drop=F]
  }

  plot_data <- plot_all_heatmap(plot_exp = plot_exp,choose_cancer = choose_cancer)
  return(list(all_data=all_data,plot_data=plot_data))
}
###################output all heatmap
plot_all_heatmap <- function(plot_exp,choose_cancer = TRUE){
  ueg <- unique(plot_exp[,"GeneName"])
  uegy <- cbind(ueg,c(0:(length(ueg)-1)))
  
  if(choose_cancer == TRUE){
    exp_cancers <- unique(plot_exp[,"Cancer"])
    ecx <- cbind(exp_cancers,c(0:(length(exp_cancers)-1)))
  }else{
    exp_cancers <- "User cancer"
    ecx <- cbind("User cancer",c(0:(length(exp_cancers)-1)))
  }
  
  plot_exp2 <- cbind(plot_exp,x=0,y=0)
  ####
  for(i in 1:nrow(uegy)){
    plot_exp2[which(plot_exp2[,"GeneName"] %in% uegy[i,1]),"y"] <- uegy[i,2]
  }
  
  for(i in 1:nrow(ecx)){
    plot_exp2[which(plot_exp2[,"Cancer"] %in% ecx[i,1]),"x"] <- ecx[i,2]
  }
  ####[y,x,value]
  xyAxis0 <- paste("[",plot_exp2[,"y"],",",plot_exp2[,"x"],",",plot_exp2[,"Exp_mean"],"]",sep='')
  xyAxis1 <- paste(xyAxis0,collapse =',')
  xyAxis <- paste("[",xyAxis1,"]",sep = "")
  y0<- paste(ueg,collapse="','")
  y <- paste("['",y0,"']",sep='')
  x0<- paste(exp_cancers,collapse="','")
  x <- paste("['",x0,"']",sep='')
  
  exp_min <- min(as.numeric(plot_exp2[,"Exp_mean"]))
  exp_max <- max(as.numeric(plot_exp2[,"Exp_mean"]))
  
  return(t(c(x,y,xyAxis,exp_min,exp_max)))
}
###########boxplot
plot_box_data <- function(all,choose_cancer = TRUE){
  
  if(choose_cancer == TRUE){
  
    all_cancer <- unique(all[,"Cancer"])
    all_test <- lapply(all_cancer,function(x){
      data <- all[all[,"Cancer"] == x,,drop=F]
      
      if(length(data[,4])==1|length(data[,5])==1){
        p <- 1
      }else if(length(unique(data[,4]))==1&length(unique(data[,5]))==1){
        p <- 1
      }else{
        p <- t.test(as.numeric(data[,4]),as.numeric(data[,5]),alternative = "two.sided",paired = T)$p.value
      }
    
      real_box0 <- paste(round(as.numeric(data[,4]),4),collapse=',')
      real_box <- paste("[",real_box0,"]",sep = '')
      
      random_box0 <- paste(round(as.numeric(data[,5]),4),collapse=',')
      random_box <- paste("[",random_box0,"]",sep = '')
      
      return(list(p_value = p,real_box = real_box,random_box = random_box))
    })
    
    p_value <- lapply(all_test,function(x) x$p_value)
    p_value <- do.call(rbind,p_value)
    
    xx <- paste(all_cancer,"(p=",round(p_value,4),")",sep='')
    all_cancer0 <- paste(xx,collapse = "','")
    all_cancer_s <- paste("'",all_cancer0,"'",sep='')

    real_box <- lapply(all_test,function(x) x$real_box)
    real_box <- do.call(rbind,real_box)
    real_box_s0 <- paste(real_box,collapse=',')
    real_box_s <- paste("[",real_box_s0,"]",sep='')
    
    random_box <- lapply(all_test,function(x) x$random_box)
    random_box <- do.call(rbind,random_box)
    random_box_s0 <- paste(random_box,collapse=',')
    random_box_s <- paste("[",random_box_s0,"]",sep='')
    
    return(t(c(all_cancer_s,real_box_s,random_box_s)))
  }else{
    
    if(length(all[,4])==1|length(all[,5])==1){
      p <- 1
    }else if(length(unique(all[,4]))==1&length(unique(all[,5]))==1){
      p <- 1
    }else{
      p <- t.test(as.numeric(all[,4]),as.numeric(all[,5]),alternative = "two.sided",paired = T)$p.value
    }
    
    cancer_s0 <-"User cancer"
    cancer_s0 <- paste(cancer_s0,"(p=",round(p,4),")",sep='')
    cancer_s <- paste("'",cancer_s0,"'",sep='')
    
    real_box_s0 <- paste(round(as.numeric(all[,4]),4),collapse=',')
    real_box_s <- paste("[[",real_box_s0,"]]",sep='')
    
    random_box_s0 <- paste(round(as.numeric(all[,5]),4),collapse=',')
    random_box_s <- paste("[[",random_box_s0,"]]",sep='')
    return(t(c(cancer_s,real_box_s,random_box_s)))
  }
}
