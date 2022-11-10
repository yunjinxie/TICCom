#coexpression analysis of tumor-immune interactions based TCGA RNA expression profiles 
#or other bulk data
#user can input expression data or gene list
#referencePath 
TCGA_coexpression <- function(IDPath=NULL,expPath=NULL,gene_listPath=NULL,cancer=NULL,
                              referencePath=NULL,pvalue=0.05,pearson=0.3,outpath){
  
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
      if(sum(cancer %in% "empty")==0){
        cancer <- unlist(strsplit(cancer,","))
        #run the method
        TCGAFile <- paste(referencePath,"/",cancer,"_cancer_genes_expression_30_05_log2.RData",sep='')
        exp <- lapply(TCGAFile,function(x){
          expl <- load(x)
          exp_data <- eval(parse(text = expl))
          return(exp_data)
        })
        
        corr_result <- lapply(1:length(cancer),function(x) Coexpression_analysis(exp=exp[[x]],pairs=pairs2,pvalue = pvalue,pearson = pearson,cancer=cancer[x],outpath = outpath))
        
        sig <- lapply(corr_result,function(x) x$sig)
        sig <- do.call(rbind,sig)
        
        all <- lapply(corr_result,function(x) x$all)
        all <- do.call(rbind,all)
        all <- na.omit(all)
        
        if(is.null(sig)){
          stopMessage <- "no results."
        }else{
          if(nrow(all)==0){
            stopMessage <- "no results."
          }else{
            ####=====================================================boxplot
            box_res <- plot_box_data(all = all,choose_cancer = TRUE)
            write.table(box_res,paste(outpath,"TCGA_coexpression_boxplot.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
            ####=====================================================heatmap
            sig[,"Pearson"] <- round(as.numeric(sig[,"Pearson"]),4)
            sig[,"Pvalue"] <- signif(as.numeric(sig[,"Pvalue"]),4)
            write.table(sig,paste(outpath,"TCGA_coexpression_heatmap_download.txt",sep='/'),sep='\t',quote = F,row.names = F)
            heatmap_res <- plot_heatmap_data(corr_result = sig,choose_cancer = TRUE)
            all_data <- heatmap_res[[1]]
            plot_data <- heatmap_res[[2]]
            write.table(all_data,paste(outpath,"TCGA_coexpression_heatmap_all.txt",sep='/'),sep='\t',row.names = F,col.names = F,quote = F)
            write.table(plot_data,paste(outpath,"TCGA_coexpression_heatmap.txt",sep='/'),sep='\t',row.names = F,col.names = F,quote = F)
            ####=====================================================output table
            sig_df <- as.data.frame(sig)
            rownames(sig_df) <- NULL
            output <- list(data = sig_df)
            jsoncars <- toJSON(output, pretty=TRUE)
            cat(jsoncars, file = paste(outpath,'TCGA_coexpression_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
          }
        }
      }else if(expPath!="empty"){
        exp0 <- read.table(expPath,sep='\t',header=T,as.is = T,fill = T,strip.white = T,check.names=F)
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
          gene_inter <- intersect(rownames(exp),union(pairs2$Gene1Id,pairs2$Gene2Id))
          if(length(gene_inter)==0){
            stopMessage <- "on interactions in expression."
          }else{
            #run the method 
            corr_result <- Coexpression_analysis(exp=exp,pairs=pairs2,pvalue = pvalue,pearson = pearson,outpath = outpath)
            sig <- corr_result$sig
            all <- corr_result$all
            all <- na.omit(all)
            
            if(is.null(sig)){
              stopMessage <- "no results."
            }else{
              if(nrow(all)==0){
                stopMessage <- "no results."
              }else{
                ####=========================================boxplot
                box_res <- plot_box_data(all = all,choose_cancer = FALSE)
                write.table(box_res,paste(outpath,"TCGA_coexpression_boxplot.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
                ####==========================================heatmap
                sig[,"Pearson"] <- round(as.numeric(sig[,"Pearson"]),4)
                sig[,"Pvalue"] <- signif(as.numeric(sig[,"Pvalue"]),4)
                write.table(sig,paste(outpath,"TCGA_coexpression_heatmap_download.txt",sep='/'),sep='\t',quote = F,row.names = F)
                
                heatmap_res <- plot_heatmap_data(corr_result = sig,choose_cancer = FALSE)
                all_data <- heatmap_res[[1]]
                plot_data <- heatmap_res[[2]]
                write.table(all_data,paste(outpath,"TCGA_coexpression_heatmap_all.txt",sep='/'),sep='\t',row.names = F,col.names = F,quote = F)
                write.table(plot_data,paste(outpath,"TCGA_coexpression_heatmap.txt",sep='/'),sep='\t',row.names = F,col.names = F,quote = F)
                ####========================================output table
                sig_df <- as.data.frame(sig)
                output <- list(data = sig_df)
                jsoncars <- toJSON(output, pretty=TRUE)
                cat(jsoncars, file = paste(outpath,'TCGA_coexpression_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
              }
            }
          }
        }
      }else{
        stopMessage <- "wrong choiese."
      }
    }
    if(is.null(stopMessage)){
      download_table="TCGA_coexpression_re.txt"
      heat_download_table="TCGA_coexpression_heatmap_download.txt"
      Res_heat_all_file="TCGA_coexpression_heatmap_all.txt"
      Res_heat_file="TCGA_coexpression_heatmap.txt"
      Res_boxplot_file="TCGA_coexpression_boxplot.txt"
      resut_merge=paste0("{",
                         '"download_table" :','"',download_table,'",',
                         '"heat_download_table" :','"',heat_download_table,'",',
                         '"Res_heat_all_file" :','"',Res_heat_all_file,'",',
                         '"Res_heat_file" :','"',Res_heat_file,'",',
                         '"Res_boxplot_file" :','"',Res_boxplot_file,'",',
                         '"error_attention" :','"no',
                         '"}')
      write.table(resut_merge,paste(outpath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
      
      Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
      Task_ID_new <- Task_ID
      Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "TCGA_coexpression+success"
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
#######################################
Coexpression_analysis <- function(exp,pairs,cancer=NULL,pvalue,pearson,outpath){
  gene1 <- rownames(exp)[rownames(exp) %in% pairs$Gene1Id]
  gene2 <- rownames(exp)[rownames(exp) %in% pairs$Gene2Id]
  
  index1 <- which(pairs$Gene1Id %in% gene1)
  index2 <- which(pairs$Gene2Id %in% gene2)
  index <- intersect(index1,index2)
  
  if(length(index)==0){
    return(list(sig = NULL,all = NULL))  
  }
  
  pairs <- unique(pairs[index,c(1:4),drop=F])
  nlen <- nrow(pairs)
####real
  result <- lapply(1:nlen,function(x){

    gene1_exp <- as.numeric(exp[pairs[x,,drop=F]$Gene1Id,])
    gene2_exp <- as.numeric(exp[pairs[x,,drop=F]$Gene2Id,])
    
    if(length(gene1_exp)==1|length(gene2_exp)==1){
      p_value <- NA
      evaluate <- NA
    }else if(length(unique(gene1_exp))==1|length(unique(gene2_exp))==1){
      p_value <- NA
      evaluate <- NA
    }else{
      p_value <- cor.test(gene1_exp,gene2_exp)$p.value
      evaluate <- cor(gene1_exp,gene2_exp)
    }

    if(!is.null(cancer)==T){
      res <- cbind(as.matrix(pairs[x,c(1:4),drop=F]),cancer,evaluate,p_value)
      colnames(res) <- c("Gene1Name","Gene1Id","Gene2Name","Gene2Id","Cancer","Pearson","Pvalue")
      return(res)
    }else{
      res <- cbind(as.matrix(pairs[x,c(1:4),drop=F]),"-",evaluate,p_value)
      colnames(res) <- c("Gene1Name","Gene1Id","Gene2Name","Gene2Id","Cancer","Pearson","Pvalue")
      return(res)
    }
  })
  ########
  result2 <- do.call(rbind,result)
  result2 <- na.omit(result2)
  
  if(nrow(result2)==0){
    return(list(sig = NULL,all = NULL)) 
  }
  #random coexpression
  r1 <- sample(1:nrow(exp),nrow(result2))
  r2 <- sample(1:nrow(exp),nrow(result2))
  
  r12 <- cbind(r1,r2)
  
  r_result <- t(apply(r12,1,function(x){
    r1_exp <- as.numeric(exp[x[1],])
    r2_exp <- as.numeric(exp[x[2],])
    
    if(length(r1_exp)==1|length(r2_exp)==1){
      p_value <- NA
      evaluate <- NA
    }else if(length(unique(r1_exp))==1|length(unique(r2_exp))==1){
      p_value <- NA
      evaluate <- NA
    }else{
      p_value <- cor.test(r1_exp,r2_exp)$p.value
      evaluate <- cor(r1_exp,r2_exp)
    }
  
    r_res <- cbind(evaluate,p_value)
    return(r_res)
  }))
  
  result2_sig <- result2[as.numeric(result2[,"Pvalue"]) < as.numeric(pvalue)&
                           abs(as.numeric(result2[,"Pearson"])) > as.numeric(pearson),,drop = F]
  
  if(nrow(result2_sig)==0){
    return(list(sig = NULL,all = NULL)) 
  }else{
    p <- apply(result2_sig,1,function(x){
      gene1 <- x[2]
      gene2 <- x[4]
      gene1_exp <- as.numeric(exp[gene1,])
      gene2_exp <- as.numeric(exp[gene2,])
      
      if(!is.null(cancer)==T){
        ll <- lm(gene1_exp~gene2_exp)
        file_name <- paste(cancer,x[1],x[3],sep="_")
        
        png(file = paste(outpath,"/coexpression_",file_name,".png",sep=''),units= "in",width = 5,height = 5,res = 300)
        pp <- plot(gene2_exp,gene1_exp,xlab = x[1],ylab=x[3],
             main = paste("Pearson coefficient of the expression in ",cancer,
                          "\n R-square: ",round(summary(ll)$r.squared,4),
                          "\n Pearson: ",round(as.numeric(x[6]),4)," p-value: ",signif(as.numeric(x[7]),4),sep=''),pch=15,col="orange")
        lines(gene2_exp,fitted(ll),col="red")
        print(pp)
        dev.off()
      }else{
        ll <- lm(gene1_exp~gene2_exp)
        file_name <- paste("-",x[1],x[3],sep="_")
        
        png (file = paste(outpath,"/coexpression_",file_name,".png",sep=''),units= "in",width = 5,height = 5,res = 300)
        pp <- plot(gene2_exp,gene1_exp,xlab = x[1],ylab = x[3],
             main = paste("Pearson coefficient of the expression\n R-square: ",round(summary(ll)$r.squared,4),
                          "\n Pearson: ",round(as.numeric(x[6]),4)," p-value: ",signif(as.numeric(x[7]),4),sep=''),pch=15,col="orange")
        lines(gene2_exp,fitted(ll),col="red")
        print(pp)
        dev.off()
      }
    })
    all_result <- cbind(result2[,c(5:7),drop=F],r_result)
    return(list(sig = result2_sig,all = all_result))
  }
}
#######################################
ID_transfer <- function(gene_list,ID_transfer_Path=NULL){
  
  l2 <- load(paste(ID_transfer_Path,"ID_transfer_file.RData",sep='/'))
  ID_transfer_data <- eval(parse(text = l2))
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
########################plot heatmap
plot_heatmap_data <- function(corr_result,choose_cancer = TRUE){
  corr_pairs <- paste(corr_result[,"Gene1Name"],corr_result[,"Gene2Name"],sep='_')
  corr_result2 <- cbind(corr_result,corr_pairs)
  colnames(corr_result2)[8] <- "pairs"

  all_data <- plot_all_heatmap(plot_corr = corr_result2,choose_cancer = choose_cancer)

  if(length(unique(as.matrix(corr_result2[,"Cancer"])))==1){
    corr_result2 <- corr_result2[order(as.numeric(corr_result2[,"Pearson"]),decreasing = T),,drop=F]
    
    if(nrow(corr_result2)>= 25){
      plot_corr <- corr_result2[1:25,]
    }else{
      plot_corr <- corr_result2
    }
  }else{
    pairs_count <- aggregate(corr_result2[,"Cancer"], list(pairs = corr_result2[,"pairs"]), function(x) length(unique(x)))
    pairs_count_order <- pairs_count[order(pairs_count[,2],decreasing = T),,drop=F]
    
    if(nrow(pairs_count_order)>=25){
      plot_pairs <- pairs_count_order[1:25,]
    }else{
      plot_pairs <- pairs_count_order
    }
    plot_corr <- corr_result2[corr_result2[,"pairs"] %in% plot_pairs[,1],]
  }
  plot_data <- plot_all_heatmap(plot_corr = plot_corr,choose_cancer = choose_cancer)
  return(list(all_data=all_data,plot_data=plot_data))
}
###################output all heatmap
plot_all_heatmap <- function(plot_corr,choose_cancer = TRUE){
  ucp <- unique(plot_corr[,"pairs"])
  ucpy <- cbind(ucp,c(0:(length(ucp)-1)))
  
  if(choose_cancer == TRUE){
    corr_cancers <- unique(plot_corr[,"Cancer"])
    ccx <- cbind(corr_cancers,c(0:(length(corr_cancers)-1)))
  }else{
    corr_cancers <- "User cancer"
    ccx <- cbind("User cancer",c(0:(length(corr_cancers)-1)))
  }
  plot_corr2 <- cbind(plot_corr,x=0,y=0)
  for(i in 1:nrow(ucpy)){
    plot_corr2[which(plot_corr2[,"pairs"] %in% ucpy[i,1]),"y"] <- ucpy[i,2]
  }
  
  for(i in 1:nrow(ccx)){
    plot_corr2[which(plot_corr2[,"Cancer"] %in% ccx[i,1]),"x"] <- ccx[i,2]
  }
  ####[y,x,value]
  xyAxis0 <- paste("[",plot_corr2[,"y"],",",plot_corr2[,"x"],",",plot_corr2[,"Pearson"],"]",sep='')
  xyAxis1 <- paste(xyAxis0,collapse =',')
  xyAxis <- paste("[",xyAxis1,"]",sep = "")
  y0<- paste(ucp,collapse="','")
  y <- paste("['",y0,"']",sep='')
  x0<- paste(corr_cancers,collapse="','")
  x <- paste("['",x0,"']",sep='')
  
  return(t(c(x,y,xyAxis)))
}
######################boxplot
plot_box_data <- function(all,choose_cancer = TRUE){
  if(choose_cancer == TRUE){
    
    all_cancer <- unique(all[,1])
    
    all_test <- lapply(all_cancer,function(x){
      data <- all[all[,1] == x,,drop=F]
      
      if(length(data[,2])==1|length(data[,4])==1){
        p <- 1
      }else if(length(unique(data[,2]))==1&length(unique(data[,4]))==1){
        p <- 1
      }else{
        p <- t.test(as.numeric(data[,2]),as.numeric(data[,4]),alternative = "two.sided",paired = T)$p.value
        p <- round(p,4)
      }
      real_box0 <- paste(round(as.numeric(data[,2]),4),collapse=',')
      real_box <- paste("[",real_box0,"]",sep = '')
      
      random_box0 <- paste(round(as.numeric(data[,4]),4),collapse=',')
      random_box <- paste("[",random_box0,"]",sep = '')

      return(list(p_value = p,real_box = real_box,random_box = random_box))
    })

    p_value <- lapply(all_test,function(x) x$p_value)
    p_value <- do.call(rbind,p_value)
   
    xx <- paste(all_cancer,"(p=",p_value,")",sep='')
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
    
    if(length(all[,2])==1|length(all[,4])==1){
      p <- 1
    }else if(length(unique(all[,2]))==1&length(unique(all[,4]))==1){
      p <- 1
    }else{
      p <- t.test(as.numeric(all[,2]),as.numeric(all[,4]),alternative = "two.sided",paired = T)$p.value
    }
      xx <- paste("User cancer(p=",signif(p,4),")",sep='')
      cancer_s <- paste("'",xx,"'",sep='')
      
      real_box_s0 <- paste(round(as.numeric(all[,2]),4),collapse=',')
      real_box_s <- paste("[[",real_box_s0,"]]",sep='')
      
      random_box_s0 <- paste(round(as.numeric(all[,4]),4),collapse=',')
      random_box_s <- paste("[[",random_box_s0,"]]",sep='')
      return(t(c(cancer_s,real_box_s,random_box_s)))
  }
}