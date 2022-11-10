TCGA_survival <- function(IDPath=NULL,expPath=NULL,clinicPath=NULL,gene_listPath=NULL,referencePath=NULL,
                          cancer=NULL,pvalue=0.05,outpath){
  suppressMessages({
    library(survival)
    library(jsonlite)
  })

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
        exp_data <- lapply(TCGAFile,function(x){
          el <- load(x)
          e_data <- eval(parse(text = el))
          return(e_data)
        })
        names(exp_data) <- cancer
        
        TCGAClinicFile <- paste(referencePath,"/",cancer,"_clinic.RData",sep='')
        clinic_data <- lapply(TCGAClinicFile,function(x){
          clil <- load(x)
          cli_data <- eval(parse(text = clil))
          return(cli_data)
        })
        names(clinic_data) <- cancer
        sur_res <- lapply(1:length(cancer),function(x) Survival_analysis(exp = exp_data[[x]],clinic = clinic_data[[x]],
                                                                         cancer = cancer[x],pairs = pairs2,
                                                                         pvalue = pvalue,outpath = outpath))
        pairs_res <- lapply(sur_res,function(x) x[[1]])
        pairs_res <- do.call(rbind,pairs_res)
        cox_res <- lapply(sur_res,function(x) x[[2]])
        cox_res <- do.call(rbind,cox_res)
        if(is.null(pairs_res)){
          stopMessage <- "no results."
        }else{
          ####====================================================barplot
          bar_res <- plot_bar_data(data = cox_res,choose_cancer =TRUE)
          write.table(bar_res,paste(outpath,"TCGA_survival_barplot_data.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
          write.table(cox_res,paste(outpath,"TCGA_survival_barplot_download.txt",sep='/'),sep='\t',quote = F,row.names = F)
          ####==================================================output table
          pairs_res[,"RScore(median)"] <- round(as.numeric(pairs_res[,"RScore(median)"]),4)
          pairs_res[,"k-m pvalue"] <- signif(as.numeric(pairs_res[,"k-m pvalue"]),4)
          
          pairs_res_df <- as.data.frame(pairs_res)
          rownames(pairs_res_df) <- NULL
          output <- list(data = pairs_res_df)
          jsoncars <- toJSON(output, pretty=TRUE)
          cat(jsoncars, file = paste(outpath,'TCGA_survival_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
          ####====================================================plot heatmap
          heatmap_res<- plot_heatmap_data(sur_result = pairs_res,choose_cancer = TRUE)
          all_data <- heatmap_res[[1]]
          plot_data <- heatmap_res[[2]]
          write.table(all_data,paste(outpath,"TCGA_survival_heatmap_all.txt",sep='/'),sep='\t',row.names = F,col.names = F,quote = F)
          write.table(plot_data,paste(outpath,"TCGA_survival_heatmap.txt",sep='/'),sep='\t',row.names = F,col.names = F,quote = F)
          write.table(pairs_res,paste(outpath,"TCGA_survival_heatmap_download.txt",sep='/'),sep='\t',quote = F,row.names = F)
        }
      }else if(expPath!="empty"&clinicPath!="empty"){
        exp0 <- read.table(expPath,sep='\t',header=T,as.is = T,fill = T,strip.white = T,check.names=F)
        if((class(exp0[,1]) == "character") | (class(exp0[,1]) == "factor")){
          rownames(exp0) <- as.matrix(exp0[,1])
          exp0 <- exp0[,-1]
        }
        clinic <- read.table(clinicPath,sep='\t',header=T,as.is = T,fill = T,strip.white = T,check.names=F)
        if(sum("submitter_id" %in% colnames(clinic))==0|sum("vital_status" %in% colnames(clinic))==0|sum("days" %in% colnames(clinic))==0){
          stopMessage <- "error setting colnames in clinic file."
        }else{
          rownames(clinic) <- as.matrix(clinic[,"submitter_id"])
          sam_inter <- intersect(colnames(exp0),rownames(clinic))
          if(length(sam_inter)==0){
            stopMessage <- "no matched samples."
          }else{
            exp0 <- exp0[,sam_inter]
            clinic <- clinic[sam_inter,]
            #change to ensembl ID
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
                sur_res <- Survival_analysis(exp = exp,clinic = clinic,pairs = pairs2,pvalue = pvalue,outpath = outpath)
                pairs_res <- sur_res[[1]]
                cox_res <- sur_res[[2]]
                if(is.null(pairs_res)){
                  stopMessage <- "no results."
                }else{
                  ####========================================plot bar
                  bar_res <- plot_bar_data(data = cox_res,choose_cancer =FALSE)
                  write.table(bar_res,paste(outpath,"TCGA_survival_barplot_data.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
                  write.table(cox_res,paste(outpath,"TCGA_survival_barplot_download.txt",sep='/'),sep='\t',quote = F,row.names = F)
                  ####========================================output table
                  pairs_res[,"RScore(median)"] <- round(as.numeric(pairs_res[,"RScore(median)"]),4)
                  pairs_res[,"k-m pvalue"] <- signif(as.numeric(pairs_res[,"k-m pvalue"]),4)
                  pairs_res_df <- as.data.frame(pairs_res)
                  rownames(pairs_res_df) <- NULL
                  output <- list(data = pairs_res_df)
                  jsoncars <- toJSON(output, pretty=TRUE)
                  cat(jsoncars, file = paste(outpath,'TCGA_survival_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
                  ####=======================================plot heatmap
                  heatmap_res <- plot_heatmap_data(sur_result = pairs_res,choose_cancer = FALSE)
                  all_data <- heatmap_res[[1]]
                  plot_data <- heatmap_res[[2]]
                  write.table(all_data,paste(outpath,"TCGA_survival_heatmap_all.txt",sep='/'),sep='\t',row.names = F,col.names = F,quote = F)
                  write.table(plot_data,paste(outpath,"TCGA_survival_heatmap.txt",sep='/'),sep='\t',row.names = F,col.names = F,quote = F)
                  write.table(pairs_res,paste(outpath,"TCGA_survival_heatmap_download.txt",sep='/'),sep='\t',quote = F,row.names = F)
                }
              }
            }
          }
        }
      }else{
        stopMessage <- "wrong choices."
      }
    }
    if(is.null(stopMessage)){
      download_table="TCGA_survival_re.txt"
      heat_download_table="TCGA_survival_heatmap_download.txt"
      Res_heat_all_file="TCGA_survival_heatmap_all.txt"
      Res_heat_file="TCGA_survival_heatmap.txt"
      Res_barplot_file="TCGA_survival_barplot_data.txt"
      barplot_download_table="TCGA_survival_barplot_download.txt"
      resut_merge=paste0("{",
                         '"download_table" :','"',download_table,'",',
                         '"heat_download_table" :','"',heat_download_table,'",',
                         '"Res_heat_all_file" :','"',Res_heat_all_file,'",',
                         '"Res_heat_file" :','"',Res_heat_file,'",',
                         '"Res_barplot_file" :','"',Res_barplot_file,'",',
                         '"barplot_download_table" :','"',barplot_download_table,'",',
                         '"error_attention" :','"no',
                         '"}')
      write.table(resut_merge,paste(outpath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
      
      Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
      Task_ID_new <- Task_ID
      Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "TCGA_survival+success"
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
  },error =function(e){
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
################
Survival_analysis <- function(exp,clinic,cancer=NULL,pairs,pvalue = pvalue,outpath = outpath){
  
  pairs <- as.matrix(pairs)
  gene1 <- rownames(exp)[rownames(exp) %in% pairs[,"Gene1Id"]]
  gene2 <- rownames(exp)[rownames(exp) %in% pairs[,"Gene2Id"]]
  index1 <- which(pairs[,"Gene1Id"] %in% gene1)
  index2 <- which(pairs[,"Gene2Id"] %in% gene2)
  index <- intersect(index1,index2)
  pairs <- unique(pairs[index,c(1:4),drop=F])
  nlen <- nrow(pairs)
  if(nlen==0){
    return(list(pair = NULL,single = NULL))
  }
  
  result <- lapply(1:nlen,function(x){
    
    gene1 <- pairs[x,"Gene1Name"]
    gene2 <- pairs[x,"Gene2Name"]
    gene1ID <- pairs[x,"Gene1Id"]
    gene2ID <- pairs[x,"Gene2Id"]
    gene1_exp <- as.numeric(exp[gene1ID,])
    gene2_exp <- as.numeric(exp[gene2ID,])
    #==================survival a pair
    a <- coxph(Surv(days,vital_status)~as.numeric(gene1_exp)+as.numeric(gene2_exp),clinic)
    aa<-summary(a)$coef
    sam <- colnames(exp)
    pairs_risk_score <- as.numeric(aa[1,1])*gene1_exp+as.numeric(aa[2,1])*gene2_exp
    cutoff <- median(pairs_risk_score)
    if(is.na(cutoff)){
      return(list(pair = NULL,single = NULL))
    }
    if(cutoff == min(pairs_risk_score)){
      group_high <- sam[pairs_risk_score > cutoff]
      group_low <- sam[pairs_risk_score <= cutoff]
    }else{
      group_high <- sam[pairs_risk_score >= cutoff]
      group_low <- sam[pairs_risk_score < cutoff]
    }
    if(length(group_high)==0|length(group_low)==0){
      return(list(pair = NULL,single = NULL))
    }
    group <- rbind(cbind(group_high,"high"),cbind(group_low,"low"))
    colnames(group) <- c("submitter_id","group")
    clinic_group <- merge(clinic,group,by="submitter_id")
    #===================================================================================
    kmfit<-survfit(Surv(days,vital_status)~group,clinic_group)
    km_p<-1-pchisq(survdiff(Surv(days,vital_status)~group,clinic_group)$chisq,1)
    
    if(km_p < as.numeric(pvalue)){
      #==================survival the first gene
      data1 <- list(time=clinic$days,status=clinic$vital_status,exp=gene1_exp)
      cox_result1 <- summary(coxph(Surv(time,status)~exp,data1))
      cox_p1 <- cox_result1$coefficients[,5]
      cox_HR1 <- cox_result1$coefficients[2]
      cutoff1 <- median(gene1_exp)
      if(is.na(cutoff1)){
        return(list(pair = NULL,single = NULL))
      }
      if(cutoff1 == min(gene1_exp)){
        gene1_exp_high <- sam[gene1_exp > cutoff1]
        gene1_exp_low <- sam[gene1_exp <= cutoff1]
      }else{
        gene1_exp_high <- sam[gene1_exp >= cutoff1]
        gene1_exp_low <- sam[gene1_exp < cutoff1]
      }
      
      if(length(gene1_exp_high)==0|length(gene1_exp_low)==0){
        return(list(pair = NULL,single = NULL))
      }
      gene1_exp_group <- rbind(cbind(gene1_exp_high,"high"),cbind(gene1_exp_low,"low"))
      colnames(gene1_exp_group) <- c("submitter_id","group")
      gene1_exp_group2 <- merge(clinic,gene1_exp_group,by="submitter_id")
      
      kmfit1<-survfit(Surv(days,vital_status)~group,gene1_exp_group2)
      km_p1<-1-pchisq(survdiff(Surv(days,vital_status)~group,gene1_exp_group2)$chisq,1)
      #====================survival the second gene
      data2 <- list(time=clinic$days,status=clinic$vital_status,exp=gene2_exp)
      cox_result2 <- summary(coxph(Surv(time,status)~exp,data2))
      cox_p2 <- cox_result2$coefficients[,5]
      cox_HR2 <- cox_result2$coefficients[2]
      cutoff2 <- median(gene2_exp)
      if(is.na(cutoff2)){
        return(list(pair = NULL,single = NULL))
      }
      if(cutoff2 == min(gene2_exp)){
        gene2_exp_high <- sam[gene2_exp > cutoff2]
        gene2_exp_low <- sam[gene2_exp <= cutoff2]
      }else{
        gene2_exp_high <- sam[gene2_exp >= cutoff2]
        gene2_exp_low <- sam[gene2_exp < cutoff2]
      }
      
      if(length(gene2_exp_high)==0|length(gene2_exp_low)==0){
        return(list(pair = NULL,single = NULL))
      }
      gene2_exp_group <- rbind(cbind(gene2_exp_high,"high"),cbind(gene2_exp_low,"low"))
      colnames(gene2_exp_group) <- c("submitter_id","group")
      gene2_exp_group2 <- merge(clinic,gene2_exp_group,by="submitter_id")
      
      kmfit2<-survfit(Surv(days,vital_status)~group,gene2_exp_group2)
      km_p2<-1-pchisq(survdiff(Surv(days,vital_status)~group,gene2_exp_group2)$chisq,1)
      #=============================boxplot data
      pair_exp <- rbind(gene1_exp,gene2_exp)
      
      colnames(pair_exp) <- sam
      rownames(pair_exp) <- c(gene1,gene2)
      
      pair_group <- clinic_group[match(sam,clinic_group[,"submitter_id"]),"group",]
      box_data <- data.frame(exp = as.numeric(pair_exp),
                             sam = rep(sam,each = 2),
                             group = rep(pair_group,each = 2),
                             gene = rep(c(gene1,gene2),times = ncol(pair_exp)))
      
      gene1_high_risk <- as.numeric(box_data[box_data[,3]=="high"&box_data[,4]==gene1,1])
      gene1_low_risk <- as.numeric(box_data[box_data[,3]=="low"&box_data[,4]==gene1,1])
      if(length(gene1_high_risk)==1|length(gene1_low_risk)==1){
        p01 <- 1
      }else if(length(unique(gene1_high_risk))==1&length(unique(gene1_low_risk))==1){
        p01 <- 1
      }else{
        p01 <- t.test(gene1_high_risk,gene1_low_risk)$p.value
      }
      
      gene2_high_risk <- as.numeric(box_data[box_data[,3]=="high"&box_data[,4]==gene2,1])
      gene2_low_risk <- as.numeric(box_data[box_data[,3]=="low"&box_data[,4]==gene2,1])
      if(length(gene2_high_risk)==1|length(gene2_low_risk)==1){
        p02 <- 1
      }else if(length(unique(gene2_high_risk))==1&length(unique(gene2_low_risk))==1){
        p02 <- 1
      }else{
        p02 <- t.test(gene2_high_risk,gene2_low_risk)$p.value
      }

      gene1_group <- gene1_exp_group2[match(sam,gene1_exp_group2[,"submitter_id"]),"group",]
      box_data1 <- data.frame(exp = gene1_exp,
                              sam = sam,
                              group = gene1_group,
                              gene = rep(gene1,times = length(gene1_exp)))
      
      gene1_high <- as.numeric(box_data1[box_data1[,3]=="high"&box_data1[,4]==gene1,1])
      gene1_low <- as.numeric(box_data1[box_data1[,3]=="low"&box_data1[,4]==gene1,1])
      if(length(gene1_high)==1|length(gene1_low)==1){
        p1 <- 1
      }else if(length(unique(gene1_high))==1&length(unique(gene1_low))==1){
        p1 <- 1
      }else{
        p1 <- t.test(gene1_high,gene1_low)$p.value
      }
      
      gene2_group <- gene2_exp_group2[match(sam,gene2_exp_group2[,"submitter_id"]),"group",]
      box_data2 <- data.frame(exp = gene2_exp,
                              sam = sam,
                              group = gene2_group,
                              gene = rep(gene2,times = length(gene2_exp)))
      
      gene2_high <- as.numeric(box_data2[box_data2[,3]=="high"&box_data2[,4]==gene2,1])
      gene2_low <- as.numeric(box_data2[box_data2[,3]=="low"&box_data2[,4]==gene2,1])
      if(length(gene2_high)==1|length(gene2_low)==1){
        p2 <- 1
      }else if(length(unique(gene2_high))==1&length(unique(gene2_low))==1){
        p2 <- 1
      }else{
        p2 <- t.test(gene2_high,gene2_low)$p.value
      }
      #===============================plot=============================================
      if(!is.null(cancer)==T){
        ##用survminer画
        file_name <- paste(cancer,gene1,gene2,sep='_')
        png(paste(outpath,"/survival_",file_name,".png",sep=''),units= "in",width = 12,height = 7,res = 300)
        layout(matrix(1:6,byrow=T,nrow=2))
       #======================plot a pair 
        plot(kmfit, conf.int=F, mark.time=TRUE,
             pch=3, col=c("#b00011","#00b033"), lty=1, lwd=2, cex.main=1.5,cex.lab=1.2,
             xlab="Time (days)", ylab="Survival probability",
             conf.type = "log",
             main=paste(gene1," and ",gene2," in ",cancer,"\np = ",signif(km_p,4),sep=''))
        legend("bottomleft",legend = c("high","low"),ncol=1,col=c("#b00011","#00b033"),lwd=2,xpd=T,cex=1.2)
        #legend("center",legend = paste0("p = ",signif(km_p,4)),cex=1.2)
        #===============plot the first gene
        plot(kmfit1, conf.int=F, mark.time=TRUE,
             pch=3, col=c("#b00011","#00b033"), lty=1, lwd=2, cex.main=1.5,cex.lab=1.2,
             xlab="Time (days)", ylab="Survival probability",
             conf.type = "log",
             main=paste("K-M plot of ",gene1," in ",cancer,"\np = ",signif(km_p1,4),sep=''))
        legend("bottomleft",legend = c("high","low"),ncol=1,col=c("#b00011","#00b033"),lwd=2,xpd=T,cex=1.2)
        #===============plot the second gene
        plot(kmfit2, conf.int=F, mark.time=TRUE,
             pch=3, col=c("#b00011","#00b033"), lty=1, lwd=2, cex.main=1.5,cex.lab=1.2,
             xlab="Time (days)", ylab="Survival probability",
             conf.type = "log",
             main=paste("K-M plot of ",gene2," in ",cancer,"\np = ",signif(km_p2,4),sep=''))
        legend("bottomleft",legend = c("high","low"),ncol=1,col=c("#b00011","#00b033"),lwd=2,xpd=T,cex=1.2)
        ########boxplot
        boxplot(exp~(group+gene),box_data,col = c("#b00011","#00b033"),xlab = "",ylab = "Expression",las=2,xaxt = 'n',cex.lab=1.2,main= paste("(p1=",signif(p01,4),", p2=",signif(p02,4),")",sep = ""))
        axis(1,at = seq(0.5,by=1,length.out = 4),labels = c("",gene1,"",gene2),tick = T,lwd.ticks = 0.1,cex.axis = 1.2)
         ##################gene1
        boxplot(exp~(group+gene),box_data1,col = c("#b00011","#00b033"),xlab = "",ylab = "Expression",las=2,xaxt = 'n',cex.lab=1.2,main = paste("(p=",signif(p1,4),")",sep = ""))
        axis(1,at = seq(0.5,by=1,length.out = 2),labels = c("",gene1),tick = T,lwd.ticks = 0.1,cex.axis = 1.2)
       ###################gene2
        boxplot(exp~(group+gene),box_data2,col = c("#b00011","#00b033"),xlab = "",ylab = "Expression",las=2,xaxt = 'n',cex.lab=1.2,main = paste("(p=",signif(p2,4),")",sep = ""))
        axis(1,at = seq(0.5,by=1,length.out = 2),labels = c("",gene2),tick = T,lwd.ticks = 0.1,cex.axis = 1.2)
   
        dev.off()
        
        pair_res <- cbind(pairs[x,c(1:4),drop=F],cancer,round(cutoff,4),signif(km_p,4))
        colnames(pair_res) <- c("Gene1Name","Gene1Id","Gene2Name","Gene2Id","Cancer","RScore(median)","k-m pvalue")
        
        single_res <- rbind(cbind(gene1,gene1ID,cancer,signif(cox_p1,4),round(cox_HR1,4)),
                            cbind(gene2,gene2ID,cancer,signif(cox_p2,4),round(cox_HR2,4)))
        colnames(single_res) <- c("GeneName","GeneID","Cancer","Cox pvalue","HR")
        return(list(pair = pair_res,single = single_res))
      }else{
        file_name <- paste("-",gene1,gene2,sep='_')
        png(paste(outpath,"/survival_",file_name,".png",sep=''),units= "in",width = 12,height = 7,res = 300)
        layout(matrix(1:6,byrow=T,nrow=2))
        #======================plot a pair 
        plot(kmfit, conf.int=F, mark.time=TRUE,
             pch=3, col=c("#b00011","#00b033"), lty=1, lwd=2, cex.main=1.5,cex.lab=1.2,
             xlab="Time (days)", ylab="Survival probability",
             conf.type = "log",
             main=paste("K-M plot of ",gene1," and ",gene2,"\n","p = ",signif(km_p,4),sep=''))
        legend("bottomleft",legend = c("high","low"),ncol=1,col=c("#b00011","#00b033"),lwd=2,xpd=T,cex=1.2)
        #===============plot the first gene
        plot(kmfit1, conf.int=F, mark.time=TRUE,
             pch=3, col=c("#b00011","#00b033"), lty=1, lwd=2, cex.main=1.5,cex.lab=1.2,
             xlab="Time (days)", ylab="Survival probability",
             conf.type = "log",
             main=paste("K-M plot of ",gene1,"\np = ",signif(km_p1,4),sep=''))
        legend("bottomleft",legend = c("high","low"),ncol=1,col=c("#b00011","#00b033"),lwd=2,xpd=T,cex=1.2)
        #===============plot the second gene
        plot(kmfit2, conf.int=F, mark.time=TRUE,
             pch=3, col=c("#b00011","#00b033"), lty=1, lwd=2, cex.main=1.5,cex.lab=1.2,
             xlab="Time (days)", ylab="Survival probability",
             conf.type = "log",
             main=paste("K-M plot of ",gene2,"\np = ",signif(km_p2,4),sep=''))
        legend("bottomleft",legend = c("high","low"),ncol=1,col=c("#b00011","#00b033"),lwd=2,xpd=T,cex=1.2)
        ########boxplot
        boxplot(exp~(group+gene),box_data,col = c("#b00011","#00b033"),xlab = "",ylab = "Expression",las=2,xaxt = 'n',cex.lab=1.2,main= paste("(p1=",signif(p01,4),", p2=",signif(p02,4),")",sep = ""))
        axis(1,at = seq(0.5,by=1,length.out = 4),labels = c("",gene1,"",gene2),tick = T,lwd.ticks = 0.1,cex.axis = 1.2)
        ##################gene1
        boxplot(exp~(group+gene),box_data1,col = c("#b00011","#00b033"),xlab = "",ylab = "Expression",las=2,xaxt = 'n',cex.lab=1.2,main= paste("(p=",signif(p1,4),")",sep = ""))
        axis(1,at = seq(0.5,by=1,length.out = 2),labels = c("",gene1),tick = T,lwd.ticks = 0.1,cex.axis = 1.2)
        ###################gene2
        boxplot(exp~(group+gene),box_data2,col = c("#b00011","#00b033"),xlab = "",ylab = "Expression",las=2,xaxt = 'n',cex.lab=1.2,main= paste("(p=",signif(p2,4),")",sep = ""))
        axis(1,at = seq(0.5,by=1,length.out = 2),labels = c("",gene2),tick = T,lwd.ticks = 0.1,cex.axis = 1.2)
        
        dev.off()
        
        pair_res <- cbind(pairs[x,c(1:4),drop=F],"-",round(cutoff,4),signif(km_p,4))
        colnames(pair_res) <- c("Gene1Name","Gene1Id","Gene2Name","Gene2Id","Cancer","RScore(median)","k-m pvalue")
        
        single_res <- rbind(cbind(gene1,gene1ID,"-",signif(cox_p1,4),round(cox_HR1,4)),
                            cbind(gene2,gene2ID,"-",signif(cox_p2,4),round(cox_HR2,4)))
        colnames(single_res) <- c("GeneName","GeneID","Cancer","Cox pvalue","HR")
        return(list(pair = pair_res,single = single_res))
      }
    }else{
      return(list(pair = NULL,single = NULL))
    }
  })
  
  pairs_result <- lapply(result,function(x) x[[1]])
  pairs_result2 <- do.call(rbind,pairs_result)
  
  single_result <- lapply(result,function(x) x[[2]])
  single_result2 <- unique(do.call(rbind,single_result))
  
  return(list(pairs_result = pairs_result2,single_result = single_result2))
}
################ID transfer
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
plot_heatmap_data <- function(sur_result,choose_cancer = TRUE){
  sur_pairs <- paste(sur_result[,"Gene1Name"],sur_result[,"Gene2Name"],sep='_')
  sur_pairs2 <- cbind(sur_result,sur_pairs)
  colnames(sur_pairs2)[8] <- "pairs"
  
  all_data <- plot_all_heatmap(plot_sur = sur_pairs2,choose_cancer = choose_cancer)

  if(length(unique(as.matrix(sur_pairs2[,"Cancer"])))==1){
    plot_sur <- sur_pairs2[order(as.numeric(sur_pairs2[,"RScore(median)"]),decreasing = T),,drop=F]
    
    if(nrow(plot_sur)>=25){
      plot_sur <- plot_sur[1:25,]
    }else{
      plot_sur <- plot_sur
    }
  }else{
    pairs_count <- aggregate(sur_pairs2[,"Cancer"], list(pairs = sur_pairs2[,"pairs"]), function(x) length(unique(x)))
    pairs_count_order <- pairs_count[order(pairs_count[,2],decreasing = T),,drop=F]
    
    if(nrow(pairs_count_order)>=25){
      plot_pairs <- pairs_count_order[1:25,]
    }else{
      plot_pairs <- pairs_count_order
    }
    plot_sur <- sur_pairs2[sur_pairs2[,"pairs"] %in% plot_pairs[,1],]
  }
  plot_data <- plot_all_heatmap(plot_sur = plot_sur,choose_cancer = choose_cancer)
  return(list(all_data=all_data,plot_data=plot_data))
}
###################output all heatmap
plot_all_heatmap <- function(plot_sur,choose_cancer = TRUE){
  usp <- unique(plot_sur[,"pairs"])
  uspy <- cbind(usp,c(0:(length(usp)-1)))
  
  if(choose_cancer == TRUE){
    sur_cancers <- unique(plot_sur[,"Cancer"])
    scx <- cbind(sur_cancers,c(0:(length(sur_cancers)-1)))
  }else{
    sur_cancers <- "User cancer"
    scx <- cbind("User cancer",c(0:(length(sur_cancers)-1)))
  }
  
  plot_sur2 <- cbind(plot_sur,x=0,y=0)
  ####
  for(i in 1:nrow(uspy)){
    plot_sur2[which(plot_sur2[,"pairs"] %in% uspy[i,1]),"y"] <- uspy[i,2]
  }
  
  for(i in 1:nrow(scx)){
    plot_sur2[which(plot_sur2[,"Cancer"] %in% scx[i,1]),"x"] <- scx[i,2]
  }
  ####[y,x,value]
  xyAxis0 <- paste("[",plot_sur2[,"y"],",",plot_sur2[,"x"],",",plot_sur2[,"RScore(median)"],"]",sep='')
  xyAxis1 <- paste(xyAxis0,collapse =',')
  xyAxis <- paste("[",xyAxis1,"]",sep = "")
  y0<- paste(usp,collapse="','")
  y <- paste("['",y0,"']",sep='')
  x0<- paste(sur_cancers,collapse="','")
  x <- paste("['",x0,"']",sep='')
  
  sur_min <- min(as.numeric(plot_sur2[,"RScore(median)"]))
  sur_max <- max(as.numeric(plot_sur2[,"RScore(median)"]))
  
  return(t(c(x,y,xyAxis,sur_min,sur_max)))
}
#######################plot barplot data
plot_bar_data <- function(data,choose_cancer = TRUE){
  
  data <- as.matrix(data)
  data <- data[as.numeric(data[,"HR"]) < 1|as.numeric(data[,"HR"]) > 1,]
  HR_type <- rep("HR > 1",length.out = nrow(data))
  data2 <- cbind(data,HR_type)
  data2[as.numeric(data[,"HR"]) < 1,"HR_type"] <- "HR < 1"
  count <- aggregate(data2[,"GeneName"],list(data2[,"Cancer"],data2[,"HR_type"]),length)
  colnames(count) <- c("Cancer","HR_type","number")
  
  if(choose_cancer == TRUE){
    
    cancer <- unique(as.matrix(count$Cancer))
    cancer2 <- paste(cancer,collapse = "','")
    cancer3 <- paste("['",cancer2,"']",sep='')
    
    HR1 <- count[count$HR_type=="HR < 1",]
    rownames(HR1) <- as.matrix(HR1$Cancer)
    HR1_num <- as.numeric(HR1[cancer,"number"])
    HR1_num[is.na(HR1_num)] <- 0
    
    HR1_min <- -max(HR1_num)
    
    HR1_num <- paste("-",HR1_num,sep='')
    HR1_num2 <- paste(HR1_num,collapse=',')
    HR1_num3 <- paste("[",HR1_num2,"]",sep = '')
    
    HR2 <- count[count$HR_type=="HR > 1",]
    rownames(HR2) <- as.matrix(HR2$Cancer)
    HR2_num <- as.numeric(HR2[cancer,"number"])
    HR2_num[is.na(HR2_num)] <- 0
    
    HR2_max <- max(HR2_num) 
    
    HR2_num2 <- paste(HR2_num,collapse=',')
    HR2_num3 <- paste("[",HR2_num2,"]",sep = '')
    
  }else{
    cancer <- "User cancer"
    cancer2 <- paste(cancer,collapse = "','")
    cancer3 <- paste("['",cancer2,"']",sep='')
    
    HR1 <- count[count$HR_type=="HR < 1",,drop=F]
    
    HR1_min <- -max(as.numeric(as.matrix(HR1[,"number"])))
    
    HR1_num <- paste("-",HR1[,"number"],sep='')
    HR1_num2 <- paste(HR1_num,collapse=',')
    HR1_num3 <- paste("[",HR1_num2,"]",sep = '')
    
    HR2 <- count[count$HR_type=="HR > 1",,drop=F]
    
    HR2_max <- max(as.numeric(as.matrix(HR2[,"number"])))
    
    HR2_num <- HR2[,"number"]
    HR2_num2 <- paste(HR2_num,collapse=',')
    HR2_num3 <- paste("[",HR2_num2,"]",sep = '')
  }
  return(t(c(cancer3,HR1_num3,HR2_num3,HR1_min,HR2_max)))
}
