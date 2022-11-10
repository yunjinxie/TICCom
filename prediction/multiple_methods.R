multiple_methods <- function(path,method_names,outpath){  
  
  suppressMessages({library(jsonlite)})
  
  stopMessage <- NULL
  
  ID_new <- paste(tail(unlist(strsplit(path,"/")),n=1),collapse = "/")
  ## updata Job_ID.txt, append a new job and state is run
  son_path <- paste(tail(unlist(strsplit(path,"/")),n=2),collapse = "/")
  referencePath <- gsub(son_path,"Rfunction/function_relay",path)
  
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
    method <- data.frame(real_name = c("iTALK(top genes)","iTALK(DEG)","CellTalker","ICELLNET","NicheNet"),
                         result_name = c("iTALK_method_top_genes","iTALK_method_DEG","celltalker_method",
                                         "ICELLNET_method","NicheNet_method"),
                         method_name = c("italk_top","italk_deg","celltalker","icellnet","nichenet"))
    
    method_names <- unlist(strsplit(method_names,','))
    method_real <- method[method[,"method_name"] %in% method_names,,drop=F]
    file <- paste(as.matrix(method_real[,"result_name"]),"_r_re.txt",sep='')
    
    pairs_data_all <- c()
    for(i in 1:nrow(method_real)){
      file_path <- file.path(path,file[i])
      if(file.exists(file_path)){
        if(file.info(file_path)$size != 0){
          json <- do.call(rbind, 
                          lapply(paste(readLines(file_path, warn=FALSE),
                                       collapse=""), 
                                 jsonlite::fromJSON))
          pairs_data <- json[[1]]
          pairs_data <- as.matrix(pairs_data)
          
          pairs_data <- pairs_data[,c("ligand","cell_from","receptor","cell_to","Function","curated")]
          pairs_data <- cbind(pairs_data,method = as.character(method_real[i,"real_name"]))
          pairs_data_all <- rbind(pairs_data_all,pairs_data)
        }
      }else{
        pairs_data<- NULL
        pairs_data_all <- rbind(pairs_data_all,pairs_data)
      }
    }
    
    if(is.null(pairs_data_all)){
      stopMessage <- "no results."
    }else{
      ####==============================================output table
      union_pair <- unique(pairs_data_all[,c("ligand","cell_from","receptor","cell_to","Function","curated"),drop=F])
      method_r <- unique(pairs_data_all[,"method"])
      union_pair_mat <- matrix("-",nrow=nrow(union_pair),ncol=nrow(method))
      colnames(union_pair_mat) <- as.matrix(method[,1])
      union_pair_mat2 <- cbind(union_pair,union_pair_mat)
      
      pairs <- apply(union_pair,1,function(x) paste(x[1:5],collapse='_'))
      rownames(union_pair_mat2) <- pairs
      
      for(j in 1:length(method_r)){
        pairs_data_mm <- pairs_data_all[pairs_data_all[,"method"] %in% method_r[j],,drop=F]
        if(nrow(pairs_data_mm)!=0){
          pairs1 <- apply(pairs_data_mm,1,function(x) paste(x[1:5],collapse='_'))
          union_pair_mat2[pairs1,method_r[j]] <- "Yes"
        }
      }
      rownames(union_pair_mat2) <- NULL
      plot_data_df <- as.data.frame(union_pair_mat2)
      output <- list(data = plot_data_df)
      jsoncars <- toJSON(output, pretty=TRUE)
      cat(jsoncars, file = paste(outpath,'multiple_methods_gene_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
      ####========================================Venn plot combn(n,m-1)
      rs1 <- c()
      for(k in 1:length(method_r)){
        group <- combn(length(method_r),k)
        for(n in 1:ncol(group)){
          index <- group[,n]
          col_method <- method_r[index]
          
          union_pair_mat3 <- union_pair_mat2[,col_method,drop=F]
          rs <- apply(union_pair_mat3,1,function(x) length(which(x=="Yes")))
          len <- length(which(rs==length(col_method)))
          
          yes_name1 <- paste(col_method,collapse="','")
          sets <- paste("sets: ['",yes_name1,"']",sep='')
          value <- paste("value: ",len,sep = "")
          rs0 <- paste("{",sets,",",value,"}",sep='')
          rs1 <- c(rs1,rs0)
        }
      } 
      rs2_s <- paste(rs1,collapse = ",")
      write.table(rs2_s,paste(outpath,"multiple_methods_gene_re_venn.txt",sep='/'),sep='\t',quote=F,row.names=F,col.names=F)
      ####========================================riverplot
      count  <-  apply(union_pair_mat2,1,function(x) {length(which(x=="Yes"))})
      c_max <- max(count)
      plot_data <- union_pair_mat2[count==c_max,,drop=F]
      
      if(nrow(plot_data) >= 25){
        plot_data <- plot_data[1:25,]
      }else{
        plot_data <- plot_data
      }
      
      plot_data[,"cell_from"] <- paste(plot_data[,"cell_from"],"(from)",sep='')
      plot_data[,"cell_to"] <- paste(plot_data[,"cell_to"],"(to)",sep='')
      plot_data[,"ligand"] <- paste(plot_data[,"ligand"],"(ligand)",sep='')
      plot_data[,"receptor"] <- paste(plot_data[,"receptor"],"(receptor)",sep='')

      node <- matrix(unique(c(plot_data[,"cell_from"],plot_data[,"ligand"],plot_data[,"Function"],plot_data[,"receptor"],plot_data[,"cell_to"])),ncol=1)
      colnames(node) <- "name"
      node <- as.data.frame(node)
      
      index <- rep(1,times= nrow(plot_data))
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
      cat(jsoncars, file = paste(outpath,'multiple_methods_riverplot.json',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
      ####=======================================circle
      plot_circle <- union_pair_mat2[,c(1,2,3,4,5),drop=F]
      gene_pair <- paste(plot_circle[,"ligand"],plot_circle[,"receptor"],sep='_')
      plot_circle2 <- cbind(plot_circle,gene_pair)
      
      cell_count <- aggregate(plot_circle2[,"gene_pair"],list(plot_circle2[,"cell_from"],plot_circle2[,"cell_to"]),length)
      cell_count_sort <- cell_count[order(cell_count[,3],decreasing = T),,drop=F]
      
      plot_circle_final <- plot_circle2[plot_circle2[,"cell_from"] == cell_count_sort[1,1]&plot_circle2[,"cell_to"] == cell_count_sort[1,2],,drop=F]
      if(nrow(plot_circle_final)>=40){
        plot_circle_final2 <-plot_circle_final[1:40,]
      }else{
        plot_circle_final2 <- plot_circle_final
      }
      
      node <- union(plot_circle_final2[,"ligand"],plot_circle_final2[,"receptor"])
      degree <- table(c(plot_circle_final2[,"ligand"],plot_circle_final2[,"receptor"]))
      degree2 <- as.matrix(degree)
      degree3 <- degree2[node,]
      
      ID <- 0:(length(node)-1)
      
      symbolSize <- rep(30,times = length(node))
      
      ligand <- unique(plot_circle_final2[,"ligand"])
      receptor <- unique(plot_circle_final2[,"receptor"])
      cell_from <- unname(plot_circle_final2[1,"cell_from"])
      cell_to <- unname(plot_circle_final2[1,"cell_to"])
      
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
      
      links_final <- cbind(node2[match(plot_circle_final2[,"ligand"],node2[,"name"]),"id"],
                           node2[match(plot_circle_final2[,"receptor"],node2[,"name"]),"id"])
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
      cat(jsoncars, file = paste(outpath,"multiple_methods_circle.json",sep='/'), fill = FALSE, labels = NULL, append = FALSE)
    }
    if(is.null(stopMessage)){
      download_table="multiple_methods_gene_re.txt"
      Res_venn_file="multiple_methods_gene_re_venn.txt"
      Res_riverplot_file="multiple_methods_riverplot.json"
      Res_circle_file="multiple_methods_circle.json"
      resut_merge=paste0("{",
                         '"download_table" :','"',download_table,'",',
                         '"Res_venn_file" :','"',Res_venn_file,'",',
                         '"Res_riverplot_file" :','"',Res_riverplot_file,'",',
                         '"Res_circle_file" :','"',Res_circle_file,'",',
                         '"error_attention" :','"no',
                         '"}')
      write.table(resut_merge,paste(outpath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
      
      Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
      Task_ID_new <- Task_ID
      Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "multiple_methods+success"
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
  