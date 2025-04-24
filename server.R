suppressPackageStartupMessages({ ##loading all necessary libraries 
library(readxl, warn=FALSE)
library(openxlsx, warn=FALSE)
library(writexl, warn=FALSE)
library(rlist, warn=FALSE)
library(seqinr, warn=FALSE)
library(gridExtra, warn=FALSE)
library(grid, warn=FALSE)
library(devtools,warn=FALSE)
library(shiny,warn=FALSE)
library(shinydashboard,warn=FALSE)
library(shinyjs,warn=FALSE)
library(NGLVieweR,warn=FALSE)
library(shinyWidgets,warn=FALSE)
library(shinydisconnect,warn=FALSE)
library(shinyBS,warn=FALSE)
library(hdxstats,warn=FALSE)
library(RColorBrewer,warn=FALSE)
library(pheatmap,warn=FALSE)
library(scales,warn=FALSE)
library(viridis,warn=FALSE)
library(patchwork,warn=FALSE)
library(Biostrings,warn=FALSE)
library(tidyverse, warn=FALSE)
library(qpdf, warn=FALSE)
library(ggpubr, warn=FALSE)
library(plyr, warn=FALSE)
})

server <- function(input, output,session) { 
  
 cantdeut=reactive({ #number of residues that cannot withhold deuteration as input by user
    as.numeric(input$cantdeutbutton)
  })
 timepoints=reactive({#number of timepoints as input by user
    input$timepoints 
  })
  replicates=reactive({ #number of replicates as input by user
    input$replicates
  })
  states=reactive({ #number of states as input by user
    input$states
  })
  
  significancelevel=reactive({ #significance level as input by user
    input$significancelevel
  })
  
  labelingtimepoints=reactive({ #labeling time points as input by user
    as.numeric(unlist(strsplit(input$labelingtimepoints,",")))
  })
  
  output$state1box<-renderUI({ #drop down menu to select reference state
    req(initialdata())
    selectInput(inputId="state1",label = "Select reference state:",choices ="", selected=NULL)
  })
  
  output$state2box<-renderUI({ #drop down menu to select protein state 2
    req(initialdata())
    selectInput(inputId="state2",label = "Select state 2:",choices ="", selected=NULL)
  })
  
  observeEvent(input$refresh, #Clear all data button
               session$reload()
               )
  
  output$analyze<-renderUI({ #add Calculate button when all files are uploaded
    req(input$hdexaminerfile)
    req(input$fastafile)
    req(input$timepoints)
    req(input$replicates)
    req(input$states)
    req(input$significancelevel)
    req(input$labelingtimepoints)
    
    actionButton("start","Calculate")
  })
  
  initialdata=eventReactive(input$start,{ #Load initial data when the 'start' button is pressed
    req(input$hdexaminerfile) # Ensure an HDExaminer file is uploaded
    req(input$fastafile) # Ensure a FASTA file is uploaded
    req(input$timepoints) # Ensure timepoints are provided
    req(input$replicates) # Ensure replicates are provided
    req(input$states) # Ensure states are provided
    req(input$significancelevel) # Ensure a significance level is provided
    req(input$labelingtimepoints) # Ensure labeling timepoints are provided
    
    if(length(input$datacheck)==0){ # If no specific data format is selected, assume HDExaminer format
      rawdata <- read.csv(input$hdexaminerfile$datapath, # Import HDExaminer output file
                          header=FALSE,
                          stringsAsFactors=FALSE)
      
      # Delete unnecessary columns from the output file and format the data properly for analysis
      HDXdata = rawdata[-1,] # Remove first row
      colnames(HDXdata)=HDXdata[1,] # Assign column names
      HDXdata = HDXdata[-1,] # Remove now redundant first row
      colnames(HDXdata)[2]="pepnumber" # Rename second column to 'pepnumber'
      HDXdata=HDXdata[-c(6,8)] # Remove specific unnecessary columns
      
      # Format numeric data and remove redundant columns dynamically based on timepoints and replicates
      for(x in 1:(timepoints()*replicates())){
        HDXdata=HDXdata[-c(x+6,x+7,x+9,x+10,x+11,x+12,x+13)]
        HDXdata[, x+6] = sapply(HDXdata[,x+6], as.numeric)
      }
    } else {
      if(input$datacheck=="WATERS/Dynamix output file"){ # Import data if using WATERS/Dynamix output file
        watersdata <- read.csv(input$hdexaminerfile$datapath, 
                               header=TRUE,
                               stringsAsFactors=FALSE)
        
        # Remove unnecessary columns from WATERS output
        watersdata=watersdata[-c(5,6,7,8,11,13,14)]
        watersdata['Exposure']=watersdata['Exposure']*60 # Convert exposure time to seconds
        watersdata$Exposure=as.integer(watersdata$Exposure)
        
        watersdata = watersdata %>% group_by(Sequence,State,z,Exposure) %>% dplyr::mutate(Exposure=paste0(Exposure,"s_",row_number()))
        watersdata = watersdata %>% pivot_wider(names_from=Exposure,values_from=Center) # Reshape data
        watersdata = watersdata %>% relocate(State) # Move State column to front
        
        # Split data
        wdata=subset(watersdata, select=c(1:6))
        wdata2=subset(watersdata, select=c(7:length(watersdata)))
        
        # Create column names for replicates and merge with data
        names=c()
        for(i in 1:replicates()){
          names=append(names, paste0(labelingtimepoints(),"s_",i))
        }
        df=data.frame(matrix(ncol=timepoints()*replicates(),nrow=0))
        names(df)=names
        wdata2[setdiff(names(df),names(wdata2))]<-NA
        wdata2=wdata2[,order(as.numeric(sub("\\D*(\\d+).*","\\1",colnames(wdata2))))]
        watersdata=cbind(wdata,wdata2)
        
        # Column names and remove unnecessary data
        names(watersdata)[2]="pepnumber"
        names(watersdata)[6]="Charge"
        watersdata=data.frame(watersdata, check.names = FALSE)
        nonddata=watersdata %>% select(starts_with("0s_"))
        nonddata$mean=rowMeans(nonddata,na.rm=TRUE) 
        
        
        watersdata=cbind(watersdata,mean=nonddata$mean)
        initialdata=watersdata %>% select(-starts_with("0s_"))
        initialdata = initialdata %>% drop_na(mean)
        for (i in 7:(timepoints()*replicates()+6)){
          initialdata[,i]=(initialdata[,i]-initialdata$mean)*initialdata[,6]
        }
        HDXdata= initialdata %>% select(-mean)
        HDXdata$pepnumber="" 
      }else{
        
        HDXdata=read.csv(input$hdexaminerfile$datapath, 
                         header=TRUE,
                         stringsAsFactors=FALSE)
        HDXdata[1:6]=lapply(HDXdata[1:6],as.character)
        HDXdata$pepnumber=""
      }
    }
    
    # Remove rows with excessive missing values
    HDXdata=HDXdata[rowSums(is.na(HDXdata[ , 6:((timepoints()*replicates())+5)])) <(timepoints()*replicates())-(replicates()), ]
    
    # Filter out inconsistent data
    HDXdata=HDXdata%>%group_by(Sequence, Charge)%>%filter(n_distinct(State)==states())
    pepcharge=HDXdata[,6] # Store charge information separately
    HDXdata=HDXdata[,-6] # Remove charge column from main dataset
    HDXdata=data.frame(HDXdata)
    
    # Rename columns based on timepoints and replicates
    for(i in 0:(timepoints()-1)){
      for(j in 1:replicates()){
        colnames(HDXdata)[(i*replicates())+5+j]=paste(labelingtimepoints()[i+1],"s_",j,sep="")
      }
    }
    
    shinyjs::show("functional") 
    return(list(HDXdata=HDXdata,pepcharge=pepcharge)) # Return processed data
  })
  
  
  
  observe({ #alert message if two check boxes are selected at the same time
    if(length(input$datacheck)==2){
      showModal(modalDialog(title="Please select only ONE option",footer=NULL,easyClose=TRUE))
      updateCheckboxGroupInput(session, "datacheck",choices=c("WATERS/Dynamix output file", "Custom .csv output file"),selected = NULL)
      }
  })
  
  observe({ #drop down input updated with protein states names in the csv file
    req(initialdata())
    HDXdata=initialdata()$HDXdata
    states=unique(HDXdata$State)
    updateSelectInput(session,"state1",label = "Select reference state:",choices=states,selected=states[1])
  })
  
  observe({ #drop down input updated with protein states names in the csv file
    req(initialdata())
    HDXdata=initialdata()$HDXdata
    states=unique(HDXdata$State)
    updateSelectInput(session,"state2",label = "Select state 2:",choices = states[!states==input$state1])
  })
  
  variables=reactive({ #calculate all variables to be used later
    
    req(initialdata(), input$state1, input$state2, input$state1 != input$state2)
    
    sequence = read.fasta(input$fastafile$datapath, seqtype = "AA",whole.header=FALSE,seqonly=FALSE,as.string = TRUE) #Import FASTA file 
    sequence=gsub(" ","",sequence) #delete blank spaces in the sequence
    withProgress(message = 'Analyzing data', value=0, { 
      incProgress(1/10)
      
      HDXdata=initialdata()$HDXdata
      HDXdata$State = factor(HDXdata$State, levels=unique(HDXdata$State))
      HDXdatastates = list()
      
      otherstates = unique(HDXdata$State)
      otherstates = otherstates[!otherstates==input$state1]
      otherstates = otherstates[!otherstates==input$state2]
      
      HDXdatastates[[1]]=filter(HDXdata, State==input$state1) #filter data by state
      HDXdatastates[[2]]=filter(HDXdata, State==input$state2)
      
      if(states()>2){
      for(i in 3:states()){
        HDXdatastates[[i]]=filter(HDXdata, State==otherstates[i-2])
          }
        }
      
      
      data=add_column(initialdata()$HDXdata,initialdata()$pepcharge,.after="Sequence")
      data <- data %>% dplyr::relocate (Charge, .after = Sequence)
      TF1 <- data
      
      
      #data transformations for functional data fitting
      TF2 <- TF1 %>%
        subset(select = -c(`pepnumber`))
      
     TF2 <- TF2 %>% 
        pivot_longer(
          cols = contains("_"),
          names_to = c("Time", "Replicate"),
          values_transform = list( Replicate = as.double, Time = as.double),
          names_pattern = "X?(.*)s_(.*)",
          values_to = "D"
        )
      
      incProgress(2/10)
      TF2 <- TF2 %>% 
        dplyr::rename("pep_sequence"="Sequence") %>%
        dplyr::rename("pep_charge"="Charge") %>%
        dplyr::rename("d"="D") %>%
        dplyr::rename("hx_time"="Time") %>%
        dplyr::rename("replicate_cnt"="Replicate") %>%
        dplyr::rename("hx_sample"="State")
      
      TF2 <- TF2 %>% mutate_at(c('hx_time', 'replicate_cnt'), as.numeric)
      
      #removal of any additional charge states
      TF2 <- TF2 %>% dplyr::group_by(pep_sequence, pep_charge, hx_time, replicate_cnt, hx_sample) %>%
        dplyr::mutate(n = dplyr::n())
      TF2 <- TF2 %>% filter(n == 1)
      
      #generation of tidy format for qfeatures
      HDX_wide <- pivot_wider(data.frame(TF2),
                              values_from = d,
                              names_from = c("hx_time", "replicate_cnt", "hx_sample"),
                              id_cols = c("pep_sequence","pep_charge"))
      
      #generation of colnames for qFeatures
      new.colnames <- gsub("0_", "0rep", paste0("X", colnames(HDX_wide)[-c(1,2)]))
      new.colnames <- gsub("_", "cond", new.colnames)
      new.colnames <- gsub(" ", "", new.colnames)
      new.colnames <-  gsub(".", "", new.colnames, fixed = TRUE)
      
      #generation of qFeatures
      HDXqf <- parseDeutData(object = DataFrame(HDX_wide),
                             design = new.colnames,
                             quantcol = 3:((replicates()*timepoints()*states())+2)) # Check had 3:50
      incProgress(3/10)
      subHDXqf=list()
      if(states()>2){
        subTF2=TF1 %>% subset(select = -c(`pepnumber`))
        sub1 = subTF2 %>% filter(State==input$state1)
        sub2 = subTF2 %>% filter(State==input$state2)
        subTF2 = rbind(sub1,sub2)
        subTF2 <- subTF2 %>% 
          pivot_longer(
            cols = contains("_"),
            names_to = c("Time", "Replicate"),
            values_transform = list( Replicate = as.double, Time = as.double),
            names_pattern = "X?(.*)s_(.*)",
            values_to = "D"
          )
        subTF2 <- subTF2 %>% 
          dplyr::rename("pep_sequence"="Sequence") %>%
          dplyr::rename("pep_charge"="Charge") %>%
          dplyr::rename("d"="D") %>%
          dplyr::rename("hx_time"="Time") %>%
          dplyr::rename("replicate_cnt"="Replicate") %>%
          dplyr::rename("hx_sample"="State")
        
        subTF2 <- subTF2 %>% mutate_at(c('hx_time', 'replicate_cnt'), as.numeric)
        
        #removal of any additional charge states
        subTF2 <- subTF2 %>% dplyr::group_by(pep_sequence, pep_charge, hx_time, replicate_cnt, hx_sample) %>%
          dplyr::mutate(n = dplyr::n())
        subTF2 <- subTF2 %>% filter(n == 1)
        
        #generation of tidy format for qfeatures
        subHDX_wide <- pivot_wider(data.frame(subTF2),
                                values_from = d,
                                names_from = c("hx_time", "replicate_cnt", "hx_sample"),
                                id_cols = c("pep_sequence","pep_charge"))
        
        #generation of colnames for qFeatures
        subnew.colnames <- gsub("0_", "0rep", paste0("X", colnames(subHDX_wide)[-c(1,2)]))
        subnew.colnames <- gsub("_", "cond", subnew.colnames)
        subnew.colnames <- gsub(" ", "", subnew.colnames)
        subnew.colnames <-  gsub(".", "", subnew.colnames, fixed = TRUE)
        
        #generation of qFeatures
        subHDXqf <- parseDeutData(object = DataFrame(subHDX_wide),
                               design = subnew.colnames,
                               quantcol = 3:((replicates()*timepoints()*2)+2)) # Check had 3:50
      }
      
      heatmap=pheatmap(t(assay(HDXqf)), #generation of functional data analysis heatmap
                       cluster_rows = FALSE, 
                       cluster_cols = FALSE,
                       color = brewer.pal(n = 9, name = "BuPu"),
                       main = "Protein: Deuterium Incorporation", 
                       fontsize = 8,
                       legend_breaks = c(0, 1, 2, 3, 4, 5, 6, max(assay(HDXqf))),
                       legend_labels = c("0", "1", "2", "3", "4", "5", "6", "Incorporation")) 
      
      var=c()
      mean1=c()
      sd1=c()
      sd2=c()
      npooledSD=0
      dpooledSD=0
      n=0
      
      for (z in 1:2){ #Calculate the average and SD for each labeling time point and add it to the data frame
        HDXdatastates[[z]][["pepnumber"]]=1:nrow(HDXdatastates[[z]])
        for(i in 1:timepoints()){            
          for(h in 1:nrow(HDXdatastates[[z]])) { 
            for(k in 1:replicates()){
              var=append(var,HDXdatastates[[z]][h,5+(replicates()*(i-1))+k])
                }
            mean1=append(mean1, mean(as.numeric(var), na.rm= TRUE))
            sd1=append(sd1, sd(as.numeric(var), na.rm =TRUE ))
            sd2=sd(as.numeric(var), na.rm =TRUE )
            n=length(na.omit(var))
            var=c()
            if (is.na(sd2)){
              sd2=0
                }
            npooledSD=npooledSD + ((sd2)^2)*(n-1)
            dpooledSD=dpooledSD +(n-1)
            } 
          HDXdatastates[[z]]=cbind(HDXdatastates[[z]],mean1)
          colnames(HDXdatastates[[z]])[ncol(HDXdatastates[[z]])]=paste('Average D at ',labelingtimepoints()[i],'sec')
          HDXdatastates[[z]]=cbind(HDXdatastates[[z]],sd1)
          colnames(HDXdatastates[[z]])[ncol(HDXdatastates[[z]])]=paste('SD for ',labelingtimepoints()[i],'sec')
          mean1=c()
          sd1=c()
          }
      }
      if(input$showallstates==TRUE){
        for (z in 3:states()){ #Calculate the average and SD for each labeling time point and add it to the data frame
          HDXdatastates[[z]][["pepnumber"]]=1:nrow(HDXdatastates[[z]])
          for(i in 1:timepoints()){            
            for(h in 1:nrow(HDXdatastates[[z]])) { 
              for(k in 1:replicates()){
                var=append(var,HDXdatastates[[z]][h,5+(replicates()*(i-1))+k])
                 }
              mean1=append(mean1, mean(as.numeric(var), na.rm= TRUE))
              sd1=append(sd1, sd(as.numeric(var), na.rm =TRUE ))
              sd2=sd(as.numeric(var), na.rm =TRUE )
              n=length(na.omit(var))
              var=c()
              if (is.na(sd2)){
                sd2=0
                }
               } 
            HDXdatastates[[z]]=cbind(HDXdatastates[[z]],mean1)
            colnames(HDXdatastates[[z]])[ncol(HDXdatastates[[z]])]=paste('Average D at ',labelingtimepoints()[i],'sec')
            HDXdatastates[[z]]=cbind(HDXdatastates[[z]],sd1)
            colnames(HDXdatastates[[z]])[ncol(HDXdatastates[[z]])]=paste('SD for ',labelingtimepoints()[i],'sec')
            mean1=c()
            sd1=c()
          }
        }
      }
      incProgress(4/10)
      pooledSD=sqrt(npooledSD/dpooledSD) #calculate pooled SD
      SEM=sqrt(2*(pooledSD^2)/replicates()) # calculate standard error of the mean
      threshold=qt(p=significancelevel()/2, df=(2*replicates())-2,lower.tail = FALSE)*SEM #Calculate significance threshold based on SEM and significance level
      
      parameters=data.frame(PooledSD= pooledSD, StrdErrorMean=SEM, alpha=significancelevel(),tValUsed=qt(p=significancelevel()/2, df=(2*replicates())-2,lower.tail = FALSE),StatisticalThreshold=threshold)
      parameters[5,1]="Peptides presenting protection and deprotection at different labeling times are colored yellow (Check data in HDExaminer)"
      
      listtestresults=data.frame(HDXdatastates[[1]][2],HDXdatastates[[1]][3],HDXdatastates[[1]][4],HDXdatastates[[1]][5]) #create data frame to save statistics results. Include pepnumber and sequence
      clusteringresults=listtestresults
      
      dataset1=c()
      dataset2=c()
      difference=c()
      ttest=c()
      normvalues=c()
      
      #t-test between two states 
      for(m in 1:timepoints()){           
        for(n in 1:nrow(HDXdatastates[[1]])){ 
          for(o in 1:replicates()){
            dataset1=append(dataset1,HDXdatastates[[1]][n,5+((m-1)*replicates())+o])
            dataset2=append(dataset2,HDXdatastates[[2]][n,5+((m-1)*replicates())+o])
             }
          
          if(sum(!is.na(dataset1))<2||sum(!is.na(dataset2))<2){
            ttest=append(ttest, NA)
            difference=append(difference, NA)
            normvalues=append(normvalues,NA)
             } else {
            
            ttest=append(ttest, t.test(na.omit(dataset2),na.omit(dataset1),var.equal = FALSE,conf.level = significancelevel()/2,alternative = "two.sided")$p.value)
            difference=append(difference, mean(dataset2,na.rm = TRUE)-mean(dataset1,na.rm=TRUE))
            normvalues=append(normvalues,((mean(dataset2,na.rm = TRUE)-mean(dataset1,na.rm=TRUE))/(str_length(HDXdatastates[[1]][n,5])-str_count(HDXdatastates[[1]][n,5],"P")-cantdeut())))
            
            }
          dataset1=c()
          dataset2=c()
        } 
        #save all results in dataframes 
        listtestresults=cbind(listtestresults, difference)
        colnames(listtestresults)[ncol(listtestresults)]=paste0('Delta D at ',labelingtimepoints()[m],'sec')
        
        listtestresults=cbind(listtestresults, ttest)
        colnames(listtestresults)[ncol(listtestresults)]=paste0('t test ',labelingtimepoints()[m],'sec')
        
        clusteringresults=cbind(clusteringresults, difference)
        colnames(clusteringresults)[ncol(clusteringresults)]=paste0('Delta D at ',labelingtimepoints()[m],'sec')
        
        clusteringresults=cbind(clusteringresults, normvalues)
        colnames(clusteringresults)[ncol(clusteringresults)]=paste0('Cluster',labelingtimepoints()[m],'sec')
        
        clusteringresults=cbind(clusteringresults, ttest)
        colnames(clusteringresults)[ncol(clusteringresults)]=paste0('t test',labelingtimepoints()[m],'sec')
        
        difference=c()
        ttest=c()
        normvalues=c()
      }  
      
      incProgress(5/10)
      significant=c()
      var2=c()
      
      for(i in 1:nrow(listtestresults)){
        for(j in 1:timepoints()){
          if (is.na(listtestresults[i,(2*j)+3]) || is.na(listtestresults[i,(2*j)+4])){
            var2=append(var2,0)
              } else{
            if (((listtestresults[i,(2*j)+3])>=threshold) && (listtestresults[i,(2*j)+4]<significancelevel())){
              var2=append(var2,1)}else{
                if (((listtestresults[i,(2*j)+3])<=-threshold) && (listtestresults[i,(2*j)+4]<significancelevel())){
                  var2=append(var2,-1)}else{
                    var2=append(var2,0)
                  }
              }
          }
        }
        if ((length(var2[var2>0]))>0 && (length(var2[var2<0])==0)){
          significant=append(significant,1)
            }else{
            if ((length(var2[var2<0]))>0 && (length(var2[var2>0])==0)){
            significant=append(significant,-1)
              }else{
              if ((length(var2[var2>0]))>0 && (length(var2[var2>0])>0)){
              significant=append(significant, 5)}else{
                significant=append(significant,0)
                }
             }
          }
        var2=c()
      }
      listtestresults=cbind(listtestresults, significant)
      colnames(listtestresults)[ncol(listtestresults)]="SignificantResults"
      significant=c()
      
      exchangableamides=c()
      normalization=c()
      allvalues=c()
      for (i in 1:nrow(HDXdatastates[[1]])){ #normalization of differences for k-means clustering
        peptideseq=listtestresults[i,4]
        exchangableamides[i]=str_length(peptideseq)-str_count(peptideseq,"P")-cantdeut()
        for(h in 1:timepoints()){
          if (is.na(listtestresults[i,(h*2)+3])){}else{
            normalization[i]=(listtestresults[i,(h*2)+3]/exchangableamides[i])
            }
          }
        allvalues=append(allvalues, normalization)
      }
      kmeansclust=kmeans(abs(na.omit(allvalues)),3,iter.max=200) #k-means clustering of normalized differences
      centers=sort(kmeansclust$centers)
      allvalues=as.data.frame(allvalues) 
      
      limit1=-(((centers[3]-centers[2])/2)+centers[2]) #k-means cluster limits
      limit2=-(((centers[2]-centers[1])/2)+centers[1])
      limit3=-limit2
      limit4=-limit1 
      incProgress(8/10)
      for (l in 1:nrow(HDXdatastates[[1]])){
        for(m in 1:timepoints()){
          if(is.na(clusteringresults[l,(m*3)+3])){clusteringresults[l,(m*3)+3]=1}else{
            if(clusteringresults[l,(m*3)+3]>limit4){clusteringresults[l,(m*3)+3]=3}else{
              if(clusteringresults[l,(m*3)+3]>limit3){clusteringresults[l,(m*3)+3]=2}else{
                if(clusteringresults[l,(m*3)+3]>limit2){clusteringresults[l,(m*3)+3]=1}else{
                  if(clusteringresults[l,(m*3)+3]>limit1){clusteringresults[l,(m*3)+3]=2}else{clusteringresults[l,(m*3)+3]=3}}}}}
          
          if(is.na(abs(clusteringresults[l,(m*3)+2]))){clusteringresults[l,(m*3)+3]=1}else{
            if(abs(clusteringresults[l,(m*3)+2])<threshold){clusteringresults[l,(m*3)+3]=1}}
          if(is.na(clusteringresults[l,(m*3)+4])){clusteringresults[l,(m*3)+3]=1}else{
            if(-log10(clusteringresults[l,(m*3)+4])<(-log10(significancelevel()))){clusteringresults[l,(m*3)+3]=1}
          }
        }
      }
      
      clustvec=c()
      for (i in 1:nrow(clusteringresults)){
        for(j in 1:timepoints()){
          if(is.na(clusteringresults[i,(j*3)+2])){
            clusteringresults[i,(j*3)+2]=0
            }
          if(is.na(clusteringresults[i,(j*3)+4])){
            clusteringresults[i,(j*3)+4]=0
            }
          if(abs(clusteringresults[i,(j*3)+2])>=threshold  && (clusteringresults[i,(j*3)+4]<=significancelevel())){
          clustvec=append(clustvec,clusteringresults[i,(3*j)+3])
            }else{
            clustvec=append(clustvec,1)
          }
        }
        clusteringresults[i,(3*timepoints())+5]=max(clustvec)
        clustvec=c()
      }
      colnames(clusteringresults)[(3*timepoints())+5]="Highest significant cluster in peptide"
      a=ncol(clusteringresults)
      clusteringresults[1,a+1]="Cluster sizes"
      colnames(clusteringresults)[a+1]=""
      clusteringresults[2,(a+1):(a+3)]=kmeansclust$size
      colnames(clusteringresults)[a+2]=""
      colnames(clusteringresults)[a+3]=""
      clusteringresults[3,a+1]="Cluster 1 (abs values)"
      clusteringresults[4,a+1]="Cluster 2 (abs values)"
      clusteringresults[5,a+1]="Cluster 3 (abs values)"
      clusteringresults[3,a+2]=paste("From",0,"to",signif(limit3,3))
      clusteringresults[4,a+2]=paste("From",signif(limit3,3),"to",signif(limit4,3))
      clusteringresults[5,a+2]=paste("From",signif(limit4,3),"to max value")
      clusteringresults[6,a+1]="k-means clustering iterations"
      clusteringresults[6,a+2]=kmeansclust$iter
      clusteringresults[7,a+1]="Datapoints in cluster 1 are considered insignificant or negligible within your dataset"
      clusteringresults[8,a+1]="Datapoints in cluster 2 are intermediate effects within your dataset"
      clusteringresults[9,a+1]="Datapoints in cluster 3 are strong effects within your dataset"
      incProgress(9/10)
    })
    
    if(nchar(sequence) < max(as.numeric(HDXdatastates[[1]]$End))){
      shinyalert(
        title = "",
        text = "You FASTA file seems to have less residues than your experiment. Please check your FASTA file to avoid errors",
        size = "xs", 
        closeOnEsc = FALSE,
        closeOnClickOutside = TRUE,
        html = FALSE,
        type = "warning",
        showConfirmButton = FALSE,
        showCancelButton = FALSE,
        timer = 0,
        imageUrl = "",
        animation = TRUE
      )
    }else{}
    return(list(HDXdatastates=HDXdatastates,listtestresults=listtestresults,parameters=parameters,threshold=threshold,clusteringresults=clusteringresults,sequence=sequence,centers=centers,allvalues=allvalues,HDXqf=HDXqf,heatmap=heatmap,subHDXqf=subHDXqf))
  })
  
  output$heatmap= renderPlot({ #show heatmap
    req(variables())
    return(variables()$heatmap)
  })
  
  output$heatbutton<-downloadHandler( #show download button
    filename=function(){
      paste0("Heat map",input$heattype)
    },
    content=function(file){
      ggsave(file,variables()$heatmap,height=7,width=15)
    }
  )
  realfit<-reactiveVal(value=1)
  
  observe({ #functional data analysis parameters boxes update after selection of different equation
    if(input$selectedeq=="equationwithoutd"){
      shinyjs::disable("d_input")
      shinyjs::enable("p_input")
      updateNumericInput(session,"d_input",value="NULL")
      updateNumericInput(session,"p_input",value=1)
    }
    if(input$selectedeq=="equationwithd"){
      shinyjs::enable("d_input")
      shinyjs::enable("p_input")
      updateNumericInput(session,"d_input",value=0)
      updateNumericInput(session,"p_input",value=1)
    }
    if(input$selectedeq=="equationwithdwithoutP"){
      shinyjs::enable("d_input")
      shinyjs::disable("p_input")
      updateNumericInput(session,"d_input",value=0)
      updateNumericInput(session,"p_input",value="NULL")
    }
    if(input$selectedeq=="equationwithoutdorP"){
      shinyjs::disable("d_input")
      shinyjs::disable("p_input")
      updateNumericInput(session,"d_input",value="NULL")
      updateNumericInput(session,"p_input",value="NULL")
    }
  })
  
 counter<-reactiveVal(value=0)
  
 observeEvent(input$reset ,{
   counter(counter()+1)
   return()
   })
  
 event_trigger <- reactive({
   list(input$state1,input$state2)
 })
 
 observeEvent(ignoreInit = TRUE, event_trigger(),{
   click("reset")
 })
 
  fitting=eventReactive(input$startfitting,{ #functional data fitting saved in a reactive variable
    
    req(variables())
    all_peptides <- rownames(variables()$HDXqf)[[1]]
    
    if(states()>2){
    suball_peptides <- rownames(variables()$subHDXqf)[[1]]
    }
    
    if(input$selectedeq=="equationwithd"){
      starting_parameters <- list(a = input$a_input, b = input$b_input, d=input$d_input , p=input$p_input )
      } else if (input$selectedeq=="equationwithdwithoutP"){
      starting_parameters <- list(a = input$a_input, b = input$b_input, d=input$d_input)
      }else if(input$selectedeq=="equationwithoutd"){
      starting_parameters <- list(a = input$a_input, b = input$b_input , p=input$p_input )
      } else{
      starting_parameters <- list(b = input$b_input, a = input$a_input)
      }
    subresults=list()
    withProgress(message = 'Fitting data', min=0,max=1, {
      incProgress(1/10)
      if(input$selectedeq=="equationwithd"){ #fit the data using the equation selected by user
        results <- analyse_kinetics(data = variables()$HDXqf, #functional data fitting
                                    method = "fit",
                                    formula = value ~ a * (1-exp(-b*(timepoint)^p))+d,
                                    peptide_selection = all_peptides,
                                    start = starting_parameters,
                                    maxAttempts = 20)
        if(states()>2){
          subresults <- analyse_kinetics(data = variables()$subHDXqf, #functional data fitting
                                      method = "fit",
                                      formula = value ~ a * (1-exp(-b*(timepoint)^p))+d,
                                      peptide_selection = suball_peptides,
                                      start = starting_parameters,
                                      maxAttempts = 20)
        }
        
       }else{ if(input$selectedeq=="equationwithoutd"){
         results <- analyse_kinetics(data = variables()$HDXqf, #functional data fitting
                                    method = "fit",
                                    formula = value ~ a * (1-exp(-b*(timepoint)^p)),
                                    peptide_selection = all_peptides,
                                    start = starting_parameters,
                                    maxAttempts = 20) 
         if(states()>2){
           subresults <- analyse_kinetics(data = variables()$subHDXqf, #functional data fitting
                                       method = "fit",
                                       formula = value ~ a * (1-exp(-b*(timepoint)^p)),
                                       peptide_selection = suball_peptides,
                                       start = starting_parameters,
                                       maxAttempts = 20) 
         }
         
        }else{ if(input$selectedeq=="equationwithdwithoutP"){
          results <- analyse_kinetics(data = variables()$HDXqf, #functional data fitting
                                    method = "fit",
                                    formula = value ~ a * (1-exp(-b*(timepoint)))+d,
                                    peptide_selection = all_peptides,
                                    start = starting_parameters,
                                    maxAttempts = 20)
          if(states()>2){
            subresults <- analyse_kinetics(data = variables()$subHDXqf, #functional data fitting
                                        method = "fit",
                                        formula = value ~ a * (1-exp(-b*(timepoint)))+d,
                                        peptide_selection = suball_peptides,
                                        start = starting_parameters,
                                        maxAttempts = 20) 
          }
          
        }else{
         results <- analyse_kinetics(data = variables()$HDXqf, #functional data fitting
                                    method = "fit",
                                    formula = value ~ a * (1-exp(-b*(timepoint))),
                                    peptide_selection = all_peptides,
                                    start = starting_parameters,
                                    maxAttempts = 20) 
         if(states()>2){
           subresults <- analyse_kinetics(data = variables()$subHDXqf, #functional data fitting
                                       method = "fit",
                                       formula = value ~ a * (1-exp(-b*(timepoint))),
                                       peptide_selection = suball_peptides,
                                       start = starting_parameters,
                                       maxAttempts = 20) 
         }
         
          }
        }
      }
      incProgress(6/10)
      graphics_kinetics <- visualise_hdx_data(results, type = "kinetics") #fitted kinetic plots
      if(states()==2){graphics_forest<-visualise_hdx_data(results, type = "forest") #forest plots
      }else{graphics_forest<-visualise_hdx_data(subresults, type = "forest")}
      summary=BiocGenerics::lapply(results[["fitted_models"]]@statmodels,summary) #summary of models
      incProgress(9/10)
    })
    shinyjs::show("fdabox")
    realfit(2)
    counter(0)
    
    return(list(kinetics=graphics_kinetics,forest=graphics_forest,fittingparameters=results$functional_analysis@results,summary=summary,results=results,subresults=subresults))
  })

  singlefittedplot=reactive({ #single fitted plot for selected peptide
    req(fitting())
    graphics_kinetics=fitting()$kinetics
    graphics_forest=fitting()$forest
    summary=fitting()$summary
    wanted=as.numeric(str_split(input$selpeptide2,":")[[1]][1])
    null=data.frame(coef(summary[[wanted]]$null),check.names=F)
    null$variable <- substr(rownames(null),1,1)
    null= null %>% select(variable, everything())
    null[,6]="Null"
    colnames(null)[6]="State"
    
    
    for (i in 1:states()){
      df=data.frame(coef(summary[[wanted]][[i+1]]),check.names=F)
      df$variable <- substr(rownames(df),1,1)
      df= df %>% select(variable, everything())
      df[,6]=paste0("State",i)
      colnames(df)[6]="State"
      null=rbind(null,df)
      } 
   if(is.null(graphics_forest[[wanted]])){
     plot2=NULL
   }else{
     plot2=(graphics_forest[[wanted]]+ggtitle(paste0("Forest plot ",input$state1," vs. ", input$state2))+theme(plot.title = element_text(hjust = 0.5)))
   }
    return(list(plot1=graphics_kinetics[[wanted]],plot2=plot2,summary=null))
  })
  
  output$kineticplot=renderPlot({
    print(singlefittedplot()$plot1)
  })
  
  output$forestplot=renderPlot({
    test=counter()
    text=paste("\n Your reference/comparison states have changed. Please refit the data.\n")
    if(test != 0){ggplot() + 
        annotate("text", x = 4, y = 25, size=7, label = text) + 
        theme_void()}else{
    singlefittedplot()$plot2}
  })  
  
  output$summary=renderTable({
    singlefittedplot()$summary
    },digits=3)
  
  observe({
    req(fitting())
    labels=c()
    for(i in 1:length(fitting()$kinetics)){
      labels=append(labels,paste0(i,": ",fitting()$kinetics[[i]]$data$rowname[1]))
      }
    updateSelectInput(session,"selpeptide2",label = "Select peptide:",choices = labels)
  })
  
  output$kineticbutton<-downloadHandler(
    filename=function(){
      paste0("Fitted plot",input$kinetictype)
      },
    content=function(file){
      pdf(file,width=7,height=7)
      print(singlefittedplot()$plot1)
      print(singlefittedplot()$plot2)
      dev.off()
    }
  )
  kineticdocument=reactive({
    graphs=list()
    a=1
    withProgress(message = 'Creating document', min=0.5,max=1, { 
      for(i in 1:length(fitting()$kinetics)){
        graphs[[a]]=fitting()$kinetics[[i]]
        graphs[[a+1]]=fitting()$forest[[i]]
        a=a+2
       }
   return(marrangeGrob(graphs,ncol=2,nrow=2,top=NULL,layout_matrix = matrix(1:4, 2, 2, TRUE)))
    })   
  })
  observeEvent(input$downloadallkinetic,{
    
    showModal(modalDialog(title="Please wait. Creating your document...",easyClose=FALSE,
                            modalButton("Cancel"),
                            downloadButton("allkineticbutton", "Download all kinetic plots"),
                          footer=NULL
                     )
                )
    shinyjs::disable(id="allkineticbutton")
    a=kineticdocument()
    shinyjs::enable(id="allkineticbutton")
   })
  
  output$allkineticbutton<-downloadHandler(
    filename="Allfittedplots.pdf",
    content=function(file){
      withProgress(message = 'Exporting', min=0,max=1, { 
        ggsave(file,kineticdocument(),width=14,height=8)
      })
    }
  )
  
  observe({
    req(initialdata())
    updateSelectInput(session,"seltimemanhattan",label = "Select labeling time:",choices = as.numeric(unlist(strsplit(input$labelingtimepoints,","))))
  })
  
  manhattan=reactive({ #calculation of Manhattan plot
    req (fitting())
    if(states()==2){
    results=fitting()$results
    HDXqf=variables()$HDXqf
    }else{
      results=fitting()$subresults
      HDXqf=variables()$subHDXqf
    }
    
    data=add_column(initialdata()$HDXdata,initialdata()$pepcharge, .after="Sequence")
    labelingtimes=as.numeric(unlist(strsplit(input$labelingtimepoints,",")))
    time=1
    for(i in 1:input$timepoints){
      if(input$seltimemanhattan==labelingtimes[i]){time=i}
    }
    
    diffdata <- rowMeans(assay(HDXqf)[,(((replicates()*(time-1))+1):(replicates()*time))],na.rm=TRUE) - rowMeans(assay(HDXqf)[,(((replicates()*timepoints())+(replicates()*(time-1))+1):((replicates()*timepoints())+(replicates()*time)))],na.rm=TRUE)
    out<-processFunctional(object=HDXqf, params=results$fitted_models)
    data$names=paste0(data$Sequence,"_",data$Charge)
    data$peptide=paste0("[",data$Start,":",data$End,"]")
    data=data %>% distinct(names, .keep_all=TRUE)
    sequences = rownames(out@results)
    
    df=data.frame(out@results)
    df$names=rownames(df)
    diffdata <-data.frame(diffdata)
    diffdata$names=rownames(diffdata)
    diffdata=drop_na(diffdata)
    signs <- inner_join(diffdata,df)
    names <-inner_join(data,df)
    all <- inner_join(signs,names)
    all <- all %>% distinct(names, .keep_all=TRUE)
    
    manhplot=ggplot(all, aes(x=1:nrow(all), y=-log10(ebayes.fdr),color=factor(as.vector(sign(diffdata))))) + geom_point(size=3.5) +
      scale_color_manual(breaks=c(-1,1,0),values=c(input$colorflex,input$colorprot,"gray"),labels=c("deprotected","protected","no difference")) + scale_x_continuous(expand = c(0, 1.2),guide = guide_axis(n.dodge=3),breaks=1:length(all$fitcomplete),labels=all$peptide) + theme_classic()+
      geom_hline(yintercept = -log10(input$significancelevel), linetype = "dashed", colour = "red",linewidth=2)+ xlab("Peptide")+ ylab("-log(ebayes FDR)")+
      theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + labs(color="Effect") +
      ggtitle(paste0("Manhattan plot ",input$state1," vs. ", input$state2))+theme(plot.title = element_text(hjust = 0.5))
    
    return(manhplot)
  })
  
  output$manhattanplot <- renderPlot({
    test=counter()
    text=paste("\n Your reference/comparison states have changed. Please refit the data.\n")
    if(test != 0){ggplot() + 
        annotate("text", x = 4, y = 25, size=10, label = text) + 
        theme_void()}else{
    print(manhattan())}
  })
  
  output$manhattanbutton<-downloadHandler(
    filename=function(){
      paste0("ManhattanPlot",input$seltimemanhattan,"s",input$manhattantype)
    },
    content=function(file){
      ggsave(file,(manhattan()+ggtitle("Manhattan plot")+theme(plot.title = element_text(hjust = 0.5,size=20))),height=7,width=14)
    }
  )
#-------------------## First tab statistical analysis##----------------  
observe({
    if(input$clusteringsel==TRUE){
      shinyjs::show("panel")
      }  
  if(input$clusteringsel==FALSE){
    shinyjs::hide("panel")
      }
  })
  
observeEvent(input$start,{
  shinyjs::show("panel2")
  })  

volcanoplot=reactive({ #calculate volcano plot (hybrid significance test)
    
    req(variables())
    threshold=variables()$threshold
    listtestresults=variables()$listtestresults
    clusteringresults=variables()$clusteringresults
    shinyjs::toggle()
    maxvalue=c()
    maxpvalue=c()
    ymaxlim=0
    for(i in 1:timepoints()){
      maxvalue=append(maxvalue,listtestresults[,(2*i)+3])
      maxpvalue=append(maxpvalue, listtestresults[,(2*i)+4])
       }
    volcanoplot= ggplot(data=clusteringresults)+lapply(seq(timepoints()),function(f){
      if(input$clusteringsel==TRUE){
        geom_point(aes(x=clusteringresults[,(f*3)+2], y=-log10(clusteringresults[,(f*3)+4]),color=factor(clusteringresults[,(f*3)+3])),na.rm=TRUE)}else{
          geom_point(aes(x=clusteringresults[,(f*3)+2], y=-log10(clusteringresults[,(f*3)+4])),na.rm=TRUE,color="black")}
      }) 
    if(-log10(min(abs(maxpvalue),na.rm=TRUE))<=-log10(significancelevel())){
      ymaxlim=-log10(significancelevel())+0.1}else{ymaxlim=-log10(min(abs(maxpvalue),na.rm=TRUE))
       }
    volcanoplot =volcanoplot+ 
      geom_segment(x=threshold, xend=100,y=-log10(significancelevel()),yend=-log10(significancelevel()),color = "red",linewidth=1.1, linetype="dashed")+
      geom_segment(x=-threshold, xend=-100,y=-log10(significancelevel()),yend=-log10(significancelevel()),color = "red",linewidth=1.1, linetype="dashed")+
      geom_segment(x=threshold, xend=threshold,y=-log10(significancelevel()),yend=100,color = "red",linewidth=1.1, linetype="dashed")+
      geom_segment(x=-threshold, xend=-threshold,y=-log10(significancelevel()),yend=100,color = "red",linewidth=1.1, linetype="dashed")+
      xlab("mass difference(Da)")+ylab("-log(p value)")+ xlim((0-max(abs(maxvalue),na.rm=TRUE)),(0+max(abs(maxvalue),na.rm=TRUE))) +
      ylim(0,ymaxlim)+theme_bw()+
      theme(axis.title=element_text(size=12), panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    if(input$clusteringsel==TRUE){
      volcanoplot=volcanoplot+scale_colour_manual(values=c("gray50",input$colorint,input$colorstr))+theme(legend.position="none")
      }else{}
    
    return(volcanoplot)
  })
  
  output$volcanoplot=renderPlot({
    print(volcanoplot())
  })
  
  output$binsize<- renderUI({
    req(input$clusteringsel==TRUE)
    selectInput(inputId = "selbinsize", label = "Select number of bins:",c("little","few","many"),selected="little")
  })
  
  histogram=reactive({ #histogram of normalized differences
    req(input$selbinsize)
    centers=variables()$centers
    allvalues=variables()$allvalues
    
    limit1=-(((centers[3]-centers[2])/2)+centers[2])
    limit2=-(((centers[2]-centers[1])/2)+centers[1])
    limit3=-limit2
    limit4=-limit1  
    bins=switch(input$selbinsize, little=15L,few=50L,many=100L)
    bins=seq(min(allvalues),max(allvalues),length.out=bins)
    hist=ggplot(allvalues,aes(x=allvalues))+
      geom_histogram(data=subset(allvalues, allvalues<=limit1),color="black",fill=input$colorstr,breaks=bins)+
      geom_histogram(data=subset(allvalues, allvalues<=limit2 & allvalues>limit1),color="black",fill=input$colorint,breaks=bins)+
      geom_histogram(data=subset(allvalues, allvalues<=limit3 & allvalues>limit2),color="black",fill="gray50",breaks=bins)+
      geom_histogram(data=subset(allvalues, allvalues<=limit4 & allvalues>limit3),color="black",fill=input$colorint,breaks=bins)+
      geom_histogram(data=subset(allvalues, allvalues>limit4),color="black",fill=input$colorstr,breaks=bins)+
      geom_vline(aes(xintercept = limit1),color = "red", linetype="dashed")+
      geom_vline(aes(xintercept = limit2),color = "red", linetype="dashed")+
      geom_vline(aes(xintercept = limit3),color = "red", linetype="dashed")+
      geom_vline(aes(xintercept = limit4),color = "red", linetype="dashed")+
      theme_bw()+xlab("normalized mass difference(Da)")+ylab("frequency")
    return(hist)
  })
  
  output$volcanotext<-renderUI({
    req(volcanoplot())
    tags$iframe(src='./volcano.html',
                width = '100%', height = '310px',
                frameborder = 0, scrolling = "no",align="center")
  })
  
  output$histogramclustplot=renderPlot({
    req(input$clusteringsel)
    if(input$clusteringsel==TRUE){
      print(histogram())}
  })
  
  output$volcanotitle <- renderUI({
    req(volcanoplot())
    h2("Volcano plots (Hybrid significance test)",align="center")
  })
  
  output$volcanotype <- renderUI({
    req(volcanoplot())
    selectInput("volcanotypes", 
                label = NULL,
                choices = c(".pdf", 
                            ".svg",
                            ".eps"
                ),
                selected = ".pdf")
  })
  
  output$volcanodownload <- renderUI({
    req(volcanoplot())
    downloadButton("volcanobutton", "Download Volcano plot")
  })
  
  output$histdownload<-renderUI({
    req(histogram())
    req(input$clusteringsel==TRUE)
    downloadButton("histbutton", "Download Histogram")
  })
  
  output$histbutton<-downloadHandler(
    filename=function(){
      paste0("Histogram",input$volcanotypes)
    },
    content=function(file){
      ggsave(file,(
        histogram()+ggtitle("Histogram of normalized differences")+theme(plot.title=element_text(hjust=0.5,size=20))
       )
      )
    }
  )
  
  output$volcanobutton<-downloadHandler(
    filename=function(){
      paste0("VolcanoPlot",input$volcanotypes)
    },
    content=function(file){
      ggsave(file,(
        volcanoplot()+ggtitle("Hybrid significance test")+theme(plot.title=element_text(hjust=0.5,size=20))
        )
      )
    })
  
  #----------------------------------## Second tab Global Woods plot###-------------------------------------------------
  woodsplots=reactive({ #calculation of global woods plot
    
    req(variables())
    listtestresults=variables()$listtestresults
    sequence=variables()$sequence
    
    listtestresults$SignificantResults=as.numeric(listtestresults$SignificantResults)
    listtestresults$Start=as.numeric(listtestresults$Start)
    listtestresults$End=as.numeric(listtestresults$End)
    
    colors=c("0"="gray43","-1"=input$colorprot,"1"=input$colorflex,"5"="yellow") #colors used for significant peptides
    allwoodsplot =  ggplot()+ geom_rect(data=listtestresults,color="black",linewidth=0.1,mapping=aes(xmin=pepnumber-0.40,xmax=pepnumber+0.40,ymin=Start,ymax=End,fill=factor(SignificantResults)))+
      scale_fill_manual(values=colors,guide="none")+
      scale_x_continuous(expand=c(0,2),breaks=seq(0,1000,by=20))+
      theme_bw()+
      labs(x="Peptide number",y="Protein residue")+
      theme(legend.position = "",panel.grid.major = element_blank())+
      scale_y_continuous(expand=c(0,2),breaks=seq(0,1000,by=25),limits=c(0,nchar(sequence)+5))
    
    return(allwoodsplot)
  })
  
  output$globwoodsplot=renderPlot({
    print(woodsplots())
  })
  
  output$globwoodstitle <- renderUI({
    req(variables())
    h2("Global woods plot",align="center")
  })
  
  output$globwoodstype <- renderUI({
    req(woodsplots())
    selectInput("woodstypes", 
                label = NULL,
                choices = c(".pdf", 
                            ".svg",
                            ".eps"
                ),
                selected = ".pdf")
  })
  
  output$woodsplottext<-renderUI({
    req(woodsplots())
    tags$iframe(src='./woodstext.html',
                width = '100%', height = '310px',
                frameborder = 0, scrolling = 'auto',align="center")
  })
  
  output$globwoodsdownload <- renderUI({
    req(woodsplots())
    downloadButton("woodsbutton", "Download")
  })
  
  output$woodsbutton<-downloadHandler(
    filename=function(){
      paste0("Global woods plots",input$woodstypes)
    },
    content=function(file){
      ggsave(file,(
        woodsplots()+ggtitle("Global woods plot")+theme(plot.title=element_text(hjust=0.5,size=20))
        )
      )
    })
  
  #------------------------------------- ## 4th tab Wood's plot by time-point ###-------------------------------------
  observe({
    req(input$hdexaminerfile)
    req(input$fastafile)
    req(input$timepoints)
    req(input$replicates)
    req(input$states)
    req(input$significancelevel)
    req(input$labelingtimepoints)
    
    updateSelectInput(session,"timepoint",label = "Select labeling time:",choices = labelingtimepoints())
    
  })
  
  woodsbytimepoint=reactive({ #wood's plots by time point
    req(variables())
    
    listtestresults=variables()$listtestresults
    HDXdatastates=variables()$HDXdatastates
    threshold=variables()$threshold
    
    plotsbytimepoint=list()
    
    listtestresults$Start=as.numeric(HDXdatastates[[1]]$Start)
    listtestresults$End=as.numeric(HDXdatastates[[1]]$End)
    
    maxvalue2=c()
    for (b in 1:timepoints()){
      maxvalue2=append(maxvalue2,listtestresults[,(b*2)+3])
      }
    maxvalue2=as.numeric(maxvalue2)
    
    for (a in 1:timepoints()){
      color1=c()
      plotsbytimepoint[[a]]=local({
        a=a
        for(i in 1:nrow(listtestresults)){
          if (is.na(listtestresults[i,(a*2)+3]) || is.na(listtestresults[i,(a*2)+4])){
            color1=append(color1,0)}else{
              if (((listtestresults[i,(a*2)+3])>=threshold) && (listtestresults[i,(a*2)+4]<significancelevel())){
                color1=append(color1,1)}else{
                  if (((listtestresults[i,(a*2)+3])<=(-threshold)) && (listtestresults[i,(a*2)+4]<significancelevel())){
                    color1=append(color1,-1)}else{
                      color1=append(color1,0)
                    }
                }
            }
        }
        color1=as.character(color1)
        
        ggplot()+ geom_rect(data=listtestresults,color="black",linewidth=0.15,aes(xmin=Start,xmax=End,ymin=listtestresults[,((2*a)+3)]-0.05,ymax=listtestresults[,((2*a)+3)]+0.05,fill=color1))+
          scale_x_continuous(expand=c(0,2),breaks=seq(-1000,1000,by=35))+
          scale_y_continuous(breaks=seq(-1000,1000,by=0.5))+theme_bw()+ ylim(-max(abs(maxvalue2),na.rm=TRUE)-0.05,max(abs(maxvalue2),na.rm=TRUE)+0.05)+
          labs(x="Residue",y=paste("delta HX (Da) ",labelingtimepoints()[a],"s",sep=""))+
          geom_hline(yintercept=threshold,color = "red",linewidth=1.0, linetype="dashed")+
          geom_hline(yintercept=-threshold,color = "red",linewidth=1.0, linetype="dashed")+
          theme(legend.position = "",panel.grid.major = element_blank())+
          scale_fill_manual(values=c("0"="gray43","-1"=input$colorprot,"1"=input$colorflex))
        
      })
    }
    
    return(plotsbytimepoint)
  })
  
  output$woodsbytimepointplot=renderPlot({  
    req(woodsbytimepoint())
    wanted=(input$timepoint)
    
    plots=woodsbytimepoint()
    for(i in 1:timepoints()){
      if(wanted==labelingtimepoints()[[i]]){number=i}
      }
    print(plots[[number]])
  })
  
  output$woodsbytimepointtype <- renderUI({
    req(woodsbytimepoint())
    selectInput("woodsbytimepointtypes", 
                label = NULL,
                choices = c(".pdf", 
                            ".svg",
                            ".eps"
                ),
                selected = ".pdf")
  })
  
  output$woodsplotbytimetext<-renderUI({
    req(woodsbytimepoint())
    tags$iframe(src='./woodsbytimetext.html',
                width = '100%', height = '310px',
                frameborder = 0, scrolling = 'auto',align="center")
  })
  
  output$woodsbytimepointtitle <- renderUI({
    req(woodsbytimepoint())
    h2("Woods-plot by timepoint",align="center")
  })
  
  output$woodsbytimepointdownload <- renderUI({
    req(woodsbytimepoint())
    downloadButton("woodsbytimepointbutton", "Download all")
  })
  
  output$woodsbytimepointbutton<-downloadHandler(
    filename=function(){
      paste0("Woodsplot by timepoint",input$woodsbytimepointtypes)
    },
    content=function(file){
      ggsave(file,marrangeGrob(woodsbytimepoint(), nrow=timepoints(), ncol=1,top="Woods-plot by timepoint"),width = 8.5, height = 11)
    })
  
  #----------------------------------## 5th tab digestion eff###-------------------------------------------------
  digestion=reactive({ #calculate digestion parameters 
    
    req(variables())
    
    HDXdatastates=variables()$HDXdatastates
    sequence=variables()$sequence
    
    coverage =as.character(sequence)
    for (i in 1:2){
      HDXdatastates[[i]]$Start=as.numeric(HDXdatastates[[i]]$Start)
      HDXdatastates[[i]]$End=as.numeric(HDXdatastates[[i]]$End)
       }
    for (i in 1:nrow(HDXdatastates[[1]])){
      substr(coverage,HDXdatastates[[1]][i,3],HDXdatastates[[1]][i,4])=strrep("Z",HDXdatastates[[1]][i,4]-HDXdatastates[[1]][i,3]+1)
       }
    
    digestion=data.frame()
    digestion[1,1]=nrow(HDXdatastates[[1]])
    digestion[1,2]=signif((str_count(coverage,"Z")/str_length(coverage))*100,4)
    digestion[1,3]=signif(mean(HDXdatastates[[1]]$End - HDXdatastates[[1]]$Start),4)
    
    redundancy=c()  
    redundancy=rep(0,str_length(sequence))
    for(i in 1:nrow(HDXdatastates[[1]])){
      for(j in HDXdatastates[[1]][i,3]:HDXdatastates[[1]][i,4]){
        redundancy[j]=redundancy[j]+1
      }
    }
    digestion[1,4]=signif(mean(redundancy),4)
    names(digestion)=c("Number of peptides","% coverage","Avg peptide length","Redundancy")
    
    avglength=ggplot(data=HDXdatastates[[1]],aes(x=End-Start))+geom_histogram(binwidth = 1,color="black",fill="gray")+
      theme_bw()+xlab("Peptide length")+ylab("frequency")+
      geom_vline(aes(xintercept = (mean(HDXdatastates[[1]]$End - HDXdatastates[[1]]$Start))),color = "red", linewidth=1.5,linetype="dashed")
    
    efftable = tableGrob(digestion, rows=NULL)
    
    return(list(plot=grid.arrange(avglength,efftable,nrow=2,ncol=1,top=""),digestion=digestion))
    
  })
  
  output$digestionplot=renderPlot({
    print(digestion()$plot)
  })
  
  output$digestiontitle <- renderUI({
    req(digestion())
    h2("Digestion performance",align="center")
  })
  
  output$digestiontype <- renderUI({
    req(digestion())
    selectInput("digestiontypes", 
                label = NULL,
                choices = c(".pdf", 
                            ".svg",
                            ".eps"
                ),
                selected = ".pdf")
  })
  
  output$digestiondownload <- renderUI({
    req(digestion())
    downloadButton("digestionbutton", "Download")
  })
  
  output$digestionbutton<-downloadHandler(
    filename=function(){
      paste0("Digestion",input$digestiontypes)
    },
    content=function(file){
      ggsave(file,digestion(),height=7,width=7)
    })
  
  output$digestiontext<-renderUI({
    req(digestion())
    tags$iframe(src='./digestiontext.html',
                width = '100%', height = '310px',
                frameborder = 0, scrolling = 'auto',align="center")
  })
  #----------------------------------## 6th tab peptide map###-------------------------------------------------
  
  peptidemap=reactive({ #calculation of peptide map
    
    req(variables())
    req(digestion())
    digestion=digestion()$digestion
    shinyjs::show("pep2")
    HDXdatastates=variables()$HDXdatastates
    sequence=variables()$sequence
    
    HDXdatastates[[1]]=HDXdatastates[[1]][,-1]
    pepmap=HDXdatastates[[1]][,1:3]
    pepmap[,4]=rep(1,nrow(pepmap))
    pepmap[,1]=as.numeric(pepmap[,1])
    pepmap[,2]=as.numeric(pepmap[,2])
    pepmap[,3]=as.numeric(pepmap[,3])               
    
    height=list()
    height[[1]]=rep(0,str_length(sequence))
    for(b in 1:nrow(pepmap)){
      c=1
      if(all(height[[c]][pepmap[b,2]:pepmap[b,3]]==0)){
        pepmap[b,4]=c
        height[[c]][pepmap[b,2]:pepmap[b,3]]=rep(1,pepmap[b,3]-pepmap[b,2]+1)
        }else{
          while (is.element(1,height[[c]][pepmap[b,2]:pepmap[b,3]])){
           c=c+1
           if((c-1)==length(height)){
            height[[c]]=rep(0,str_length(sequence))
             }
           if(all(height[[c]][pepmap[b,2]:pepmap[b,3]]==0)){
            pepmap[b,4]=c
            height[[c]][pepmap[b,2]:pepmap[b,3]]=rep(1,pepmap[b,3]-pepmap[b,2]+1)
            break
          }
        }
      }
    }
    
    coveragemap= ggplot(data=pepmap)+geom_rect(aes(xmin=Start,xmax=End,ymin=pepmap[,4]-0.5,ymax=pepmap[,4]+0.5),fill=input$colorpep,color="black",linewidth=0.07)+
      theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank(),axis.text=element_text(size=8),panel.background = element_rect(fill = 'white', color = 'white'),panel.grid.major = element_line(color = 'white'),panel.grid.minor = element_line(color = 'white'),axis.line.x.top = element_line(color="white"),axis.text.x.top = element_text(face="bold"))
    
    labels <- data.frame(
      x = c(
        seq(1, str_length(sequence), length.out = str_length(sequence)),
        c(1, (1:(str_length(sequence)/25))*25)
      ),
      label = c(
        format(s2c(sequence), length.out = str_length(sequence)),
        paste0("\n", c(1, (1:(str_length(sequence)/25))*25)))
    )
    
    plots=list()
    for(d in 1:(ceiling(str_length(sequence)/100))){
      plots[[d]]=coveragemap+
        scale_x_continuous(breaks = labels$x,labels = labels$label,minor_breaks = NULL,expand = c(0, 0), limits = c(0,NA))+
        coord_cartesian(xlim=c(((d-1)*100)+1,(d*100)))
    }
    
    coveragemap2= ggplot(data=pepmap)+geom_rect(aes(xmin=Start,xmax=End,ymin=0,ymax=1),fill=input$colorpep)+theme_void()+ theme(panel.background = element_rect(fill = 'gray'))+
      scale_x_continuous(expand = c(0, 0),limits=c(1,nchar(sequence))) + scale_y_continuous(expand = c(0, 0))+ggtitle(paste0(digestion[1,2],"% coverage with ", digestion[1,1]," peptides"))+ theme(plot.title = element_text(face = "bold"))
    
    
    coveragepep= ggplot(data=pepmap)+geom_rect(aes(xmin=Start,xmax=End,ymin=-pepmap[,4]-0.44,ymax=-pepmap[,4]+0.44),fill=input$colorpep,color="black",linewidth=0.03)+
      scale_x_continuous(expand=c(0,0),breaks=seq(0,nchar(sequence),by=50),limits=c(1,nchar(sequence))) + scale_y_continuous(expand = c(0, 0))+ xlab("residue number")+
      theme(axis.text.x = element_text(face="bold"),axis.line.x = element_line(color="black", linewidth = 1),axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
    
    summary=coveragemap2/coveragepep +plot_layout(heights=c(1,3))
    
   return(list(peptidemap= marrangeGrob(plots, nrow=ceiling(str_length(sequence)/100),ncol=1, top=NULL),summary=summary))
  })
  
  output$summaryplot=renderPlot({
    peptidemap()$summary}
    )
  
  output$summarypepbutton<-downloadHandler(
    filename=function(){
      paste0("Summary Peptide map",input$peptidetypes)
    },
    content=function(file){
      ggsave(file,peptidemap()$summary, height=3.5,width=14)
    })
  
  output$peptideplot=renderPlot({
    req(peptidemap())
    shinyjs::show("pep")
    print(peptidemap()$peptidemap)
  })
  
  output$peptidetitle <- renderUI({
    req(peptidemap())
    h2("Detailed peptide coverage map",align="center")
  })
  
  output$peptidetype <- renderUI({
    req(peptidemap())
    selectInput("peptidetypes", 
                label = NULL,
                choices = c(".pdf", 
                            ".svg",
                            ".eps"
                ),
                selected = ".pdf")
  })
  
  output$peptidedownload <- renderUI({
    req(peptidemap())
    downloadButton("peptidebutton", "Download")
  })
  
  output$peptidebutton<-downloadHandler(
    filename=function(){
      paste0("Peptide map",input$peptidetypes)
    },
    content=function(file){
      ggsave(file,peptidemap())
    })
  
  output$peptidemaptext<-renderUI({
    req(peptidemap())
    tags$iframe(src='./peptidemaptext.html',
                width = '100%', height = '95px',
                frameborder = 0, scrolling = 'no',align="center")
  })
  #----------------------------------## 7th tab uptake plots###-------------------------------------------------
  
  observe({
    req(variables())
    HDXdatastates=variables()$HDXdatastates
    updateSelectInput(session,"selpeptide",label = "Select peptide:",choices = paste0(HDXdatastates[[1]]$pepnum,": (", HDXdatastates[[1]]$Start,"-",HDXdatastates[[1]]$End,") ",HDXdatastates[[1]]$Sequence))
  })
  
  output$uptaketype <- renderUI({
    req(variables())
    selectInput("uptaketypes", 
                label = NULL,
                choices = c(".pdf", 
                            ".svg",
                            ".eps"
                ),
                selected = ".pdf")
  })
  
  singleplot = reactive({  # Uptake plot of selected peptide
    req(variables())
    shinyjs::show("plotssettings")
    shinyjs::show("limits")
    req(input$coloring, input$thickness, input$fontsize)
    
    if (states() > 2) {
      shinyjs::show("allplotscheck")
    }
    
    wanted = input$selpeptide
    selected = as.numeric(str_split(wanted, ":")[[1]][1])
    
    HDXdatastates = variables()$HDXdatastates
    listtestresults = variables()$listtestresults
    
    library(ggthemr)
    ggthemr(input$coloring)
    
    deut = c()
    sd = c()
    state = c()
    times = c()
    title = c()
    
    if (input$showallstates) {
      for (h in 1:states()) {
        deut = append(deut, as.vector(as.matrix(HDXdatastates[[h]] %>% filter(pepnumber == selected) %>% select(matches("Average")))))
        sd = append(sd, as.vector(as.matrix(HDXdatastates[[h]] %>% filter(pepnumber == selected) %>% select(matches("SD")))))
        state = append(state, rep(as.vector(as.matrix(HDXdatastates[[h]] %>% filter(pepnumber == selected) %>% select(matches("State")))), each = timepoints()))
        times = append(times, labelingtimepoints())
      }
    } else {
      for (h in 1:2) {
        deut = append(deut, as.vector(as.matrix(HDXdatastates[[h]] %>% filter(pepnumber == selected) %>% select(matches("Average")))))
        sd = append(sd, as.vector(as.matrix(HDXdatastates[[h]] %>% filter(pepnumber == selected) %>% select(matches("SD")))))
        state = append(state, rep(as.vector(as.matrix(HDXdatastates[[h]] %>% filter(pepnumber == selected) %>% select(matches("State")))), each = timepoints()))
        times = append(times, labelingtimepoints())
      }
    }
    
    uptake = data.frame(Prostate = state, deuteration = deut, stddev = sd, time = times)
    
    sig_result = listtestresults %>% filter(pepnumber == selected) %>% select('SignificantResults')
    title_suffix = ifelse(sig_result == 5, "**", ifelse(abs(sig_result) == 1, "*", ""))
    title = paste(HDXdatastates[[1]][selected, 3], "-", HDXdatastates[[1]][selected, 4], ":", HDXdatastates[[1]][selected, 5], title_suffix, sep = "")
    
    plot = ggplot(data = uptake, aes(x = time, y = deuteration, color = Prostate)) +
      scale_x_continuous(trans = 'log10') +
      scale_color_discrete(name = "Protein State") +
      xlab("Labeling time (sec)") +
      ylab("Average deuteration (Da)") +
      ggtitle(title) +
      theme_bw() +
      expand_limits(y = 0) +
      geom_line(linewidth = input$thickness) +
      geom_point(na.rm = TRUE, size = input$thickness * 1.5) +
      geom_errorbar(aes(ymin = deuteration - stddev, ymax = deuteration + stddev),
                    width = 0.05, linewidth = input$thickness * 0.5, color = "black") +
      theme(
        legend.title = element_text(size = input$fontsize + 2),
        axis.text = element_text(size = input$fontsize + 2),
        axis.title = element_text(size = input$fontsize + 2),
        legend.text = element_text(size = input$fontsize),
        plot.title = element_text(hjust = 0.5, size = input$fontsize + 8)
      )
    if(is.na(input$yminvalue)  | is.na(input$ymaxvalue)){
      return(plot)
    }else {
    return(plot+ylim(input$yminvalue,input$ymaxvalue))
    }
    
  })
  
  output$uptakeplot=renderPlot({  
    singleplot()
  })
  
  output$uptakedownload <- renderUI({
    req(variables())
    actionButton("uptakebutton", "Download all",icon=icon("download"))
  })
  
  observeEvent(input$selpeptide,{
    updateNumericInput(session, "yminvalue", value="",min=0,max=25)
    updateNumericInput(session, "ymaxvalue", value="",min=0,max=25)
  })
  
  output$singledownload<-renderUI({
    req(variables())
    downloadButton("singlebutton", "Download",icon=icon("download"))
  })
  
  x=3
  y=6
  
  observeEvent(input$uptakebutton,{
    showModal(modalDialog(title="Please select number of rows and columns to export your uptakeplots",easyClose=TRUE,
                          column(5,numericInput("columns", 
                                                "Select number of columns:", 
                                                value = 2,min=1,max=6)),
                          column(4,numericInput("rows", 
                                                "Select number of rows:", 
                                                value = 4,min=1,max=6)),br(),
                          footer = tagList (
                            column(2,modalButton("Cancel"),offset=6),
                            downloadButton("downloadall", "Download uptake plots")
                                          )
                        )
           )
    x=input$columns
    y=input$rows  
  })
  

  output$singlebutton<-downloadHandler(
    filename=function(){
      paste0("Uptakeplot",input$uptaketypes)
    },
    content=function(file){
      withProgress(message = 'Creating plot', min=0,max=1, { 
        ggsave(file,plot=singleplot(),width=8.5)
      })
    })
  
  allplots = reactive({  # Creation of all uptake plots
    req(variables())
    
    HDXdatastates = variables()$HDXdatastates
    listtestresults = variables()$listtestresults
    
    withProgress(message = 'Creating plots', min = 0, max = 1, { 
      
      library(ggthemr)
      ggthemr(input$coloring)
      
      uptakeplots = list()
      deut = c()
      sd = c()
      state = c()
      times = c()
      title = c()
      
      end = if (input$showallstates) states() else 2
      
      for (i in 1:nrow(HDXdatastates[[1]])) {
        for (h in 1:end) {
          deut = append(deut, as.vector(as.matrix(HDXdatastates[[h]] %>% filter(pepnumber == i) %>% select(matches("Average")))))
          sd = append(sd, as.vector(as.matrix(HDXdatastates[[h]] %>% filter(pepnumber == i) %>% select(matches("SD")))))
          state = append(state, rep(as.vector(as.matrix(HDXdatastates[[h]] %>% filter(pepnumber == i) %>% select(matches("State")))), each = timepoints()))
          times = append(times, labelingtimepoints())
        }
        
        uptake = data.frame(Prostate = state, deuteration = deut, stddev = sd, time = times)
        title = paste(HDXdatastates[[1]][i, 3], "-", HDXdatastates[[1]][i, 4], ":", HDXdatastates[[1]][i, 5], sep = "")
        
        uptakeplots[[i]] = local({
          i = i
          ggplot(data = uptake, aes(x = time, y = deuteration, color = Prostate)) +
            geom_point(na.rm = TRUE, size = 1.2) +
            geom_line(linewidth = 0.7) +
            geom_errorbar(aes(ymin = deuteration - stddev, ymax = deuteration + stddev),
                          width = 0.05, linewidth = 0.1, color = "black") +
            scale_x_continuous(trans = 'log10') +
            scale_color_discrete(name = "Protein State") +
            xlab("Labeling time (sec)") +
            ylab("Average deuteration (Da)") +
            ggtitle(title) +
            theme_bw() +
            expand_limits(y = 0) +
            theme(
              legend.title = element_text(size = 40 / (x + y)),
              axis.text = element_text(size = 64 / (x + y)),
              axis.title = element_text(size = 40 / (x + y / 2)),
              legend.text = element_text(size = 32 / (x + y)),
              plot.title = element_text(hjust = 0.5, size = 70 / (x + y)),
              legend.position = c(0.87, 0.26),
              legend.background = element_rect(fill = "white", color = "gray")
            )
        })
        
        deut = c()
        sd = c()
        state = c()
        times = c()
        title = c()
        incProgress(i / nrow(HDXdatastates[[1]]))
      }
      
      return(uptakeplots)
    })
  })
  
  
  output$downloadall<-downloadHandler(
    filename=function(){
      paste0("All uptakeplots.pdf")
    },
    content=function(file){
      withProgress(message = 'Exporting', min=0,max=1, { 
        on.exit(removeModal())
        ggsave(file,marrangeGrob(allplots(), nrow=input$rows, ncol=input$columns),width = 8.5, height = 11)
      })
    })
  
  output$uptaketext<-renderUI({
    req(variables())
    tags$iframe(src='./uptaketext.html',
                width = '100%', height = '310px',
                frameborder = 0, scrolling = 'auto',align="center")
  })
  
  #--------## 8th tab 3D Structure ###-------------------------------------
  output$pdbselector <-renderUI({
    req(variables())
    
    radioButtons("pdbfiletype","Please select how to import your PDB file:",
                 choices=c("Using PDB file","Using PDB number"),selected="Using PDB file")
  })
  
  output$pdbtext<-renderUI({
    req(input$pdbfiletype)
    req(variables())
    if(input$pdbfiletype=="Using PDB number"){
      textInput("pdbinputtext", "Type PDB code here:", value="")}else{NULL}
  })
  
  output$pdbinput<-renderUI({
    req(input$pdbfiletype)
    req(variables())
    if(input$pdbfiletype=="Using PDB file"){
      fileInput("pdbfile", "Select PDB file:",
                multiple = FALSE,
                accept = c(".pdb"))}else{NULL}
  })
  
  output$showactionbutton<-renderUI({
    req(variables())
    actionButton("showplot","Calculate")
  })
  
  pdbstructure=eventReactive(input$showplot,{ #visualize pdb file 
    
    if(input$pdbfiletype=="Using PDB file"){
      req(input$pdbfile)
      pdbstructure=input$pdbfile$datapath}else{
        if(input$pdbfiletype=="Using PDB number"){
          req(input$pdbinputtext)
          pdbstructure=input$pdbinputtext
        }else{pdbstructure=NULL}
      }
    return(pdbstructure)
  })
  
  significantresidues<-reactive({ #calculation of significant peptides/residues according to the hybrid significant test  
    req(variables())
    sequence=variables()$sequence
    
    listtestresults=variables()$listtestresults
    listtestresults$Start=as.numeric(listtestresults$Start)
    listtestresults$End=as.numeric(listtestresults$End)
    clusteringresults=variables()$clusteringresults
    
    if(input$clustercoloring==TRUE){
      sigresidues=c()
      sigresidues=as.numeric(rep(0,str_length(sequence)))
      for(i in 1:nrow(clusteringresults)){
        for(j in clusteringresults[i,2]:clusteringresults[i,3]){
          if(sigresidues[j]==0 ){
            sigresidues[j]=clusteringresults[i,(timepoints()*3)+5]}else{
              if(sigresidues[j]!=0 && (clusteringresults[i,(timepoints()*3)+5]>sigresidues[j])){
                sigresidues[j]=clusteringresults[i,(timepoints()*3)+5]
              }
            }
         }
      }
    }else{
      sigresidues=c()
      sigresidues=as.numeric(rep(0,str_length(sequence)))
    for(i in 1:nrow(listtestresults)){
      for(j in listtestresults[i,2]:listtestresults[i,3]){
        if(sigresidues[j]==0){
          sigresidues[j]=listtestresults[i,(timepoints()*2)+5]}else{
            if(sigresidues[j]==1 && (listtestresults[i,(timepoints()*2)+5]==-1|listtestresults[i,(timepoints()*2)+5]==5)){
              sigresidues[j]=5
            }else{
              if(sigresidues[j]==-1 && (listtestresults[i,(timepoints()*2)+5]==1|listtestresults[i,(timepoints()*2)+5]==5)){
                sigresidues[j]=5
              }
            }
          }
        
      }
    }
    }
     if (input$clustercoloring == FALSE) {
       return(list(
        all= paste((input$aaoffset + which(sigresidues == 0)), collapse = " or "),
        a = paste((input$aaoffset + which(sigresidues == 5)), collapse = " or "),
        b = paste((input$aaoffset + which(sigresidues == 1)), collapse = " or "),
        c = paste((input$aaoffset + which(sigresidues == -1)), collapse = " or ")
      ))
    } else {
      return(list(
        all= paste((input$aaoffset + which(sigresidues == 0)), collapse = " or "),
        a = paste((input$aaoffset + which(sigresidues == 1)), collapse = " or "),
        b = paste((input$aaoffset + which(sigresidues == 2)), collapse = " or "),
        c = paste((input$aaoffset + which(sigresidues == 3)), collapse = " or ")
      ))
    }
    
  })
  
  output$structure<-renderNGLVieweR({ #3d structure with coloring of significant peptides
    
    if(input$colorflex=="firebrick2"){flex="#EE2C2C"} else if(input$colorflex=="cyan"){flex="#00FFFF"}else if(input$colorflex=="deeppink"){flex="#FF1493"} else {flex="#EE2C2C"}
    if(input$colorprot=="blue"){prot="#0000FF"}else if(input$colorprot=="darkcyan"){prot="#008B8B"}else if(input$colorprot=="chocolate"){prot="#D2691E"} else {prot="#0000FF"}
   
    if(significantresidues()$a==""){a="none"}else{a=significantresidues()$a}
    if(significantresidues()$b==""){b="none"}else{b=significantresidues()$b}
    if(significantresidues()$c==""){c="none"}else{c=significantresidues()$c}
    
    if(input$clustercoloring==FALSE){NGLVieweR(pdbstructure(),width="50%") %>%
      stageParameters(backgroundColor ="white", zoomSpeed = 1) %>%
      addRepresentation(input$structuretype, param = list(sele=significantresidues()$all,color="gray"))%>%
      addRepresentation(input$structuretype,param =list(sele=a,color="yellow"))%>%
      addRepresentation(input$structuretype,param =list(sele=b,color=flex))%>%
      addRepresentation(input$structuretype,param =list(sele=c,color=prot))%>%
      setQuality("high")%>%setSpin(input$spinstructure)}else{
        NGLVieweR(pdbstructure(),width="50%") %>%
          stageParameters(backgroundColor ="white", zoomSpeed = 1) %>%
          addRepresentation(input$structuretype, param = list(sele=significantresidues()$all,color="gray"))%>%
          addRepresentation(input$structuretype,param=list(sele=a,color="gray"))%>%
          addRepresentation(input$structuretype,param=list(sele=b,color=input$colorint))%>%
          addRepresentation(input$structuretype,param=list(sele=c,color=input$colorstr))%>%
          setQuality("high")%>%setSpin(input$spinstructure)
      }
    
  })
  
  output$hrline<-renderUI({
    req(pdbstructure())
    hr(style = "border-top: 1px solid #000000;")
  })
  
  output$representation<-renderUI({
    req(pdbstructure())
    shinyjs::show("3dtoggles")
    selectInput("structuretype", 
                label = "Select Representation",
                choices = c("ball+stick", 
                            "surface",
                            "cartoon"
                ),
                selected = "cartoon")
  })
  
  output$spinbox<-renderUI({
    req(pdbstructure())
    materialSwitch("spinstructure", label = "Spin", value = TRUE, status = "primary")
  })
  
  observeEvent(input$structuretype,{
    NGLVieweR_proxy("structure", session=session)%>%removeSelection("allwaters")
  })
  
  output$structureoffset<-renderUI({
    req(pdbstructure())
    numericInput("aaoffset", 
                 label="Amino acid offset", 
                 value = 0,min=-50,max=50,step=1)
  })
  
  observeEvent(input$spinstructure,{
    NGLVieweR_proxy("structure", session=session)%>%updateSpin(input$spinstructure)
  })
  
  output$snapshot<-renderUI({
    req(pdbstructure())
    actionButton("screenshot", "Save screenshot")
  })
  
  observeEvent(input$screenshot,{
    NGLVieweR_proxy("structure",session=session)%>%setQuality("high")%>%
      snapShot("Structurescreenshot.PNG",param = list(
        antialias = TRUE,
        trim = TRUE,
        transparent = TRUE,
        scale = 1))
  })
  
  output$structuretitle<-renderUI({
    req(significantresidues())
    h2("3D Structure", align="left")
  })
  
  output$pymolbutton<-renderUI({
    req(pdbstructure())
    downloadButton("pymolscript", "Download pymol script")
  })
  
  observe({
    if(input$clustercoloring==TRUE){
    showModal(modalDialog(title="IMPORTANT!",
    "When coloring using clustering results, colors indicate if the peptides have an intermediate or strong effect. Colors DO NOT indicate protection or deprotection",footer=NULL,easyClose=TRUE))
    }
  })
  
  pymoloutput=reactive({ #calculate of pymol script with significant residues
    req(pdbstructure)
    req(significantresidues())
    if(input$clustercoloring==TRUE){
      colorstr=col2rgb(input$colorstr)
      colorstr=paste0("[",colorstr[1],",",colorstr[2],",",colorstr[3],"]")
      colorint=col2rgb(input$colorint)
      colorint=paste0("[",colorint[1],",",colorint[2],",",colorint[3],"]")
    }
    significantresidues=significantresidues()
    export=data.frame()
    mixed=gsub(" or ","+",significantresidues$a)
    deprotected=gsub(" or ","+",significantresidues$b)
    protected=gsub(" or ","+",significantresidues$c)
    export[1,1]="hide everything"
    export[2,1]="show cartoon"
    export[3,1]="color gray,all"
    if(input$clustercoloring==TRUE){
      export[4,1]=paste0("set_color strong, ",colorstr)
      export[5,1]=paste0("select strongeffects, (i. ",protected,")")
      export[6,1]="color strong, strongeffects" 
    }else{
    export[4,1]=paste0("select mixedpep, (i. ",mixed,")")
    export[5,1]="color yellow, mixedpep"}
    if(input$clustercoloring==TRUE){
      export[7,1]=paste0("set_color inter, ",colorint)
      export[8,1]=paste0("select intermediateeffects, (i. ",deprotected,")")
      export[9,1]="color inter, intermediateeffects"
    }else{
    if(input$colorflex=="firebrick2"){flexcode="[238,44,44]"} else if(input$colorflex=="cyan"){flexcode="[0,255,255]"}else if(input$colorflex=="deeppink"){flexcode="[255,20,147]"}
    export[6,1]=paste0("set_color deprotected, ",flexcode)
    export[7,1]=paste0("select deprotectedpep, (i. ",deprotected,")")
    export[8,1]="color deprotected, deprotectedpep"}
    if(input$clustercoloring==TRUE){
      export[10,1]=paste0("select negligibleeffects, (i. ",mixed,")")
      export[11,1]="color gray, negligibleeffects"
      export[12,1]="deselect"
    }else{
    if(input$colorprot=="blue"){protcode="[0,0,255]"}else if(input$colorprot=="darkcyan"){protcode="[0,139,139]"}else if(input$colorprot=="chocolate"){protcode="[210,105,30]"}
    export[9,1]=paste0("set_color protected, ",protcode)
    export[10,1]=paste0("select protectedpep, (i. ",protected,")")
    export[11,1]="color protected, protectedpep"
    export[12,1]="deselect"}
    return(export)
  })
  
  output$pymolscript<-downloadHandler(
    filename=function(){
      paste0("Pymol script.pml")
    },
    content=function(file){
      write.table(pymoloutput(),file, na='',row.names=FALSE,col.names=FALSE,quote=FALSE)
    })
  #------------## 9th tab Export results ###------------------------------------
  output$datadownload <- renderUI({
    req(variables())
    downloadButton("exportbutton", "Export results")
  })
  
  output$exportheader<-renderUI({
    req(variables())
    h2("Export Results")
  })
  
  output$exporttext<-renderUI({
    req(variables())
    p("Click the buttons below to download the results and all the generated plots")
  })
  
  exportdata=reactive({ #export excel spreadsheet with all data/results
    req(variables())
    req(realfit())
    withProgress(message = 'Exporting data', min=0,max=1, {
      library(openxlsx, warn=FALSE)
      
      pepcharge=initialdata()$pepcharge
      HDXdatastates=variables()$HDXdatastates
      clusteringresults=variables()$clusteringresults
      listtestresults=variables()$listtestresults
      parameters=variables()$parameters
      parameters[3,1]=paste0("Positive deltaD values (State2-State ref) means increase in flexibility by state 2 compared to the reference (Peptides colored ",input$colorflex,")")
      parameters[4,1]=paste0("Negative deltaD values (State2-State ref) means increase in protection by state 2 (Peptides colored ",input$colorprot,")")
      threshold=variables()$threshold
      incProgress(1/4)
      for (i in 1:2){
        HDXdatastates[[i]]=cbind(HDXdatastates[[i]][,1:5], pepcharge, HDXdatastates[[i]][,6:((timepoints()*replicates())+(timepoints()*2)+5)])
      }
      
      wb <<- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb,paste("State 1", HDXdatastates[[1]][1,1],"(Ref)"))
      openxlsx::writeData(wb,paste("State 1",HDXdatastates[[1]][1,1],"(Ref)"),HDXdatastates[[1]])
      openxlsx::addWorksheet(wb,paste("State 2",as.character(HDXdatastates[[2]][1,1])))
      openxlsx::writeData(wb,paste("State 2",as.character(HDXdatastates[[2]][1,1])),HDXdatastates[[2]])
      
      openxlsx::addWorksheet(wb,"Statistical parameters")
      openxlsx::writeData(wb,"Statistical parameters",parameters)
      openxlsx::addWorksheet(wb,"Statistics 2 vs ref")
      listtestresults$SignificantResults[listtestresults$SignificantResults==5]="check"
      openxlsx::writeData(wb,"Statistics 2 vs ref",listtestresults)
      incProgress(2/4)
      openxlsx::addWorksheet(wb,paste('Clustering', 2,'vs ref'))
      openxlsx::writeData(wb,paste('Clustering', 2, 'vs ref'),clusteringresults)
      
      style1 = openxlsx::createStyle(fontColour = input$colorflex)
      style2 = openxlsx::createStyle(fontColour = input$colorprot)
      style3 = openxlsx::createStyle(fontColour = "yellow3")
      style4 = openxlsx::createStyle(fontColour=input$colorprot,textDecoration = "bold")
      style5 = openxlsx::createStyle(fontColour=input$colorflex,textDecoration = "bold")
      
      if(realfit()==2){
        peptide=c()
        for(i in 1:length(fitting()$kinetics)){
          peptide=append(peptide,fitting()$kinetics[[i]]$data$rowname[1])
        }
        
        openxlsx::addWorksheet(wb, "Fitting parameters")
        openxlsx::writeData(wb,"Fitting parameters",cbind(fitting()$fittingparameters,peptide))
        table=data.frame()
        summary=fitting()$summary
        
        for(k in 1:length(fitting()$kinetics)){
          null=data.frame(coef(summary[[k]]$null),check.names=F)
          null$variable <- substr(rownames(null),1,1)
          null= null %>% select(variable, everything())
          null[,6]="Null"
          colnames(null)[6]="State"
        
          for (l in 1:states()){
            df=data.frame(coef(summary[[k]][[l+1]]),check.names=F)
            df$variable <- substr(rownames(df),1,1)
            df= df %>% select(variable, everything())
            df[,6]=paste0("State",l)
            colnames(df)[6]="State"
            null<-rbind.fill(null,df)
           }
          
          null=null %>% unite('variable',c(variable,State),remove=TRUE)
          null=pivot_wider(null, names_from = variable,values_from = c(Estimate, 'Std. Error','t value','Pr(>|t|)'))
          
          table<-rbind.fill(table,null)
        }
        
        openxlsx::addWorksheet(wb, "Calculated parameters")
        openxlsx::writeData(wb,"Calculated parameters",cbind(peptide,table))
        
      }else if (realfit()==1){}
      
      openxlsx::addStyle(wb, sheet="Statistical parameters", style=style1, rows=4, cols=1, gridExpand=TRUE)
      openxlsx::addStyle(wb, sheet="Statistical parameters", style=style2, rows=5, cols=1, gridExpand=TRUE)
      openxlsx::addStyle(wb, sheet="Statistical parameters", style=style3, rows=6, cols=1, gridExpand=TRUE)
      
      incProgress(3/4)
      for(j in 1:nrow(listtestresults)){
        if (listtestresults$SignificantResults[j]=="-1"){
          openxlsx::addStyle(wb, sheet=paste('Statistics', 2,'vs ref'), style=style2, rows=j+1, cols=1:ncol(listtestresults), gridExpand=TRUE)
        } else{
          if (listtestresults$SignificantResults[j]=="1"){
            openxlsx::addStyle(wb, sheet=paste('Statistics', 2,'vs ref'), style=style1, rows=j+1, cols=1:ncol(listtestresults), gridExpand=TRUE)
          } else{
            if (listtestresults$SignificantResults[j]=="check"){
              openxlsx::addStyle(wb, sheet=paste('Statistics', 2,'vs ref'), style=style3, rows=j+1, cols=1:ncol(listtestresults), gridExpand=TRUE)
            }else{}
          }
        }
      }
      
      for(i in 1:nrow(listtestresults)){
        for(j in 1:timepoints()){
          if (is.na(listtestresults[i,(2*j)+3]) || is.na(listtestresults[i,(2*j)+4])){
          }else{
            if (((listtestresults[i,(2*j)+3])>=threshold) && (listtestresults[i,(2*j)+4]<significancelevel())){
              openxlsx::addStyle(wb, sheet=paste('Statistics', 2,'vs ref'), style=style5, rows=i+1, cols=(((2*j)+3):((2*j)+4)))}else{
                if (((listtestresults[i,(2*j)+3])<=-threshold) && (listtestresults[i,(2*j)+4]<significancelevel())){
                  openxlsx::addStyle(wb, sheet=paste('Statistics', 2,'vs ref'), style=style4, rows=i+1, cols=(((2*j)+3):((2*j)+4)))}else{
                    
                  }
              }
          }  
        }
      }
      incProgress(4/4)
    return(wb)
    })
  })
  
  output$tableresults<-renderUI({
    req(variables())
    tableOutput("outputresults")
  })
  
  output$outputresults<-renderTable({
    req(variables())
    listtestresults=variables()$listtestresults
    listtestresults$SignificantResults[listtestresults$SignificantResults==5]="check"
    listtestresults
  }, digits = 3)
  
  

  output$exportbutton<-downloadHandler(
    filename=function(){
      paste0("HDX Analysis results.xlsx")
    },
    content=function(file){
      openxlsx::saveWorkbook(exportdata(),file)
    })
  
  output$allplotsdownload <-renderUI({
    req(variables())
    actionButton("exportallaction","Export all plots", icon=icon("download"))
  })
  
  document2=reactive({
    volcano=volcanoplot()
    if(input$clusteringsel==TRUE){
      histogram=histogram()
    }
    woodsplots=woodsplots()
    woods=woodsbytimepoint()
    digestion=digestion()
    peptide=peptidemap()
    return(list(volcano=volcano,histogram=histogram,woodsplots=woodsplots,digestion=digestion,peptide=peptide,woods=woods))
  })
  
  observeEvent(input$exportallaction,{
    showModal(modalDialog(title="Please wait. Creating your document...",easyClose=FALSE,
                          modalButton("Cancel"),
                          downloadButton("exportallplotsbutton", "Export all plots"),
                          footer=NULL
                 )
             )
    shinyjs::disable(id="exportallplotsbutton")
    if(realfit()==2){
    a=kineticdocument()}
    b=allplots()
    c=document2()
    shinyjs::enable(id="exportallplotsbutton")
  })
  
  output$exportallplotsbutton<-downloadHandler( #export of all plots in pdf format
    filename="All plots.pdf",
    content=function(file){
      withProgress(message = 'Exporting', min=0,max=1, { 
        pdf(paste0(file, ".7x7"), width=7, height=7)
        print(document2()$volcano)
        incProgress(1/7)
        if(input$clusteringsel==TRUE){print(document2()$histogram)}
        print(document2()$woodsplots)
        print(marrangeGrob(document2()$woods, nrow=timepoints(), ncol=1,top="Woods-plot by timepoint"))
        incProgress(2/7)
        print(document2()$peptide)
        incProgress(3/7)
        print(digestion())
        dev.off()
        pdf(paste0(file, ".8x11"), width=8, height=11)
        incProgress(4/7)
        print(marrangeGrob(allplots(), nrow=4, ncol=2, top=NULL))
        dev.off()
        if(realfit()==2){
        pdf(paste0(file, ".14x8"), width=14, height=8)
        print(kineticdocument())
        dev.off()
        }else{}
        incProgress(5/7)
        qpdf::pdf_subset(paste0(file, ".8x11"), pages=2:(qpdf::pdf_length(paste0(file, ".8x11"))),output=paste0(file, ".8"))
        if(realfit()==1){
          qpdf::pdf_combine(paste0(file, c(".7x7", ".8")), output=file)
          file.remove(paste0(file, c(".7x7", ".8")))
        }else{
            qpdf::pdf_combine(paste0(file, c(".7x7", ".8",".14x8")), output=file)
            file.remove(paste0(file, c(".7x7", ".8",".14x8")))
          }
        incProgress(6/7)
      
      })
     }
    )
}

