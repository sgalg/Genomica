#' Genomica - analysis with only one predictor
#'
#' This function will perform a linear mixed model on all the features in your data set, whilst also taking into account of your experimental design (i.e., fixed and random effects). Please see the vignette for more information on how to prepare your data frames for the analysis.
#' @param Data A data frame with samples as columns and features as rows.
#' @param Metadata A data frame with samples as rows.
#' @param Predictor example: c('Treatment') - This is the main predictor (independent variable) according to your study design.
#' @param P_levels example: c('Control','Treatment1','Treatment2) - If your main predictor is not a numerical variable, this will set the levels for your main predictor (the model will test the comparisons towards the first term in this vector, such as Control Vs Treatment1 and Control Vs Treatment2 in the example).
#' @param R_Effects example: c('Room') - This vector defines the random effects according to your study design.
#' @param R1_Levels example: c('Room1,'Room2,Room3,Room4,'Room5,Room6,Room7) - This vector will specify the levels for the random effect, according to your study design.
#' @param Already_Log10_transformed Default = c('NO'). If c('YES'), the log10 (n+1) transformation of the features in the object Data will not be performed.
#' @param Folder_Name Example: c('Test') - This string will be attached to the folder name in which the results will be generated (e.g.,"Genomica_Output_Test/")
#' @import readxl
#' @import writexl
#' @import lmerTest
#' @import dplyr
#' @import fuzzySim
#' @export
#' @examples
#' #Load Data data frame (Data_Demo)
#' Data<-Data_Demo
#' #In this example Data is a data frame of 3 taxonomical features (rows) across 42 samples (columns)
#' print(Data[,1:5])
#' Metadata<-Metadata_Demo
#' #in this example Metadata is a data frame with 42 rows (sample IDs) and 5 columns (Sample ID, Treatment, Sex, Room and Block). Treatment contains the treatment layout (Control, T1 or T2) and will be the fixed effect of the model, whilst room describes the stratification of this experimental layout according to the random effect "Room".
#' print(Metadata[1:5,])
#' genomica_one_predictor(Data = Data,Metadata = Metadata,
#' Predictor = c('Treatment'),P_Levels=c('Control','T1','T2'),
#' R_Effects=c('Room'),R1_Levels=c('Room1','Room2','Room3','Room4','Room5','Room6','Room7'),
#' Already_Log10_transformed = c('NO'),Folder_Name=c('Test'))


genomica_one_predictor<-function(Data,Metadata,Predictor,
                                 P_Levels=c(0),R_Effects,
                                 R1_Levels,Already_Log10_transformed=c('NO'),
                                 Folder_Name){


  FOLDER<-paste0('Genomica_Output_',Folder_Name)
  FILE_COMBINED<-'Combined_All_Features'
  PREVIOUS<-list.files(path = paste0(FOLDER,'/'),pattern='*.xlsx')

  if (length(PREVIOUS)>0){
    for (i in PREVIOUS[1:length(PREVIOUS)]){
      unlink(paste0(FOLDER,'/',i))
    }}


  ifelse(!dir.exists(FOLDER), dir.create(FOLDER), "A folder with the same name already exists")
  Path<-paste0(FOLDER,'/',FILE_COMBINED,'.xlsx')

  Rosetta_Sone<-as.data.frame(cbind(rownames(Data),seq(1,nrow(Data),1)))
  colnames(Rosetta_Sone)<-c('Feature','Feature_Code')
  rownames(Data)<-Rosetta_Sone$Feature_Code
  Data<-as.data.frame(t(Data))
  SampleID<-rownames(Data)
  Data<-as.data.frame(apply(Data[,1:length(Data)], 2, function(x) as.numeric(as.character(x))))
  if (Already_Log10_transformed=='NO'){
  Data<-log10(Data+1)
  }
  rownames(Data)<-SampleID
  Data$ID<-rownames(Data)
  Metadata$ID<-rownames(Metadata)
  DF<-Metadata
  DF<-dplyr::left_join(Metadata,Data)
  Random1<-R_Effects[1]
  Random2<-R_Effects[2]

  if(length(R_Effects)==1){
  DF[,Random1]<-factor((DF[,Random1]),levels=c(R1_Levels))
  Random1<-DF[,Random1]
  } else if (length(R_Effects)==2){
    DF[,Random1]<-factor((DF[,Random1]),levels=c(R1_Levels))
    DF[,Random2]<-factor((DF[,Random2]),levels=c(R2_Levels))
    Random1<-DF[,Random1]
    Random2<-DF[,Random2]
  }

  # Random2<-DF[,Random2]

  if (length(P_Levels)==1){
    DF[,Predictor]<-as.numeric(as.character(DF[,Predictor]))
  } else if(length(P_Levels)>1){
    DF[,Predictor]<-factor(DF[,Predictor],levels = P_Levels)
  }

  y<-colnames(DF)
  y<-y[(length(Metadata)+1):length(y)]



  for (i in y[1:length(y)]){

    Results<-matrix(NA,ncol=6,1,dimnames = list(c(i),
                                                c("Sum Sq","Mean Sq","NumDF","DenDF","F value","Pr(>F)")))




    if ((sum(DF[,i]))>0){
      model<-lmerTest::lmer(DF[,i]~DF[,Predictor]+(1|Random1),data=DF)
      res<-stats::anova(model)
      Results[i,]<-as.numeric(c(res[,1:length(res)]))
      Results<-as.data.frame(Results)
      Results$Feature_Code<-rownames(Results)
      Results<-Results[,c(7,1:6)]
      filepath<-paste(FOLDER,'/',i,'.xlsx',sep='')
      writexl::write_xlsx(Results,path =filepath)
  }}

  # Combined_Results<-matrix(NA,ncol=7,length(y),dimnames = list(c(as.vector(as.character(y))),
  #                                                              c('Feature_Code',"Sum Sq","Mean Sq","NumDF","DenDF","F value","Pr(>F)")))
  # Combined_Results[,1]<-as.vector(as.character(y))
  #
  # for (i in y[1:length(y)]){
  #   filename<-paste(FOLDER,'/',i,'.xlsx',sep='')
  #   Combined_Results[i,2:7]<-as.numeric(as.vector(readxl::read_xlsx(filename)))
  #   Combined_Results<-as.data.frame(Combined_Results)
  #   writexl::write_xlsx(Combined_Results,path = Path)
  #   unlink(filename)
  # }

  file.list1 <- list.files(path = paste0(FOLDER,'/'),pattern='*.xlsx')
  df.list1 <- lapply(paste0(FOLDER,'/',file.list1), readxl::read_xlsx)
  for (i in y[1:length(y)]){
    filen<-paste(FOLDER,'/',i,'.xlsx',sep='')
    unlink(filen)}

  Combined_Results<-do.call(rbind.data.frame, df.list1)

  Combined_Results$`Pr(>F)`<-as.numeric(as.character(Combined_Results$`Pr(>F)`))
  p.values<-data.frame(var=Combined_Results$Feature_Code,pval = Combined_Results$`Pr(>F)`)
  p.adj<-fuzzySim::FDR(pvalues = p.values)

  FDR_DF<-p.adj$exclude
  FDR_DF<-rbind(p.adj$exclude,p.adj$select)
  FDR_DF$Feature_Code<-rownames(FDR_DF)
  FDR_DF<-FDR_DF[,-1]
  Combined_Results<-dplyr::left_join(Combined_Results,FDR_DF)
  Col<-seq(2,8,1)
  Combined_Results[, Col] <- apply(Combined_Results[, Col], 2, function(x) as.numeric(as.character(x)))
  Combined_Results$P_Significant<-ifelse(Combined_Results$`Pr(>F)`<0.05,'Significant','No')
  Combined_Results$Q_Significant<-ifelse(Combined_Results$p.adjusted<0.05,'Significant','No')
  Combined_Results2<-dplyr::left_join(Rosetta_Sone,Combined_Results,by = 'Feature_Code')
  Combined_Results2<-Combined_Results2[,-2]
  writexl::write_xlsx(Combined_Results2,path = Path)

  Significant<-subset(Combined_Results,Q_Significant=='Significant')

  x<-Significant$Feature_Code
  rows<-c(paste0('(Intercept)',as.character(levels(DF[,Predictor]))))

  if (as.numeric(length(x))<1){
    return(print(paste0('No significant correlations found, the complete list of results can be found in ',FOLDER)))
    } else{
  for (i in x[1:length(x)]){

    model2<-lmerTest::lmer(DF[,i]~DF[,Predictor]+(1|Random1),data=DF)
    sum<-summary(model2)
    summ<-as.data.frame(sum$coefficients)
    summ$Comparison<-rownames(summ)
    summ$compared_to<-c(rep(i,nrow(summ)))
    Results<-as.data.frame(cbind(summ$compared_to,summ$Comparison,summ[,1:(length(summ)-2)]))
    filepath<-paste(FOLDER,'/',i,'_Significant.xlsx',sep='')
    writexl::write_xlsx(Results,path =filepath)
  }
  }

  file.list <- list.files(path = paste0(FOLDER,'/'),pattern='*_Significant.xlsx')
  df.list <- lapply(paste0(FOLDER,'/',file.list), readxl::read_xlsx)

  for (i in x[1:length(x)]){
    filename2<-paste(FOLDER,'/',i,'_Significant.xlsx',sep='')
    unlink(filename2)}

  Combined_Significant<-do.call(rbind.data.frame, df.list)
  Combined_Significant<-subset(Combined_Significant,Combined_Significant$`summ$Comparison`!='(Intercept)')
  Combined_Significant$`Pr(>|t|)`<-as.numeric(as.character(Combined_Significant$`Pr(>|t|)`))
  Combined_Significant$Significant<-ifelse(Combined_Significant$`Pr(>|t|)`<0.05,
                                           'Significant','No')
  Combined_Significant<-subset(Combined_Significant,Significant=='Significant')
  Combined_Significant<-dplyr::rename(Combined_Significant,'Feature_Code'=`summ$compared_to`)
  Combined_Significant<-dplyr::rename(Combined_Significant,'Predictor'=`summ$Comparison`)
  Combined_Significant<-Combined_Significant[-length(Combined_Significant)]

  Combined_Significant$Predictor<-sub(pattern = 'DF',replacement = '',Combined_Significant$Predictor)
  Combined_Significant$Predictor<-sub(pattern = 'Predictor',replacement = '',Combined_Significant$Predictor)
  Combined_Significant$Predictor<-sub(pattern = '[, ]',replacement = '',Combined_Significant$Predictor)
  Combined_Significant$Predictor<-sub(pattern = '[^[:alnum: ]]',replacement = '',Combined_Significant$Predictor)
  Combined_Significant$Predictor<-gsub(pattern = '[[]',replacement = '`',Combined_Significant$Predictor)
  Combined_Significant$Predictor<-sub(pattern = '`',replacement = '',Combined_Significant$Predictor)
  Combined_Significant$Status<-ifelse(Combined_Significant$Estimate<0,'Depleted','Enriched')
  Combined_Significant<-dplyr::left_join(Combined_Significant,Rosetta_Sone,by='Feature_Code')
  Combined_Significant<-Combined_Significant[,-1]
  Combined_Significant<-Combined_Significant[,c(8,1:7)]
  writexl::write_xlsx(Combined_Significant,
             path = paste0(FOLDER,'/Significant_Comparisons.xlsx'))
}
