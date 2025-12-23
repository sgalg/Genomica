#' Genomica - LMM driven (FDR adjusted) KEGG Orthologs enrichment analysis
#'
#' This function will perform a linear mixed model on all the features in your data set, whilst also taking into account of your experimental design (i.e., fixed and random effects). Thus, Genomica will carry out an enrichment analysis of the orthologs based on the given predictor list. Please see the vignette for more information on how to prepare your data frames for the analysis.
#' @param Data A data frame with samples as columns and features as rows.
#' @param Metadata A data frame with samples as rows.
#' @param Predictors example: c('Treatment','Growing_Phase') - This vector will be used to specify the predictors independent variables, according to your study design. Currently, Genomica accepts maximum two predictors, in which case it will automatically calculate the interaction between the two according to the formula predictor 1 * predictor 2 within the linear mixed model.
#' @param P1_levels example: c('Control','Treatment1','Treatment2) - If your main predictor is not a numerical variable, this will set the levels for your main predictor (the model will test the comparisons towards the first term in this vector, such as Treatment1 Vs Control and Treatment2 Vs Control in the example).
#' @param P2_levels example: c('Starter','Grower') - If your second predictor is not a numerical variable, this will set the levels for your main predictor (the model will test the comparisons towards the first term in this vector, such as Grower Vs Starter in the example).
#' @param R_Effects example: c('Block') - This vector defines the random effects according to your study design.
#' @param R1_Levels example: c('Block1,'Block2,Block3,Block4,'Block5,Block6,Block7) - This vector will specify the levels for the random effect, according to your study design.
#' @param Already_Log10_transformed Default = c('NO'). If c('YES'), the log10 (n+1) transformation of the features in the object Data will not be performed.
#' @param Folder_Name Example: c('Test') - This string will be attached to the folder name in which the results will be generated (e.g.,"Genomica_Output_Test/")
#' @param FDR_Level Default: 0.05 - Target False Discovery Rate (FDR). Benjamini-Hochberg multiple testing correction. Results with BH-adjusted P Values ≤ FDR_Level are considered significant.
#' @param ImgRes Default: 300 - Dots Per Inch (DPI) resolution of the .tiff output files.
#' @import readxl
#' @import writexl
#' @import lmerTest
#' @import dplyr
#' @import fuzzySim
#' @import MicrobiomeProfiler
#' @import enrichplot
#' @import clusterProfiler
#' @import tictoc
#' @import tidyr
#' @export
#' @examples
#' #Load Data data frame (Data_Demo)
#' Data<-Data_Demo
#' #In this example Data is a data frame of 500 KEGG orthologs (KOs) in rows, across 35 samples (columns)
#' print(Data[1:5,1:5])
#' #Make sure that Data data frame is formatted with the feature names as rows, and that it only contains the feature abundance as columns (i.e., delete eventual column characters):
#' rownames(Data)<-Data$KO
#' Data<-Data[,-1]
#' Metadata<-Metadata_Demo
#' #in this example Metadata is a data frame with 35 rows (sample IDs) and 5 columns (SampleID, Treatment and Block). Treatment contains the treatment layout (1,2,3,4,5) and will be the fixed effect of the model, whilst block describes the stratification of this experimental layout according to the random effect "Room".
#' print(Metadata[1:5,])
#' #Make sure that Metadata data frame is formatted with the sample names as rows:
#' rownames(Metadata)<-Metadata$SampleID
#' genomica(Data = Data,Metadata = Metadata,
#' Predictors = c('Treatment'),P1_Levels=c('1','2','3','4','5'),
#' R_Effects=c('Block'),R1_Levels=c('1','2','3','5','6','7'),
#' Already_Log10_transformed = c('NO'),Folder_Name=c('Test'))


genomica<-function(Data,Metadata,
                   Predictors,P1_Levels=c(0),P2_Levels=c(0),
                   R_Effects,R1_Levels,
                   Already_Log10_transformed=c('NO'),Folder_Name,
                   FDR_Level=0.05,ImgRes=300){


  if (length(Predictors)>2){
    print('Genomica currently works with maximum two predictors, please amend the `Predictors` vector in order to include maximum 2 predictors')
  }

Check_Predictors<-setdiff(Predictors,colnames(Metadata))

if (length(Check_Predictors)>0){
  stop(
    paste0('Genomica has detected an error in the data frame setup. The Predictor(s) "', Check_Predictors, '" could not be found in Metadata. Please ensure that the Metadata column names are identical to the Predictors specified here.')
  )
}


if (length(P1_Levels)==1){
  if(!is.numeric(P1_Levels)){
    stop(
      paste0('Genomica has detected an error in the data frame setup. ', Predictors[1],' appears to be a numerical predictor, thus please specify P1_Levels = c(0). OR please set the appropriate number of levels (e.g., P1_Levels = c("1","2")), if ',
             Predictors[1], ' is a categorical predictor.')
    )
  }
}

if (length(P1_Levels)>1){
if (!is.character(P1_Levels)){
  stop(
    paste0('Genomica has detected an error in the data frame setup. P1_Levels for ', Predictors[1],' must be a character vector. Please transform to character before running Genomica, OR please specify P1_Levels = c(0) if ',
           Predictors[1], ' is a numerical predictor.')
  )
}

  Predictor1_Duplicate<-P1_Levels[duplicated(P1_Levels)]
  if (length(Predictor1_Duplicate)>0){
    stop(
      paste0('WARNING: Genomica found that one or more levels in P1_Levels is duplicated: ',
             Predictor1_Duplicate,' can be found twice in P1_Level - Please rename levels to avoid duplicates.')
    )
  }

Predictor1_L<-unique(as.character(Metadata[,Predictors[1]]))
Predictor1_Missing<-setdiff(P1_Levels,Predictor1_L)
Predictor1_Extra<-setdiff(Predictor1_L,P1_Levels)
if (length(Predictor1_Missing)>0){
  stop(
    paste0('Genomica has detected an error in the data frame setup. The levels specified for Predictor 1 are not the same as the ones in Metadata. Level(s) "',
           Predictor1_Missing,
           '": not found in Metadata.')

  )
}
if (length(Predictor1_Extra)>0){
  stop(
    paste0('Genomica has detected an error in the data frame setup. The levels specified for Predictor 1 are not the same as the ones in Metadata. Level(s) "',
           Predictor1_Extra,
           '": not found in P1_Levels')
  )
}
}

if (length(Predictors)==2){
  if (length(P2_Levels)==1){
    if(!is.numeric(P2_Levels)){
      stop(
        paste0('Genomica has detected an error in the data frame setup. ', Predictors[2],' appears to be a numerical predictor, thus please specify P2_Levels = c(0). OR please set the appropriate number of levels (e.g., P2_Levels = c("1","2")), if ',
               Predictors[2], ' is a categorical predictor.')
      )
    }
  }

  if (length(P2_Levels)>1){
    if (!is.character(P2_Levels)){
      stop(
        paste0('Genomica has detected an error in the data frame setup. P2_Levels for ', Predictors[2],' must be a character vector. Please transform to character before running Genomica, OR please specify P2_Levels = c(0) if ',
               Predictors[2], ' is a numerical predictor.')
      )
    }

    Predictor2_Duplicate<-P2_Levels[duplicated(P2_Levels)]
    if (length(Predictor2_Duplicate)>0){
      stop(
        paste0('WARNING: Genomica found that one or more levels in P2_Levels is duplicated: ',
               Predictor2_Duplicate,' can be found twice in P2_Level - Please rename levels to avoid duplicates.')
      )
    }

    Predictor2_L<-unique(as.character(Metadata[,Predictors[2]]))
    Predictor2_Missing<-setdiff(P2_Levels,Predictor2_L)
    Predictor2_Extra<-setdiff(Predictor2_L,P2_Levels)
    if (length(Predictor2_Missing)>0){
      stop(
        paste0('Genomica has detected an error in the data frame setup. The levels specified for Predictor 2 are not the same as the ones in Metadata. Level(s) "',
               Predictor2_Missing,
               '": not found in Metadata.')

      )
    }
    if (length(Predictor2_Extra)>0){
      stop(
        paste0('Genomica has detected an error in the data frame setup. The levels specified for Predictor 2 are not the same as the ones in Metadata. Level(s) "',
               Predictor2_Extra,
               '": not found in P2_Levels')
      )
    }
  }
}



ALL_NUMERIC<-all(sapply(Data, is.numeric))

if (!ALL_NUMERIC){
  WHICH_NON_NUMERIC<-names(Data)[!sapply(Data,is.numeric)]

  stop(
    'Genomica has detected an error in the data frame setup. The column(s) named: "',
    WHICH_NON_NUMERIC,
    '" contains non-numeric values. Please make sure that Data only contains numerical values (i.e., Data must be formatted with features as rows and samples as columns).'
    )
}


SAMPLE_MATCH_Data<-setdiff(as.character(colnames(Data)),as.character(rownames(Metadata)))
SAMPLE_MATCH_Metadata<-setdiff(as.character(rownames(Metadata)),as.character(colnames(Data)))

if (length(SAMPLE_MATCH_Data)>0){
  stop(
    paste0('Genomica has detected an error in the data frame setup. The column(s) named "',
           SAMPLE_MATCH_Data,
           '" can be found in Data but not in Metadata. Is "', SAMPLE_MATCH_Data,'" a sample?')
  )
}

if (length(SAMPLE_MATCH_Metadata)>0){
  stop(
    paste0('Genomica has detected an error in the data frame setup. The row(s) named "',
           SAMPLE_MATCH_Metadata,
           '" can be found in Metadata but not in Data. Is "', SAMPLE_MATCH_Metadata,'" a sample, and does it appear in one of the columns in Data?')
  )
}



Check_Random<-setdiff(R_Effects,colnames(Metadata))

if (length(Check_Random)>0){
  stop(
    paste0('Genomica has detected an error in the data frame setup. The Random effect "', Check_Random, '" could not be found in Metadata.')
  )
}



if (length(R1_Levels)==1){
    stop(
      paste0('Genomica has detected an error in the data frame setup. Only one level appears to be imput for ', R_Effects,'. Please specify at least 2 levels in R1_Levels.')
    )
}

if (length(R1_Levels)>1){
  if (!is.character(R1_Levels)){
    stop(
      paste0('Genomica has detected an error in the data frame setup. R1_Levels for ', R_Effects,' must be a character vector. Please transform to character before running Genomica.')
    )
  }

  Random1_Duplicate<-R1_Levels[duplicated(R1_Levels)]
  if (length(Random1_Duplicate)>0){
    stop(
      paste0('WARNING: Genomica found that one or more levels in R1_Levels is duplicated: ',
             Random1_Duplicate,' can be found twice in R1_Level - Please rename levels to avoid duplicates.')
    )
  }

  Random1_L<-unique(as.character(Metadata[,R_Effects]))
  Random1_Missing<-setdiff(R1_Levels,Random1_L)
  Random1_Extra<-setdiff(Random1_L,R1_Levels)
  if (length(Random1_Missing)>0){
    stop(
      paste0('Genomica has detected an error in the data frame setup. One or more levels specified for the random effect ',
      R_Effects,' could not be found in Metadata. Level(s) "',
             Random1_Missing,
             '": not found in Metadata.')

    )
  }
  if (length(Random1_Extra)>0){
    stop(
      paste0('Genomica has detected an error in the data frame setup. One or more levels specified for the random effect ',
      R_Effects,' does not match the levels found in Metadata. Level(s) "',
      Random1_Extra,
             '": not found in R1_Levels')
    )
  }
}



  if (length(Predictors)<3){
  print('Genomica is running')
  tictoc::tic()
  shh<- function(...)
  {
    m <- match.call()[-1]
    c <- capture.output(
      tryCatch(
        suppressMessages(suppressWarnings(
          eval(as.list(m)[[1]])
        )), error = function(e) ""))
  }

  shh({
      FOLDER<-paste0('Genomica_Output_',Folder_Name)
      FILE_COMBINED<-'Combined_All_Features'
      unlink(FOLDER,recursive=TRUE,force = TRUE)
      dir.create(FOLDER)

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

      SIGNIFICANCE_THRESHOLD<-FDR_Level

      RESOLUTION<-ImgRes


      if(length(R_Effects)==1){
      DF[,Random1]<-factor((DF[,Random1]),levels=c(R1_Levels))
      Random1<-DF[,Random1]
      } else if (length(R_Effects)==2){
        DF[,Random1]<-factor((DF[,Random1]),levels=c(R1_Levels))
        DF[,Random2]<-factor((DF[,Random2]),levels=c(R2_Levels))
        Random1<-DF[,Random1]
        Random2<-DF[,Random2]
      }

      if (length(Predictors)==1){
        Predictor<-Predictors
      if (length(P1_Levels)==1){
        DF[,Predictor]<-as.numeric(as.character(DF[,Predictor]))
      } else if(length(P1_Levels)>1){
        DF[,Predictor]<-factor(DF[,Predictor],levels = P1_Levels)
      }}

      if (length(Predictors)==2){
        Predictor<-Predictors[1]
        Predictor2<-Predictors[2]

        if (length(P1_Levels)==1){
          DF[,Predictor]<-as.numeric(as.character(DF[,Predictor]))
        } else if(length(P1_Levels)>1){
          DF[,Predictor]<-factor(DF[,Predictor],levels = P1_Levels)
        }
        if (length(P2_Levels)==1){
          DF[,Predictor2]<-as.numeric(as.character(DF[,Predictor2]))
        } else if(length(P2_Levels)>1){
          DF[,Predictor2]<-factor(DF[,Predictor2],levels = P2_Levels)
        }
        }



      y<-colnames(DF)
      y<-y[(length(Metadata)+1):length(y)]



      Combined_Results<-data.frame()
      Model_Diagnosis<-data.frame()


      if (length(Predictors)==1){


        for(i in y)  {


          if (sum(DF[,i],na.rm=TRUE)>0){


            LMM<-as.formula(paste0('`',i,'`', "~", Predictor, '+(1|',R_Effects,')'))
            model<-lmerTest::lmer(LMM,data=DF,na.action='na.pass')

            res<-stats::anova(model)
            Results<-as.data.frame(res)
            Results$Feature_Code<-i
            Results<-Results[,c(7,1:6)]
            Combined_Results<-rbind(Combined_Results,Results)


            Shapiro<-shapiro.test(residuals(model))

            Heteroscedasticity<-lmtest::bptest(residuals(model)~fitted(model))

            Test<-DHARMa::simulateResiduals(model)
            Test_Values<-DHARMa::testResiduals(Test,plot=FALSE)

            Model_Diagnosis<-rbind(Model_Diagnosis,
                                   data.frame(Feature_Code=i,
                                              Normality_Shapiro_Wilk=Shapiro$p.value,
                                              Heteroscedasticity_Breusch_Pagan=Heteroscedasticity$p.value,
                                              Outlier_P_value=Test_Values$outliers$p.value))
          }
        }

      } else if (length(Predictors)==2){

        for(i in y)  {


          if (sum(DF[,i],na.rm=TRUE)>0){

            LMM<-as.formula(paste0('`',i,'`', "~", Predictor, '*', Predictor2, '+(1|',R_Effects,')'))
            model<-lmerTest::lmer(LMM,data=DF,na.action='na.pass')


            res<-stats::anova(model)
            Results<-as.data.frame(res)
            Results$Feature_Code<-i
            Results<-Results[,c(7,1:6)]

            Combined_Results<-rbind(Combined_Results,Results)

            Shapiro<-shapiro.test(residuals(model))

            Heteroscedasticity<-lmtest::bptest(residuals(model)~fitted(model))

            Test<-DHARMa::simulateResiduals(model)
            Test_Values<-DHARMa::testResiduals(Test,plot=FALSE)

            Model_Diagnosis<-rbind(Model_Diagnosis,
                                   data.frame(Feature_Code=i,
                                              Normality_Shapiro_Wilk=Shapiro$p.value,
                                              Heteroscedasticity_Breusch_Pagan=Heteroscedasticity$p.value,
                                              Outlier_P_value=Test_Values$outliers$p.value))

          }}
      }

      rownames(Model_Diagnosis)<-Model_Diagnosis$Feature_Code
      Model_Diagnosis$Normality<-ifelse(Model_Diagnosis$Normality_Shapiro_Wilk>0.05,'Yes','No')
      Model_Diagnosis$Heteroscedasticity<-ifelse(Model_Diagnosis$Heteroscedasticity_Breusch_Pagan<0.05,'Yes','No')
      Model_Diagnosis$Outliers<-ifelse(Model_Diagnosis$Outlier_P_value<0.05,'Yes','No')






      Combined_Results$`Pr(>F)`<-as.numeric(as.character(Combined_Results$`Pr(>F)`))
      if (length(Predictors)==1){
        Combined_Results$Predictor<-rep(c(Predictor),nrow(Combined_Results))
      } else if (length(Predictors)==2){
        Combined_Results$Predictor<-rep(c(Predictor,Predictor2,paste0(Predictor,':',Predictor2,' (interaction)')),nrow(Combined_Results)/3)
      }
      Combined_Results<-Combined_Results[,c(1,8,2:7)]

      Combined_Results$ID<-paste0(Combined_Results$Feature_Code,'_',Combined_Results$Predictor)
      p.values<-data.frame(var=Combined_Results$ID,pval = Combined_Results$`Pr(>F)`)
      p.adj<-fuzzySim::FDR(pvalues = p.values)

      FDR_DF<-p.adj$exclude
      FDR_DF<-rbind(p.adj$exclude,p.adj$select)
      FDR_DF$ID<-rownames(FDR_DF)
      FDR_DF<-FDR_DF[,-1]
      Combined_Results<-dplyr::left_join(Combined_Results,FDR_DF,by='ID')
      Combined_Results<-Combined_Results[,-9]
      Col<-seq(3,9,1)
      Combined_Results[, Col] <- apply(Combined_Results[, Col], 2, function(x) as.numeric(as.character(x)))


      Combined_Results$P_Significant<-ifelse(Combined_Results$`Pr(>F)`<0.05,'Significant','No')


      Combined_Results$P_Adjusted_Significant<-ifelse(Combined_Results$p.adjusted<SIGNIFICANCE_THRESHOLD,'Significant','No')

      Combined_Results2<-dplyr::left_join(Rosetta_Sone,Combined_Results,by = 'Feature_Code')
      Combined_Results2<-Combined_Results2[,-2]

      names(Combined_Results2)[names(Combined_Results2)=='Pr(>F)']<-'P Value'
      names(Combined_Results2)[names(Combined_Results2)=='p.adjusted']<-'P Adjusted (BH)'
      names(Combined_Results2)[names(Combined_Results2)=='P_Significant']<-'Significant (P Value)'
      names(Combined_Results2)[names(Combined_Results2)=='P_Adjusted_Significant']<-'Significant (P Adjusted)'
      writexl::write_xlsx(Combined_Results2,path = Path)
      write.table(Combined_Results2,file = paste0(FOLDER,'/Combined_All_Features.txt'),
                  row.names = FALSE,sep='\t')

      Significant<-subset(Combined_Results,P_Adjusted_Significant=='Significant')

      x<-unique(Significant$Feature_Code)



      DIAGNOSIS_FOLDER<-paste0(FOLDER,'/Model_Diagnosis')
      ifelse(!dir.exists(DIAGNOSIS_FOLDER), dir.create(DIAGNOSIS_FOLDER), "A model diagnosis folder already exists")
      DIAGNOSIS_FOLDER_LMM<-paste0(FOLDER,'/Model_Diagnosis/LMM')
      ifelse(!dir.exists(DIAGNOSIS_FOLDER_LMM), dir.create(DIAGNOSIS_FOLDER_LMM), "A model diagnosis folder already exists")

      Model_Diagnosis<-Model_Diagnosis[Model_Diagnosis$Feature_Code %in% x,]
      Model_Diagnosis<-dplyr::left_join(Rosetta_Sone,Model_Diagnosis,by = 'Feature_Code')
      Model_Diagnosis<-Model_Diagnosis[!is.na(Model_Diagnosis$Normality_Shapiro_Wilk),]
      Model_Diagnosis$Feature_Code<-NULL
      Model_Diagnosis$Concern_Towards_LMM_Appropriateness<-ifelse(Model_Diagnosis$Normality=='No'&Model_Diagnosis$Heteroscedasticity=='Yes'&Model_Diagnosis$Outliers=='Yes',
                                                                  'Yes, Concerning','No')



      Percentage_Non_Normal<-sum(Model_Diagnosis$Normality=='No',na.rm=TRUE)/nrow(Model_Diagnosis)*100
      Percentage_Heteroscedasticity<-sum(Model_Diagnosis$Heteroscedasticity=='Yes',na.rm=TRUE)/nrow(Model_Diagnosis)*100
      Percentage_Significant_Outliers<-sum(Model_Diagnosis$Outliers=='Yes',na.rm=TRUE)/nrow(Model_Diagnosis)*100
      Percentage_Concerning_Cases<-sum(Model_Diagnosis$Concern_Towards_LMM_Appropriateness=='Yes',na.rm=TRUE)/nrow(Model_Diagnosis)*100


      Model_Diagnosis_Print<-Model_Diagnosis


      Description_diagnosis<-Model_Diagnosis[1:11,]
      Description_diagnosis[]<-NA
      Description_diagnosis[1,1]<-'Model (LMM) diagnosis on significant correlations (main effects)'
      Description_diagnosis[2,1]<-paste0('Percentage of LMM outputs with non normal residuals: ',round(Percentage_Non_Normal,2),'%')
      Description_diagnosis[3,1]<-paste0('Percentage of significant residual heteroscedasticity across the significant LMM outputs: ',round(Percentage_Heteroscedasticity,2),'%')
      Description_diagnosis[4,1]<-paste0('Percentage of significant residual outliers across the significant LMM outputs: ',round(Percentage_Significant_Outliers,2),'%')
      Description_diagnosis[5,1]<-paste0('Percentage of LMM outputs with non-normal and significantly heteroscedastic residuals, and significant residual outliers: ',round(Percentage_Concerning_Cases,2),'%')


      THRESHOLD<-10
      Description_diagnosis[7,1]<-
        ifelse(Percentage_Concerning_Cases<=THRESHOLD,paste0(round(Percentage_Concerning_Cases,2),
                                                             '% of the outputs violate all the LMM conditions, thus in support of the validity of the current approach (No warnings or problems detected)'),
               paste0(round(Percentage_Concerning_Cases,2),'% of the outputs violate all the LMM conditions: PROBLEMS DETECTED WITH THE CURRENT MODELLING APPROACH, PLEASE CONSIDER USING DIFFERENT MODELLING STRATEGIES (E.G., GLMM)'))


      Description_diagnosis[9,1]<-'The following table depicts the model diagnosis (i.e., normality, heteroscedasticity and outliers) for the residuals of each significant LMM'

      Description_diagnosis[11,]<-colnames(Model_Diagnosis)

      Model_Diagnosis_Print<-rbind(Description_diagnosis,Model_Diagnosis)
      names(Model_Diagnosis_Print)<-c('Genomica Version 2.0.0',rep("",(length(Model_Diagnosis_Print)-1)))



      writexl::write_xlsx(Model_Diagnosis_Print,path = paste0(DIAGNOSIS_FOLDER_LMM,'/LMM_Model_Diagnosis.xlsx'))
      write.table(Model_Diagnosis_Print,file = paste0(DIAGNOSIS_FOLDER_LMM,'/LMM_Model_Diagnosis.txt'),
                  row.names = FALSE,sep='\t')



      Combined_Significant<-data.frame()


      if (length(Predictors)==1){
      if (length(x)<1){
        return(print(paste0('No significant correlations found, the complete list of results can be found in ',FOLDER)))
        } else if (length(x)>0){
      for (i in x[1:length(x)]){

        model2<-lmerTest::lmer(DF[,i]~DF[,Predictor]+(1|Random1),data=DF,na.action='na.pass')
        sum<-summary(model2)
        summ<-as.data.frame(sum$coefficients)
        summ$Predictor<-rownames(summ)
        summ$Feature_Code<-c(rep(i,nrow(summ)))
        Results2<-as.data.frame(cbind(summ['Feature_Code'],summ['Predictor'],summ[,1:(length(summ)-2)]))
        Combined_Significant<-rbind(Combined_Significant,Results2)

      }}}else if (length(Predictors)==2){
        if (as.numeric(length(x))<1){
          return(print(paste0('No significant correlations found, the complete list of results can be found in ',FOLDER)))
        } else{
          for (i in x[1:length(x)]){

            model2<-lmerTest::lmer(DF[,i]~DF[,Predictor]*DF[,Predictor2]+(1|Random1),data=DF,na.action='na.pass')
            sum<-summary(model2)
            summ<-as.data.frame(sum$coefficients)
            summ$Predictor<-rownames(summ)
            summ$Feature_Code<-c(rep(i,nrow(summ)))
            Results2<-as.data.frame(cbind(summ['Feature_Code'],summ['Predictor'],summ[,1:(length(summ)-2)]))
            Combined_Significant<-rbind(Combined_Significant,Results2)

          }}}


      Combined_Significant<-subset(Combined_Significant,Combined_Significant$Predictor!='(Intercept)')



      if (length(Predictors)==1){
        if (length(P1_Levels)>1){
        Combined_Significant$Predictor<-gsub('DF\\[, Predictor\\]',paste0(Predictor,'_'),Combined_Significant$Predictor)
      }else if (length(P1_Levels)==1){
        Combined_Significant$Predictor<-gsub('DF\\[, Predictor\\]',Predictor,Combined_Significant$Predictor)
      }}else if (length(Predictors)==2){
        if (length(P1_Levels)>1){
          Combined_Significant$Predictor<-gsub('DF\\[, Predictor\\]',paste0(Predictor,'_'),Combined_Significant$Predictor)
        }else if (length(P1_Levels)==1){
          Combined_Significant$Predictor<-gsub('DF\\[, Predictor\\]',Predictor,Combined_Significant$Predictor)
        }
        if (length(P2_Levels)>1){
          Combined_Significant$Predictor<-gsub('DF\\[, Predictor2\\]',paste0(Predictor2,'_'),Combined_Significant$Predictor)
        }else if (length(P2_Levels)==1){
          Combined_Significant$Predictor<-gsub('DF\\[, Predictor2\\]',Predictor2,Combined_Significant$Predictor)
        }}


      Combined_Significant$ID<-paste0(Combined_Significant$Feature_Code,'_',Combined_Significant$Predictor)
      p.values2<-data.frame(var=Combined_Significant$ID,pval = Combined_Significant$`Pr(>|t|)`)

      p.adj2<-fuzzySim::FDR(pvalues = p.values2)
      FDR_DF2<-p.adj2$exclude
      FDR_DF2<-rbind(p.adj2$exclude,p.adj2$select)
      FDR_DF2$ID<-rownames(FDR_DF2)
      FDR_DF2<-FDR_DF2[,-1]
      Combined_Significant<-dplyr::left_join(Combined_Significant,FDR_DF2,by = 'ID')
      Combined_Significant<-Combined_Significant[,-8]
      Combined_Significant$p.adjusted<-as.numeric(as.character(Combined_Significant$p.adjusted))

      Combined_Significant$Significant<-ifelse(Combined_Significant$p.adjusted<SIGNIFICANCE_THRESHOLD,
                                               'Significant','No')
      Combined_Significant<-subset(Combined_Significant,Significant=='Significant')
      Combined_Significant<-Combined_Significant[-length(Combined_Significant)]


      Combined_Significant$Status<-ifelse(Combined_Significant$Estimate<0,'Depleted','Enriched')
      Combined_Significant<-dplyr::left_join(Combined_Significant,Rosetta_Sone,by='Feature_Code')
      Combined_Significant<-Combined_Significant[,-1]
      Combined_Significant<-Combined_Significant[,c(9,1:8)]

      Combined_Significant2<-Combined_Significant



      names(Combined_Significant2)[names(Combined_Significant2)=='Pr(>|t|)']<-'P Value'
      names(Combined_Significant2)[names(Combined_Significant2)=='p.adjusted']<-'P Adjusted (BH)'

      write.table(Combined_Significant2,file = paste0(FOLDER,'/Significant_Comparisons.txt'),
                  row.names = FALSE,sep='\t')
      writexl::write_xlsx(Combined_Significant2,
                 path = paste0(FOLDER,'/Significant_Comparisons.xlsx'))



      ##ENRICHMENT

      #Cumulative enrichment


      Enrichment_List<-as.data.frame(Combined_Significant)
      Enrichment_List$Feature<-sub("\\:.*","",Enrichment_List$Feature)
      test<-Enrichment_List[1,1]
      testk<-substr(test,1,2)
      testk2<-rep('K0',each=1)
      test_do_KO<-identical(testk,testk2)

    if(test_do_KO==TRUE){
      ENRICHMENT_FOLDER<-paste0(FOLDER,'/Enrichment')
      ENRICHED_FOLDER<-paste0(ENRICHMENT_FOLDER,'/Enriched')
      DEPLETED_FOLDER<-paste0(ENRICHMENT_FOLDER,'/Depleted')
      ifelse(!dir.exists(ENRICHMENT_FOLDER), dir.create(ENRICHMENT_FOLDER), "A folder with the same name already exists")
      ifelse(!dir.exists(ENRICHED_FOLDER), dir.create(ENRICHED_FOLDER), "A folder with the same name already exists")
      ifelse(!dir.exists(DEPLETED_FOLDER), dir.create(DEPLETED_FOLDER), "A folder with the same name already exists")




        EnLis<-Combined_Significant
        EnLis$Feature<-sub("\\:.*","",EnLis$Feature)
        EnLis$ID<-EnLis$Predictor
        EnLis<-EnLis[!grepl(':', EnLis$ID),]
        EnLis<-dplyr::mutate(
            tidyr::separate(
              EnLis,
              col=Predictor,
              into=c('Left','Level'),
              sep='_'
            ),
            Predictor1=ifelse(Left==Predictor,Left,NA_character_)
            )
        if (length(Predictors)==2){
             EnLis<-dplyr::mutate(
               EnLis,
               Predictor2=ifelse(Left==Predictor2,Left,NA_character_)
             )}

        En_Pred1<-subset(EnLis,Left==Predictor)
        # En_Pred1<-En_Pred1[,-(length(En_Pred1)-2):-length(En_Pred1)]
        names(En_Pred1)[names(En_Pred1)=='Left']<-'Predictor'
        if (length(Predictors)==2){
          En_Pred2<-subset(EnLis,Left==Predictor2)
          # En_Pred2<-En_Pred2[,-(length(En_Pred2)-2):-length(En_Pred2)]
          names(En_Pred2)[names(En_Pred2)=='Left']<-'Predictor'
        }



        dir.create(paste0(ENRICHED_FOLDER,'/',Predictor))
        dir.create(paste0(DEPLETED_FOLDER,'/',Predictor))

        if (length(Predictors)==2){
          dir.create(paste0(ENRICHED_FOLDER,'/',Predictor2))
          dir.create(paste0(DEPLETED_FOLDER,'/',Predictor2))
        }

          #######Predictor 1 - Cumulative

          Predictor1_Enriched<-as.data.frame(subset(En_Pred1,Status=='Enriched'))
          Predictor1_Enriched_KO<-as.vector(unique(Predictor1_Enriched$Feature))
          Predictor1_Enriched_KO_Enrichment<-MicrobiomeProfiler::enrichKO(gene = Predictor1_Enriched_KO)
          Predictor1_Enriched_KO_Enrichment2<-as.data.frame(Predictor1_Enriched_KO_Enrichment)
          if (nrow(Predictor1_Enriched_KO_Enrichment2)==0){
            Predictor1_Enriched_KO_Enrichment2[2,1]<-paste0(length(Predictor1_Enriched_KO),' over-abundant KOs were found for this predictor = No significant enrichment found')
          }
          write.table(Predictor1_Enriched_KO_Enrichment2,
                      file = paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,'_Cumulative_Vs_Control_Enriched.txt'),
                      row.names = FALSE,sep='\t')
          writexl::write_xlsx(Predictor1_Enriched_KO_Enrichment2,
                              path = paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,'_Cumulative_Vs_Control_Enriched.xlsx'))


          Predictor1_Depleted<-as.data.frame(subset(En_Pred1,Status=='Depleted'))
          Predictor1_Depleted_KO<-as.vector(unique(Predictor1_Depleted$Feature))
          Predictor1_Depleted_KO_Enrichment<-MicrobiomeProfiler::enrichKO(gene = Predictor1_Depleted_KO)
          Predictor1_Depleted_KO_Enrichment2<-as.data.frame(Predictor1_Depleted_KO_Enrichment)

          if (nrow(Predictor1_Depleted_KO_Enrichment2)==0){
            Predictor1_Depleted_KO_Enrichment2[2,1]<-paste0(length(Predictor1_Depleted_KO),' under-abundant KOs were found for this predictor = No significant enrichment found')
          }
          write.table(Predictor1_Depleted_KO_Enrichment2,
                      file = paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,'_Cumulative_Vs_Control_Depleted.txt'),
                      row.names = FALSE,sep='\t')
          writexl::write_xlsx(Predictor1_Depleted_KO_Enrichment2,
                              path = paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,'_Cumulative_Vs_Control_Depleted.xlsx'))


          if (nrow(Predictor1_Enriched_KO_Enrichment2)>5){
            Predictor1_Enriched_KO_Enrichment_pairwise<-enrichplot::pairwise_termsim(Predictor1_Enriched_KO_Enrichment)
            Predictor1_Enriched_KO_Enrichment_Tree<-enrichplot::treeplot(Predictor1_Enriched_KO_Enrichment_pairwise,cladelab_offset=5)
            ggplot2::ggsave(plot = Predictor1_Enriched_KO_Enrichment_Tree,
                            filename = paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,'_Cumulative_Vs_Control_Enriched_tree.tiff'),
                            width = 35, height = 17.5, units = 'cm',dpi=RESOLUTION)}

          if (nrow(Predictor1_Depleted_KO_Enrichment2)>5){
            Predictor1_Depleted_KO_Enrichment_pairwise<-enrichplot::pairwise_termsim(Predictor1_Depleted_KO_Enrichment)
            Predictor1_Depleted_KO_Enrichment_Tree<-enrichplot::treeplot(Predictor1_Depleted_KO_Enrichment_pairwise,cladelab_offset=5)
              ggplot2::ggsave(plot = Predictor1_Depleted_KO_Enrichment_Tree,
                            filename = paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,'_Cumulative_Vs_Control_Depleted_tree.tiff'),
                             width = 35, height = 17.5, units = 'cm',dpi=RESOLUTION)
            }



          #######Predictor 1 - Treatment-wise

          if (length(P1_Levels)>1){
            for (i in P1_Levels[2:length(P1_Levels)]){
              P1_pathtoenr<-paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,'_',i,'_Vs_Control')
              ifelse(!dir.exists(P1_pathtoenr), dir.create(P1_pathtoenr), "A folder with the same name already exists")
              P1_pathtodep<-paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,'_',i,'_Vs_Control')
              ifelse(!dir.exists(P1_pathtodep), dir.create(P1_pathtodep), "A folder with the same name already exists")
            }

            for (i in P1_Levels[2:length(P1_Levels)]){
              P1_DF_LIST_EN<-as.data.frame(subset(Predictor1_Enriched,Level==i))
              P1_Enriched<-MicrobiomeProfiler::enrichKO(gene = P1_DF_LIST_EN$Feature)
              P1_ForTree_Enriched<-P1_Enriched
              P1_Enriched<-as.data.frame(P1_Enriched)
              if (nrow(P1_Enriched)>0){
              P1_E<-clusterProfiler::dotplot(P1_ForTree_Enriched)
              ggplot2::ggsave(plot = P1_E,
                              filename = paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,'_',i,'_Vs_Control/',i,'_Vs_Control_Enrichment.tiff'),
                              width = 25, height = 25, units = 'cm',dpi=RESOLUTION)
              }

              if (nrow(P1_Enriched)==0){
                P1_Enriched[2,1]<-paste0(length(P1_DF_LIST_EN$Feature),' over-abundant KOs were found for this predictor and for this contrast= No significant enrichment found')
              }

              write.table(P1_Enriched,file = paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,'_',i,'_Vs_Control/',i,'_Vs_Control_Enriched.txt'),
                          row.names = FALSE,sep='\t')
              writexl::write_xlsx(P1_Enriched,
                                  path = paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,'_',i,'_Vs_Control/',i,'_Vs_Control_Enriched.xlsx'))
            }

            for (i in P1_Levels[2:length(P1_Levels)]){
              P1_DF_LIST_DEP<-as.data.frame(subset(Predictor1_Depleted,Level==i))
              P1_Depleted<-MicrobiomeProfiler::enrichKO(gene = P1_DF_LIST_DEP$Feature)
              P1_ForTree_Depleted<-P1_Depleted
              P1_Depleted<-as.data.frame(P1_Depleted)
              if (nrow(P1_Depleted)>0){
              P1_DP<-clusterProfiler::dotplot(P1_ForTree_Depleted)
              ggplot2::ggsave(plot = P1_DP,
                              filename = paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,'_',i,'_Vs_Control/',i,'_Vs_Control_Enrichment.tiff'),
                              width = 25, height = 25, units = 'cm',dpi=RESOLUTION)
              }

              if (nrow(P1_Depleted)==0){
                P1_Depleted[2,1]<-paste0(length(P1_DF_LIST_DEP$Feature),' under-abundant KOs were found for this predictor and for this contrast= No significant enrichment found')
              }
              write.table(P1_Depleted,
                          file = paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,'_',i,'_Vs_Control/',i,'_Vs_Control_Depleted.txt'),
                          row.names = FALSE,sep='\t')
              writexl::write_xlsx(P1_Depleted,
                                  path = paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,'_',i,'_Vs_Control/',i,'_Vs_Control_Depleted.xlsx'))

              }
              }



          #######Predictor 2

          if (length(Predictors)==2){

            #######Predictor 2 - Cumulative

            Predictor2_Enriched<-as.data.frame(subset(En_Pred2,Status=='Enriched'))
            Predictor2_Enriched_KO<-as.vector(unique(Predictor2_Enriched$Feature))
            Predictor2_Enriched_KO_Enrichment<-MicrobiomeProfiler::enrichKO(gene = Predictor2_Enriched_KO)
            Predictor2_Enriched_KO_Enrichment2<-as.data.frame(Predictor2_Enriched_KO_Enrichment)

            if (nrow(Predictor2_Enriched_KO_Enrichment2)==0){
              Predictor2_Enriched_KO_Enrichment2[2,1]<-paste0(length(Predictor2_Enriched_KO),' over-abundant KOs were found for this predictor = No significant enrichment found')
            }
            write.table(Predictor2_Enriched_KO_Enrichment2,
                        file = paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,'_Cumulative_Vs_Control_Enriched.txt'),
                        row.names = FALSE,sep='\t')
            writexl::write_xlsx(Predictor2_Enriched_KO_Enrichment2,
                                path = paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,'_Cumulative_Vs_Control_Enriched.xlsx'))


            Predictor2_Depleted<-as.data.frame(subset(En_Pred2,Status=='Depleted'))
            Predictor2_Depleted_KO<-as.vector(unique(Predictor2_Depleted$Feature))
            Predictor2_Depleted_KO_Enrichment<-MicrobiomeProfiler::enrichKO(gene = Predictor2_Depleted_KO)
            Predictor2_Depleted_KO_Enrichment2<-as.data.frame(Predictor2_Depleted_KO_Enrichment)

            if (nrow(Predictor2_Depleted_KO_Enrichment2)==0){
              Predictor2_Depleted_KO_Enrichment2[2,1]<-paste0(length(Predictor2_Depleted_KO),' under-abundant KOs were found for this predictor = No significant enrichment found')
            }
            write.table(Predictor2_Depleted_KO_Enrichment2,
                        file = paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,'_Cumulative_Vs_Control_Depleted.txt'),
                        row.names = FALSE,sep='\t')
            writexl::write_xlsx(Predictor2_Depleted_KO_Enrichment2,
                                path = paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,'_Cumulative_Vs_Control_Depleted.xlsx'))


            if (nrow(Predictor2_Enriched_KO_Enrichment2)>5){
              Predictor2_Enriched_KO_Enrichment_pairwise<-enrichplot::pairwise_termsim(Predictor2_Enriched_KO_Enrichment)
              Predictor2_Enriched_KO_Enrichment_Tree<-enrichplot::treeplot(Predictor2_Enriched_KO_Enrichment_pairwise,cladelab_offset=5)
              ggplot2::ggsave(plot = Predictor2_Enriched_KO_Enrichment_Tree,
                              filename = paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,'_Cumulative_Vs_Control_Enriched_tree.tiff'),
                              width = 35, height = 17.5, units = 'cm',dpi=RESOLUTION)}

            if (nrow(Predictor2_Depleted_KO_Enrichment2)>5){
              Predictor2_Depleted_KO_Enrichment_pairwise<-enrichplot::pairwise_termsim(Predictor2_Depleted_KO_Enrichment)
              Predictor2_Depleted_KO_Enrichment_Tree<-enrichplot::treeplot(Predictor2_Depleted_KO_Enrichment_pairwise,cladelab_offset=5)
              ggplot2::ggsave(plot = Predictor2_Depleted_KO_Enrichment_Tree,
                              filename = paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,'_Cumulative_Vs_Control_Depleted_tree.tiff'),
                              width = 35, height = 17.5, units = 'cm',dpi=RESOLUTION)}


        #######Predictor 2 - Treatment-wise

            ###########################################

            if (length(P2_Levels)>1){
              for (i in P2_Levels[2:length(P2_Levels)]){
                P2_pathtoenr<-paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,'_',i,'_Vs_Control')
                ifelse(!dir.exists(P2_pathtoenr), dir.create(P2_pathtoenr), "A folder with the same name already exists")
                P2_pathtodep<-paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,'_',i,'_Vs_Control')
                ifelse(!dir.exists(P2_pathtodep), dir.create(P2_pathtodep), "A folder with the same name already exists")
              }

              for (i in P2_Levels[2:length(P2_Levels)]){
                P2_DF_LIST_EN<-as.data.frame(subset(Predictor2_Enriched,Level==i))
                P2_Enriched<-MicrobiomeProfiler::enrichKO(gene = P2_DF_LIST_EN$Feature)
                P2_ForTree_Enriched<-P2_Enriched
                P2_Enriched<-as.data.frame(P2_Enriched)
                if (nrow(P2_Enriched)>0){
                  P2_E<-clusterProfiler::dotplot(P2_ForTree_Enriched)
                  ggplot2::ggsave(plot = P2_E,
                                  filename = paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,'_',i,'_Vs_Control/',i,'_Vs_Control_Enrichment.tiff'),
                                  width = 25, height = 25, units = 'cm',dpi=RESOLUTION)
                  }

                if (nrow(P2_Enriched)==0){
                  P2_Enriched[2,1]<-paste0(length(P2_DF_LIST_EN$Feature),' over-abundant KOs were found for this predictor and for this contrast= No significant enrichment found')
                }

                write.table(P2_Enriched,file = paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,'_',i,'_Vs_Control/',i,'_Vs_Control_Enriched.txt'),
                            row.names = FALSE,sep='\t')
                writexl::write_xlsx(P2_Enriched,
                                    path = paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,'_',i,'_Vs_Control/',i,'_Vs_Control_Enriched.xlsx'))

              }


              for (i in P2_Levels[2:length(P2_Levels)]){
                P2_DF_LIST_DEP<-as.data.frame(subset(Predictor2_Depleted,Level==i))
                P2_Depleted<-MicrobiomeProfiler::enrichKO(gene = P2_DF_LIST_DEP$Feature)
                P2_ForTree_Depleted<-P2_Depleted
                P2_Depleted<-as.data.frame(P2_Depleted)
                if (nrow(P2_Depleted)>0){
                  P2_DP<-clusterProfiler::dotplot(P2_ForTree_Depleted)
                  ggplot2::ggsave(plot = P2_DP,
                                  filename = paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,'_',i,'_Vs_Control/',i,'_Vs_Control_Enrichment.tiff'),
                                  width = 25, height = 25, units = 'cm',dpi=RESOLUTION)
                }

                if (nrow(P2_Depleted)==0){
                  P2_Depleted[2,1]<-paste0(length(P2_DF_LIST_DEP$Feature),' under-abundant KOs were found for this predictor and for this contrast= No significant enrichment found')
                }
                write.table(P2_Depleted,
                            file = paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,'_',i,'_Vs_Control/',i,'_Vs_Control_Depleted.txt'),
                            row.names = FALSE,sep='\t')
                writexl::write_xlsx(P2_Depleted,
                                    path = paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,'_',i,'_Vs_Control/',i,'_Vs_Control_Depleted.xlsx'))

              }
            }

        }
    }




  })

  FOLDER<-paste0('Genomica_Output_',Folder_Name)
  FILE_COMBINED<-'Combined_All_Features'
  Path<-paste0(FOLDER,'/',FILE_COMBINED,'.xlsx')

  tictoc::toc()

  SignificantDF<-as.data.frame(readxl::read_xlsx(Path))
  SignificantDF<-subset(SignificantDF,`Significant (P Adjusted)`=='Significant')
  if (nrow(SignificantDF)==0){
    return(print(paste0('No significant correlations found, the complete list of results can be found in ',FOLDER)))
  }

  Enrichment_List_DF<-as.data.frame(readxl::read_xlsx(Path))
  Enrichment_List_DF$Feature<-sub("\\:.*","",Enrichment_List_DF$Feature)
  test_Post<-Enrichment_List_DF[1,1]
  testk_Post<-substr(test_Post,1,2)
  testk2_Post<-rep('K0',each=1)
  test_do_KO_Post<-identical(testk_Post,testk2_Post)
  if(test_do_KO_Post==FALSE){
    print('Features are not KOs; enrichment analysis was not performed')
  }

  print(paste0('Done! Outputs saved in ',FOLDER))
  }}

