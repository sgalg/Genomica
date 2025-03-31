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
#' @import readxl
#' @import writexl
#' @import lmerTest
#' @import dplyr
#' @import fuzzySim
#' @import MicrobiomeProfiler
#' @import enrichplot
#' @import clusterProfiler
#' @import tictoc
#' @export
#' @examples
#' #Load Data data frame (Data_Demo)
#' Data<-Data_Demo
#' #In this example Data is a data frame of 500 KEGG orthologs (KOs) in rows, across 35 samples (columns)
#' print(Data[1:5,1:5])
#' Metadata<-Metadata_Demo
#' #in this example Metadata is a data frame with 35 rows (sample IDs) and 5 columns (SampleID, Treatment and Block). Treatment contains the treatment layout (1,2,3,4,5) and will be the fixed effect of the model, whilst block describes the stratification of this experimental layout according to the random effect "Room".
#' print(Metadata[1:5,])
#' genomica(Data = Data,Metadata = Metadata,
#' Predictors = c('Treatment'),P1_Levels=c('1','2','3','4','5'),
#' R_Effects=c('Block'),R1_Levels=c('1','2','3','5','6','7'),
#' Already_Log10_transformed = c('NO'),Folder_Name=c('Test'))


genomica<-function(Data,Metadata,
                   Predictors,P1_Levels=c(0),P2_Levels=c(0),
                   R_Effects,R1_Levels,
                   Already_Log10_transformed=c('NO'),Folder_Name){


  if (length(Predictors)>2){
    print('Genomica currently works with maximum two predictors, please amend the `Predictors` vector in order to include maximum 2 predictors')
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



      if (length(Predictors)==1){
      for(i in y[1:length(y)])  {


        if ((sum(DF[,i]))>0){

          model<-lmerTest::lmer(DF[,i]~DF[,Predictor]+(1|Random1),data=DF,na.action='na.pass')
          res<-stats::anova(model)
          Results<-as.data.frame(res)
          Results$Feature_Code<-c(i)
          Results<-Results[,c(7,1:6)]
          Combined_Results<-rbind(Combined_Results,Results)

        }}}else if (length(Predictors)==2){

          for(i in y[1:length(y)])  {


            if ((sum(DF[,i]))>0){

              model<-lmerTest::lmer(DF[,i]~DF[,Predictor]*DF[,Predictor2]+(1|Random1),data=DF,na.action='na.pass')
              res<-stats::anova(model)
              Results<-as.data.frame(res)
              Results$Feature_Code<-c(i)
              Results<-Results[,c(7,1:6)]
              Combined_Results<-rbind(Combined_Results,Results)
        }}}





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
      Combined_Results$Q_Significant<-ifelse(Combined_Results$p.adjusted<0.05,'Significant','No')

      Combined_Results2<-dplyr::left_join(Rosetta_Sone,Combined_Results,by = 'Feature_Code')
      Combined_Results2<-Combined_Results2[,-2]
      writexl::write_xlsx(Combined_Results2,path = Path)
      write.table(Combined_Results2,file = paste0(FOLDER,'/Combined_All_Features.txt'),
                  row.names = FALSE,sep='\t')



      Significant<-subset(Combined_Results,Q_Significant=='Significant')

      x<-unique(Significant$Feature_Code)

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
      Combined_Significant$Significant<-ifelse(Combined_Significant$p.adjusted<0.05,
                                               'Significant','No')
      Combined_Significant<-subset(Combined_Significant,Significant=='Significant')
      Combined_Significant<-Combined_Significant[-length(Combined_Significant)]

      Combined_Significant$Predictor<-sub(pattern = '[[]',replacement = '`',Combined_Significant$Predictor)
      Combined_Significant$Predictor<-sub(pattern = '[]]',replacement = '`',Combined_Significant$Predictor)
      Combined_Significant$Predictor<-sub(pattern = '[[]',replacement = '`',Combined_Significant$Predictor)
      Combined_Significant$Predictor<-sub(pattern = '[]]',replacement = '`',Combined_Significant$Predictor)
      Combined_Significant$Predictor<-sub(pattern = 'DF`,',replacement = '',Combined_Significant$Predictor)
      Combined_Significant$Predictor<-sub(pattern = 'DF`,',replacement = '',Combined_Significant$Predictor)

      Combined_Significant$Status<-ifelse(Combined_Significant$Estimate<0,'Depleted','Enriched')
      Combined_Significant<-dplyr::left_join(Combined_Significant,Rosetta_Sone,by='Feature_Code')
      Combined_Significant<-Combined_Significant[,-1]
      Combined_Significant<-Combined_Significant[,c(9,1:8)]
      Combined_Significant2<-Combined_Significant
      Combined_Significant2$Predictor<-sub(pattern = ' Predictor`',replacement = paste0(Predictor,'_'),Combined_Significant2$Predictor)
      if (length(Predictors)==2){
        Combined_Significant2$Predictor<-sub(pattern = ' Predictor2`',replacement = paste0(Predictor2,'_'),Combined_Significant2$Predictor)
      }
      write.table(Combined_Significant2,file = paste0(FOLDER,'/Significant_Comparisons.txt'),
                  row.names = FALSE,sep='\t')
      writexl::write_xlsx(Combined_Significant2,
                 path = paste0(FOLDER,'/Significant_Comparisons.xlsx'))



      ##ENRICHMENT

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
        EnLis$ID<-sub(' Predictor`',Predictor,EnLis$ID)
        if (length(Predictors)==2){
        EnLis$ID<-sub(' Predictor2`',Predictor2,EnLis$ID)}
        EnLis$ID2<-EnLis$ID
        if (length(P1_Levels)>1){
        for (i in P1_Levels[2:length(P1_Levels)]){
          EnLis$ID2<-sub(paste0(Predictor,i),Predictor,EnLis$ID2)}
        }
        if (length(Predictors)==2){
          if (length(P2_Levels)>1){
          for (i in P2_Levels[2:length(P2_Levels)]){
          EnLis$ID2<-sub(paste0(Predictor2,i),Predictor2,EnLis$ID2)}
        }}



        Predictor1_Enrichment_List<-EnLis
        SFAS<-Predictor
        Predictor1_Enrichment_List<-subset(Predictor1_Enrichment_List,ID2==SFAS)
        if (length(Predictors)==2){
          Predictor2_Enrichment_List<-EnLis
          SFAS2<-Predictor2
          Predictor2_Enrichment_List<-subset(Predictor2_Enrichment_List,ID2==SFAS2)
        }


        dir.create(paste0(ENRICHED_FOLDER,'/',Predictor))
        dir.create(paste0(DEPLETED_FOLDER,'/',Predictor))

        if (length(Predictors)==2){
          dir.create(paste0(ENRICHED_FOLDER,'/',Predictor2))
          dir.create(paste0(DEPLETED_FOLDER,'/',Predictor2))
        }

          Predictor1_Enriched<-as.data.frame(subset(Predictor1_Enrichment_List,Status=='Enriched'))
          Predictor1_Enriched_KO<-as.vector(unique(Predictor1_Enriched$Feature))
          Predictor1_Enriched_KO_Enrichment<-MicrobiomeProfiler::enrichKO(gene = Predictor1_Enriched_KO)
          Predictor1_Enriched_KO_Enrichment2<-as.data.frame(Predictor1_Enriched_KO_Enrichment)
          write.table(Predictor1_Enriched_KO_Enrichment2,
                      file = paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,'_Cumulative_Vs_Control_Enriched.txt'),
                      row.names = FALSE,sep='\t')
          writexl::write_xlsx(Predictor1_Enriched_KO_Enrichment2,
                              path = paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,'_Cumulative_Vs_Control_Enriched.xlsx'))


          Predictor1_Depleted<-as.data.frame(subset(Predictor1_Enrichment_List,Status=='Depleted'))
          Predictor1_Depleted_KO<-as.vector(unique(Predictor1_Depleted$Feature))
          Predictor1_Depleted_KO_Enrichment<-MicrobiomeProfiler::enrichKO(gene = Predictor1_Depleted_KO)
          Predictor1_Depleted_KO_Enrichment2<-as.data.frame(Predictor1_Depleted_KO_Enrichment)
          write.table(Predictor1_Depleted_KO_Enrichment2,
                      file = paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,'_Cumulative_Vs_Control_Depleted.txt'),
                      row.names = FALSE,sep='\t')
          writexl::write_xlsx(Predictor1_Depleted_KO_Enrichment2,
                              path = paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,'_Cumulative_Vs_Control_Depleted.xlsx'))

          if (length(Predictor1_Enriched_KO_Enrichment2$ID)>5){
            Predictor1_Enriched_KO_Enrichment_pairwise<-enrichplot::pairwise_termsim(Predictor1_Enriched_KO_Enrichment)
            Predictor1_Enriched_KO_Enrichment_Tree<-enrichplot::treeplot(Predictor1_Enriched_KO_Enrichment_pairwise)
            ggplot2::ggsave(plot = Predictor1_Enriched_KO_Enrichment_Tree,
                            filename = paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,'_Cumulative_Vs_Control_Enriched_tree.tiff'),
                            width = 35, height = 17.5, units = 'cm',dpi=1200)}

          if (length(Predictor1_Depleted_KO_Enrichment2$ID)>5){
            Predictor1_Depleted_KO_Enrichment_pairwise<-enrichplot::pairwise_termsim(Predictor1_Depleted_KO_Enrichment)
            Predictor1_Depleted_KO_Enrichment_Tree<-enrichplot::treeplot(Predictor1_Depleted_KO_Enrichment_pairwise)
            ggplot2::ggsave(plot = Predictor1_Depleted_KO_Enrichment_Tree,
                            filename = paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,'_Cumulative_Vs_Control_Depleted_tree.tiff'),
                            width = 35, height = 17.5, units = 'cm',dpi=1200)}


          if (length(P1_Levels)>1){
            for (i in P1_Levels[2:length(P1_Levels)]){
              P1_pathtoenr<-paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,i,'_Vs_Control')
              ifelse(!dir.exists(P1_pathtoenr), dir.create(P1_pathtoenr), "A folder with the same name already exists")
              P1_pathtodep<-paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,i,'_Vs_Control')
              ifelse(!dir.exists(P1_pathtodep), dir.create(P1_pathtodep), "A folder with the same name already exists")
            }
            for (i in P1_Levels[2:length(P1_Levels)]){
              P1_DF_LIST_EN<-c(paste0(i,'_Enriched'))
              P1_DF_LIST_EN<-as.data.frame(subset(Predictor1_Enrichment_List,Predictor==paste0(' Predictor`',i)&Status=='Enriched'))
              P1_Enriched<-as.data.frame(MicrobiomeProfiler::enrichKO(gene = P1_DF_LIST_EN$Feature))
              write.table(P1_Enriched,file = paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,i,'_Vs_Control/',i,'_Vs_Control_Enriched.txt'),
                          row.names = FALSE,sep='\t')
              writexl::write_xlsx(P1_Enriched,
                                  path = paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,i,'_Vs_Control/',i,'_Vs_Control_Enriched.xlsx'))
            }
            for (i in P1_Levels[2:length(P1_Levels)]){
              P1_DF_LIST_DEP<-c(paste0(i,'_Depleted'))
              P1_DF_LIST_DEP<-as.data.frame(subset(Predictor1_Enrichment_List,Predictor==paste0(' Predictor`',i)&Status=='Depleted'))
              P1_Depleted<-as.data.frame(MicrobiomeProfiler::enrichKO(gene = P1_DF_LIST_DEP$Feature))
              write.table(P1_Depleted,
                          file = paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,i,'_Vs_Control/',i,'_Vs_Control_Depleted.txt'),
                          row.names = FALSE,sep='\t')
              writexl::write_xlsx(P1_Depleted,
                                  path = paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,i,'_Vs_Control/',i,'_Vs_Control_Depleted.xlsx'))
            }
            for (i in P1_Levels[2:length(P1_Levels)]){
              P1_DF_LIST_EN<-c(paste0(i,'_Enriched'))
              P1_DF_LIST_EN<-as.data.frame(subset(Predictor1_Enrichment_List,Predictor==paste0(' Predictor`',i)&Status=='Enriched'))
              P1_ForTree_Enriched<-MicrobiomeProfiler::enrichKO(gene = P1_DF_LIST_EN$Feature)
              P1_DP<-clusterProfiler::dotplot(P1_ForTree_Enriched)
              ggplot2::ggsave(plot = P1_DP,
                              filename = paste0(ENRICHED_FOLDER,'/',Predictor,'/',Predictor,i,'_Vs_Control/',i,'_Vs_Control_Enrichment.tiff'),
                              width = 25, height = 25, units = 'cm',dpi=1200)
              }
            for (i in P1_Levels[2:length(P1_Levels)]){
              P1_DF_LIST_DEP<-c(paste0(i,'_Depleted'))
              P1_DF_LIST_DEP<-as.data.frame(subset(Predictor1_Enrichment_List,Predictor==paste0(' Predictor`',i)&Status=='Depleted'))
              P1_ForTree_Depleted<-MicrobiomeProfiler::enrichKO(gene = P1_DF_LIST_DEP$Feature)
              P1_DP<-clusterProfiler::dotplot(P1_ForTree_Depleted)
              ggplot2::ggsave(plot = P1_DP,
                              filename = paste0(DEPLETED_FOLDER,'/',Predictor,'/',Predictor,i,'_Vs_Control/',i,'_Vs_Control_Enrichment.tiff'),
                              width = 25, height = 25, units = 'cm',dpi=1200)
            }}



        if (length(Predictors)==2){

            Predictor2_Enriched<-as.data.frame(subset(Predictor2_Enrichment_List,Status=='Enriched'))
            Predictor2_Enriched_KO<-as.vector(unique(Predictor2_Enriched$Feature))
            Predictor2_Enriched_KO_Enrichment<-MicrobiomeProfiler::enrichKO(gene = Predictor2_Enriched_KO)
            Predictor2_Enriched_KO_Enrichment2<-as.data.frame(Predictor2_Enriched_KO_Enrichment)
            write.table(Predictor2_Enriched_KO_Enrichment2,
                        file = paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,'_Cumulative_Vs_Control_Enriched.txt'),
                        row.names = FALSE,sep='\t')
            writexl::write_xlsx(Predictor2_Enriched_KO_Enrichment2,
                                path = paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,'_Cumulative_Vs_Control_Enriched.xlsx'))


            Predictor2_Depleted<-as.data.frame(subset(Predictor2_Enrichment_List,Status=='Depleted'))
            Predictor2_Depleted_KO<-as.vector(unique(Predictor2_Depleted$Feature))
            Predictor2_Depleted_KO_Enrichment<-MicrobiomeProfiler::enrichKO(gene = Predictor2_Depleted_KO)
            Predictor2_Depleted_KO_Enrichment2<-as.data.frame(Predictor2_Depleted_KO_Enrichment)
            write.table(Predictor2_Depleted_KO_Enrichment2,
                        file = paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,'_Cumulative_Vs_Control_Depleted.txt'),
                        row.names = FALSE,sep='\t')
            writexl::write_xlsx(Predictor2_Depleted_KO_Enrichment2,
                                path = paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,'_Cumulative_Vs_Control_Depleted.xlsx'))
            if (length(Predictor2_Enriched_KO_Enrichment2$ID)>5){
              Predictor2_Enriched_KO_Enrichment_pairwise<-enrichplot::pairwise_termsim(Predictor2_Enriched_KO_Enrichment)
              Predictor2_Enriched_KO_Enrichment_Tree<-enrichplot::treeplot(Predictor2_Enriched_KO_Enrichment_pairwise)
              ggplot2::ggsave(plot = Predictor2_Enriched_KO_Enrichment_Tree,
                              filename = paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,'_Cumulative_Vs_Control_Enriched_tree.tiff'),
                              width = 35, height = 17.5, units = 'cm',dpi=1200)}
            if (length(Predictor2_Depleted_KO_Enrichment2$ID)>5){
              Predictor2_Depleted_KO_Enrichment_pairwise<-enrichplot::pairwise_termsim(Predictor2_Depleted_KO_Enrichment)
              Predictor2_Depleted_KO_Enrichment_Tree<-enrichplot::treeplot(Predictor2_Depleted_KO_Enrichment_pairwise)
              ggplot2::ggsave(plot = Predictor2_Depleted_KO_Enrichment_Tree,
                              filename = paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,'_Cumulative_Vs_Control_Depleted_tree.tiff'),
                              width = 35, height = 17.5, units = 'cm',dpi=1200)}

      if (length(P2_Levels)>1){
              for (i in P2_Levels[2:length(P2_Levels)]){
                P2_pathtoenr<-paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,i,'_Vs_Control')
                ifelse(!dir.exists(P2_pathtoenr), dir.create(P2_pathtoenr), "A folder with the same name already exists")
                P2_pathtodep<-paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,i,'_Vs_Control')
                ifelse(!dir.exists(P2_pathtodep), dir.create(P2_pathtodep), "A folder with the same name already exists")
              }
        for (i in P2_Levels[2:length(P2_Levels)]){
                P2_DF_LIST_EN<-c(paste0(i,'_Enriched'))
                P2_DF_LIST_EN<-as.data.frame(subset(Predictor2_Enrichment_List,Predictor==paste0(' Predictor2`',i)&Status=='Enriched'))
                P2_Enriched<-as.data.frame(MicrobiomeProfiler::enrichKO(gene = P2_DF_LIST_EN$Feature))
                write.table(P2_Enriched,file = paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,i,'_Vs_Control/',i,'_Vs_Control_Enriched.txt'),
                            row.names = FALSE,sep='\t')
                writexl::write_xlsx(P2_Enriched,
                                    path = paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,i,'_Vs_Control/',i,'_Vs_Control_Enriched.xlsx'))
        }
        for (i in P2_Levels[2:length(P2_Levels)]){
            P2_DF_LIST_DEP<-c(paste0(i,'_Depleted'))
                P2_DF_LIST_DEP<-as.data.frame(subset(Predictor2_Enrichment_List,Predictor==paste0(' Predictor2`',i)&Status=='Depleted'))
                P2_Depleted<-as.data.frame(MicrobiomeProfiler::enrichKO(gene = P2_DF_LIST_DEP$Feature))
                write.table(P2_Depleted,
                            file = paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,i,'_Vs_Control/',i,'_Vs_Control_Depleted.txt'),
                            row.names = FALSE,sep='\t')
                writexl::write_xlsx(P2_Depleted,
                                    path = paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,i,'_Vs_Control/',i,'_Vs_Control_Depleted.xlsx'))
        }
        for (i in P2_Levels[2:length(P2_Levels)]){
                P2_DF_LIST_EN<-c(paste0(i,'_Enriched'))
                P2_DF_LIST_EN<-as.data.frame(subset(Predictor2_Enrichment_List,Predictor==paste0(' Predictor2`',i)&Status=='Enriched'))
                P2_ForTree_Enriched<-MicrobiomeProfiler::enrichKO(gene = P2_DF_LIST_EN$Feature)
                P2_DP<-clusterProfiler::dotplot(P2_ForTree_Enriched)
                ggplot2::ggsave(plot = P2_DP,
                                filename = paste0(ENRICHED_FOLDER,'/',Predictor2,'/',Predictor2,i,'_Vs_Control/',i,'_Vs_Control_Enrichment.tiff'),
                                width = 25, height = 25, units = 'cm',dpi=1200)
        }
        for (i in P2_Levels[2:length(P2_Levels)]){
                P2_DF_LIST_DEP<-c(paste0(i,'_Depleted'))
                P2_DF_LIST_DEP<-as.data.frame(subset(Predictor2_Enrichment_List,Predictor==paste0(' Predictor2`',i)&Status=='Depleted'))
                P2_ForTree_Depleted<-MicrobiomeProfiler::enrichKO(gene = P2_DF_LIST_DEP$Feature)
                P2_DP<-clusterProfiler::dotplot(P2_ForTree_Depleted)
                ggplot2::ggsave(plot = P2_DP,
                                filename = paste0(DEPLETED_FOLDER,'/',Predictor2,'/',Predictor2,i,'_Vs_Control/',i,'_Vs_Control_Enrichment.tiff'),
                                width = 25, height = 25, units = 'cm',dpi=1200)
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
  SignificantDF<-subset(SignificantDF,Q_Significant=='Significant')
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

