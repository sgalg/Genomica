# Introduction to Genomica

## Welcome to Genomica

Genomica provides a suite for the analysis of functional orthologs from the KEGG Orthology (KO) database. Genomica is based on a linear mixed model approach, whilst also allowing to test for interactions between maximum two predictors, therefore outputting both a list of significant, FDR adjusted KOs across the treatment layout and the results of the enrichment analysis based on the LMM-informed enriched or depleted orthologs, both in the different treatment groups and cumulatively.

Only two data frames are needed as input to run Genomica, one containing the data (KO relative abundance) and one containing metadata information (e.g., treatment layout), moreover if required, Genomica will carry out the log10 transformation of pre-normalised data.
The outputs from Genomica are summarised both in (tab-delimited) .txt and in .xlsx files describing the whole list of statistical results per feature and a list of significant comparisons, respectively. Moreover, an Enrichment directory is created within the output directory, which contains the results (.txt, .xlsx and .tiff) from the enrichment analysis.

The tests carried out in Genomica develop through linear mixed models (LMM) via calling the package lmer4 (Bates et al., 2015) and calculating the P-value via type III ANOVA using the Satterthwaite's method through the package “lmerTest” (Kuznetsova et al., 2017).
Moreover, the tests, and ultimately the significance levels are based on the calculation of the false discovery rate (type I error) under repeated testing. This is carried out via calling the package fuzzySim (A. Marcia Barbosal, 2015), which uses the Benjamini & Hochberg correction to generate p.adjust values based on the p values generated by the LMM.
Finally, the enrichment analysis is performed via calling MicrobiomeProfiler (Yu G., Chen M, 2024), which, through the analyses performed in Genomica also requires enrichplot (Yu G., 2025) and clusterProfiler (Xu, S. et al., 2024, Yu G, 2012).

If used to analyse other types of complex datasets instead of KOs (e.g., AMR tables or quantitative gene tables from qPCR experiments), Genomica will only return the lists of significant comparisons across the predictors, as established by the user, without performing the enrichment analysis.

#### If you use Genomica, please cite it as below:

Salvatore, G. Genomica; Linear mixed model based, multiple hypothesis testing corrected, ortholog enrichment analysis (Version 1.2.0). https://doi.org/https://doi.org/10.58073/SRUC.28695350.v1

## Installation

The enrichment analysis carried out in Genomica depends on MicrobiomeProfiler (Yu G., Chen M, 2024) and devtools, thus you would need to install these packages first:

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MicrobiomeProfiler")

install.packages("devtools")
```
Thus, you could proceed installing Genomica directly from GitHub:

```{r}
devtools::install_github('sgalg/Genomica')
```

In order to carry out the analyses, you will need to load two data frames (i.e., Data and Metadata). In this section, some suggestions are provided in order to prepare the data frames for the analyses.

In Genomica, you will find two pre-loaded demo data frames, i.e., Data_Demo (500 cpm-normalised KOs across 35 samples) and Metadata_Demo (informing on the study layout, such as treatment and random factor allocation relative to these 35 samples). Hereafter, we will use these demo data frames for tutorial purposes. To call these data frames simply load the package via using:

```{r}
library(Genomica)
```
and then you can store the demo data frames into variables, by typing, for example:
```{r}
Data<-Data_Demo
Metadata<-Metadata_Demo
```

### Data

IMPORTANT, the data frame "Data" must be formatted with features as rows and samples as columns, however, currently we can see that the list of features (KOs) is embedded within the first column:

```{r}
print(Data[1:5,1:5])
```

In order to proceed with the analysis, you must assign the row names and delete the feature-embedded column. This can be done as follows:
```{r}
rownames(Data)<-Data$KO
Data<-Data[,-1]
print(Data[1:5,1:5])
```
At this point, your Data is ready for downstream analyses with Genomica.

IMPORTANT: Please delete the rows containing the abundance of UNMAPPED and UNGROUPED KOs, as these would interfere with the enrichment.


### Metadata
Your metadata must be formatted with samples as rows and features as columns (opposite to what seen for Data).
IMPORTANT: the sample names must be assigned as rows (row names):

```{r}
print(Metadata[1:5,])
rownames(Metadata)<-Metadata$SampleID
```

## Performing an analysis with Genomica

Genomica first carries out a linear mix model on all the features in your Data, thus it will further analyse the significant comparisons (according to the false discovery rate) and provide information on the eventual significant comparisons within the specified predictor. All the significant comparisons found at this stage will be sorted according to the model estimate in enriched or depleted KOs, which will then used for the enrichment analysis.

IMPORTANT, Genomica will not carry out the linear mixed model of the feature whose quantity is 0 throughout all the samples in the data frame.

Currently, Genomica allows to perform the analysis either with a single or with two predictors, and in the latter case, it allows to test for the interaction between the two predictors. Everything is done by calling the function genomica():

```{r}
genomica(Data = Data, Metadata = Metadata,
         Predictors = c('Treatment'),P1_Levels = c('1','2','3','4','5'),
         R_Effects = c('Block'),R1_Levels = c('1','2','3','4','5','6','7'),
         Already_Log10_transformed = c('NO'),Folder_Name = c('Test'))
```

So, in this particular case, the metadata only contains one predictor (Treatment), whose levels (i.e., different treatment groups) are Treatment 1 to 5, with Treatment 1 being the control group.

Therefore, when setting the levels for this predictor through P1_Levels, we would organise the vector in such a way that "1" (i.e., with levels as characters) is the first number.

#### IMPORTANT: if your predictor is a numerical variable, you can set the levels via: "P1_Levels=c(0)" (i.e., with 0 as numerical) which indicates 0 levels in your predictor.

#### IMPORTANT: if you had two predictors, you could add these in the “Predictors” vector (e.g., Predictors=c(‘Treatment’, ‘Phase’)), and specify the levels for the second predictor in P2_Levels (e.g., P2_Levels=c(‘Starter’, ‘Grower’, ‘Finisher’)).

The random factor in this metadata is "Block", indeed in this example animal study, the different statistical units were organised in a total of 7 blocks (which are summarised in the levels for the random effect with R1_Levels).

Folder name is used to label your Genomica_Output directory. In this case the output directory name will be "Genomica_Output_Test".

## Results

The results for the analysis are all organised in the the "Genomica_Output" directory:
* Genomica_Output
  * Combined_All_Features (This file, saved both as .txt and .xslx summarises the LMM results for all the feature in Data)
  * Significant_Comparison (This file, saved both as .txt and .xslx summarises the significant comparisons, via also including a pair-wise analysis for all the levels in the predictor)
  * Enrichment (This directory will store the results of the enrichment analysis)
    * Enriched (Directory storing the enrichment analysis results for the enriched orhtologs)
      * Predictor 1 (Genomica will create a folder for each predictor)
        * Cumulative_Vs_Control_Enriched (This file, saved both as .txt and .xslx summarises the p.adjusted enriched functions across the orthologs)
        * If there are more than five function a publication-rady 1,200 dpi tree.tiff figure will be generated.
        * Level 1 (a directory will be created for every predictor level, in which the p.adjusted enriched functions are stored together with a publication-ready 1,200 dpi dot plot.tiff file)
    * Depleted (Directory storing the enrichment analysis results for the depleted orhtologs)
      * Predictor 1 (Genomica will create a folder for each predictor)
        * Cumulative_Vs_Control_Enriched (This file, saved both as .txt and .xslx summarises the p.adjusted enriched functions across the orthologs)
        * If there are more than five function a publication-rady 1,200 dpi tree.tiff figure will be generated.
        * Level 1 (a directory will be created for every predictor level, in which the p.adjusted enriched functions are stored together with a publication-ready 1,200 dpi dot plot.tiff file)


## References

Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1–48. https://doi.org/10.18637/jss.v067.i01

Kuznetsova, A., Brockhoff, P. B., & Christensen, R. H. B. (2017). lmerTest Package: Tests in Linear Mixed Effects Models. Journal of Statistical Software, 82(13), 1–26. https://doi.org/10.18637/jss.v082.i13

Barbosa, A.M. (2015), fuzzySim: applying fuzzy logic to binary similarity indices in ecology. Methods Ecol Evol, 6: 853-858. https://doi.org/10.1111/2041-210X.12372

Yu G, Chen M (2024). MicrobiomeProfiler: An R/shiny package for microbiome functional enrichment analysis. R package version 1.12.0, https://yulab-smu.top/contribution-knowledge-mining/, https://github.com/YuLab-SMU/MicrobiomeProfiler/.

Yu G (2025). enrichplot: Visualization of Functional Enrichment Result.
doi:10.18129/B9.bioc.enrichplot <https://doi.org/10.18129/B9.bioc.enrichplot>, R package

Xu, S., Hu, E., Cai, Y. et al. Using clusterProfiler to characterize multiomics data. Nat Protoc 19, 3292–3320 (2024). https://doi.org/10.1038/s41596-024-01020-z

Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012 May;16(5):284-7. doi: 10.1089/omi.2011.0118. Epub 2012 Mar 28. PMID: 22455463; PMCID: PMC3339379.


