#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("timecourse")

## 15-20 needed for adding to the cluster
options(repos = BiocManager::repositories())

list.of.packages <- c("shiny", "shinythemes","shinyWidgets","reactable","tidyverse","plotly","cowplot","GGally","ggrepel","stringi","FactoMineR","factoextra","igraph","ggraph","tidygraph","pheatmap","readr", "ComplexUpset","eulerr","RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#"BiocManager","timecourse",

library(shiny)
library(shinythemes)
library(shinyWidgets)
library(reactable)
library(tidyverse)
library(plotly)
library(cowplot)
library(GGally)
library(ggrepel)
library(stringi)
library(FactoMineR)
library(factoextra)
library(igraph)
library(ggraph)
library(tidygraph)
library(pheatmap)
library(timecourse)
library(readr)
library(ComplexUpset)
library(eulerr)
library(RColorBrewer)



# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("timecourse")
#library(timecourse)

"%!in%" <- function(x,y)!('%in%'(x,y))

locale_input <- read.csv("org_shortcut_prot_ribo_withMouse2_YEAST.csv")
# go.df<-read.csv("human_go_referencefile.csv")%>%
#   dplyr::mutate(GO=if_else(stringr::str_sub(GO,1,1)==" ",stringr::str_sub(GO,1,-1),GO))


goinput_start<-read.csv("human_mouse_yeast_go.csv")%>%
  dplyr::mutate(GO=if_else(stringr::str_sub(GO,1,1)==" ",stringr::str_sub(GO,2,-1),GO))%>%
  dplyr::filter(GO != "")

#go.df<-read.csv("C:/Harper/Side_analysis/go_parser/human_go_referencefile.csv")
int_df<-read.delim2("BioPlex_293T_Network_10K_Dec_2019.tsv")
bioplex_binary<-read.csv("bioplex_binary_interactions_final.csv")
#C:/Harper/Experiments/Exp0001_CellLines_Chaperone_Choice/BioPlex_293T_Network_10K_Dec_2019.tsv
corum_initial<-read.csv("corum_ref3.csv")
#C:/Harper/Experiments/Exp0003_MLN_Ctrl_Proteome/tmt_proteome_2way/
locale2<- read.csv("MH_organelle list_fromAlban_editedERproteins.csv")

userinput_example <- read.csv("userinput3.csv")

# locale<-locale%>%
#   dplyr::select(Gene.Symbol,Compartment)%>%
#   dplyr::distinct()

locale2<-locale2%>%
  dplyr::select(Gene.Symbol,Compartment)%>%
  dplyr::distinct()


expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  y <- unique(y)
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

ContrastMatrix.list_FIX = function(contrasts, Conditions, labels) {
  num_Conditions = length(Conditions)
  contrast_matrix = matrix(0, nrow = length(contrasts),
                           ncol = num_Conditions)
  for (contrast_id in seq_along(contrasts)) {
    contrast = contrasts[[contrast_id]]
    contrast_vector = rep(0, num_Conditions)
    positive = Conditions %in% contrast[[2]]
    negative = Conditions %in% contrast[[1]]
    contrast_vector[positive] = 1 / sum(positive)
    contrast_vector[negative] = -1 / sum(negative)
    contrast_matrix[contrast_id, ] = contrast_vector
  }

  row.names(contrast_matrix) = labels

  colnames(contrast_matrix) = Conditions
  contrast_matrix
}

anova_test<-function(data){
  df<-data.frame(data)
  test<-aov(Abundance ~ Condition, data= df)
  return(summary(test)[[1]][["Pr(>F)"]][1])
}


# Define UI for application that draws a histogram
ui <- fluidPage(theme=shinytheme("united"),
  
  # Application title
  titlePanel("TMT Protein Abundance Analysis"),
  navbarPage(
    "ShinyApp Multiple Condition Comparisons",
    tabPanel("Data Import",
             sidebarPanel(
               tags$h3("Protein TMT Results Import:"),
               # textInput("uniqueID","Give analysis uniqueID","exp1analysis1"),
               selectInput("txthumanmouse", "Is this human or mouse data?",
                           choices=list("human"="human","mouse"="mouse","yeast"="yeast"),selected = "human"),
               selectInput("inputtype", "PAbeta or MSstats?",
                           choices=list("MSstats"="MSstats","PAbeta"="PAbeta"),selected = "MSstats"),
               fileInput('file1', 'Select TMT stats results file',
                         accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
               fileInput('file2', 'Select the TMT channel design .csv file',
                         accept=c('text/csv','text/comma-separated-values,text/plain','.csv'))

             ),
             mainPanel(
               h1("TMT Data Overview"),
               h3("Following import, explore other panels for analysis")
              
             )
             
    ),
    tabPanel("Correlations",
             sidebarPanel(
               tags$h3("Number of proteins post-filtering:"),
               textOutput("totalproteins"),
               numericInput("mincorrval","Correlation plot minimum value",-1,min=-1),
               downloadButton('CorrDownloadPlot', 'Download Replicate Correlation Plot')
             ),
             mainPanel(
               h1("Replicate Correlation Plot"),
               h4("Plot"),
               plotOutput("plot4")
             )),
    tabPanel("PCA",
             sidebarPanel(
               downloadButton('PCADownloadPlot', 'Download PCA Plot'),
               selectInput("txtXpca", "PCA x-axis component:",
                           choices=list("PC1"="PC1","PC2"="PC2","PC3"="PC3","PC4"="PC4","PC5"="PC5",
                                        "PC6"="PC6","PC7"="PC7","PC8"="PC8","PC9"="PC9","PC10"="PC10"),selected = "PC1"
               ),
               selectInput("txtYpca", "PCA y-axis component:",
                           choices=list("PC1"="PC1","PC2"="PC2","PC3"="PC3","PC4"="PC4","PC5"="PC5",
                                        "PC6"="PC6","PC7"="PC7","PC8"="PC8","PC9"="PC9","PC10"="PC10"),selected = "PC2"
               ),
               selectInput("txtZpca", "Proteins that Drive this PC:",
                                  choices=list("PC1"="PC1","PC2"="PC2","PC3"="PC3","PC4"="PC4","PC5"="PC5",
                                               "PC6"="PC6","PC7"="PC7","PC8"="PC8","PC9"="PC9","PC10"="PC10"),selected = "PC1"
               ),
               plotOutput("plotPCA"),
               plotOutput("plotPCAadd")
             ),
             mainPanel(
               h1("PCA PC Drivers"),
               downloadButton('driversDownload', 'Download All PC Weights'),
               reactableOutput("PCAdrivers")
             )),
    tabPanel("Heatmap/Clustering",
             sidebarPanel(
               tags$h3("Individual Protein Abundance Plots"),
               selectInput("txtScale", "Row Scaling options",
                           choices=list("RowMean"="RowMean","Zscore"="Zscore"),selected = "Zscore"),
               selectInput("anovaorall", "ANOVA Significant or All Proteins for Heatmap",
                           choices=list("ANOVA"="ANOVA","All"="All"),selected = "ANOVA"),
               numericInput("num77","q-value cut-off",0.05,min=0,max=1),
               numericInput("num78","fold change cut-off (|abs|)",1.5,min=0),
               selectInput("txtclustercols", "Cluster columns?",
                           choices=list("Yes"="Yes","No"="No"),selected = "No"),
               numericInput("numclusters","Number of Clusters=",2,min=2),
               reactableOutput("numbercluster"),
               plotOutput("numclusterplot"),
               width=3
             ),
             mainPanel(
               tags$h3("If heatmap does not render, resize the window (will take a few minutes)"),
               downloadButton('plotheatmapprint', 'Download Heatmap'),
               downloadButton('printclustersplot', 'Download Clusters Plot'),
               plotOutput("heatmapplot"),
               plotOutput("plotclustertrace"),
               width=9
             )

    ),
    tabPanel("Results Export",
             mainPanel(
               h1("Protein-level export"),
               h3("Enables table regeneration for alternative cluster designation."),
               textInput("uniqueID","Give analysis uniqueID","exp1analysis1"),
               downloadButton("alldownload", label="Download All Results"),
               # downloadButton("parameter", label="Download Analysis Parameters"),
               reactableOutput("table16")
               #plotlyOutput("plot6"),
               #reactableOutput("brush")
               
             )
             
    )
    ,
    tabPanel("Timecourse",
             sidebarPanel(
               h1("Timecourse analysis: Hotelling T2"),
               h3("Requires temporal data ordered by priority (earliest time point priority 1)"),
               actionButton("runTimecourse",label="Run Timecourse Analysis"),
               # selectInput("txttimecourse", "Is this timecourse data?",
               #             choices=list("Yes"="Yes","No"="No"),selected = "No"),
               reactableOutput("priorityTable5"),
               numericInput("numrepstimecourse","Number of replicates for each time point:",3,min=2,max=4),
               numericInput("numtimecourse1","Time point 0 for Hotelling:",1,min=1,max=8),
               numericInput("numtimecourse2","Comparison time point for Hotelling (based on priority):",2,min=1,max=8),
               textInput("text113","Gene Name (case sensitive)","GAPDH"),
               plotOutput("plottcprot")
               # ,
               # plotOutput("plotclustertrace")
             ),
             mainPanel(
               downloadButton("printhotellingplot", label="Download Hotelling Plot"),
               downloadButton("plottcprotprint2", label="Download Protein Timecourse"),

               #numericInput("alphanum1","Alpha for volcano and correlation plot points",0.2,min=0, max=1),
               plotlyOutput("plothotelling"),
               reactableOutput("hotellingclick")

             )

    ),
    tabPanel("Timecourse Export",
             sidebarPanel(
               h1("Timecourse analysis: Hotelling T2"),
               h3("Requires temporal data ordered by priority (earliest time point priority 1)")
               # ,
               # plotOutput("plotclustertrace")
             ),
             mainPanel(
               downloadButton("hotellingdownload", label="Download Hotelling Results"),
               #numericInput("alphanum1","Alpha for volcano and correlation plot points",0.2,min=0, max=1),
               reactableOutput("hotellingdf2")
               #plotlyOutput("plot6"),
               #reactableOutput("brush")

             )

    ),
    tabPanel("Volcano Plot",
             sidebarPanel(
               h1("Volcano Comparisons across different Conditions: pairwise"),
               h3("Compares two Conditions at a time"),
               reactableOutput("priorityTable2"),
               numericInput("numtimecourse1vol","Baseline Condition for Volcano:",1,min=1),
               numericInput("numtimecourse2vol","Comparison Condition (based on priority):",2,min=1),
               selectInput("ratioratioindicatorRank", "Rank Plot Ratio or Ratio of Ratio (RoR) Indicator:",
                           choices=list("Single Ratio"="Single Ratio", "RoR" = "RoR"),selected = "Single Ratio"),
               numericInput("denomRank1","RoR Denominator Baseline Condition for Rank Plot:",2,min=1),
               numericInput("denomRank2","RoR Denominator Experimental Condition for Rank Plot:",3,min=1),
               selectInput("txtinversion", "Reverese Condition Order (main volcano only):",
                           choices=list("Yes"="Yes","No"="No"),selected = "No"),
               numericInput("num55","Bottom and Top Ranked Proteins by Fold Change",5,min=1,max=40),
               plotOutput("rankfoldchange"),

               # numericInput("num3","q-value cut-off",0.05,min=0,max=1),
               # numericInput("num4","fold change cut-off (|abs|)",1.5,min=0),
               textInput("text11","Gene Name (case sensitive)","GAPDH"),
               plotOutput("plot5"),
               h3("Multiple Condition Correlation Plot inputs:"),
               numericInput("numtimecourse3vol","x-axis baseline Condition (RofR:Numerator):",1,min=1),
               numericInput("numtimecourse4vol","x-axis experimental Condition (RofR:Numerator):",2,min=1),
               numericInput("numtimecourse5vol","y-axis baseline Condition (RofR:Numerator):",1,min=1),
               numericInput("numtimecourse6vol","y-axis experimental Condition (RofR:Numerator):",3,min=1),
               selectInput("ratioratioindicatorVolcano", "Ratio or Ratio of Ratio (RoR) Indicator:",
                           choices=list("Single Ratio"="Single Ratio", "X axis RoR" = "X axis RoR",
                                        "Y axis RoR" = "Y axis RoR","RoR same Denominator" = "RoR same Denominator",
                                        "RoR different Denominator" = "RoR different Denominator"),selected = "Single Ratio"),
               numericInput("denom1","Secondary x-axis/both denom baseline (RofR:Denominator):",2,min=1),
               numericInput("denom2","Secondary x-axis/both denom experimental  (RofR:Denominator):",3,min=1),
               numericInput("denom3","Secondary y-axis baseline Condition (RofR:Denominator):",2,min=1),
               numericInput("denom4","Secondary y-axis experimental Condition(RofR:Denominator):",4,min=1),
               plotlyOutput("volcanoplotnorm"),
               width=5
               # ,
               # plotOutput("plotclustertrace")
             ),
             mainPanel(
               downloadButton('printvolcano', 'Download Volcano Plot'),
               downloadButton('rankfoldchangeprint2', 'Download Top Up/Down FC Plot'),
               downloadButton('proteinBarplotDownloadPlot','Download Single Protein Barplot'),
               downloadButton('printcorrtwocond', 'Download Ratio Correlation Plot'),

               #numericInput("alphanum1","Alpha for volcano and correlation plot points",0.2,min=0, max=1),
               plotlyOutput("volcanoplot"),
               reactableOutput("volcanoclick"),

               width=7
               #plotlyOutput("plot6"),
               #reactableOutput("brush")

             )

    ),

    tabPanel("Volcano (labels)",
             sidebarPanel(
               h1("Labeled Volcano Plot"),
               reactableOutput("priorityTable10"),
               numericInput("numtimecourse7vol","Baseline Condition for Volcano:",1,min=1,max=8),
               numericInput("numtimecourse8vol","Comparison Experimental Condition (based on priority):",2,min=1,max=8),
               numericInput("alphanum2","Alpha for volcano and correlation plot points",0.5,min=0, max=1),
               selectInput("volcanoLabelOptions", "Options Volcano Labeling:",
                           choices=list("CustomN qvalue"="CustomN qvalue","CustomN Log2FC"="CustomN Log2FC","All"="All"),selected = "CustomN Log2FC"),
               numericInput("maxNlabelVolcano","Maximum number of labeled proteins:",25,min=1,max=100)
               # ,
               # plotOutput("plotclustertrace")
             ),
             mainPanel(


               downloadButton('volcanoDownloadPlot', 'Download Volcano Plot (labels)'),
               plotOutput("plotVL")
             )

    ),

    tabPanel("Venn/Upset",
             sidebarPanel(
               h1("Venn Diagram"),
               selectInput("customall", "Custom condition comparisons or All comparisons:",
                           choices=list("Custom"="Custom","All"="All"),selected = "Custom"),
               reactableOutput("priorityTable1"),
               textInput("txtlistsig","Custom Setting based on Priority list: 1,2; 1,3; 2,3"),
               selectInput("txtupdownoverlap", "Which regulated proteins for Venn:",
                           choices=list("up"="up","down"="down","both"="both"),selected = "both"),
               downloadButton('printploteuler', 'Download Venn Diagram'),
               plotOutput("ploteuler")
             ),
             mainPanel(
               h1("Upset Plot"),
               downloadButton('printplotupset', 'Download Upset Plot'),
               plotOutput("plotupset"), 
               downloadButton("venndiagramtableexport2", label="Download Venn Diagram Table"),
               reactableOutput("venndiagramtable")
             )

    ),
    tabPanel("Multiple Condition Comparisons",
             sidebarPanel(
               reactableOutput("priorityTable100"),
               # selectInput("performlineartest", "Start performing linear model multi-comparison test?",
               #             choices=list("Do Stats"="Do Stats","Hold"="Hold"),selected = "Hold"),
               h3("Ensure parameters and contrast matrix below are correct before running analysis"),
               actionButton("runLinearModel",label="Run Linear Model Analysis"),
               
               # ,
               textInput("variableTest1","Variable 1 Name", "var1"),
               textInput("variableTest2","Variable 2 Name","var2"),
               textInput("multiplecomparisontest","Generate Priority Comparison List: 1,2; 3,4"),
               selectInput("betachoice", "Volcano Beta to Plot",
                           choices=list("Var1"="Var1","Var2"="Var2","Var1Var2 interaction" = "Var1Var2 interaction"),selected = "Var1Var2 interaction"),
               numericInput("qcutoffLinear","q-value cut-off",0.05,min=0,max=1),
               numericInput("fccutoffLinear","fold change cut-off (|abs|)",1.5,min=0),
               
               reactableOutput("tableContrastJoiner"),
               textInput("linearIndivid","Gene Name (case sensitive)","GAPDH"),
               plotOutput("linearModelPlotIndividual")
             ),
             mainPanel(
               h2("Analysis takes several minutes."),
               h3("Volcano Plot based on each beta term (var1, var2, or var1var2 interaction term):"),
               downloadButton("LMStatsDownloadTable", label="Download All Linear Stats Results"),
               downloadButton("DownloadVolcanoLinear", label="Download Volcano PDF"),
               plotlyOutput("volcanoplotLinear"),
               reactableOutput("LinearStatsTable")
               
             )
             
    ),

    tabPanel("BioPlex Top",
             sidebarPanel(
               h1("BioPlex Interactome Top Hits"),
               reactableOutput("priorityTable3"),
               numericInput("numconbioplex3","Baseline Condition:",1,min=1,max=8),
               numericInput("numconbioplex4","Comparison condtion (based on priority):",2,min=1,max=8),
               numericInput("num10001","minimium number of interactors",4,min=2),
               downloadButton("alldownloadbioplex", label="Download Bioplex Results"),
               reactableOutput("tableBinaryBioplex"),
               width = 5
             ),
             mainPanel(

               numericInput("numbioplex1","TopN most up regulated interactomes:",10,min=1),
               downloadButton("printbiolplextop", label="Download Top Up Bioplex Plot"),
               plotOutput("bioplextopN"),
               numericInput("numbioplex2","TopN most down regulated interactomes:",10,min=1),
               downloadButton("printbiolplexbottom", label="Download Top Down Bioplex Plot"),


               plotOutput("bioplexbottomN"),
               width = 7

             )

    ),
    tabPanel("BioPlex User",
             sidebarPanel(
               h1("BioPlex Interactome User Defined"),
               reactableOutput("priorityTable4"),
               textInput("txtbioplex","Gene Name for Bioplex","GAPDH"),
               numericInput("numconbioplex1","Baseline Condition:",1,min=1,max=8),
               numericInput("numconbioplex2","Comparison condtion (based on priority):",2,min=1,max=8),
               downloadButton("printbioplexExp", label="Download Single Condition Bioplex Plot"),
               plotOutput("bioplexExp"),
               downloadButton("printbioplexExpAll", label="Download Multiple Condition Bioplex Plot"),
               plotOutput("bioplexExpAll"),
               reactableOutput("tablebioplex"),
               h3("Protein Interactome Level: (baseline:Condition above are different than last tab)"),
               numericInput("num10002","minimium number of interactors",4,min=2),
               reactableOutput("tableBinaryBioplex2"),
               width = 6
             ),
             mainPanel(
               # downloadButton("printbioplexNetwork", label="Download Network Plot"),
               plotOutput("bioplexNetwork"),
               width = 6

             )

    ),
    tabPanel("Corum",
             sidebarPanel(
               h1("Corum Complex Analysis"),
               downloadButton("alldownloadCorum",label = "Download All Corum Results"),
               h3("Complex Selection Parameter"),
               reactableOutput("priorityTable6"),
               numericInput("numcon1","Baseline Condition:",1,min=1,max=8),
               numericInput("numcon2","Comparison condtion (based on priority):",2,min=1,max=8),
               numericInput("numcomplexmin","Minimum number of identified complex members",4,min=2),
               numericInput("complexranknum","Complex Rank Number for Selection",1,min=1),
               h3("Complex-level Table"),
               reactableOutput("tablecorum"),
               h3("Top Regulated Complexes:"),
               downloadButton("printcorumExp", label="Download Top Up Corum Plot"),
               numericInput("numcorum","TopN most up regulated complexes:",10,min=1),
               plotOutput("corumExp"),
               downloadButton("printcorumExp2", label="Download Top Down Corum Plot"),
               numericInput("numcorum2","TopN most down regulated complexes:",10,min=1),
               plotOutput("corumExp2"),
               width = 6
             ),
             mainPanel(
               h3("Individual Complex Selection"),
               plotOutput("corumNetwork"),
               downloadButton("printcorumExp3", label="Download Protein (Single Comparison) Corum Plot"),
               plotOutput("corumExp3"),
               reactableOutput("tablecorumOneComplex"),
               downloadButton("printcorumExpallcond", label="Download Protein (Mulit-comparison) Corum Plot"),
               plotOutput("corumExpallcond"),
               width = 6

             )


    ),
    tabPanel("GO Clusters",
             sidebarPanel(


               numericInput("numcluster","cluster for GO enrichment",1,min=1),
               actionButton("runGOcluster",label="Run GO enrichment (Cluster)"),
               numericInput("numgo","q-value cut-off",0.05,min=0,max=1),
               numericInput("numtopn","topN number of GO terms (q-value ranked):",10,min=1,max=100),
               #downloadButton('plotGOall', 'Download GO scatter'),
               #downloadButton('TopGODownloadPlot', 'Download Top GO'),
               downloadButton("GOclusterdownload", label="Download Cluster GO Results"),
               downloadButton("printplotGO", label="Download All GO Plot (Cluster)"),
               plotlyOutput("plotGO")
             ),
             mainPanel(
               h3("Click on point to explore Proteins in GO term:"),
               downloadButton("printplotGOtopn2", label="Download TopN GO Plot (Cluster)"),
               plotlyOutput("plotGOtopn2"),
               reactableOutput("goclick")

             )

    ),
    tabPanel("GO Clusters 2",
             sidebarPanel(
               h3("Significance designations derived from user-defined cut-offs from the GO Clusters enrichment tab"),
               h3("GO terms with associated proteins:"),
               # downloadButton("alldownloadGO", label="Download All Results"),
               reactableOutput("tableGO"),
               width=5
             ),
             mainPanel(
               textInput("txtlistnum","List of GO ranks separated by commas ex: 1,14,24,16"),
               numericInput("numgo2","q-value cut-off",0.05,min=0,max=1),
               downloadButton("printplotGOuser", label="Download User-defined GO Plot (Cluster)"),
               plotOutput("plotGOuser"),
               h3("GO terms Ranked:"),
               reactableOutput("tableGO2"),
               width=7

             )


    ),
    tabPanel("GO Regulated",
             sidebarPanel(
               h3("Signficance designations from statistics section"),
               selectInput("txtupdown", "Which regulated proteins for GO:",
                           choices=list("up"="up","down"="down","both"="both"),selected = "up"),
               reactableOutput("priorityTable7"),
               numericInput("numpriority1","Baseline Condition for Volcano:",1,min=1,max=8),
               numericInput("numpriority2","Comparison condtion (based on priority):",2,min=1,max=8),
               actionButton("runGOreg",label="Run GO enrichment (Regulated)"),
               numericInput("numgo1","q-value cut-off",0.05,min=0,max=1),
               numericInput("numtopn1","topN number of GO terms (q-value ranked):",10,min=1,max=100),
               downloadButton("GOregulateddownload", label="Download Regulated GO Results"),
               downloadButton("printplotGO2", label="Download All GO Plot (Regulated)"),
               plotlyOutput("plotGO2")
             ),
             mainPanel(
               h3("Click on point to explore Proteins in GO term:"),
               downloadButton("printplotGOtopn3", label="Download TopN GO Plot (Regulated)"),
               plotlyOutput("plotGOtopn3"),
               reactableOutput("goclick3")

             )

    ),
    tabPanel("GO Regulated 2",
             sidebarPanel(
               h3("Significance designations derived from user-defined cut-offs from the GO Regulated enrichment tab"),
               h3("GO terms with associated proteins:"),
               # downloadButton("alldownloadGO", label="Download All Results"),
               reactableOutput("tableGOreg"),
               width=5
             ),
             mainPanel(
               textInput("txtlistnum3","List of GO ranks separated by commas ex: 1,14,24,16"),
               numericInput("numgo3","q-value cut-off",0.05,min=0,max=1),
               downloadButton("printplotGOuserReg", label="Download User-defined GO Plot (Regulated)"),
               plotOutput("plotGOuserReg"),
               h3("GO terms Ranked:"),
               reactableOutput("tableGO2reg"),
               width=7

             )


    ),

    tabPanel("Localization v1",
             sidebarPanel(
               tags$h3("Localization Volcano"),
               selectInput("localization", "Localization for Volcano Plot (yeast vacuole = lysosome):",
                           choices=list("Actin binding proteins"="Actin binding proteins", "Cytoplasm"="Cytoplasm",
                                        "Endosome"="Endosome",   "ER"= "ER", "ER_high_curvature"= "ER_high_curvature",
                                        "ER-Lumen"="ER-Lumen", "ER-Membrane"="ER-Membrane", "Ergic/cisGolgi"="Ergic/cisGolgi",
                                        "Golgi"="Golgi" ,"Large Protein Complex"="Large Protein Complex", "Lysosome"="Lysosome",
                                        "Mitochondria"="Mitochondria", "Mitochondrion-IM"="Mitochondrion-IM",
                                        "Mitochondrion-Mtx"="Mitochondrion-Mtx", "Mitochondrion-OM"="Mitochondrion-OM",
                                        "Nuclear pore complex"="Nuclear pore complex", "Nucleus"="Nucleus"  ,
                                        "Peroxisome"="Peroxisome", "Plasma membrane"="Plasma membrane","Proteasome"="Proteasome",
                                        "Ribosome"="Ribosome",
                                        "Golgi.apparatus..KH."="Golgi.apparatus..KH." ,"Golgi.membrane..KH."="Golgi.membrane..KH.",
                                        "Golgi.membrane.associated..KH."="Golgi.membrane.associated..KH.",
                                        "Cis.Golgi..KH."="Cis.Golgi..KH.","Trans.Golgi..KH."="Trans.Golgi..KH."),selected = "Lysosome"),
               reactableOutput("priorityTable8"),
               numericInput("numtimecourse1loc","Baseline Condition (based on priority):",1,min=1,max=8),
               numericInput("numtimecourse2loc","Comparison Condition (based on priority):",2,min=1,max=8),
               multiInput("locales", "Fold change across localizations (Select all that apply):",
                           choices=list("Actin binding proteins"="Actin binding proteins", "Cytoplasm"="Cytoplasm",
                                        "Endosome"="Endosome",   "ER"= "ER", "ER_high_curvature"= "ER_high_curvature",
                                        "ER-Lumen"="ER-Lumen", "ER-Membrane"="ER-Membrane", "Ergic/cisGolgi"="Ergic/cisGolgi",
                                        "Golgi"="Golgi" ,"Large Protein Complex"="Large Protein Complex", "Lysosome"="Lysosome",
                                        "Mitochondria"="Mitochondria", "Mitochondrion-IM"="Mitochondrion-IM",
                                        "Mitochondrion-Mtx"="Mitochondrion-Mtx", "Mitochondrion-OM"="Mitochondrion-OM",
                                        "Nuclear pore complex"="Nuclear pore complex", "Nucleus"="Nucleus"  ,
                                        "Peroxisome"="Peroxisome", "Plasma membrane"="Plasma membrane","Proteasome"="Proteasome",
                                        "Ribosome"="Ribosome",
                                        "Golgi.apparatus..KH."="Golgi.apparatus..KH." ,"Golgi.membrane..KH."="Golgi.membrane..KH.",
                                        "Golgi.membrane.associated..KH."="Golgi.membrane.associated..KH.",
                                        "Cis.Golgi..KH."="Cis.Golgi..KH.","Trans.Golgi..KH."="Trans.Golgi..KH."),selected = "Lysosome"),
               selectInput("txtinversionloc", "Reverese Condition Order (main volcano only):",
                           choices=list("Yes"="Yes","No"="No"),selected = "No"),
               # ,
               downloadButton('downloadlocalizationdata', 'Download Localization Data'),
               downloadButton('printcorrcolorcompartment', 'Download Correlation Localization'),
               plotOutput("corrcolorcompartment"),
               numericInput("numtimecourse1loc2","x-axis baseline Condition (RofR:Numerator):",1,min=1),
               numericInput("numtimecourse2loc2","x-axis experimental Condition (RofR:Numerator):",2,min=1),
               numericInput("numtimecourse3loc2","y-axis baseline Condition (RofR:Numerator):",1,min=1),
               numericInput("numtimecourse4loc2","y-axis experimental Condition (RofR:Numerator):",3,min=1),
               selectInput("ratioratioindicator", "Ratio or Ratio of Ratio (RoR) Indicator:",
                           choices=list("Single Ratio"="Single Ratio", "X axis RoR" = "X axis RoR",
                                        "Y axis RoR" = "Y axis RoR","RoR same Denominator" = "RoR same Denominator",
                                        "RoR different Denominator" = "RoR different Denominator"),selected = "Single Ratio"),
               numericInput("numtimecourse5loc2","Secondary x-axis/both denom baseline (RofR:Denominator):",2,min=1),
               numericInput("numtimecourse6loc2","Secondary x-axis/both denom experimental  (RofR:Denominator):",3,min=1),
               numericInput("numtimecourse7loc2","Secondary y-axis baseline Condition (RofR:Denominator):",2,min=1),
               numericInput("numtimecourse8loc2","Secondary y-axis experimental Condition(RofR:Denominator):",4,min=1),
               downloadButton('printcorrcolorcompartmentmulti', 'Download Correlation Localization Ratio'),
               plotOutput("corrcolorcompartmentmulti")
             ),
             mainPanel(
               downloadButton('printplot7', 'Download Volcano Localization'),
               plotOutput("plot7"),
               numericInput("numlocale1","x-axis minimium",-2),
               numericInput("numlocale2","x-axis maximium",2),
               downloadButton('printplot8', 'Download Single Condition Foldchange for Localization'),
               plotOutput("plot8"),
               downloadButton('printplotmultilocale', 'Download Multi-Condition Foldchange for Localization'),
               plotOutput("plotmultilocale"),
               downloadButton('printpointPlot', 'Localization Correlation Plot'),
               plotOutput("pointPlot")
             )

    ),
    tabPanel("Localization (labels) v1",
             sidebarPanel(
               tags$h3("Localization Volcano with Labels"),
               reactableOutput("priorityTable9"),
               numericInput("numtimecourse1loc3","Baseline Condition (based on priority):",1,min=1,max=8),
               numericInput("numtimecourse2loc3","Comparison Condition (based on priority):",2,min=1,max=8),
               selectInput("localization2", "Localization for Volcano Plot (yeast vacuole = lysosome):",
                           choices=list("Actin binding proteins"="Actin binding proteins", "Cytoplasm"="Cytoplasm",
                                        "Endosome"="Endosome",   "ER"= "ER", "ER_high_curvature"= "ER_high_curvature",
                                        "ER-Lumen"="ER-Lumen", "ER-Membrane"="ER-Membrane", "Ergic/cisGolgi"="Ergic/cisGolgi",
                                        "Golgi"="Golgi" ,"Large Protein Complex"="Large Protein Complex", "Lysosome"="Lysosome",
                                        "Mitochondria"="Mitochondria", "Mitochondrion-IM"="Mitochondrion-IM",
                                        "Mitochondrion-Mtx"="Mitochondrion-Mtx", "Mitochondrion-OM"="Mitochondrion-OM",
                                        "Nuclear pore complex"="Nuclear pore complex", "Nucleus"="Nucleus"  ,
                                        "Peroxisome"="Peroxisome", "Plasma membrane"="Plasma membrane","Proteasome"="Proteasome",
                                        "Ribosome"="Ribosome",
                                        "Golgi.apparatus..KH."="Golgi.apparatus..KH." ,"Golgi.membrane..KH."="Golgi.membrane..KH.",
                                        "Golgi.membrane.associated..KH."="Golgi.membrane.associated..KH.",
                                        "Cis.Golgi..KH."="Cis.Golgi..KH.","Trans.Golgi..KH."="Trans.Golgi..KH."),selected = "Lysosome"),
               selectInput("volcanoLabelOptions2", "Options Volcano Labeling:",
                           choices=list("CustomN qvalue"="CustomN qvalue","CustomN Log2FC"="CustomN Log2FC","All"="All"),selected = "CustomN Log2FC"),
               numericInput("maxNlabelVolcano2","Maximum number of labeled proteins:",10,min=1,max=100)
             ),
             mainPanel(
               downloadButton('printplot44', 'Download Volcano Localization Labels'),
               plotOutput("plot44")
             )

    ),
    
    
  tabPanel("User Select",
             sidebarPanel(
               tags$h3("Select proteins to plot"),
               fileInput('file4', 'Select Proteins of Interest Dataframe',
                         accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
               downloadButton('exampleUserInput', 'Download Example User Annotation File'),
               numericInput("numconduser1","Select baseline condition for Volcano (by priority):",1,min=1,max=8),
               numericInput("numconduser2","Comparison condition (by priority):",2,min=1,max=8),
               selectInput("txtinversionuser", "Reverese Condition Order (for volcano):",
                           choices=list("Yes"="Yes","No"="No"),selected = "No"),
               selectInput("txtlabeluser", "Label User-selected Proteins:",
                           choices=list("Yes"="Yes","No"="No"),selected = "Yes"),
               numericInput("alpha5","Alpha background points",0.2,min=0,max=1),
               numericInput("alpha6","Alpha user-defined points",0.9,min=0,max=1),
               numericInput("size1","Background point size",1,min=0,max=20),
               numericInput("size2","User-defined point size",3,min=0,max=20),
               selectInput("txtclustercolsUser", "Cluster columns and or rows?",
                           choices=list("Both"="Both","Rows Only"="Rows Only", "Columns Only"="Columns Only", "Neither"="Neither"), selected = "Rows Only"),
               selectInput("heatmapLabelRow", "Label Heatmap Rows?",
                           choices=list("Yes"="Yes","No"="No"),selected = "No")
             ),
             mainPanel(
               downloadButton('printplotuser', 'Download User-defined Volcano'),
               plotOutput("plotuser"),
               downloadButton('printcategviolin', 'Download User-defined Violin Plot'),
               plotOutput("categviolin"),
               # plotOutput("corruserplot"),
               
               downloadButton('printheatmapUser', 'Download User-defined Violin Plot'),
               plotOutput("heatmapUser")
             )

    )

  )

)

# Define server logic required to draw a histogram
options(shiny.maxRequestSize=1000*1024^2)
server <- function(input, output, session) {
  
  
  go.df <- reactive({
    go_input <- goinput_start%>%
      dplyr::filter(organism == input$txthumanmouse)%>%
      dplyr::select(-organism)
    return(go_input)
    # print(go_input[1:50,])
    # print(unique(go_input$organism))
    # print(unique(go_input$GO_cat))
      
  })
  
  corum <- reactive({

    corum_input <- corum_initial%>%
      dplyr::filter(Organism == input$txthumanmouse)%>%
      dplyr::select(-Organism)
    return(corum_input)
    
  })
  
  locale <- reactive({
    locale33<-locale_input %>%
      dplyr::select(Gene.Symbol,Compartment,Organism)%>%
      dplyr::filter(Organism == input$txthumanmouse)%>%
      dplyr::ungroup()%>%
      dplyr::select(-Organism)%>%
      dplyr::distinct()
    return(locale33)
  })
    
    
  
  df <- reactive({
    req(input$file1)
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    tbl <- read.csv(inFile$datapath, check.names=FALSE)
    
    return(tbl)
  })

  exp_design <- reactive({
    req(input$file2)
    inFile1 <- input$file2
    if (is.null(inFile1))
      return(NULL)
    tbl2 <- read.csv(inFile1$datapath)
    return(tbl2)
  })

  # exp_design_frac <- reactive({
  #   req(input$file3)
  #   inFile1 <- input$file3
  #   if (is.null(inFile1))
  #     return(NULL)
  #   tbl2 <- read.csv(inFile1$datapath)
  #   return(tbl2)
  # })
  
  # annotation_file <- reactive({
  #   annotation_file<-dplyr::full_join(exp_design(),exp_design_frac(),by="Mixture")%>%
  #     dplyr::select(Run,Fraction,TechRepMixture,Channel,Condition, Replicate,Mixture,priority)%>%
  #     dplyr::mutate(Replicate = paste(Condition, Replicate,sep="_"))%>%
  #     dplyr::rename(BioReplicate = Replicate)
  #   return(annotation_file)
  # })

  
  ### start to remove this code and re-create from the output files
  # exp_design_factors<- reactive({
  #   factor_order <- exp_design()%>%
  #     mutate(BioReplicate = paste(Condition,Replicate,sep="_"))%>%
  #     select(Condition,BioReplicate,priority)%>%
  #     distinct()
  #   exp_design_factors$BioReplicate<-reorder(exp_design_factors$BioReplicate,exp_design_factors$priority,min)
  #   exp_design_factors$Condition<-reorder(exp_design_factors$Condition,exp_design_factors$priority,min)
  #   return(exp_design_factors)
  # })
  
  norm_sum_final2 <-reactive({
    if (input$inputtype == "MSstats") {
      factor_order <- exp_design()%>%
        dplyr::filter(Condition %!in% c("Empty","empty"))%>%
        dplyr::mutate(BioReplicate = paste(Condition,Replicate,sep="_"))%>%
        dplyr::select(Condition,BioReplicate,priority)%>%
        dplyr::distinct()%>%
        dplyr::arrange(priority,Condition,BioReplicate )%>%
        dplyr::mutate(Condition = factor(Condition, unique(Condition)),
                      BioReplicate = factor(BioReplicate, unique(BioReplicate)))
      
        
    }
    if (input$inputtype == "PAbeta") {
      factor_order <- exp_design()%>%
        dplyr::rename(Condition = condition,
                      Replicate = replicate)%>%
        dplyr::filter(Condition %!in% c("Empty","empty"))%>%
        dplyr::mutate(BioReplicate = paste(Condition,Replicate,sep="_"))%>%
        dplyr::select(Condition,BioReplicate,priority)%>%
        dplyr::distinct()%>%
        dplyr::arrange(priority,Condition,BioReplicate )%>%
        dplyr::mutate(Condition = factor(Condition, unique(Condition)),
                      BioReplicate = factor(BioReplicate, unique(BioReplicate)))
    }
    
    
    # factor_order$BioReplicate<-reorder(factor_order$BioReplicate,factor_order$priority,min)
    # factor_order$Condition<-reorder(factor_order$Condition,factor_order$priority,min)
    # 
    # factor_order$Condition<-forcats::fct_reorder(factor_order$Condition, factor_order$priority, min)
    # factor_order$BioReplicate<-forcats::fct_reorder(factor_order$BioReplicate, factor_order$priority, min)
    return(factor_order)
  })
  
  multipleSampleLinearDesignation <- reactive({
    req(input$multiplecomparisontest)
    factor_order <- norm_sum_final2()
    factor_order$BioReplicate<-reorder(factor_order$BioReplicate,factor_order$priority,min)
    factor_order$Condition<-reorder(factor_order$Condition,factor_order$priority,min)

    # string_val <- "1,2 ; 3,4"
    value<-stringr::str_remove_all(input$multiplecomparisontest," ")
    
    list_index_ranks<-data.frame("contrasts"=unlist(strsplit(value, "\\;")))%>%
      tidyr::separate(contrasts, into=c("first","second"),sep="\\,")%>%
      dplyr::mutate(first=as.numeric(as.character(first)),
                    second=as.numeric(as.character(second)))
    ###variable 1 generation
    list_index_ranks_1<-list_index_ranks%>%
      dplyr::select(first)%>%
      dplyr::rename(priority= first)%>%
      dplyr::mutate(var1=0)
    
    list_index_ranks_2<-list_index_ranks%>%
      dplyr::select(second)%>%
      dplyr::rename(priority= second)%>%
      dplyr::mutate(var1=1)%>%
      dplyr::bind_rows(list_index_ranks_1,.)
    
    var2_df<-data.frame(t(list_index_ranks))
    
    colnames(var2_df)<-c("first","second")
    
    list_index_ranks_3<-var2_df%>%
      dplyr::select(first)%>%
      dplyr::rename(priority= first)%>%
      dplyr::mutate(var2=0)
    
    list_index_ranks_4<-var2_df%>%
      dplyr::select(second)%>%
      dplyr::rename(priority= second)%>%
      dplyr::mutate(var2=1)%>%
      dplyr::bind_rows(list_index_ranks_3,.)
    
    rownames(list_index_ranks_4) <- 1:4
    
    contrast_joiner<-dplyr::left_join(list_index_ranks_2,list_index_ranks_4,by="priority")%>%
      dplyr::mutate(var1= as.numeric(as.character(var1)),
                    var2 = as.numeric(as.character(var2)))%>%
      dplyr::right_join(factor_order%>%dplyr::mutate(priority = as.numeric(as.character(priority))),.,by="priority")
    # names(contrast_joiner) <- gsub(x = names(contrast_joiner), pattern = "var1", replacement = "#") 
    return(contrast_joiner)
  })
  
  output$tableContrastJoiner<-renderReactable({
    newDF <- multipleSampleLinearDesignation()
    names(newDF) <- gsub(x = names(newDF), pattern = "var1", replacement = input$variableTest1) 
    names(newDF) <- gsub(x = names(newDF), pattern = "var2", replacement = input$variableTest2) 
    
    reactable(newDF, filterable = TRUE, defaultPageSize  = dim(multipleSampleLinearDesignation())[[1]])
  })
  
  
  
  
  ### Start linear model test
  observeEvent(input$runLinearModel, {
    linearmodelDF <- reactive({
      req(input$multiplecomparisontest)
      # if (input$performlineartest == "Do Stats") {
      prot_norm_data_reps_tidy_stat <- prot_biorep_final()%>%
        dplyr::group_by(ProtID)%>%
        tidyr::gather("BioReplicate", "Abundance", 2:dim(prot_biorep_final())[[2]])%>%
        dplyr::ungroup()

      stat_df<-dplyr::inner_join(prot_norm_data_reps_tidy_stat,multipleSampleLinearDesignation(),by="BioReplicate" )
      
      stat_df$BioReplicate<-reorder(stat_df$BioReplicate,stat_df$priority,min)
      stat_df$Condition<-reorder(stat_df$Condition,stat_df$priority,min)
      
      
      ## nest data for protein-level linear model fit
      nest_df<-stat_df%>%
        dplyr::group_by(ProtID)%>%
        tidyr::nest()
      
      ## linear model as function of the classifiers defined above
      linear_fit_fun <- function(data){
        input<-data.frame(data)
        model<-lm(data=input,Abundance ~   var1*var2)
        #summary(model)
        ## generate dataframe that appends the linear model beta estimates and p-values
        lm_results<-data.frame(coeff = summary(model)$coefficients[1:4,1],p.value = summary(model)$coefficients[1:4,4])%>%
          tibble::rownames_to_column("Contrast")%>%
          dplyr::mutate(Contrast = as.character(Contrast))
        return(lm_results) 
      }
      
      #perform linear model using the above function
      nest_df2<-nest_df%>%
        dplyr::mutate(lm_results = data %>% map(linear_fit_fun))
      
      ## pivot data to untidy data type
      convert_long<-function(lm_results){
        df<-data.frame(lm_results)
        df_all<-df%>%
          tidyr::pivot_wider( 
            names_from = Contrast, 
            values_from = c(coeff,p.value ))
        return(df_all)
      }
      
      ### convert linear results to untidy format which we unnest
      nest_df3<-nest_df2%>%
        dplyr::mutate(results= lm_results %>% map(convert_long))%>%
        tidyr::unnest(results)%>%
        dplyr::select(-data, -lm_results)%>%
        dplyr::distinct()
      
      
      nest_df3$q.value_var1 <- p.adjust(nest_df3$p.value_var1,method="fdr")
      nest_df3$q.value_var2 <- p.adjust(nest_df3$p.value_var2,method="fdr")
      
      
      colnames(nest_df3)[5]<-"coeff_var1_var2"
      colnames(nest_df3)[9]<-"p.value_var1_var2"
      
      nest_df3$q.value_var1_var2 <- p.adjust(nest_df3$p.value_var1_var2,method="fdr")
      
      nest_df4<-nest_df3%>%
        tidyr::separate("ProtID", into = c("Reference","Gene.Symbol","Annotation"),sep="__X__")
      # }
      # if (input$performlineartest == "Hold") {
      #   nest_df4 <- data.frame("NotRun"=c(0))
      # }
      return(nest_df4)
    })
    
    output$linearModelPlotIndividual <- renderPlot({
      req(input$multiplecomparisontest)
      # if (input$performlineartest == "Do Stats") {
      prot_norm_data_reps_tidy_stat <- prot_biorep_final()%>%
        dplyr::group_by(ProtID)%>%
        tidyr::gather("BioReplicate", "Abundance", 2:dim(prot_biorep_final())[[2]])%>%
        dplyr::ungroup()
      
      stat_df<-dplyr::inner_join(prot_norm_data_reps_tidy_stat,multipleSampleLinearDesignation(),by="BioReplicate" )
      
      stat_df$BioReplicate<-reorder(stat_df$BioReplicate,stat_df$priority,min)
      stat_df$Condition<-reorder(stat_df$Condition,stat_df$priority,min)
      
      
      stat_df_prot<-stat_df%>%
        tidyr::separate("ProtID", into = c("Reference","Gene.Symbol","Annotation"),sep="__X__")%>%
        tidyr::separate("Reference",into=c("type","Reference","extra"),sep="\\|")%>%
        dplyr::mutate(uniqueID = paste(Gene.Symbol,Reference,sep="_"))%>%
        dplyr::filter(Gene.Symbol == input$linearIndivid)%>%
        dplyr::mutate(var3 = dplyr::if_else(var2 == 0, paste("- ",input$variableTest2), paste("+ ",input$variableTest2)))
      
      stat_df_prot$var3 <- reorder(stat_df_prot$var3,stat_df_prot$var2,min)
      
      ggplot(stat_df_prot,aes(var1, color=var3, y=Abundance))+
        geom_point(size=3,alpha=0.5)+
        theme_classic()+
        scale_color_viridis_d(end=0.8)+facet_wrap(.~uniqueID,scales = "free")+
        xlim(c(-0.25,1.25))+geom_smooth(se=F,method="lm")+
        scale_x_continuous(labels= c(paste("- ",input$variableTest1),paste("+ ",input$variableTest1)), breaks=0:1, limits=c(-0.5,1.5))+
        # scale_color_manual(labels = c(paste("- ",input$variableTest2), paste("+ ",input$variableTest2)))+
        labs(color=paste(input$variableTest2), x=paste(input$variableTest1),y="Log2(MS Intensity)")+
        theme(axis.title = element_text(size = 18),
              axis.text = element_text(size=16),
              strip.text.x = element_text(size = 16),
              legend.text=element_text(size=14),
              legend.title=element_text(size=16))
      # }
    })
    
    linearmodelDFSigCall <- reactive({
      linearmodelDFSigCall<-linearmodelDF()%>%
        dplyr::mutate(sig_var1=  dplyr::if_else(q.value_var1 < input$qcutoffLinear , dplyr::if_else(abs(coeff_var1)>log2(input$fccutoffLinear), dplyr::if_else(coeff_var1>0,"up","down"),"n.s."),"n.s."),
                      sig_var2=  dplyr::if_else(q.value_var2 < input$qcutoffLinear , dplyr::if_else(abs(coeff_var2)>log2(input$fccutoffLinear), dplyr::if_else(coeff_var2>0,"up","down"),"n.s."),"n.s."),
                      sig_var1_var2=  dplyr::if_else(q.value_var1_var2 < input$qcutoffLinear , dplyr::if_else(abs(coeff_var1_var2)>log2(input$fccutoffLinear), dplyr::if_else(coeff_var1_var2>0,"up","down"),"n.s."),"n.s."))%>%
        # tidyr::separate("ProtID", into = c("Reference","Gene.Symbol","Annotation"),sep="__X__")%>%
        tidyr::separate("Reference",into=c("type","Reference","GeneNameExtra"),sep="\\|")%>%
        dplyr::mutate(uniqueID = paste(Gene.Symbol,Reference,sep="_"))
      return(linearmodelDFSigCall)
    })
    
    output$LinearStatsTable<-renderReactable({
      d <- event_data("plotly_selected")
      req(d)
      newDF<-linearmodelDFSigCall()%>%
        dplyr::filter(uniqueID %in% d$key )%>%
        dplyr::select(-GeneNameExtra,-type)%>%
        dplyr::ungroup()
      names(newDF) <- gsub(x = names(newDF), pattern = "var1", replacement = input$variableTest1) 
      names(newDF) <- gsub(x = names(newDF), pattern = "var2", replacement = input$variableTest2) 
      reactable(newDF, filterable = TRUE)
    })
    
    
    LinearModelStatsDownloadTable<-reactive({
      # req(input$performlineartest == "Do Stats")
      newDF<-linearmodelDFSigCall()%>%
        dplyr::ungroup()
      names(newDF) <- gsub(x = names(newDF), pattern = "var1", replacement = input$variableTest1)
      names(newDF) <- gsub(x = names(newDF), pattern = "var2", replacement = input$variableTest2)
      return(newDF)
    })
    
    output$LMStatsDownloadTable <- downloadHandler(
      filename = function() {
        paste("LinearModel_Stats",".csv", sep = "")
      },
      content = function(file) {
        write.csv(LinearModelStatsDownloadTable(), file, row.names = FALSE)
      }
    )
    
    
    
    output$volcanoplotLinear <- renderPlotly({
      # req(input$performlineartest == "Do Stats")
      if (input$betachoice == "Var1Var2 interaction" ) {
        volcano<-ggplot(linearmodelDFSigCall(), aes( coeff_var1_var2, -log10(q.value_var1_var2),color=sig_var1_var2,label=Gene.Symbol,label1=coeff_var1, label2=coeff_var2, key=uniqueID))+
          geom_point(alpha=0.75,size=0.75)+theme_classic()+
          scale_color_viridis_d(end=0.8)+
          geom_hline(yintercept = -log10(input$qcutoffLinear),linetype="dashed")+
          geom_vline(xintercept=log2(input$fccutoffLinear),linetype="dashed")+
          geom_vline(xintercept=-log2(input$fccutoffLinear),linetype="dashed")+
          labs(x= paste("Beta Coef. ",input$variableTest1,"_",input$variableTest2," Interaction Term",sep=""),
               y= "-Log10(q-value)",
               color=paste("sig_",input$variableTest1,"_",input$variableTest2,sep=""))+
          theme(axis.title = element_text(size = 18),
                axis.text = element_text(size=16),
                legend.text=element_text(size=14),
                legend.title=element_text(size=16))
          
      }
      if (input$betachoice == "Var1" ) {
        volcano<-ggplot(linearmodelDFSigCall(), aes( coeff_var1, -log10(q.value_var1),color=sig_var1,label1=coeff_var2, label2=coeff_var1_var2,label=Gene.Symbol, key=uniqueID))+
          geom_point(alpha=0.75,size=0.75)+theme_classic()+
          scale_color_viridis_d(end=0.8)+
          geom_hline(yintercept = -log10(input$qcutoffLinear),linetype="dashed")+
          geom_vline(xintercept=log2(input$fccutoffLinear),linetype="dashed")+
          geom_vline(xintercept=-log2(input$fccutoffLinear),linetype="dashed")+
          labs(x= paste("Beta Coef. ",input$variableTest1,sep=""), y= "-Log10(q-value)",
               color=paste("sig_",input$variableTest1,sep=""))+
          theme(axis.title = element_text(size = 18),
                axis.text = element_text(size=16),
                legend.text=element_text(size=14),
                legend.title=element_text(size=16))
      }
      if (input$betachoice == "Var2" ) {
        volcano<-ggplot(linearmodelDFSigCall(), aes( coeff_var2, -log10(q.value_var2),color=sig_var2,label=Gene.Symbol,label1=coeff_var1, label2=coeff_var1_var2, key=uniqueID))+
          geom_point(alpha=0.75,size=0.75)+theme_classic()+
          scale_color_viridis_d(end=0.8)+
          geom_hline(yintercept = -log10(input$qcutoffLinear),linetype="dashed")+
          geom_vline(xintercept=log2(input$fccutoffLinear),linetype="dashed")+
          geom_vline(xintercept=-log2(input$fccutoffLinear),linetype="dashed")+
          labs(x= paste("Beta Coef. ",input$variableTest2,sep=""), y= "-Log10(q-value)",
               color=paste("sig_",input$variableTest2,sep=""))+
          theme(axis.title = element_text(size = 18),
                axis.text = element_text(size=16),
                legend.text=element_text(size=14),
                legend.title=element_text(size=16))
      }
      
      ggplotly(volcano)%>% layout(dragmode = "select")
      
      
    })
    
    
    volcanoplotLinearPrint <- reactive({
      # req(input$performlineartest == "Do Stats")
      if (input$betachoice == "Var1Var2 interaction" ) {
        volcano<-ggplot(linearmodelDFSigCall(), aes( coeff_var1_var2, -log10(q.value_var1_var2),color=sig_var1_var2))+
          geom_point(alpha=0.75,size=0.75)+theme_classic()+
          scale_color_viridis_d(end=0.8)+
          geom_hline(yintercept = -log10(input$qcutoffLinear),linetype="dashed")+
          geom_vline(xintercept=log2(input$fccutoffLinear),linetype="dashed")+
          geom_vline(xintercept=-log2(input$fccutoffLinear),linetype="dashed")+
          labs(x= paste("Beta Coef. ",input$variableTest1,"_",input$variableTest2," Interaction Term",sep=""),
               y= "-Log10(q-value)",
               color=paste("sig_",input$variableTest1,"_",input$variableTest2,sep=""))+
          theme(axis.title = element_text(size = 18),
                axis.text = element_text(size=16),
                legend.text=element_text(size=14),
                legend.title=element_text(size=16))
        
      }
      if (input$betachoice == "Var1" ) {
        volcano<-ggplot(linearmodelDFSigCall(), aes( coeff_var1, -log10(q.value_var1),color=sig_var1))+
          geom_point(alpha=0.75,size=0.75)+theme_classic()+
          scale_color_viridis_d(end=0.8)+
          geom_hline(yintercept = -log10(input$qcutoffLinear),linetype="dashed")+
          geom_vline(xintercept=log2(input$fccutoffLinear),linetype="dashed")+
          geom_vline(xintercept=-log2(input$fccutoffLinear),linetype="dashed")+
          labs(x= paste("Beta Coef. ",input$variableTest1,sep=""), y= "-Log10(q-value)",
               color=paste("sig_",input$variableTest1,sep=""))+
          theme(axis.title = element_text(size = 18),
                axis.text = element_text(size=16),
                legend.text=element_text(size=14),
                legend.title=element_text(size=16))
      }
      if (input$betachoice == "Var2" ) {
        volcano<-ggplot(linearmodelDFSigCall(), aes( coeff_var2, -log10(q.value_var2),color=sig_var2))+
          geom_point(alpha=0.75,size=0.75)+theme_classic()+
          scale_color_viridis_d(end=0.8)+
          geom_hline(yintercept = -log10(input$qcutoffLinear),linetype="dashed")+
          geom_vline(xintercept=log2(input$fccutoffLinear),linetype="dashed")+
          geom_vline(xintercept=-log2(input$fccutoffLinear),linetype="dashed")+
          labs(x= paste("Beta Coef. ",input$variableTest2,sep=""), y= "-Log10(q-value)",
               color=paste("sig_",input$variableTest2,sep=""))+
          theme(axis.title = element_text(size = 18),
                axis.text = element_text(size=16),
                legend.text=element_text(size=14),
                legend.title=element_text(size=16))
      }
  
      return(volcano)
    })
  
    output$DownloadVolcanoLinear <- downloadHandler(
      filename = function() { paste("VolcanoPlot_LinearBetas_",input$variableTest1,"_",input$variableTest2, '.pdf', sep='') },
      content = function(file) {
        ggsave(file, plot = volcanoplotLinearPrint(), width = 10, height = 6)
      }
    )
    
  })

  #### END LINEAR MODEL
  
  priorityTable <- reactive({
    priorityTable <- norm_sum_final2()%>%
      dplyr::ungroup()%>%
      dplyr::distinct(Condition,priority)%>%
      dplyr::arrange(priority)
    return(priorityTable)
  })
  output$priorityTable100 <- renderReactable({
    reactable(priorityTable())
  })
  output$priorityTable1 <- renderReactable({
    reactable(priorityTable())
  })
  output$priorityTable2 <- renderReactable({
    reactable(priorityTable())
  })
  output$priorityTable3 <- renderReactable({
    reactable(priorityTable())
  })
  output$priorityTable4 <- renderReactable({
    reactable(priorityTable())
  })
  output$priorityTable5 <- renderReactable({
    reactable(priorityTable())
  })
  output$priorityTable6 <- renderReactable({
    reactable(priorityTable())
  })
  output$priorityTable7 <- renderReactable({
    reactable(priorityTable())
  })
  output$priorityTable8 <- renderReactable({
    reactable(priorityTable())
  })
  output$priorityTable9 <- renderReactable({
    reactable(priorityTable())
  })
  output$priorityTable10 <- renderReactable({
    reactable(priorityTable())
  })
  
  
  anova_table <- reactive({
    
    anova_table <- df()%>%
      dplyr::select(Reference,Gene.Symbol,Annotation, contains("anova"))%>%
      tidyr::unite("ProtID",c(Reference, Gene.Symbol, Annotation),sep="__X__")
    if (input$inputtype == "PAbeta") {
      anova_table<-anova_table%>%
        dplyr::rename(sig_anova = anova_sig)
    }
    return(anova_table)
  })
  
  
  prot_biorep_final <- reactive({
    bioreps <- df()%>%
      dplyr::select(Reference,Gene.Symbol,Annotation,unique(norm_sum_final2()$BioReplicate) )%>%
      tidyr::unite("ProtID",c(Reference, Gene.Symbol, Annotation),sep="__X__")
    return(bioreps)
  })
  
  results_msstats <- reactive({
    msstatsresults <- df()%>%
      # dplyr::select(Reference,Gene.Symbol,Annotation,!unique(norm_sum_final2()$BioReplicate) )%>%
      # dplyr::select(!unique(norm_sum_final2()$Condition))%>%
      # dplyr::select(!contains("anova"))%>%
      dplyr::select(!contains("cluster"))
    return(msstatsresults)
  })
  

  
  


  
  df_final_res <- reactive({
    
    if (input$inputtype == "MSstats") {
      exp_design2<-exp_design()%>%
        dplyr::filter(Condition %!in% c("Empty","empty"))%>%
        dplyr::mutate(BioReplicate = paste(Condition,Replicate,sep="_"))%>%
        dplyr::select(BioReplicate,Condition, priority)
    }
    if (input$inputtype == "PAbeta") {
      exp_design2<-exp_design()%>%
        dplyr::rename(Condition = condition,
                      Replicate = replicate)%>%
        dplyr::filter(Condition %!in% c("Empty","empty"))%>%
        # dplyr::filter(Condition != "empty")%>%
        dplyr::mutate(BioReplicate = paste(Condition,Replicate,sep="_"))%>%
        dplyr::select(BioReplicate,Condition, priority)
    }
    
    
    prot_norm_data_reps_tidy <- prot_biorep_final()%>%
      dplyr::group_by(ProtID)%>%
      tidyr::gather("BioReplicate", "Abundance", 2:dim(prot_biorep_final())[[2]])%>%
      dplyr::ungroup()%>%
      dplyr::left_join(.,exp_design2, by="BioReplicate")
    
    prot_norm_data_reps_tidy$BioReplicate<-reorder(prot_norm_data_reps_tidy$BioReplicate,prot_norm_data_reps_tidy$priority,min)
    prot_norm_data_reps_tidy$Condition<-reorder(prot_norm_data_reps_tidy$Condition,prot_norm_data_reps_tidy$priority,min)

    prot_norm_data_reps_tidy<-prot_norm_data_reps_tidy%>%
      dplyr::ungroup()%>%
      dplyr::select(-priority)%>%
      dplyr::select(ProtID,  Condition,BioReplicate, Abundance)%>%
      dplyr::group_by(ProtID,Condition,  BioReplicate)%>%
      dplyr::summarise(Abundance = mean(Abundance,na.rm=TRUE))%>%
      dplyr::mutate(Abundance = replace_na(Abundance, 0))
    return(prot_norm_data_reps_tidy)
  })
  

  
  untidy_stats_all <- reactive({
    untidy_stats_all <- dplyr::full_join(df_final_res()%>%
                                           tidyr::separate(ProtID, into=c("Reference","Gene.Symbol","Annotation"),sep="__X__"),
                                         results_msstats(),by=c("Reference","Gene.Symbol","Annotation"))
    return(untidy_stats_all)
  })
  
 

  
  output$totalproteins <- renderText({
    paste(dim(prot_biorep_final())[[1]]," total proteins analyzed")
  })
  
  output$plot4 <- renderPlot({
    lowerfun <- function(data,mapping){
      ggplot(data = data, mapping = mapping)+
        geom_point(size=0.1,alpha=0.1)+
        scale_x_continuous(limits = c(input$mincorrval,NA))+
        scale_y_continuous(limits = c(input$mincorrval,NA))
    } 
    
    #ggpairs(prot_biorep_final()[2:dim(prot_biorep_final())[[2]]],lower = list(continuous = wrap("points",size=0.1)))+theme_classic()
    ggpairs(prot_biorep_final()[2:dim(prot_biorep_final())[[2]]],lower = list(continuous = wrap(lowerfun)))+theme_classic()
    
  }, height = 1000, width = 1000 )
  
  
 
  
  
  output$plotPCA <- renderPlot({
    
    df_pca<-prot_biorep_final()%>%
      tibble::column_to_rownames(var="ProtID")
    example_pca_test<-prcomp(t(df_pca))
    
    df_pca_use<-data.frame(example_pca_test$x)%>%
      tibble::rownames_to_column("sample")

    PoV <- example_pca_test$sdev^2/sum(example_pca_test$sdev^2)
    pc_val1<-input$txtXpca
    pc_val2<-input$txtYpca

    pc_val1_num <- as.numeric(as.character(str_remove(pc_val1,"PC")))

    pc_val2_num <- as.numeric(as.character(str_remove(pc_val2,"PC")))
    
    pc_val1_num_plot<-pc_val1_num+1
    pc_val2_num_plot<-pc_val2_num+1

    var_pcVal1<-round(PoV[pc_val1_num]*100,4)
    var_pcVal2<-round(PoV[pc_val2_num]*100,4)
    
    df_pca_use<-df_pca_use%>%
               dplyr::mutate(loc_last = stringi::stri_locate_last_fixed(sample,"_")[2]-1,
                      Condition = stringr::str_sub(sample, start=1,end=-3))

    
    #print(df_pca_use)

    ggplot(df_pca_use,aes(df_pca_use[,pc_val1_num_plot],df_pca_use[,pc_val2_num_plot],color=Condition,label=sample))+
      geom_point(size=3,alpha=0.5)+
      geom_text_repel(max.overlaps = Inf,box.padding=1,show.legend=FALSE)+
      theme_classic()+scale_color_viridis_d(end=0.8)+
      labs(x=paste(input$txtXpca,": ", var_pcVal1,"%",sep=""),
           y=paste(input$txtYpca,": ", var_pcVal2,"%",sep=""))+
      theme(axis.title = element_text(size=16),axis.text = element_text(size=14),legend.position = "none")
  })
  
  
  output$PCAdrivers <- renderReactable({
    correlation_plot<-prot_biorep_final()%>%
      #tibble::remove_rownames %>%
      tibble::column_to_rownames(var="ProtID")
    example_pca_test<-prcomp(t(correlation_plot))

    var <- factoextra::get_pca_var(example_pca_test)
    contrib <- as.data.frame(var$contrib)

    contrib<-contrib%>%
      tibble::rownames_to_column("ProtID")%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    ##choose PC for drivers
    pc_val_drive<-input$txtZpca
    pc_val_drive1 <- as.numeric(as.character(str_remove(pc_val_drive,"PC")))
    # print(pc_val_drive1)
    pc_val_drive_sel<-pc_val_drive1 +3

    sel_column<-contrib%>%
      dplyr::select(1:3, all_of(pc_val_drive_sel))

    sel_column<-sel_column%>%
      plotly::arrange(desc(.[,4]))

    colnames(sel_column) <- sub("Dim\\.", "PC", colnames(sel_column))
    #sel_column<- sel_column[1:100,]
    reactable(sel_column, filterable = TRUE)
  })
  
  driverAll<- reactive({
    correlation_plot<-prot_biorep_final()%>%
      tibble::column_to_rownames(var="ProtID")
    example_pca_test<-prcomp(t(correlation_plot))
    
    var <- factoextra::get_pca_var(example_pca_test)
    contrib <- as.data.frame(var$contrib)
    
    contrib<-contrib%>%
      tibble::rownames_to_column("ProtID")%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    ##choose PC for drivers
    
    sel_column<-contrib%>%
      plotly::arrange(desc(.[,4]))
    
    colnames(sel_column) <- sub("Dim\\.", "PC", colnames(sel_column))
    #sel_column<- sel_column[1:100,]
    return(sel_column)
  })
  
  
  
  output$driversDownload <- downloadHandler(
    filename = function() {
      paste("PCA_protein_weights_all",".csv", sep = "")
    },
    content = function(file) {
      write.csv(driverAll(), file, row.names = FALSE)
    }
  )
  
  
  output$plotPCAadd <-renderPlot({
    correlation_plot<-prot_biorep_final()%>%
      #tibble::remove_rownames %>%
      tibble::column_to_rownames(var="ProtID")
    example_pca_test<-prcomp(t(correlation_plot))
    factoextra::fviz_eig(example_pca_test, addlabels = TRUE) + theme_classic()
  })
  
  
  pcaPlotInput <- reactive({
    correlation_plot<-prot_biorep_final()%>%
      tibble::column_to_rownames(var="ProtID")
    example_pca_test<-prcomp(t(correlation_plot))

    df_pca_use<-data.frame(example_pca_test$x)%>%
      tibble::rownames_to_column("sample")

    PoV <- example_pca_test$sdev^2/sum(example_pca_test$sdev^2)
    pc_val1<-input$txtXpca
    pc_val2<-input$txtYpca

    pc_val1_num <- as.numeric(as.character(str_remove(pc_val1,"PC")))

    pc_val2_num <- as.numeric(as.character(str_remove(pc_val2,"PC")))

    pc_val1_num_plot<-pc_val1_num+1
    pc_val2_num_plot<-pc_val2_num+1

    var_pcVal1<-round(PoV[pc_val1_num]*100,4)
    var_pcVal2<-round(PoV[pc_val2_num]*100,4)

    df_pca_use<-df_pca_use%>%
      dplyr::mutate(loc_last = stringi::stri_locate_last_fixed(sample,"_")[2]-1,
                    Condition = stringr::str_sub(sample, start=1,end=loc_last))

    p70<-ggplot(df_pca_use,aes(df_pca_use[,pc_val1_num_plot],df_pca_use[,pc_val2_num_plot],color=Condition,label=sample))+
      geom_point(size=3,alpha=0.5)+
      geom_text_repel(max.overlaps = Inf,box.padding=1,show.legend=FALSE)+
      theme_classic()+scale_color_viridis_d(end=0.8)+
      labs(x=paste(input$txtXpca,": ", var_pcVal1,"%",sep=""),
           y=paste(input$txtYpca,": ", var_pcVal2,"%",sep=""))+
      theme(axis.title = element_text(size=16),axis.text = element_text(size=14),legend.position = "none")

  })
  

  
  output$PCADownloadPlot <- downloadHandler(
    filename = function() { paste("PCA_Plot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = pcaPlotInput(), width = 10, height = 10)
    }
  )
  
  corrPlotInput <- reactive({
    # p72<-ggpairs(prot_biorep_final()[2:dim(prot_biorep_final())[[2]]],lower = list(continuous = wrap("points",size=0.1)))+theme_classic()
    
    lowerfun <- function(data,mapping){
      ggplot(data = data, mapping = mapping)+
        geom_point(size=0.1,alpha=0.1)+
        scale_x_continuous(limits = c(input$mincorrval,NA))+
        scale_y_continuous(limits = c(input$mincorrval,NA))
    } 
    
    #ggpairs(prot_biorep_final()[2:dim(prot_biorep_final())[[2]]],lower = list(continuous = wrap("points",size=0.1)))+theme_classic()
    p72<-ggpairs(prot_biorep_final()[2:dim(prot_biorep_final())[[2]]],lower = list(continuous = wrap(lowerfun)))+theme_classic()
    
  })
  
  
  output$CorrDownloadPlot <- downloadHandler(
    filename = function() { paste("Rep_Corr_Plot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corrPlotInput(), width = 12, height = 12)
    }
  )
  
  

  heatmapinput <-  reactive({
    
    if (input$inputtype == "MSstats") {
      exp_design2<-exp_design()%>%
        # dplyr::filter(Condition != "Empty")%>%
        # dplyr::rename(Condition = condition,
        #               Replicate = replicate)%>%
        dplyr::filter(Condition %!in% c("Empty","empty"))%>%
        dplyr::mutate(BioReplicate = paste(Condition,Replicate,sep="_"))%>%
        dplyr::select(BioReplicate, priority)
    }
    if (input$inputtype == "PAbeta") {
      exp_design2<-exp_design()%>%
        # dplyr::filter(Condition != "Empty")%>%
        dplyr::rename(Condition = condition,
                      Replicate = replicate)%>%
        dplyr::filter(Condition %!in% c("Empty","empty"))%>%
        dplyr::mutate(BioReplicate = paste(Condition,Replicate,sep="_"))%>%
        dplyr::select(BioReplicate, priority)
    }
    
    heatmap_input<-prot_biorep_final()%>%
      tidyr::gather("BioReplicate","Abundance",2:dim(prot_biorep_final())[[2]])%>%
      # dplyr::select(-intensity,-correct_factor,-relative_intensity,-norm_int)%>%
      dplyr::ungroup()%>%
      dplyr::group_by(ProtID)%>%
      dplyr::mutate(mean_row = mean(Abundance,na.rm=TRUE),
             sd_row = sd(Abundance,na.rm=TRUE),
             z_scale = (Abundance - mean_row)/sd_row,
             mean_scale = Abundance - mean_row)%>%
      dplyr::ungroup()%>%
      dplyr::left_join(.,exp_design2 , by="BioReplicate")
    
    heatmap_input$BioReplicate<-reorder(heatmap_input$BioReplicate,heatmap_input$priority,min)
    
    heatmap_input<-heatmap_input%>%
      dplyr::select(-priority)
    
    
    if (input$txtScale == "RowMean") {
      df_final2<-heatmap_input%>%
        dplyr::ungroup()%>%
        dplyr::select(ProtID,BioReplicate,mean_scale)%>%
        tidyr::spread(BioReplicate,mean_scale)%>%
        tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    } 
    if (input$txtScale == "Zscore") {
      df_final2<-heatmap_input%>%
        dplyr::ungroup()%>%
        dplyr::select(ProtID,BioReplicate,z_scale)%>%
        tidyr::spread(BioReplicate,z_scale)%>%
        tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    }
    
    
    if (input$anovaorall == "ANOVA") {

      sig<-anova_table()%>%dplyr::filter(sig_anova == "sig")%>%
        tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")

      
    }
    
    if (input$anovaorall == "All") {
      sig<-anova_table()%>%
        tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    }
    
    Refs<-sig$Reference
    
    
    df_final2<-df_final2%>%
      dplyr::filter(Reference %in% Refs)
    
    df_final2<-df_final2%>%
      tidyr::unite("ProtID",Reference,Gene.Symbol,Annotation, sep="__X__")

    df4<-df_final2%>%
      tibble::column_to_rownames("ProtID")
    return(df4)
  })
  
  userheatmapreference <- reactive({
    df<-heatmapinput()%>%
      tibble::rownames_to_column(var = "ProtID")%>%
      tidyr::separate(ProtID,into=c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    return(df)
  })
  
  output$numclusterplot <- renderPlot({
    
    k.means.opt <- function(my.clusters) {
      k <- stats::kmeans(heatmapinput(), centers=my.clusters, iter.max = 500)
      if (is.finite(k$betweenss/k$totss)==T) {
        value<-k$betweenss/k$totss
      }
      return(value)
    }
    
    cluster.max <- 10
    ## make empty list to put the data
    k.means.opt.list <- vector("list", length=cluster.max) 
    clus_num<-c()
    max_var<-c()
    ## replicate function where, for each cluster from 1:cluster.max, replicates 1000 times and returns the (between.ss/tot.ss) ratio. Each of the [[elements]] of the list is 1000 ratios for each number of clusters. Also I want to time this step.
    for (i in 1:cluster.max) {
      k.means.opt.list[[i]] <- replicate(100, k.means.opt(i))
      clus_num <- c(clus_num,i)
      max_var<-c(max_var,max(k.means.opt.list[[i]]))
    }

    ## unlist your k.means.opt.list object to extract Ratio.BSS.TSS. Take the maximum Ratio.BSS.TSS of each 1000 kmeans replicates for each clutser, i (in this case, 1000 kmeans replicates for each cluster from 1 cluster (should be 0%) to 20 clusters (should be close to 100%)). If the number of clusters = number of data points, the "Ratio" will be undefined, which will show on the graph as 0% even though it's undefined.
    # k.means.data <- data.frame(Number.of.clusters= clusters.x, Ratio.BSS.TSS=unlist(lapply(k.means.opt.list, max)))
    
    k.means.data<-data.frame("cluster_num"=clus_num, "VarExp" = max_var)
    
    ## now plot these two things
    ggplot(k.means.data, aes(x = cluster_num, y = VarExp)) +
      geom_point(size=5,alpha=0.5) + geom_line()+
      theme_classic()+
      labs(y="BSS/TSS : Variance Explained B/w Groups", x="Cluster Number")+
      theme(axis.text = element_text(size = 18), axis.title = element_text(size=19))+ 
      scale_x_continuous(breaks=seq(0,10,2))+
      ylim(c(0,1))
    
  
  })
  
  output$heatmapplot <- renderPlot({
    #cols <- colorRampPalette(brewer.pal(6,name="RdBu"))(12)
    # brks <- seq(-2,2,length.out=12)
    rg <- max(abs(heatmapinput()))
    
    cluster_design_df<-clusterTraceDFinitial()%>%
      dplyr::ungroup()%>%
      dplyr::select(ProtID,cluster)%>%
      dplyr::mutate(cluster=as.factor(cluster))
    
    cluster_order <- paste(sort(as.integer(levels(cluster_design_df$cluster))))
    
    cluster_design_df$cluster <- factor(cluster_design_df$cluster, levels = cluster_order)
    
    AnnoB<-data.frame(row.names=cluster_design_df$ProtID, ClusterNumber=cluster_design_df$cluster)
    
    #, cluster_cols = FALSE,color=cols,
    if(input$txtclustercols == "Yes"){
      names_vals<-heatmapinput()%>%
        tibble::rownames_to_column("ProtID")
      names_cols<-data.frame("sample" = colnames(names_vals)[2:dim(names_vals)[[2]]])%>%
        dplyr::mutate(Condition = str_sub(sample,start = 1,end = -3))%>%
        dplyr::left_join(.,norm_sum_final2()%>%dplyr::distinct(priority,Condition),by="Condition")%>%
        dplyr::arrange(priority,Condition,sample)%>%
        dplyr::mutate(sample = factor(sample, unique(sample)),
                      values = factor(Condition,unique(Condition)))
      
      annoD<-data.frame(row.names=names_cols$sample, Condition=names_cols$values)
      
      
      p1<-pheatmap::pheatmap(heatmapinput(),   cutree_rows = input$numclusters,main = "Heatmap",annotation = annoD,annotation_row = AnnoB,show_rownames=F, breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercols == "No"){

      names_vals<-heatmapinput()%>%
        tibble::rownames_to_column("ProtID")
      names_cols<-data.frame("sample" = colnames(names_vals)[2:dim(names_vals)[[2]]])%>%
        dplyr::mutate(Condition = str_sub(sample,start = 1,end = -3))%>%
        dplyr::left_join(.,norm_sum_final2()%>%dplyr::distinct(priority,Condition),by="Condition")%>%
        dplyr::arrange(priority,Condition,sample)%>%
        dplyr::mutate(sample = factor(sample, unique(sample)),
                      values = factor(Condition,unique(Condition)))
      
      annoD<-data.frame(row.names=names_cols$sample, Condition=names_cols$values)
      
      p1<-pheatmap::pheatmap(heatmapinput(),  cutree_rows = input$numclusters,main = "Heatmap",annotation = annoD,annotation_row = AnnoB,show_rownames=F, cluster_cols = FALSE, breaks = seq(-rg, rg, length.out = 100))
    }
    p1
    
    
  })
  
  
  print <- reactive({
    #cols <- colorRampPalette(brewer.pal(6,name="RdBu"))(12)
    # brks <- seq(-2,2,length.out=12)
    rg <- max(abs(heatmapinput()))
    cluster_design_df<-clusterTraceDFinitial()%>%
      dplyr::ungroup()%>%
      dplyr::select(ProtID,cluster)%>%
      dplyr::mutate(cluster=as.factor(cluster))
    
    cluster_order <- paste(sort(as.integer(levels(cluster_design_df$cluster))))
    
    cluster_design_df$cluster <- factor(cluster_design_df$cluster, levels = cluster_order)
    
    AnnoB<-data.frame(row.names=cluster_design_df$ProtID, ClusterNumber=cluster_design_df$cluster)
    
    #, cluster_cols = FALSE,color=cols,
    if(input$txtclustercols == "Yes"){
      names_vals<-heatmapinput()%>%
        tibble::rownames_to_column("ProtID")
      names_cols<-data.frame("sample" = colnames(names_vals)[2:dim(names_vals)[[2]]])%>%
        dplyr::mutate(Condition = str_sub(sample,start = 1,end = -3))%>%
        dplyr::left_join(.,norm_sum_final2()%>%dplyr::distinct(priority,Condition),by="Condition")%>%
        dplyr::arrange(priority,Condition,sample)%>%
        dplyr::mutate(sample = factor(sample, unique(sample)),
                      values = factor(Condition,unique(Condition)))
      
      annoD<-data.frame(row.names=names_cols$sample, Condition=names_cols$values)
      p1<-pheatmap::pheatmap(heatmapinput(),   cutree_rows = input$numclusters,main = "Heatmap",annotation_row = AnnoB,show_rownames=F, breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercols == "No"){
      names_vals<-heatmapinput()%>%
        tibble::rownames_to_column("ProtID")
      names_cols<-data.frame("sample" = colnames(names_vals)[2:dim(names_vals)[[2]]])%>%
        dplyr::mutate(Condition = str_sub(sample,start = 1,end = -3))%>%
        dplyr::left_join(.,norm_sum_final2()%>%dplyr::distinct(priority,Condition),by="Condition")%>%
        dplyr::arrange(priority,Condition,sample)%>%
        dplyr::mutate(sample = factor(sample, unique(sample)),
                      values = factor(Condition,unique(Condition)))
      
      annoD<-data.frame(row.names=names_cols$sample, Condition=names_cols$values)
      p1<-pheatmap::pheatmap(heatmapinput(),  cutree_rows = input$numclusters,main = "Heatmap",show_rownames=F,annotation_row = AnnoB, cluster_cols = FALSE, breaks = seq(-rg, rg, length.out = 100))
    }
    return(p1)
    
    
  })
  
  output$plotheatmapprint <- downloadHandler(
    filename = function() { paste("Heatmap", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = heatmapplotprint(), width = 10, height = 10)
    }
  )
  

  clusterTraceDFinitial <- reactive({
    #, cluster_cols = FALSE
    
    if(input$txtclustercols == "Yes"){
      out<-pheatmap::pheatmap(heatmapinput(),   cutree_rows = input$numclusters,main = "Heatmap",show_rownames=F)
      heatmapinput2<-heatmapinput()%>%
        dplyr::select(colnames(heatmapinput()[,out$tree_col[["order"]]]))
      
      order_levels<-colnames(heatmapinput()[,out$tree_col[["order"]]])
    }
    
    if(input$txtclustercols == "No"){
      out<-pheatmap::pheatmap(heatmapinput(),  cutree_rows = input$numclusters,main = "Heatmap",show_rownames=F, cluster_cols = FALSE)
      heatmapinput2<-heatmapinput()
    }
    
    df4.clust <- cbind(heatmapinput2,
                       cluster = cutree(out$tree_row,
                                        k = input$numclusters))
    df4.clust2<-df4.clust%>%
      tibble::rownames_to_column("ProtID")%>%
      dplyr::group_by(ProtID,cluster)

    return(df4.clust2)
    
  })
  
  ### START TIMECOURSE
  observeEvent(input$runTimecourse, {
    
  
    hotellingdf<- reactive({
      # if (input$txttimecourse == "Yes") {
      prot_values <- prot_biorep_final()$ProtID
  
      df_timecourse <- prot_biorep_final()%>%
        tibble::column_to_rownames(var = "ProtID")
  
      values <- c(names(df_timecourse))
      time.grp <- as.numeric(readr::parse_number(values))
      assay <- as.numeric(stringr::str_sub(values,-1,-1))
      if (input$numrepstimecourse == 3) {
        size <- rep(3, dim(df_timecourse)[[1]])
      }
  
      if (input$numrepstimecourse == 4) {
        size <- rep(4, dim(df_timecourse)[[1]])
      }
  
      if (input$numrepstimecourse == 2) {
        size <- rep(2, dim(df_timecourse)[[1]])
      }
  
  
      df_timecourse <- as.matrix(df_timecourse)
  
      hotelling <- timecourse::mb.long(df_timecourse, times=length(unique(time.grp)), reps=size, rep.grp=assay, time.grp=time.grp)
      hotelling_df<-data.frame("ProtID" = prot_values, "hotellingT2" = hotelling$HotellingT2)
      # }
  
      # if (input$txttimecourse == "No") {
      #   prot_values <- prot_biorep_final()$ProtID
      #   hotelling_df<-data.frame("ProtID" = prot_values, "hotellingT2" = 0)
      # }
  
      return(hotelling_df)
    })
  
    clustertimecourse<- reactive({
  
      clusterdesignation<-clusterTraceDFinitial()%>%
        dplyr::ungroup()%>%
        dplyr::select(ProtID,cluster)
  
      HotelClusterDF <- dplyr::left_join(hotellingdf(),clusterdesignation,by="ProtID")%>%
        dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
        # the good stuff here
        dplyr::mutate_if(is.numeric,dplyr::coalesce,0)
  
      HotelClusterDF <- HotelClusterDF%>%
        #dplyr::mutate(cluster = dplyr::if_else(is.na(cluster), "0", cluster))%>%
        tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
  
      all_table<-left_join(tableresults(),HotelClusterDF, by= c("Reference","Gene.Symbol","Annotation","cluster"))%>%
        dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))
  
  
      return(all_table%>%dplyr::arrange(desc(hotellingT2)))
    })
    
    
   
    
    
    output$hotellingdf2 <- renderReactable({
      reactable(clustertimecourse(),filterable = T)
    })
    
    
    output$hotellingdownload <- downloadHandler(
      filename = function() {
        paste("hotelling_results", ".csv", sep = "")
      },
      content = function(file) {
        write.csv(clustertimecourse(), file, row.names = FALSE)
      }
    )
  
  
    output$plothotelling <- renderPlotly({
  
      con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2])
  
      foldchangevalue2<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
  
      photelling<-ggplot(clustertimecourse(), aes( .data[[foldchangevalue2]], log10(hotellingT2), color=factor(cluster),label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        labs(x= paste(foldchangevalue2), y= "Log10(Hotelling T^2)",color="Cluster")
  
      ggplotly(photelling)%>% layout(dragmode = "select")
    })
    
    
    plothotellingprint <- reactive({
  
      con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2])
  
      foldchangevalue2<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
  
      photelling<-ggplot(clustertimecourse(), aes( .data[[foldchangevalue2]], log10(hotellingT2), color=factor(cluster)))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        labs(x= paste(foldchangevalue2), y= "Log10(Hotelling T^2)",color="Cluster")
  
      return(photelling)
  
  
    })
    
    
    output$printhotellingplot<- downloadHandler(
      filename = function() { paste("HotellingPlot", '.pdf', sep='') },
      content = function(file) {
        ggsave(file, plot = plothotellingprint(), width = 9, height = 5)
      }
    )
    
    
    
    
    output$hotellingclick <- renderReactable({
      d <- event_data("plotly_selected")
      req(d)
      df_sub3<-clustertimecourse()%>%
        dplyr::filter(ProtID %in% d$key )%>%
        dplyr::ungroup()
  
      reactable(df_sub3, filterable = TRUE)
    })
    
    
    output$plottcprot <- renderPlot({
      req(input$text113)
      gene_name_input <- input$text113
  
      clusterdesignation<-clusterTraceDFinitial()%>%
        dplyr::ungroup()%>%
        dplyr::select(ProtID,cluster)
  
      HotelClusterDF <- dplyr::left_join(hotellingdf(),clusterdesignation,by="ProtID")%>%
        dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
        # the good stuff here
        dplyr::mutate_if(is.numeric,dplyr::coalesce,0)
  
      HotelClusterDF <- HotelClusterDF%>%
        tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
  
      plot_protein<-left_join(untidy_stats_all(),HotelClusterDF, by= c("Reference","Gene.Symbol","Annotation"))%>%
        dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))
  
  
      plot_protein<-plot_protein%>%
        dplyr::filter(Gene.Symbol==gene_name_input)%>%
        dplyr::mutate(Condition = as.character(Condition),
                      timepoint = readr::parse_number(Condition))%>%
        dplyr::mutate(ProtID2 = paste(Gene.Symbol,Reference,sep="_"))
  
  
      plot_protein_v2<-plot_protein%>%
        dplyr::ungroup()%>%
        dplyr::group_by(Condition,timepoint,ProtID2,Reference)%>%
        dplyr::summarise(Abundance=median(Abundance))%>%
        dplyr::mutate(timepoint = as.numeric(as.character(timepoint)))
  
      plot_protein<-plot_protein%>%
        dplyr::mutate(timepoint=as.numeric(as.character(timepoint)))
  
      ggplot()+
        geom_line(data=plot_protein_v2,aes(x=timepoint,y=2^Abundance,color=Reference),alpha=0.5)+
        geom_point(data = plot_protein ,aes(x=timepoint,y=2^Abundance,color=Reference),size=3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=24,hjust=0.5),legend.position = "top")+
        labs(title = gene_name_input, x="Time",y="MS Intensity")
  
    })
  
    plottcprotprint <- reactive({
      req(input$text113)
      gene_name_input <- input$text113
  
      clusterdesignation<-clusterTraceDFinitial()%>%
        dplyr::ungroup()%>%
        dplyr::select(ProtID,cluster)
  
      HotelClusterDF <- dplyr::left_join(hotellingdf(),clusterdesignation,by="ProtID")%>%
        dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
        # the good stuff here
        dplyr::mutate_if(is.numeric,dplyr::coalesce,0)
  
      HotelClusterDF <- HotelClusterDF%>%
        tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
  
      plot_protein<-left_join(untidy_stats_all(),HotelClusterDF, by= c("Reference","Gene.Symbol","Annotation"))%>%
        dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))
  
      plot_protein<-plot_protein%>%
        dplyr::filter(Gene.Symbol==gene_name_input)%>%
        dplyr::mutate(Condition = as.character(Condition),
                      timepoint = readr::parse_number(Condition))%>%
        dplyr::mutate(ProtID2 = paste(Gene.Symbol,Reference,sep="_"))
  
  
  
      plot_protein_v2<-plot_protein%>%
        dplyr::ungroup()%>%
        dplyr::group_by(Condition,timepoint,ProtID2,Reference)%>%
        dplyr::summarise(Abundance=median(Abundance))%>%
        dplyr::mutate(timepoint = as.numeric(as.character(timepoint)))
  
      plot_protein<-plot_protein%>%
        dplyr::mutate(timepoint=as.numeric(as.character(timepoint)))
  
      p1<-ggplot()+
        geom_line(data=plot_protein_v2,aes(x=timepoint,y=2^Abundance,color=Reference),alpha=0.5)+
        geom_point(data = plot_protein ,aes(x=timepoint,y=2^Abundance,color=Reference),size=3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=24,hjust=0.5),legend.position = "top")+
        labs(title = gene_name_input, x="Time",y="MS Intensity")
      return(p1)
  
    })
    
  
    
    output$plottcprotprint2 <- downloadHandler(
      filename = function() { paste("Timecourse_", input$text113 ,'.pdf', sep='') },
      content = function(file) {
        ggsave(file, plot = plottcprotprint(), width = 10, height = 10)
      }
    )
  })
  
  ### END time course
  
  output$volcanoplot <- renderPlotly({

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1vol])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2vol])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))

    if (input$txtinversion == "No") {
      volcano<-ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con2,"-",con1,sep=""), y= "-Log10(q-value)")
    }

    if (input$txtinversion == "Yes") {
      volcano<-ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con1,"-",con2,sep=""), y= "-Log10(q-value)")
    }

    ggplotly(volcano)%>% layout(dragmode = "select")


  })
  

  volcanoplotprint <- reactive({

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1vol])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2vol])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))

    if (input$txtinversion == "No") {
      volcano<-ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con2,"-",con1,sep=""), y= "-Log10(q-value)")
    }

    if (input$txtinversion == "Yes") {
      # foldchangevalue2<-as.character(paste("log2FC_",con1,"-",con2,sep=""))
      volcano<-ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con1,"-",con2,sep=""), y= "-Log10(q-value)")
    }

    return(volcano)


  })
  
  
  output$printvolcano<- downloadHandler(
    filename = function() { paste("Volcano",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2vol]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1vol]) ,'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = volcanoplotprint(), width = 9, height = 5)
    }
  )
  
  
  
  
  output$volcanoclick <- renderReactable({
    d <- event_data("plotly_selected")
    req(d)
    df_sub2<-tableresults()%>%
      dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))%>%
      dplyr::filter(ProtID %in% d$key )%>%
      dplyr::ungroup()

    reactable(df_sub2, filterable = TRUE)
  })

  # output$volcanoplotnorm <- renderPlotly({
  # 
  #   con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse3vol])
  #   con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse4vol])
  #   con3<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse5vol])
  #   con4<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse6vol])
  # 
  # 
  #   pair_plot<-tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))
  # 
  #   volcano2<-ggplot(pair_plot, aes( .data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]],label=Gene.Symbol, key=ProtID))+
  #     geom_point(alpha=0.3,size=0.5)+theme_classic()+
  #     geom_abline(slope=1,linetype="dashed",size=0.5)+
  #     geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
  #     geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
  #     labs(x= paste(con2,"-",con1), y= paste(con4,"-",con3))
  # 
  #   ggplotly(volcano2)%>% layout(dragmode = "select")
  # 
  # })
  
  ### adjust to multipl
  
  output$volcanoplotnorm<-renderPlotly({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse3vol])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse4vol])
    con3<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse5vol])
    con4<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse6vol])
    
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    fc2<-as.character(paste("log2FC_",con4,"-",con3,sep=""))
    
    
    pair_plot<-tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))
    
    
    
    if (input$ratioratioindicatorVolcano == "Single Ratio") {
      
      pair_plot2 <- pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID, all_of(fc1), all_of(fc2) )
      
      volcano2<-ggplot(pair_plot2, aes( .data[[fc1]], .data[[fc2]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.3,size=0.5)+theme_classic()+
        geom_abline(slope=1,linetype="dashed",size=0.5)+
        geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
        geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
        labs(x= paste(fc1), y= paste(fc2))
      
    }
    
    if (input$ratioratioindicatorVolcano == "X axis RoR") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$denom1])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$denom2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))
      
      join1_df<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc3) )%>%
        dplyr::rename(delta_w = 5)%>%
        dplyr::ungroup()
      
      
      pair_plot2<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_num = 5,
                      delta_y = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "ProtID"))%>%
        dplyr::mutate(delta_x = delta_num - delta_w)

      
      volcano2<-ggplot(pair_plot2, aes( delta_x, delta_y,label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.3,size=0.5)+theme_classic()+
        geom_abline(slope=1,linetype="dashed",size=0.5)+
        geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
        geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
        labs(y=paste(fc2),x=paste(fc1,"-",fc3,sep=""))
    }
    
    if (input$ratioratioindicatorVolcano == "Y axis RoR") {
      con7<-as.character(unique(norm_sum_final2()$Condition)[input$denom3])
      con8<-as.character(unique(norm_sum_final2()$Condition)[input$denom4])
      fc4<-as.character(paste("log2FC_",con8,"-",con7,sep=""))

      join1_df<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc4) )%>%
        dplyr::rename(delta_w = 5)%>%
        dplyr::ungroup()
      
      
      pair_plot2<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_x = 5,
                      delta_num = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "ProtID"))%>%
        dplyr::mutate(delta_y = delta_num - delta_w)
      
      
      volcano2<-ggplot(pair_plot2, aes( delta_x, delta_y,label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.3,size=0.5)+theme_classic()+
        geom_abline(slope=1,linetype="dashed",size=0.5)+
        geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
        geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
        labs(y=paste(fc2,"-",fc4,sep=""),x=paste(fc1))
    }
    if (input$ratioratioindicatorVolcano == "RoR same Denominator") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$denom1])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$denom2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))

      join1_df<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc3) )%>%
        dplyr::rename(delta_w = 5)%>%
        dplyr::ungroup()
      
      
      pair_plot2<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_numx = 5,
                      delta_numy = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "ProtID"))%>%
        dplyr::mutate(delta_x = delta_numx - delta_w,
                      delta_y = delta_numy - delta_w)
      
      
      volcano2<-ggplot(pair_plot2, aes( delta_x, delta_y,label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.3,size=0.5)+theme_classic()+
        geom_abline(slope=1,linetype="dashed",size=0.5)+
        geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
        geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
        labs(y=paste(fc2,"-",fc3,sep=""),x=paste(fc1,"-",fc3,sep=""))
      
    }
    if (input$ratioratioindicatorVolcano == "RoR different Denominator") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$denom1])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$denom2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))
      
      con7<-as.character(unique(norm_sum_final2()$Condition)[input$denom3])
      con8<-as.character(unique(norm_sum_final2()$Condition)[input$denom4])
      fc4<-as.character(paste("log2FC_",con8,"-",con7,sep=""))
      
      join1_df<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc3),all_of(fc4) )%>%
        dplyr::rename(delta_w = 5,
                      delta_z = 6)%>%
        dplyr::ungroup()

      pair_plot2<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_numx = 5,
                      delta_numy = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "ProtID"))%>%
        dplyr::mutate(delta_x = delta_numx - delta_w,
                      delta_y = delta_numy - delta_z)
      
      volcano2<-ggplot(pair_plot2, aes( delta_x, delta_y,label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.3,size=0.5)+theme_classic()+
        geom_abline(slope=1,linetype="dashed",size=0.5)+
        geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
        geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
        labs(y=paste(fc2,"-",fc4,sep=""),x=paste(fc1,"-",fc3,sep=""))
    }
    
    ggplotly(volcano2)%>% layout(dragmode = "select")
    
    
    
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  volcanoplotnormprint <- reactive({
    
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse3vol])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse4vol])
    con3<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse5vol])
    con4<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse6vol])
    
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    fc2<-as.character(paste("log2FC_",con4,"-",con3,sep=""))
    
    
    pair_plot<-tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))
    
    
    
    if (input$ratioratioindicatorVolcano == "Single Ratio") {
      
      pair_plot2 <- pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID, all_of(fc1), all_of(fc2) )
      
      volcano2<-ggplot(pair_plot2, aes( .data[[fc1]], .data[[fc2]]))+
        geom_point(alpha=0.3,size=0.5)+theme_classic()+
        geom_abline(slope=1,linetype="dashed",size=0.5)+
        geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
        geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
        labs(x= paste(fc1), y= paste(fc2))
      
    }
    
    if (input$ratioratioindicatorVolcano == "X axis RoR") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$denom1])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$denom2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))
      
      join1_df<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc3) )%>%
        dplyr::rename(delta_w = 5)%>%
        dplyr::ungroup()
      
      
      pair_plot2<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_num = 5,
                      delta_y = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "ProtID"))%>%
        dplyr::mutate(delta_x = delta_num - delta_w)
      
      
      volcano2<-ggplot(pair_plot2, aes( delta_x, delta_y))+
        geom_point(alpha=0.3,size=0.5)+theme_classic()+
        geom_abline(slope=1,linetype="dashed",size=0.5)+
        geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
        geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
        labs(y=paste(fc2),x=paste(fc1,"-",fc3,sep=""))
    }
    
    if (input$ratioratioindicatorVolcano == "Y axis RoR") {
      con7<-as.character(unique(norm_sum_final2()$Condition)[input$denom3])
      con8<-as.character(unique(norm_sum_final2()$Condition)[input$denom4])
      fc4<-as.character(paste("log2FC_",con8,"-",con7,sep=""))
      
      join1_df<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc4) )%>%
        dplyr::rename(delta_w = 5)%>%
        dplyr::ungroup()
      
      
      pair_plot2<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_x = 5,
                      delta_num = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "ProtID"))%>%
        dplyr::mutate(delta_y = delta_num - delta_w)
      
      
      volcano2<-ggplot(pair_plot2, aes( delta_x, delta_y))+
        geom_point(alpha=0.3,size=0.5)+theme_classic()+
        geom_abline(slope=1,linetype="dashed",size=0.5)+
        geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
        geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
        labs(y=paste(fc2,"-",fc4,sep=""),x=paste(fc1))
    }
    if (input$ratioratioindicatorVolcano == "RoR same Denominator") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$denom1])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$denom2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))
      
      join1_df<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc3) )%>%
        dplyr::rename(delta_w = 5)%>%
        dplyr::ungroup()
      
      
      pair_plot2<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_numx = 5,
                      delta_numy = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "ProtID"))%>%
        dplyr::mutate(delta_x = delta_numx - delta_w,
                      delta_y = delta_numy - delta_w)
      
      
      volcano2<-ggplot(pair_plot2, aes( delta_x, delta_y))+
        geom_point(alpha=0.3,size=0.5)+theme_classic()+
        geom_abline(slope=1,linetype="dashed",size=0.5)+
        geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
        geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
        labs(y=paste(fc2,"-",fc3,sep=""),x=paste(fc1,"-",fc3,sep=""))
      
    }
    if (input$ratioratioindicatorVolcano == "RoR different Denominator") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$denom1])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$denom2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))
      
      con7<-as.character(unique(norm_sum_final2()$Condition)[input$denom3])
      con8<-as.character(unique(norm_sum_final2()$Condition)[input$denom4])
      fc4<-as.character(paste("log2FC_",con8,"-",con7,sep=""))
      
      join1_df<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc3),all_of(fc4) )%>%
        dplyr::rename(delta_w = 5,
                      delta_z = 6)%>%
        dplyr::ungroup()
      
      pair_plot2<-pair_plot%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, ProtID,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_numx = 5,
                      delta_numy = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "ProtID"))%>%
        dplyr::mutate(delta_x = delta_numx - delta_w,
                      delta_y = delta_numy - delta_z)
      
      volcano2<-ggplot(pair_plot2, aes( delta_x, delta_y))+
        geom_point(alpha=0.3,size=0.5)+theme_classic()+
        geom_abline(slope=1,linetype="dashed",size=0.5)+
        geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
        geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
        labs(y=paste(fc2,"-",fc4,sep=""),x=paste(fc1,"-",fc3,sep=""))
    }
  
    # con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse3vol])
    # con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse4vol])
    # con3<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse5vol])
    # con4<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse6vol])
    # 
    # 
    # pair_plot<-tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))
    # 
    # volcano2<-ggplot(pair_plot, aes( .data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]))+
    #   geom_point(alpha=0.3,size=0.5)+theme_classic()+
    #   geom_abline(slope=1,linetype="dashed",size=0.5)+
    #   geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
    #   geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
    #   labs(x= paste(con2,"-",con1), y= paste(con4,"-",con3))

    return(volcano2)

  })
  
  output$printcorrtwocond<- downloadHandler(
    filename = function() { paste("CorrPlot_x_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse4vol]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse3vol]),"_y_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse6vol]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse5vol]) ,'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = volcanoplotnormprint(), width = 9, height = 5)
    }
  )

  clusterTraceDF<-reactive({
    if(input$txtclustercols == "Yes"){
      out<-pheatmap::pheatmap(heatmapinput(),   cutree_rows = input$numclusters,main = "ANOVA Significant Proteins",show_rownames=F)
      heatmapinput2<-heatmapinput()%>%
        dplyr::select(colnames(heatmapinput()[,out$tree_col[["order"]]]))
      order_levels<-colnames(heatmapinput()[,out$tree_col[["order"]]])
      
      length_dim<-dim(clusterTraceDFinitial())[[2]]-1
      
      df4.clust3<-clusterTraceDFinitial()%>%
        tidyr::gather("BioReplicate","value",2:all_of(length_dim))
      
      df4.clust3$BioReplicate <- factor(df4.clust3$BioReplicate,levels = c(order_levels))
    }
    
    if(input$txtclustercols == "No"){
      out<-pheatmap::pheatmap(heatmapinput(),  cutree_rows = input$numclusters,main = "ANOVA Significant Proteins",show_rownames=F, cluster_cols = FALSE)
      heatmapinput2<-heatmapinput()
      
      length_dim<-dim(clusterTraceDFinitial())[[2]]-1
      
      df4.clust3<-clusterTraceDFinitial()%>%
        tidyr::gather("BioReplicate","value",2:all_of(length_dim))
      
      df4.clust3$BioReplicate <- factor(df4.clust3$BioReplicate,levels = names(clusterTraceDFinitial())[-1])
  
    }
    
    
    return(df4.clust3)
    
  })

  output$plotclustertrace <- renderPlot({
    
    df_cluster<-clusterTraceDF()%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")%>%
      dplyr::mutate(Condition = stringr::str_sub(BioReplicate,start=1,end=-3))
    
    df_cluster_med <- df_cluster%>%
      dplyr::group_by(BioReplicate,cluster,Condition)%>%
      dplyr::summarise(value = median(value))
      
    
    ggplot()+geom_violin(data=df_cluster,aes(BioReplicate,value,color=Condition),draw_quantiles = c(0.25,0.5,0.75))+
      geom_point(data=df_cluster_med,aes(BioReplicate,value,color=Condition))+
      theme_classic()+labs(x="",y="Scaled Intensity")+
      facet_wrap(vars(cluster),nrow=2)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      geom_hline(yintercept = 0,linetype="dashed")
  })
  
  
  plotclustertraceprint <- reactive({
    
    df_cluster<-clusterTraceDF()%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")%>%
      dplyr::mutate(Condition = stringr::str_sub(BioReplicate,start=1,end=-3))
    
    df_cluster_med <- df_cluster%>%
      dplyr::group_by(BioReplicate,cluster,Condition)%>%
      dplyr::summarise(value = median(value))
    
    
    return(ggplot()+geom_violin(data=df_cluster,aes(BioReplicate,value,color=Condition),draw_quantiles = c(0.25,0.5,0.75))+
             geom_point(data=df_cluster_med,aes(BioReplicate,value,color=Condition))+
      theme_classic()+labs(x="",y="Scaled Intensity")+facet_wrap(vars(cluster),nrow=2)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept = 0,linetype="dashed"))
  })
  
  output$printclustersplot <- downloadHandler(
    filename = function() { paste("ClusterViolins", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotclustertraceprint(), width = 16, height = 10)
    }
  )
  
  

  output$numbercluster<-renderReactable({
    df_clustercount<-clusterTraceDF()%>%
      dplyr::select(ProtID,cluster)%>%
      dplyr::distinct()%>%
      dplyr::group_by(cluster)%>%
      dplyr::summarise(Num_Proteins_perCluster = n())
    reactable(df_clustercount)
  })

  output$bioplexNetwork <-renderPlot({


    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))


    binary_slim <- bioplex_binary%>%
      dplyr::select(Gene.Symbol_target,Gene.Symbol_interactor)%>%
      dplyr::rename(SymbolA = Gene.Symbol_target ,
                    SymbolB = Gene.Symbol_interactor)

    if (input$txthumanmouse == "mouse") {
      results<-tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::rename(name = Gene.Symbol)
    }

    if (input$txthumanmouse == "human") {
      results<-tableresults()%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::rename(name = Gene.Symbol)
    }




    results2<-results%>%
      dplyr::mutate(med_two = (.data[[con1]] + .data[[con2]])/2)%>%
      dplyr::group_by(name)%>%
      dplyr::mutate(max_val = max(med_two))%>%
      dplyr::filter(med_two == max_val)%>%
      dplyr::ungroup()

    ref_name_df<-results2%>%
      dplyr::filter(name == input$txtbioplex)

    to_bind<-data.frame("SymbolA" = input$txtbioplex,"SymbolB" = input$txtbioplex)

    binary_select<-binary_slim%>%
      dplyr::filter(SymbolA == input$txtbioplex)%>%
      dplyr::bind_rows(.,to_bind)

    results3<-results2%>%
      filter(name %in% unique(binary_select$SymbolB))

    binary_select2<-binary_select%>%dplyr::filter(SymbolB %in% unique(results3$name))

    network<- igraph::graph_from_data_frame(binary_select2,directed=F)

    N<-tidygraph::as_tbl_graph(network)

    ggraph::set_graph_style()


    N%>%
      tidygraph::activate(nodes)%>%
      tidygraph::left_join(.,results3,by="name")%>%
      ggraph::ggraph(layout = "stress")+
      ggraph::geom_edge_fan(width=0.5,alpha=0.5)+
      ggraph::geom_node_point(aes(color=.data[[fc1]]),size=25,alpha=0.5)+
      ggraph::geom_node_text(aes(label=name),size=3)+ggraph::scale_color_viridis(end=0.8)


  },width=800,height=800)
  
  


  
  
  
  
  output$bioplexExp<-renderPlot({
    req(input$txtbioplex)
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))

    int_df_ex<-int_df%>%dplyr::select(SymbolA,SymbolB)%>%dplyr::filter(SymbolA == input$txtbioplex | SymbolB == input$txtbioplex )
    bioplex_interactors_target<- unique(append(int_df_ex$SymbolA,int_df_ex$SymbolB))
    plot_value <-paste(input$txtbioplex, " interactors")


    if (input$txthumanmouse == "mouse") {
      target_interactors<- tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)
    }

    if (input$txthumanmouse == "human") {
      target_interactors<- tableresults()%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)
    }

    ggplot()+geom_boxplot(data=target_interactors,aes(plot_value,.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=target_interactors,aes(plot_value,.data[[fc1]],color=.data[[sig1]]),size=5)+
      theme_classic()+labs(x="Protein",y=paste(fc1))+
      theme(axis.text = element_text(size=18),axis.title = element_text(size=20),legend.title=element_text(size=17),
            legend.text=element_text(size=15))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip()
  })
  
  bioplexExpprint <- reactive({
    req(input$txtbioplex)
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))

    int_df_ex<-int_df%>%dplyr::select(SymbolA,SymbolB)%>%dplyr::filter(SymbolA == input$txtbioplex | SymbolB == input$txtbioplex )
    bioplex_interactors_target<- unique(append(int_df_ex$SymbolA,int_df_ex$SymbolB))
    plot_value <-paste(input$txtbioplex, " interactors")


    if (input$txthumanmouse == "mouse") {
      target_interactors<- tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)
    }

    if (input$txthumanmouse == "human") {
      target_interactors<- tableresults()%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)
    }

    return(ggplot()+geom_boxplot(data=target_interactors,aes(plot_value,.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=target_interactors,aes(plot_value,.data[[fc1]],color=.data[[sig1]]),size=5)+
      theme_classic()+labs(x="Protein",y=paste(fc1))+
      theme(axis.text = element_text(size=18),axis.title = element_text(size=20),legend.title=element_text(size=17),
            legend.text=element_text(size=15))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip())
  })
  
  
  output$printbioplexExp <- downloadHandler(
    filename = function() { paste("Bioplex_",input$txtbioplex,"_",as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex2]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex1]) ,'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = bioplexExpprint(), width = 12, height = 8)
    }
  )
  
  
  output$bioplexExpAll<-renderPlot({


    int_df_ex<-int_df%>%
      dplyr::select(SymbolA,SymbolB)%>%
      dplyr::filter(SymbolA == input$txtbioplex | SymbolB == input$txtbioplex )
    bioplex_interactors_target<- unique(append(int_df_ex$SymbolA,int_df_ex$SymbolB))
    plot_value <-paste(input$txtbioplex, " interactors",sep="")

    if (input$txthumanmouse == "human") {
      target_interactors<- tableresults()%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)%>%
        dplyr::mutate(target = paste(input$txtbioplex))
    }

    if (input$txthumanmouse == "mouse") {
      target_interactors<- tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)%>%
        dplyr::mutate(target = paste(input$txtbioplex))
    }



    numcon<-length(unique(norm_sum_final2()$Condition))

    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])

      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))

      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]])
    }

    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))

      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    }

    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))

      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    }

    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))


      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    }

    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))

      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    }

    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      con8<-as.character(unique(norm_sum_final2()$Condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))

      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    }

    size_df <- dim(df)[2]

    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,target)%>%
      tidyr::gather("Condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()


    dftidy<-dftidy%>%
      dplyr::mutate(Condition = str_remove(Condition,pattern=paste("-",con1,sep="")))



    if (numcon  == 3) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3))
    }

    if (numcon  == 4) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4))
    }

    if (numcon  == 5) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5))
    }

    if (numcon  == 6) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6))
    }

    if (numcon  == 7) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7))
    }

    if (numcon  == 8) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }



    #plotname=unique(dftidy$target)

    maxvalue<-abs(max(dftidy$log2_FC))+0.2

    ggplot(dftidy, aes(Condition,log2_FC,color=Condition,fill=Condition))+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      geom_jitter(alpha=0.5,size=5)+theme_classic()+
      geom_hline(yintercept = 0,linetype="dashed")+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),subtitle="Bioplex Foldchange over Conditions",title=paste(plot_value))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=20,hjust=0.5))+
      ylim(c(-maxvalue,maxvalue))


  })
  
  
  bioplexExpAllprint <- reactive({


    int_df_ex<-int_df%>%
      dplyr::select(SymbolA,SymbolB)%>%
      dplyr::filter(SymbolA == input$txtbioplex | SymbolB == input$txtbioplex )
    bioplex_interactors_target<- unique(append(int_df_ex$SymbolA,int_df_ex$SymbolB))
    plot_value <-paste(input$txtbioplex, " interactors",sep="")

    if (input$txthumanmouse == "human") {
      target_interactors<- tableresults()%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)%>%
        dplyr::mutate(target = paste(input$txtbioplex))
    }

    if (input$txthumanmouse == "mouse") {
      target_interactors<- tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)%>%
        dplyr::mutate(target = paste(input$txtbioplex))
    }



    numcon<-length(unique(norm_sum_final2()$Condition))

    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])

      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))

      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]])
    }

    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))

      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    }

    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))

      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    }

    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))


      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    }

    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))

      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    }

    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      con8<-as.character(unique(norm_sum_final2()$Condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))

      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    }

    size_df <- dim(df)[2]

    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,target)%>%
      tidyr::gather("Condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()


    dftidy<-dftidy%>%
      dplyr::mutate(Condition = str_remove(Condition,pattern=paste("-",con1,sep="")))


    if (numcon  == 3) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3))
    }

    if (numcon  == 4) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4))
    }

    if (numcon  == 5) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5))
    }

    if (numcon  == 6) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6))
    }

    if (numcon  == 7) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7))
    }

    if (numcon  == 8) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }



    maxvalue<-abs(max(dftidy$log2_FC))+0.2

    return(ggplot(dftidy, aes(Condition,log2_FC,color=Condition,fill=Condition))+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      geom_jitter(alpha=0.5,size=5)+theme_classic()+
      geom_hline(yintercept = 0,linetype="dashed")+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),subtitle="Bioplex Foldchange over Conditions",title=paste(plot_value))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=20,hjust=0.5))+
      ylim(c(-maxvalue,maxvalue)))


  })
  
  output$printbioplexExpAll <- downloadHandler(
    filename = function() { paste("Bioplex_",input$txtbioplex,"_multicomparison",'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = bioplexExpAllprint(), width = 16, height = 8)
    }
  )
  
  
  output$tablebioplex<-renderReactable({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))

    int_df_ex<-int_df%>%dplyr::select(SymbolA,SymbolB)%>%dplyr::filter(SymbolA == input$txtbioplex | SymbolB == input$txtbioplex )
    bioplex_interactors_target<- unique(append(int_df_ex$SymbolA,int_df_ex$SymbolB))
    plot_value <-paste(input$txtbioplex, " interactors")

    if (input$txthumanmouse == "human") {
      target_interactors<- tableresults()%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)
    }
    if (input$txthumanmouse == "mouse") {
      target_interactors<- tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)
    }



    reactable(target_interactors, filterable = TRUE)
  })
  
###based on gene.symbol now because of losing info with reference
  
  binary_bioplex_join<-reactive({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qval1<-as.character(paste("q.val_",con2,"-",con1,sep=""))
    qval1target<-paste(qval1,"_target",sep="")
    qval1interactor <- paste(qval1,"_interactor",sep="")
    fc1interactor <- paste(fc1,"_interactor",sep="")


    if (input$txthumanmouse == "human") {
      simplify_df<-tableresults()%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::select(Gene.Symbol, .data[[qval1]],.data[[fc1]],.data[[sig1]])
    }
    if (input$txthumanmouse == "mouse") {
      simplify_df<-tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::select(Gene.Symbol, .data[[qval1]],.data[[fc1]],.data[[sig1]])
    }


    # %>%
    #   dplyr::rename(log2FC = log2_foldchange)
    simplify_df_target<- simplify_df%>%
      dplyr::rename(Gene.Symbol_target= Gene.Symbol)
    simplify_df_interactor<- simplify_df%>%
      dplyr::rename(Gene.Symbol_interactor= Gene.Symbol)
    bioplex_ref<-left_join(simplify_df_target,bioplex_binary,by="Gene.Symbol_target")
    bioplex_combine<- left_join(bioplex_ref,simplify_df_interactor,by="Gene.Symbol_interactor",suffix=c("_target","_interactor"))
    bioplex_combine2<-bioplex_combine%>%
      dplyr::filter(!is.na(.data[[qval1target]]),!is.na(.data[[qval1interactor]]))%>%
      dplyr::group_by(Gene.Symbol_target)%>%
      dplyr::mutate(interactomeFC=median(.data[[fc1interactor]]),
                    interactome_count=n())%>%
      dplyr::arrange(dplyr::desc(interactomeFC))
  })
  
  
  output$alldownloadbioplex <- downloadHandler(
    filename = function() {
      paste("all_biolpex_results_proteinLevel",as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex2]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex1]), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(binary_bioplex_join()%>%dplyr::filter(interactome_count >= input$num10001), file, row.names = FALSE)
    }
  )
  
  #binary_bioplex_join
  
  output$tableBinaryBioplex<-renderReactable({
    req(input$num10001)
    bioplex_output<-binary_bioplex_join2()%>%
      dplyr::ungroup()%>%
      dplyr::select(Gene.Symbol_target,interactomeFC,interactome_count)%>%
      dplyr::distinct()%>%
      dplyr::filter(interactome_count >= input$num10001)%>%
      tibble::rownames_to_column("BioplexRank")%>%
      dplyr::mutate(BioplexRank=as.numeric(as.character(BioplexRank)))%>%
      dplyr::arrange(BioplexRank)


    reactable(bioplex_output, filterable = TRUE)
  })
  
  output$tableBinaryBioplex2<-renderReactable({
    req(input$num10002)
    bioplex_output<-binary_bioplex_join()%>%
      dplyr::ungroup()%>%
      dplyr::select(Gene.Symbol_target,interactomeFC,interactome_count)%>%
      dplyr::distinct()%>%
      dplyr::filter(interactome_count >= input$num10002)%>%
      tibble::rownames_to_column("BioplexRank")%>%
      dplyr::mutate(BioplexRank=as.numeric(as.character(BioplexRank)))%>%
      dplyr::arrange(BioplexRank)


    reactable(bioplex_output, filterable = TRUE)
  })
  
  bioplex_binary_final<-reactive({
    req(input$num10001)
    bioplex_output<-binary_bioplex_join()%>%
      dplyr::ungroup()%>%
      dplyr::select(Gene.Symbol_target,interactomeFC,interactome_count)%>%
      dplyr::distinct()%>%
      dplyr::filter(interactome_count >= input$num10001)%>%
      tibble::rownames_to_column("BioplexRank")%>%
      dplyr::mutate(BioplexRank=as.numeric(as.character(BioplexRank)))%>%
      dplyr::arrange(BioplexRank)

    all_bioplex_info<-dplyr::inner_join(bioplex_output,binary_bioplex_join(),by=c("Gene.Symbol_target","interactomeFC","interactome_count"))
  })
  
  binary_bioplex_join2<-reactive({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex3])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex4])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qval1<-as.character(paste("q.val_",con2,"-",con1,sep=""))
    qval1target<-paste(qval1,"_target",sep="")
    qval1interactor <- paste(qval1,"_interactor",sep="")
    fc1interactor <- paste(fc1,"_interactor",sep="")


    if (input$txthumanmouse == "human") {
      simplify_df<-tableresults()%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::select(Gene.Symbol, .data[[qval1]],.data[[fc1]],.data[[sig1]])
    }
    if (input$txthumanmouse == "mouse") {
      simplify_df<-tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::select(Gene.Symbol, .data[[qval1]],.data[[fc1]],.data[[sig1]])
    }


    simplify_df_target<- simplify_df%>%
      dplyr::rename(Gene.Symbol_target= Gene.Symbol)
    simplify_df_interactor<- simplify_df%>%
      dplyr::rename(Gene.Symbol_interactor= Gene.Symbol)
    bioplex_ref<-left_join(simplify_df_target,bioplex_binary,by="Gene.Symbol_target")
    bioplex_combine<- left_join(bioplex_ref,simplify_df_interactor,by="Gene.Symbol_interactor",suffix=c("_target","_interactor"))
    bioplex_combine2<-bioplex_combine%>%
      dplyr::filter(!is.na(.data[[qval1target]]),!is.na(.data[[qval1interactor]]))%>%
      dplyr::group_by(Gene.Symbol_target)%>%
      dplyr::mutate(interactomeFC=median(.data[[fc1interactor]]),
                    interactome_count=n())%>%
      dplyr::arrange(dplyr::desc(interactomeFC))
  })
  
  bioplex_binary_final2<-reactive({
    req(input$num10001)
    bioplex_output<-binary_bioplex_join2()%>%
      dplyr::ungroup()%>%
      dplyr::select(Gene.Symbol_target,interactomeFC,interactome_count)%>%
      dplyr::distinct()%>%
      dplyr::filter(interactome_count >= input$num10001)%>%
      tibble::rownames_to_column("BioplexRank")%>%
      dplyr::mutate(BioplexRank=as.numeric(as.character(BioplexRank)))%>%
      dplyr::arrange(BioplexRank)

    all_bioplex_info<-dplyr::inner_join(bioplex_output,binary_bioplex_join2(),by=c("Gene.Symbol_target","interactomeFC","interactome_count"))
  })
  
    
  output$bioplextopN<-renderPlot({
    req(input$numbioplex1)

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex3])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex4])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qval1<-as.character(paste("q.val_",con2,"-",con1,sep=""))
    qval1target<-paste(qval1,"_target",sep="")
    qval1interactor <- paste(qval1,"_interactor",sep="")
    fc1interactor <- paste(fc1,"_interactor",sep="")
    sig1interactor <- paste(sig1,"_interactor",sep="")

    select_bioplexdf<-bioplex_binary_final2()%>%
      dplyr::filter(BioplexRank <= input$numbioplex1)



    ggplot()+geom_boxplot(data=select_bioplexdf,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]]),outlier.shape=NA)+
      geom_jitter(data=select_bioplexdf,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]],color=.data[[sig1interactor]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Protein Interactome",y=bquote(""~Log[2]~"(Fold Change Interactors)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip()
    })
  
  bioplextopNprint <- reactive({
    req(input$numbioplex1)

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex3])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex4])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qval1<-as.character(paste("q.val_",con2,"-",con1,sep=""))
    qval1target<-paste(qval1,"_target",sep="")
    qval1interactor <- paste(qval1,"_interactor",sep="")
    fc1interactor <- paste(fc1,"_interactor",sep="")
    sig1interactor <- paste(sig1,"_interactor",sep="")

    select_bioplexdf<-bioplex_binary_final2()%>%
      dplyr::filter(BioplexRank <= input$numbioplex1)



    return(ggplot()+geom_boxplot(data=select_bioplexdf,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]]),outlier.shape=NA)+
      geom_jitter(data=select_bioplexdf,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]],color=.data[[sig1interactor]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Protein Interactome",y=bquote(""~Log[2]~"(Fold Change Interactors)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip())

  })
  
  

  output$bioplexbottomN<-renderPlot({
    req(input$numbioplex2)
    num_limit<-input$numbioplex2

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex3])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex4])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qval1<-as.character(paste("q.val_",con2,"-",con1,sep=""))
    qval1target<-paste(qval1,"_target",sep="")
    qval1interactor <- paste(qval1,"_interactor",sep="")
    fc1interactor <- paste(fc1,"_interactor",sep="")
    sig1interactor <- paste(sig1,"_interactor",sep="")


    shrink_group<-bioplex_binary_final2()%>%
      dplyr::arrange(desc(BioplexRank))%>%
      dplyr::select(BioplexRank)%>%
      dplyr::distinct()

    list_ranks<-c(shrink_group$BioplexRank[1:num_limit])

    select_bioplexdf2<-bioplex_binary_final2()%>%
      dplyr::filter(BioplexRank %in% list_ranks)

    ggplot()+geom_boxplot(data=select_bioplexdf2,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]]),outlier.shape=NA)+
      geom_jitter(data=select_bioplexdf2,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]],color=.data[[sig1interactor]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Protein Interactome",y=bquote(""~Log[2]~"(Fold Change Interactors)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip()


  })
  
  bioplexbottomNprint <- reactive({
    req(input$numbioplex2)
    num_limit<-input$numbioplex2

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex3])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex4])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qval1<-as.character(paste("q.val_",con2,"-",con1,sep=""))
    qval1target<-paste(qval1,"_target",sep="")
    qval1interactor <- paste(qval1,"_interactor",sep="")
    fc1interactor <- paste(fc1,"_interactor",sep="")
    sig1interactor <- paste(sig1,"_interactor",sep="")


    shrink_group<-bioplex_binary_final2()%>%
      dplyr::arrange(desc(BioplexRank))%>%
      dplyr::select(BioplexRank)%>%
      dplyr::distinct()

    list_ranks<-c(shrink_group$BioplexRank[1:num_limit])

    select_bioplexdf2<-bioplex_binary_final2()%>%
      dplyr::filter(BioplexRank %in% list_ranks)

    return(ggplot()+geom_boxplot(data=select_bioplexdf2,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]]),outlier.shape=NA)+
      geom_jitter(data=select_bioplexdf2,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]],color=.data[[sig1interactor]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Protein Interactome",y=bquote(""~Log[2]~"(Fold Change Interactors)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip())


  })
  

  
  output$printbiolplextop <- downloadHandler(
    filename = function() { paste("TopN_UpBioplex",as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex4]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex3]), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = bioplextopNprint(), width = 12, height = 8)
    }
  )
  
  output$printbiolplexbottom <- downloadHandler(
    filename = function() { paste("TopN_DownBioplex",as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex4]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numconbioplex3]), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = bioplexbottomNprint(), width = 12, height = 8)
    }
  )
  
  corum_df<-reactive({

      con1<-as.character(unique(norm_sum_final2()$Condition)[input$numcon1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[input$numcon2])
      fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
      medname <- paste("median_log2FC_",con2,"-",con1,sep="")

    corum_binder<-tableresults()%>%
      tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")
    corum_binder<-dplyr::left_join(corum(),corum_binder,by="Reference")
    corum_binder<-corum_binder%>%
      dplyr::group_by(ComplexName)%>%
      dplyr::mutate(total_complexMembers = n())
    corum_binder2<-corum_binder%>%
      dplyr::select(ComplexName,.data[[fc1]],total_complexMembers)%>%
      dplyr::filter(!is.na(.data[[fc1]]))%>%
      dplyr::group_by(ComplexName,total_complexMembers)%>%
      dplyr::summarise(complexMembersIDed=n(),
                !!medname:= median(.data[[fc1]]))%>%
      dplyr::filter(complexMembersIDed>input$numcomplexmin-1)%>%
      dplyr::arrange(desc(.data[[medname]]))%>%
      tibble::rownames_to_column("ComplexRank")%>%
      dplyr::mutate(ComplexRank=as.numeric(as.character(ComplexRank)))%>%
      dplyr::arrange(ComplexRank)
    return(corum_binder2)
  })
  
  complete_corum_df2<-reactive({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    medname <- paste("median_log2FC_",con2,"-",con1,sep="")

    corum_binder<-results_msstats()%>%
      tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")
    corum_binder<-dplyr::left_join(corum(),corum_binder,by="Reference")
    corum_binder<-corum_binder%>%
      dplyr::group_by(ComplexName)%>%
      dplyr::mutate(total_complexMembers = n())
    corum_binder2<-corum_binder%>%
      dplyr::distinct()%>%
      dplyr::select(ComplexName,.data[[fc1]],total_complexMembers)%>%
      dplyr::filter(!is.na(.data[[fc1]]))%>%
      dplyr::group_by(ComplexName,total_complexMembers)%>%
      dplyr::summarise(complexMembersIDed=n(),
                       !!medname:= median(.data[[fc1]]))%>%
      dplyr::filter(complexMembersIDed>input$numcomplexmin-1)%>%
      dplyr::arrange(desc(.data[[medname]]))%>%
      tibble::rownames_to_column("ComplexRank")%>%
      dplyr::mutate(ComplexRank=as.numeric(as.character(ComplexRank)))%>%
      dplyr::arrange(ComplexRank)
    complete_corum_df<-dplyr::left_join(corum_binder,corum_binder2,by=c("ComplexName","total_complexMembers"))%>%
      dplyr::ungroup()%>%
      dplyr::filter(!is.na(Gene.Symbol),!is.na(complexMembersIDed))
    return(complete_corum_df)
  })
  
  
  output$corumNetwork <-renderPlot({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    medname <- paste("median_log2FC_",con2,"-",con1,sep="")

    df_title<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank == input$complexranknum  )
    plot_title<-unique(df_title$ComplexName)
    df_for_network<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank == input$complexranknum  )%>%
      dplyr::mutate(complex="complex")%>%
      dplyr::ungroup()%>%
      dplyr::select(Gene.Symbol,complex)
    network<- igraph::graph_from_data_frame(df_for_network,directed=F)

    N<-tidygraph::as_tbl_graph(network)

    ggraph::set_graph_style()

    df_for_network2<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank == input$complexranknum  )%>%
      dplyr::mutate(complex="complex")%>%
      dplyr::ungroup()%>%
      dplyr::rename(name = Gene.Symbol,
                    !!fc1:= .data[[fc1]])


    N%>%
      tidygraph::activate(nodes)%>%
      tidygraph::left_join(.,df_for_network2,by="name")%>%
      ggraph::ggraph(layout = "stress")+
      ggraph::geom_edge_fan(width=0.5,alpha=0.5)+
      ggraph::geom_node_point(aes(color=.data[[fc1]]),size=25,alpha=0.5)+
      ggraph::geom_node_text(aes(label=name),size=3)+ggraph::scale_color_viridis(end=0.8)+
      ggplot2::ggtitle(paste(plot_title))

  })
  
  output$corumExp<-renderPlot({
    req(input$numcorum)
    choosen_prot<-input$numcorum

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    medname <- paste("median_log2FC_",con2,"-",con1,sep="")

    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank<choosen_prot+1)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))

    total_comp<-unique(select_complex$total_complexMembers)
    id_comp <-unique(select_complex$complexMembersIDed)

    ggplot()+geom_boxplot(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]],color=.data[[sig1]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Complex",y=bquote(""~Log[2]~"(Fold Change)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip()
  })
  
  corumExpprint<-reactive({
    req(input$numcorum)
    choosen_prot<-input$numcorum

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    medname <- paste("median_log2FC_",con2,"-",con1,sep="")

    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank<choosen_prot+1)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))

    # median_delta_value <- unique(select_complex$median_log2FC)
    total_comp<-unique(select_complex$total_complexMembers)
    id_comp <-unique(select_complex$complexMembersIDed)

    return(ggplot()+geom_boxplot(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]],color=.data[[sig1]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Complex",y=bquote(""~Log[2]~"(Fold Change)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip())
  })
  
  
  
  output$corumExp2<-renderPlot({
    req(input$numcorum2)
    num_limit<-input$numcorum2

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))


    shrink_group<-complete_corum_df2()%>%
      dplyr::arrange(desc(ComplexRank))%>%
      dplyr::select(ComplexRank)%>%
      dplyr::distinct()

    list_ranks<-c(shrink_group$ComplexRank[1:num_limit])

    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank %in% list_ranks)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))

    # median_delta_value <- unique(select_complex$median_log2FC)
    total_comp<-unique(select_complex$total_complexMembers)
    id_comp <-unique(select_complex$complexMembersIDed)

    ggplot()+geom_boxplot(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]],color=.data[[sig1]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Complex",y=bquote(""~Log[2]~"(Fold Change)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip()
  })
  
  corumExp2print<-reactive({
    req(input$numcorum2)
    num_limit<-input$numcorum2

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))


    shrink_group<-complete_corum_df2()%>%
      dplyr::arrange(desc(ComplexRank))%>%
      dplyr::select(ComplexRank)%>%
      dplyr::distinct()

    list_ranks<-c(shrink_group$ComplexRank[1:num_limit])

    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank %in% list_ranks)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))

    total_comp<-unique(select_complex$total_complexMembers)
    id_comp <-unique(select_complex$complexMembersIDed)

    return(ggplot()+geom_boxplot(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]],color=.data[[sig1]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Complex",y=bquote(""~Log[2]~"(Fold Change)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip())
  })


  output$printcorumExp <- downloadHandler(
    filename = function() { paste("TopN_UpCorum",as.character(unique(norm_sum_final2()$Condition)[input$numcon2]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numcon1]), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corumExpprint(), width = 12, height = 8)
    }
  )
  
  output$printcorumExp2 <- downloadHandler(
    filename = function() { paste("TopN_DownCorum",as.character(unique(norm_sum_final2()$Condition)[input$numcon2]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numcon1]), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corumExp2print(), width = 12, height = 8)
    }
  )
  
  
  
  output$corumExp3<-renderPlot({
    choosen_prot<-input$numcorum

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))

    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank==input$complexranknum)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))

    total_comp<-unique(select_complex$total_complexMembers)
    id_comp <-unique(select_complex$complexMembersIDed)

    ggplot()+geom_boxplot(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]],color=.data[[sig1]]),size=5)+
      theme_classic()+labs(x="Complex",y=paste(fc1))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip()
  })
  
  corumExp3print <- reactive({
    choosen_prot<-input$numcorum

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))

    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank==input$complexranknum)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))

    total_comp<-unique(select_complex$total_complexMembers)
    id_comp <-unique(select_complex$complexMembersIDed)

    return(ggplot()+geom_boxplot(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]],color=.data[[sig1]]),size=5)+
      theme_classic()+labs(x="Complex",y=paste(fc1))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip())
  })
  
  nameComplex <- reactive({
    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank==input$complexranknum)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))
    return(unique(select_complex$ComplexName_slim))
  })
  
  
  output$printcorumExp3 <- downloadHandler(
    filename = function() { paste("Corum_SingleCompare_",nameComplex(),"_",as.character(unique(norm_sum_final2()$Condition)[input$numcon2]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numcon1]), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corumExp3print(), width = 12, height = 8)
    }
  )
  

  
  output$corumExpallcond<-renderPlot({
    choosen_prot<-input$numcorum

    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank==input$complexranknum)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))


    numcon<-length(unique(norm_sum_final2()$Condition))

    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])

      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))

      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]])
    }

    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))

      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    }

    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))

      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    }

    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))


      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    }

    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))

      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    }

    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      con8<-as.character(unique(norm_sum_final2()$Condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))

      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    }

    size_df <- dim(df)[2]

    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,ComplexName_slim)%>%
      tidyr::gather("Condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()


    dftidy<-dftidy%>%
      dplyr::mutate(Condition = str_remove(Condition,pattern=paste("-",con1,sep="")))



    if (numcon  == 3) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3))
    }

    if (numcon  == 4) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4))
    }

    if (numcon  == 5) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5))
    }

    if (numcon  == 6) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6))
    }

    if (numcon  == 7) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7))
    }

    if (numcon  == 8) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }




    plotname=unique(dftidy$ComplexName_slim)

    maxvalue<-abs(max(dftidy$log2_FC))+0.2

    ggplot(dftidy, aes(Condition,log2_FC,color=Condition,fill=Condition))+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      geom_jitter(alpha=0.5,size=5)+theme_classic()+
      geom_hline(yintercept = 0,linetype="dashed")+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),subtitle="Complex Foldchange over Conditions",title=paste(plotname))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=20,hjust=0.5))+
      ylim(c(-maxvalue,maxvalue))


  })
  
  output$alldownloadCorum <- downloadHandler(
    filename = function() {
      paste("all_corum_results", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(complete_corum_df2(), file, row.names = FALSE)
    }
  )
  
  
  corumExpallcondprint <- reactive({
    choosen_prot<-input$numcorum

    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank==input$complexranknum)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))


    numcon<-length(unique(norm_sum_final2()$Condition))

    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])

      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))

      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]])
    }

    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))

      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    }

    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))

      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    }

    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))


      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    }

    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))

      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    }

    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      con8<-as.character(unique(norm_sum_final2()$Condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))

      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    }

    size_df <- dim(df)[2]

    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,ComplexName_slim)%>%
      tidyr::gather("Condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()


    dftidy<-dftidy%>%
      dplyr::mutate(Condition = str_remove(Condition,pattern=paste("-",con1,sep="")))



    if (numcon  == 3) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3))
    }

    if (numcon  == 4) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4))
    }

    if (numcon  == 5) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5))
    }

    if (numcon  == 6) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6))
    }

    if (numcon  == 7) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7))
    }

    if (numcon  == 8) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }


    plotname=unique(dftidy$ComplexName_slim)

    maxvalue<-abs(max(dftidy$log2_FC))+0.2

    return(ggplot(dftidy, aes(Condition,log2_FC,color=Condition,fill=Condition))+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      geom_jitter(alpha=0.5,size=5)+theme_classic()+
      geom_hline(yintercept = 0,linetype="dashed")+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),subtitle="Complex Foldchange over Conditions",title=paste(plotname))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=20,hjust=0.5))+
      ylim(c(-maxvalue,maxvalue)))


  })
  
  output$printcorumExpallcond <- downloadHandler(
    filename = function() { paste("Corum_MultiComparison_",nameComplex(), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corumExpallcondprint(), width = 16, height = 8)
    }
  )
  
  output$tablecorum<-renderReactable({
    reactable(corum_df(), filterable = TRUE)
  })
  
  output$tablecorumOneComplex<-renderReactable({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qname1<-paste("q.val_",con2,"-",con1,sep="")
    medname <- paste("median_log2FC_",con2,"-",con1,sep="")

    df_out<-complete_corum_df2()%>%
      dplyr::select(ComplexRank,ComplexName,.data[[medname]],total_complexMembers,complexMembersIDed,Gene.Symbol,.data[[sig1]],Annotation,.data[[fc1]],.data[[qname1]])%>%
      dplyr::filter(ComplexRank == input$complexranknum)


    reactable(df_out, filterable = TRUE)
  })
  
  observeEvent(input$runGOcluster, {
  
  
    go_df<-reactive({
      filt_ttest2<-tableresults()%>%
        dplyr::ungroup()%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::select(Reference,cluster)%>%
        dplyr::mutate(cluster = as.numeric(as.character(cluster)))%>%
        dplyr::mutate(set = if_else(cluster == input$numcluster, "sig","n.s."))
  
  
  
      df.bg <- filt_ttest2 %>%
        plotly::distinct(Reference) %>%
        dplyr::left_join(go.df()) %>%
        dplyr::mutate(all = n_distinct(Reference)) %>%
        dplyr::group_by(GO,GO_cat) %>%
        dplyr::summarise("pos" = n_distinct(Reference),
                  "neg" = all - pos) %>%
        dplyr::distinct() %>%
        dplyr::mutate(selection="background") %>%
        plotly::ungroup()
      # print("first")
      # print(head(df.bg))
      go.df <- df.bg %>%
        dplyr::select(GO,GO_cat) %>%
        dplyr::left_join(go.df())
      # print("second")
      # print(head(go.df))
      df.fg <- filt_ttest2 %>%
        dplyr::mutate(selection = "foreground") %>%
        dplyr::distinct(Reference, selection, set) %>%
        dplyr::right_join(go.df) %>%
        dplyr::group_by(selection, set) %>%
        dplyr::mutate(term_all = n_distinct(Reference)) %>%
        dplyr::group_by(GO,GO_cat, selection, set) %>%
        dplyr::summarise("pos" = n_distinct(Reference),
                  "neg" = term_all - pos) %>%
        plotly::ungroup() %>%
        dplyr::filter(set == "sig")%>%
        dplyr::distinct()
      # print("third")
      # print(head(df.fg))
      df.nest <- df.bg %>% dplyr::left_join(df.fg %>% dplyr::distinct(GO,GO_cat, set)) %>%
        dplyr::distinct() %>%
        dplyr::full_join(df.fg) %>%
        tidyr::drop_na(set) %>%
        tidyr::drop_na(GO) %>%
        dplyr::distinct() %>%
        dplyr::group_by(GO,GO_cat, set) %>%
        tidyr::nest()
      # print("fourth")
      # print(head(df.nest))
  
      exact_pval_fun<-function(data){
        dataframe<-data.frame(data) %>% tibble::column_to_rownames(var = "selection")
        matrix<-as.matrix(dataframe)
        exact_test_result<-fisher.test(matrix, alternative = "two.sided")
        return(exact_test_result$p.value)
      }
  
      ## function to perform fisher test and get enrichment
      enrich_fun<-function(data){
        dataframe<-data.frame(data) %>%
          dplyr::mutate(frac = pos/neg) %>%
          dplyr::select(-c(pos,neg)) %>%
          tidyr::pivot_wider(names_from = selection, values_from = frac) %>%
          dplyr::summarise(enrichment = foreground/background)
  
        return(dataframe$enrichment[1])
      }
  
  
      fish.df <- df.nest %>%
        dplyr::mutate(p.val = data %>% purrr::map_dbl(exact_pval_fun)) %>% #calc p val
        dplyr::mutate(enrich = data %>% purrr::map_dbl(enrich_fun)) %>%   #calc enrichment
        dplyr::select(-c(data)) %>%
        dplyr::filter(GO != "")%>%
        dplyr::arrange(p.val)
      # print("fifth")
  
  
      fish.df$p.val.adj = p.adjust(fish.df$p.val,method="fdr")
  
      fish.df<-fish.df%>%
        dplyr::arrange(p.val.adj)%>%
        tibble::rownames_to_column("GORank")%>%
        dplyr::mutate(GORank=as.numeric(as.character(GORank)))
      # write.csv(fish.df,"C:/Harper/Side_analysis/go_output_test.csv")
      # print("sixth")
      # fish.df<-fish.df%>%
      #   dplyr::mutate(p.val.adj = p.adjust(p.val, method = "fdr"))
  
      fish.df.sig<-fish.df%>%
        dplyr::mutate(sig_GO = if_else(p.val.adj<input$numgo,"sig","n.s."))%>%
        dplyr::mutate(GO_id = str_sub(GO,-11,-2),
               GO_shrink1 = str_sub(GO,1,-13),
               GO_shrink = str_sub(GO_shrink1,1,70))
  
      fish.df.sig$GO_cat <- factor(fish.df.sig$GO_cat, levels = c("molecular_function","cellular_component","biological_process"))
      # print(head(fish.df.sig))
      # write.csv(fish.df,"C:/Harper/Side_analysis/go_output_test2.csv")
      return(fish.df.sig)
    })
    
    
    output$tableGO2<-renderReactable({
      reactable(goExport()%>%
                  dplyr::ungroup()%>%
                  dplyr::select(GORank,GO,GO_cat,enrich,p.val.adj)%>%
                  dplyr::distinct(), filterable = TRUE)
    })

    
    goExport<-reactive({
      filt_ttest2<-tableresults()%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")
      
      
      #dplyr::select(Reference, Gene.Symbol,q.val,log2_foldchange,Annotation,Significant)
      annotate_go<-dplyr::inner_join(go.df()%>%dplyr::select(-Gene.names...primary..),filt_ttest2,by="Reference")
      go_output<-dplyr::inner_join(go_df(),annotate_go,by=c("GO","GO_cat"),suffix=c("_go","_protein"))
      return(go_output)
    })

    
    output$tableGO<-renderReactable({
      reactable(goExport()%>%dplyr::ungroup(),filterable = TRUE)
    })
    
    output$GOclusterdownload <- downloadHandler(
      filename = function() {
        paste("All_GO_results_Cluster",input$numcluster, ".csv", sep = "")
      },
      content = function(file) {
        write.csv(goExport(), file, row.names = FALSE)
      }
    )
    
    output$plotGO <- renderPlotly({
      p10101<-ggplot(go_df(),  aes(log2(enrich), -log10(p.val.adj), color=sig_GO,label=GO_id,label1=GO_shrink,key=GORank))+geom_point(alpha=0.75)+theme_classic()+
        geom_hline(yintercept = -log10(input$numgo))+scale_color_viridis_d(end=0.8)+
        labs(x="log2(Enrichment)", y="-log10(q-value)",title=paste("Cluster ",input$numcluster, sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18))
      ggplotly(p10101)%>% layout(dragmode = "select")
      
    })
    
    plotGOprint <- reactive({
      return(ggplot(go_df(),  aes(log2(enrich), -log10(p.val.adj), color=sig_GO))+geom_point(alpha=0.75)+theme_classic()+
               geom_hline(yintercept = -log10(input$numgo))+scale_color_viridis_d(end=0.8)+
               labs(x="log2(Enrichment)", y="-log10(q-value)",title=paste("Cluster ",input$numcluster, sep=""))+
               theme(axis.title = element_text(size=20),axis.text = element_text(size = 18)))
    })
    
    output$printplotGO <- downloadHandler(
      filename = function() { paste("GO_plot_AllTerms_Cluster",input$numcluster, '.pdf', sep='') },
      content = function(file) {
        ggsave(file, plot = plotGOprint(), width =8, height = 8)
      }
    )

    output$plotGOuser <-renderPlot({
      req(input$txtlistnum)
      req(input$numgo2)
      
      list_index_ranks<-unlist(strsplit(input$txtlistnum, ","))
      
      ggplot(goExport()%>%filter(GORank %in% list_index_ranks),  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat))+
        geom_point()+theme_classic()+coord_flip()+
        scale_color_viridis_d(end=0.8)+
        labs(x="GO terms", y="-log10(q-value)",title=paste("Cluster ",input$numcluster, sep=""))+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 12))+
        guides(colour = guide_legend(override.aes = list(size=3)))+geom_hline(yintercept = -log10(input$numgo2),linetype="dashed")
      
      
    })
    
    plotGOuserprint <- reactive({
      req(input$txtlistnum)
      req(input$numgo2)
      
      list_index_ranks<-unlist(strsplit(input$txtlistnum, ","))
      
      return(ggplot(goExport()%>%filter(GORank %in% list_index_ranks),  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat))+
               geom_point()+theme_classic()+coord_flip()+
               scale_color_viridis_d(end=0.8)+
               labs(x="GO terms", y="-log10(q-value)",title=paste("Cluster ",input$numcluster, sep=""))+
               theme(axis.title = element_text(size=16),axis.text = element_text(size = 12))+
               guides(colour = guide_legend(override.aes = list(size=3)))+geom_hline(yintercept = -log10(input$numgo2),linetype="dashed"))
      
      
    })
    
    
    output$printplotGOuser <- downloadHandler(
      filename = function() { paste("GO_plot_UserTerms_Cluster",input$numcluster, '.pdf', sep='') },
      content = function(file) {
        ggsave(file, plot = plotGOuserprint(), width = 12, height = 8)
      }
    )
    
    
    
    output$plotGOtopn2 <-renderPlotly({
      
      p333<-ggplot(go_df()[1:input$numtopn,],  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat,key=GORank))+
        geom_point()+theme_classic()+coord_flip()+
        scale_color_viridis_d(end=0.8)+
        labs(x="GO terms", y="-log10(q-value)",title=paste("Cluster ",input$numcluster, sep=""))+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 12))
      # +
      #   guides(colour = guide_legend(override.aes = list(size=3)))
      ggplotly(p333)%>% layout(dragmode = "select")
      
    })
    
    plotGOtopn2print <-reactive({
      #,key=GORank
      return(ggplot(go_df()[1:input$numtopn,],  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat))+
               geom_point()+theme_classic()+coord_flip()+
               scale_color_viridis_d(end=0.8)+
               labs(x="GO terms", y="-log10(q-value)",title=paste("Cluster ",input$numcluster, sep=""))+
               theme(axis.title = element_text(size=16),axis.text = element_text(size = 12)))
      
    })
    
    output$printplotGOtopn2 <- downloadHandler(
      filename = function() { paste("GO_plot_Top",input$numtopn,"_Cluster",input$numcluster,"_GOterms", '.pdf', sep='') },
      content = function(file) {
        ggsave(file, plot = plotGOtopn2print(), width = 12, height = 8)
      }
    )
    
    
    
    
    output$goclick <- renderReactable({
      d <- event_data("plotly_click")
      req(d)
      df_sub2<-goExport()%>%dplyr::filter(GORank %in% d$key )%>%
        dplyr::ungroup()
      
      reactable(df_sub2, filterable = TRUE)
    })
    
    
  })
  
  ### ALL of GO regulation code

  observeEvent(input$runGOreg, {
    go_df2<-reactive({
      con1<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority2])
  
      sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
  
      filt_ttest<-tableresults()%>%
        dplyr::ungroup()%>%
        #tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::select(Reference,.data[[sig_name]])
  
      if (input$txtupdown == "up") {
        filt_ttest2 <- filt_ttest%>%
          dplyr::mutate(set = if_else(.data[[sig_name]] == "up", "sig","n.s."))
      }
  
      if (input$txtupdown == "down") {
        filt_ttest2 <- filt_ttest%>%
          dplyr::mutate(set = if_else(.data[[sig_name]] == "down", "sig","n.s."))
      }
  
      if (input$txtupdown == "both") {
        filt_ttest2 <- filt_ttest%>%
          dplyr::mutate(set = if_else(.data[[sig_name]] == "down", "sig",if_else(.data[[sig_name]] == "up","sig","n.s.")))
      }
  
  
  
  
  
      df.bg <- filt_ttest2 %>%
        plotly::distinct(Reference) %>%
        dplyr::left_join(go.df()) %>%
        dplyr::mutate(all = n_distinct(Reference)) %>%
        dplyr::group_by(GO,GO_cat) %>%
        dplyr::summarise("pos" = n_distinct(Reference),
                         "neg" = all - pos) %>%
        dplyr::distinct() %>%
        dplyr::mutate(selection="background") %>%
        plotly::ungroup()
      go.df <- df.bg %>%
        dplyr::select(GO,GO_cat) %>%
        dplyr::left_join(go.df())
      df.fg <- filt_ttest2 %>%
        dplyr::mutate(selection = "foreground") %>%
        dplyr::distinct(Reference, selection, set) %>%
        dplyr::right_join(go.df) %>%
        dplyr::group_by(selection, set) %>%
        dplyr::mutate(term_all = n_distinct(Reference)) %>%
        dplyr::group_by(GO,GO_cat, selection, set) %>%
        dplyr::summarise("pos" = n_distinct(Reference),
                         "neg" = term_all - pos) %>%
        plotly::ungroup() %>%
        dplyr::filter(set == "sig")%>%
        dplyr::distinct()
      df.nest <- df.bg %>% dplyr::left_join(df.fg %>% dplyr::distinct(GO,GO_cat, set)) %>%
        dplyr::distinct() %>%
        dplyr::full_join(df.fg) %>%
        tidyr::drop_na(set) %>%
        tidyr::drop_na(GO) %>%
        dplyr::distinct() %>%
        dplyr::group_by(GO,GO_cat, set) %>%
        tidyr::nest()
  
      exact_pval_fun<-function(data){
        dataframe<-data.frame(data) %>% tibble::column_to_rownames(var = "selection")
        matrix<-as.matrix(dataframe)
        exact_test_result<-fisher.test(matrix, alternative = "two.sided")
        return(exact_test_result$p.value)
      }
  
      ## function to perform fisher test and get enrichment
      enrich_fun<-function(data){
        dataframe<-data.frame(data) %>%
          dplyr::mutate(frac = pos/neg) %>%
          dplyr::select(-c(pos,neg)) %>%
          tidyr::pivot_wider(names_from = selection, values_from = frac) %>%
          dplyr::summarise(enrichment = foreground/background)
  
        return(dataframe$enrichment[1])
      }
  
  
      fish.df <- df.nest %>%
        dplyr::mutate(p.val = data %>% purrr::map_dbl(exact_pval_fun)) %>% #calc p val
        dplyr::mutate(enrich = data %>% purrr::map_dbl(enrich_fun)) %>%   #calc enrichment
        dplyr::select(-c(data)) %>%
        dplyr::filter(GO != "")%>%
        dplyr::arrange(p.val)
  
  
      fish.df$p.val.adj = p.adjust(fish.df$p.val,method="fdr")
  
      fish.df<-fish.df%>%
        dplyr::arrange(p.val.adj)%>%
        tibble::rownames_to_column("GORank")%>%
        dplyr::mutate(GORank=as.numeric(as.character(GORank)))
  
      # fish.df<-fish.df%>%
      #   dplyr::mutate(p.val.adj = p.adjust(p.val, method = "fdr"))
  
      fish.df.sig<-fish.df%>%
        dplyr::mutate(sig_GO = if_else(p.val.adj<input$numgo,"sig","n.s."))%>%
        dplyr::mutate(GO_id = str_sub(GO,-11,-2),
                      GO_shrink1 = str_sub(GO,1,-13),
                      GO_shrink = str_sub(GO_shrink1,1,70))
  
      fish.df.sig$GO_cat <- factor(fish.df.sig$GO_cat, levels = c("molecular_function","cellular_component","biological_process"))
      return(fish.df.sig)
    })
    
    output$plotGOtopn3 <-renderPlotly({
      con1<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority2])
      
      sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
      
      p333<-ggplot(go_df2()[1:input$numtopn1,],  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat,key=GORank))+
        geom_point()+theme_classic()+coord_flip()+
        scale_color_viridis_d(end=0.8)+
        labs(x="GO terms", y="-log10(q-value)", title = paste(sig_name,"_",input$txtupdown,sep = ""))+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 12))
      
      ggplotly(p333)%>% layout(dragmode = "select")
      
    })
    
    plotGOtopn3print <-reactive({
      con1<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority2])
      
      sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
      #,key=GORank
      return(ggplot(go_df2()[1:input$numtopn1,],  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat))+
               geom_point()+theme_classic()+coord_flip()+
               scale_color_viridis_d(end=0.8)+
               labs(x="GO terms", y="-log10(q-value)", title = paste(sig_name,"_",input$txtupdown,sep = ""))+
               theme(axis.title = element_text(size=16),axis.text = element_text(size = 12)))
      
    })
    
    output$printplotGOtopn3 <- downloadHandler(
      
      filename = function() { paste("GO_plot_Top",input$numtopn1,"_Regulated_",as.character(unique(norm_sum_final2()$Condition)[input$numpriority2]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numpriority1]),"_",input$txtupdown, '.pdf', sep='') },
      content = function(file) {
        ggsave(file, plot = plotGOtopn3print(), width = 12, height = 8)
      }
    )
    
    output$plotGO2 <- renderPlotly({
      con1<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority2])
      
      sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
      
      p10101<-ggplot(go_df2(),  aes(log2(enrich), -log10(p.val.adj), color=sig_GO,label=GO_id,label1=GO_shrink,key=GORank))+geom_point(alpha=0.75)+theme_classic()+
        geom_hline(yintercept = -log10(input$numgo))+scale_color_viridis_d(end=0.8)+
        labs(x="log2(Enrichment)", y="-log10(q-value)",title = paste(sig_name,"_",input$txtupdown,sep = ""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18))
      ggplotly(p10101)%>% layout(dragmode = "select")
      
    })
    
    goExport2<-reactive({
      filt_ttest2<-tableresults()%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")
      
      annotate_go<-dplyr::inner_join(go.df()%>%dplyr::select(-Gene.names...primary..),filt_ttest2,by="Reference")
      go_output<-dplyr::inner_join(go_df2(),annotate_go,by=c("GO","GO_cat"),suffix=c("_go","_protein"))
      return(go_output)
    })
    
    output$tableGO2reg<-renderReactable({
      reactable(goExport2()%>%
                  dplyr::ungroup()%>%
                  dplyr::select(GORank,GO,GO_cat,enrich,p.val.adj)%>%
                  dplyr::distinct(), filterable = TRUE)
    })
    
    output$tableGOreg<-renderReactable({
      reactable(goExport2()%>%dplyr::ungroup(),filterable = TRUE)
    })
    
    output$GOregulateddownload <- downloadHandler(
      filename = function() {
        paste("All_GO_results_Regulated_",as.character(unique(norm_sum_final2()$Condition)[input$numpriority2]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numpriority1]),"_",input$txtupdown, ".csv", sep = "")
      },
      content = function(file) {
        write.csv(goExport2(), file, row.names = FALSE)
      }
    )
    
    plotGO2print <- reactive({
      con1<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority2])
      
      sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
      return(ggplot(go_df2(),  aes(log2(enrich), -log10(p.val.adj), color=sig_GO))+geom_point(alpha=0.75)+theme_classic()+
               geom_hline(yintercept = -log10(input$numgo))+scale_color_viridis_d(end=0.8)+
               labs(x="log2(Enrichment)", y="-log10(q-value)",title = paste(sig_name,"_",input$txtupdown,sep = ""))+
               theme(axis.title = element_text(size=20),axis.text = element_text(size = 18)))
      
    })
    output$printplotGO2 <- downloadHandler(
      filename = function() { paste("GO_plot_AllTerms_Regulated_",as.character(unique(norm_sum_final2()$Condition)[input$numpriority2]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numpriority1]),"_",input$txtupdown, '.pdf', sep='') },
      content = function(file) {
        ggsave(file, plot = plotGO2print(), width = 8, height = 8)
      }
    )
    
    output$plotGOuserReg <-renderPlot({
      req(input$txtlistnum3)
      req(input$numgo3)
      
      list_index_ranks2<-unlist(strsplit(input$txtlistnum3, ","))
      
      con1<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority2])
      
      sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
      
      ggplot(goExport2()%>%filter(GORank %in% list_index_ranks2),  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat))+
        geom_point()+theme_classic()+coord_flip()+
        scale_color_viridis_d(end=0.8)+
        labs(x="GO terms", y="-log10(q-value)",title = paste(sig_name,"_",input$txtupdown,sep = ""))+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 12))+
        guides(colour = guide_legend(override.aes = list(size=3)))+geom_hline(yintercept = -log10(input$numgo3),linetype="dashed")
      
      
    })
    
    plotGOuserRegprint <-reactive({
      req(input$txtlistnum3)
      req(input$numgo3)
      
      list_index_ranks2<-unlist(strsplit(input$txtlistnum3, ","))
      
      con1<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[input$numpriority2])
      
      sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
      
      return(ggplot(goExport2()%>%filter(GORank %in% list_index_ranks2),  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat))+
               geom_point()+theme_classic()+coord_flip()+
               scale_color_viridis_d(end=0.8)+
               labs(x="GO terms", y="-log10(q-value)",title = paste(sig_name,"_",input$txtupdown,sep = ""))+
               theme(axis.title = element_text(size=16),axis.text = element_text(size = 12))+
               guides(colour = guide_legend(override.aes = list(size=3)))+geom_hline(yintercept = -log10(input$numgo3),linetype="dashed"))
      
      
    })
    
    output$printplotGOuserReg <- downloadHandler(
      
      filename = function() { paste("GO_plot_UserTerms_Regulated_",as.character(unique(norm_sum_final2()$Condition)[input$numpriority2]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numpriority1]),"_",input$txtupdown, '.pdf', sep='') },
      content = function(file) {
        ggsave(file, plot = plotGOuserRegprint(), width = 12, height = 8)
      }
    )
    
    output$goclick3 <- renderReactable({
      d <- event_data("plotly_click")
      req(d)
      df_sub2<-goExport2()%>%dplyr::filter(GORank %in% d$key )%>%
        dplyr::ungroup()
      
      reactable(df_sub2, filterable = TRUE)
    })
    
    
    
  })
  
  
  
  ## first one
 
  
  


  
  output$corrplotMeans<-renderPlotly({
    x_axis = unique(norm_sum_final2()$Condition)[1]
    y_axis = unique(norm_sum_final2()$Condition)[2]


    plot_df<-results_msstats()



    plot_df<-plot_df%>%
      dplyr::rename(!!unique(norm_sum_final2()$Condition)[1]:=cond1,!!unique(norm_sum_final2()$Condition)[2]:=cond2)%>%
      dplyr::select(-data)
    x_axis<-unique(norm_sum_final2()$Condition)[1]
    y_axis<-unique(norm_sum_final2()$Condition)[2]



    p4021<-ggplot(data=plot_df, aes(.data[[unique(norm_sum_final2()$Condition)[1]]],.data[[unique(norm_sum_final2()$Condition)[2]]],label=Gene.Symbol,color=Significant))+
      geom_point(alpha=input$alphanum1)+theme_classic()+
      geom_abline(slope = 1,color="turquoise")+
      labs(x=x_axis, y=y_axis)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18))+scale_color_viridis_d(end=0.8)
    ggplotly(p4021)
  })
  
  
  
  output$rankfoldchange<-renderPlot({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1vol])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2vol])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    
    if (input$ratioratioindicatorRank == "Single Ratio") {
      sortedDF<-tableresults()%>%
        dplyr::select(Gene.Symbol,Reference,all_of(foldchangevalue))%>%
        dplyr::distinct()%>%
        dplyr::arrange(desc(.data[[foldchangevalue]]))%>%
        tibble::rownames_to_column("colnum")%>%
        dplyr::mutate(colnum=as.numeric(as.character(colnum)))
      
      top_n<-sortedDF$Reference[1:input$num55]
      last_num<-max(sortedDF$colnum)
      last_num_min<-last_num-input$num55+1
      bottom_n<-sortedDF$Reference[last_num_min:last_num]
      
      
      pRank1<-ggplot(data=sortedDF, aes(colnum,.data[[foldchangevalue]]))+
        geom_point()+
        geom_text_repel(data=sortedDF%>%dplyr::mutate(label=dplyr::if_else(Reference %in% top_n,Gene.Symbol,dplyr::if_else(Reference %in% bottom_n,Gene.Symbol,""))),aes(label=label),max.overlaps = Inf,box.padding = 1)+
        theme_classic()+labs(x="Protein Rank",y=paste(foldchangevalue))
    }
    
    if (input$ratioratioindicatorRank == "RoR") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$denomRank1])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$denomRank2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))
      
      sortedDF<-tableresults()%>%
        dplyr::select(Gene.Symbol,Reference,all_of(foldchangevalue),all_of(fc3))%>%
        dplyr::distinct()%>%
        dplyr::rename(delta_num = 3,
                      delta_denom = 4 )%>%
        dplyr::mutate(RoR = delta_num - delta_denom)%>%
        dplyr::arrange(desc(RoR))%>%
        tibble::rownames_to_column("colnum")%>%
        dplyr::mutate(colnum=as.numeric(as.character(colnum)))
      
      top_n<-sortedDF$Reference[1:input$num55]
      last_num<-max(sortedDF$colnum)
      last_num_min<-last_num-input$num55+1
      bottom_n<-sortedDF$Reference[last_num_min:last_num]
      
      pRank1<-ggplot(data=sortedDF, aes(colnum,RoR))+
        geom_point()+
        geom_text_repel(data=sortedDF%>%dplyr::mutate(label=dplyr::if_else(Reference %in% top_n,Gene.Symbol,dplyr::if_else(Reference %in% bottom_n,Gene.Symbol,""))),aes(label=label),max.overlaps = Inf,box.padding = 1)+
        theme_classic()+labs(x="Protein Rank",y=paste(foldchangevalue," - ",fc3,sep=""))
    }
    
    pRank1
    
  })
  
  rankfoldchangeprint<-reactive({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1vol])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2vol])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    
    if (input$ratioratioindicatorRank == "Single Ratio") {
      sortedDF<-tableresults()%>%
        dplyr::select(Gene.Symbol,Reference,all_of(foldchangevalue))%>%
        dplyr::distinct()%>%
        dplyr::arrange(desc(.data[[foldchangevalue]]))%>%
        tibble::rownames_to_column("colnum")%>%
        dplyr::mutate(colnum=as.numeric(as.character(colnum)))
      
      top_n<-sortedDF$Reference[1:input$num55]
      last_num<-max(sortedDF$colnum)
      last_num_min<-last_num-input$num55+1
      bottom_n<-sortedDF$Reference[last_num_min:last_num]
      
      
      pRank1<-ggplot(data=sortedDF, aes(colnum,.data[[foldchangevalue]]))+
        geom_point()+
        geom_text_repel(data=sortedDF%>%dplyr::mutate(label=dplyr::if_else(Reference %in% top_n,Gene.Symbol,dplyr::if_else(Reference %in% bottom_n,Gene.Symbol,""))),aes(label=label),max.overlaps = Inf,box.padding = 1)+
        theme_classic()+labs(x="Protein Rank",y=paste(foldchangevalue))
    }
    
    if (input$ratioratioindicatorRank == "RoR") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$denomRank1])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$denomRank2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))
      
      sortedDF<-tableresults()%>%
        dplyr::select(Gene.Symbol,Reference,all_of(foldchangevalue),all_of(fc3))%>%
        dplyr::distinct()%>%
        dplyr::rename(delta_num = 3,
                      delta_denom = 4 )%>%
        dplyr::mutate(RoR = delta_num - delta_denom)%>%
        dplyr::arrange(desc(RoR))%>%
        tibble::rownames_to_column("colnum")%>%
        dplyr::mutate(colnum=as.numeric(as.character(colnum)))
      
      top_n<-sortedDF$Reference[1:input$num55]
      last_num<-max(sortedDF$colnum)
      last_num_min<-last_num-input$num55+1
      bottom_n<-sortedDF$Reference[last_num_min:last_num]
      
      pRank1<-ggplot(data=sortedDF, aes(colnum,RoR))+
        geom_point()+
        geom_text_repel(data=sortedDF%>%dplyr::mutate(label=dplyr::if_else(Reference %in% top_n,Gene.Symbol,dplyr::if_else(Reference %in% bottom_n,Gene.Symbol,""))),aes(label=label),max.overlaps = Inf,box.padding = 1)+
        theme_classic()+labs(x="Protein Rank",y=paste(foldchangevalue," - ",fc3,sep=""))
    }
    return(pRank1)
  })
  
  output$rankfoldchangeprint2 <- downloadHandler(
    filename = function() { paste("Top_rankedUpDownFC_",input$ratioratioindicatorRank, '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = rankfoldchangeprint(), width = 13, height = 10)
    }
  )
  
  

  output$plot6 <- renderPlotly({
    x_axis <- paste(unique(norm_sum_final2()$Condition)[2], "-", unique(norm_sum_final2()$Condition)[1],sep=" ")

    p1<-ggplot(data=results_msstats(), aes(log2_foldchange, -log2(q.val),color=Significant, label=Gene.Symbol,key=Reference))+
      geom_point(alpha=input$alphanum1)+theme_classic()+
      geom_vline(xintercept = log2(input$num4))+geom_vline(xintercept = -log2(input$num4))+
      geom_hline(yintercept = -log2(input$num3))+scale_color_viridis_d(end=0.8)+
      labs(x=paste("log2(",x_axis,")"), y="-log2(q-value)")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18))


    ggplotly(p1) %>% layout(dragmode = "select")

  })

  output$brush <- renderReactable({
    d <- event_data("plotly_selected")
    req(d)
    df_sub<-results_msstats()%>%dplyr::filter(Reference %in% d$key )%>%
      dplyr::rename(!!unique(norm_sum_final2()$Condition)[1]:=cond1,!!unique(norm_sum_final2()$Condition)[2]:=cond2)%>%
      dplyr::select(-data)

    reactable(df_sub, filterable = TRUE)
  })

  volcanoPlotInput <- reactive({
    x_axis <- paste(unique(norm_sum_final2()$Condition)[2], "-", unique(norm_sum_final2()$Condition)[1],sep=" ")

    p771<-ggplot(data=results_msstats(), aes(log2_foldchange, -log2(q.val),color=Significant))+
      geom_point(alpha=input$alphanum1)+theme_classic()+
      geom_vline(xintercept = log2(input$num4))+geom_vline(xintercept = -log2(input$num4))+
      geom_hline(yintercept = -log2(input$num3))+scale_color_viridis_d(end=0.8)+
      labs(x=paste("log2(",x_axis,")"), y="-log2(q-value)")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18))
  })
  
  inputvennsigcols <- reactive({
    
    req(input$txtlistsig)
    sig_match_input<-norm_sum_final2()%>%
      dplyr::ungroup()%>%
      dplyr::distinct(Condition,priority)
    sig_match_input$Condition<-reorder(sig_match_input$Condition,sig_match_input$priority,min)
    sig_match_input<-sig_match_input%>%
      mutate(priority=as.numeric(as.character(priority)))
    
    sig_match1<-sig_match_input%>%
      dplyr::rename(first=priority)
    sig_match2<-sig_match_input%>%
      dplyr::rename(second=priority)
    
    value<-stringr::str_remove_all(input$txtlistsig," ")
    
    list_index_ranks<-data.frame("contrasts"=unlist(strsplit(value, "\\;")))%>%
      tidyr::separate(contrasts, into=c("first","second"),sep="\\,")%>%
      dplyr::mutate(first=as.numeric(as.character(first)),
                    second=as.numeric(as.character(second)))%>%
      dplyr::left_join(.,sig_match1,by="first")%>%
      dplyr::left_join(.,sig_match2,by="second",suffix=c("1","2"))%>%
      dplyr::mutate(condlist = paste("sig_",Condition2,"-",Condition1,sep=""))
  })
  
  
  output$ploteuler <- renderPlot({

    
    
    if (input$customall == "All") {
      df<-tableresults()%>%
        dplyr::ungroup()%>%
        dplyr::select(contains("sig_"))
    }
    
    if (input$customall == "Custom") {
      df<-tableresults()%>%
        dplyr::ungroup()%>%
        dplyr::select(dplyr::one_of(inputvennsigcols()$condlist))
    }
    

    if (input$txtupdownoverlap == "both") {
      df2 <- (df == "up") | (df == "down")
    }

    if (input$txtupdownoverlap == "up") {
      df2 <- (df == "up")
    }

    if (input$txtupdownoverlap == "down") {
      df2 <- (df == "down")
    }

    a<-as.matrix(df2)

    a[is.na(a)] = F

    v = gsub("sig_","",colnames(a)) #Specific nomenclature
    fit = eulerr::euler(a)
    plot(fit,
                   quantities = list(cex = 0.5),
                   fills = list(fill = brewer.pal(ncol(a),"Set3"), alpha = 0.8, cex=0.7),
                   # lty = 1:3,
                   adjust_labels=T,
                   labels = list(labels=v, font = 12, cex=1),
                   main = paste("Venn Diagram:",input$txtupdownoverlap))


  })
  
  ploteulerprint <- reactive({

    
    if (input$customall == "All") {
      df<-tableresults()%>%
        dplyr::ungroup()%>%
        dplyr::select(contains("sig_"))
    }
    
    if (input$customall == "Custom") {
      df<-tableresults()%>%
        dplyr::ungroup()%>%
        dplyr::select(dplyr::one_of(inputvennsigcols()$condlist))
    }

    if (input$txtupdownoverlap == "both") {
      df2 <- (df == "up") | (df == "down")
    }

    if (input$txtupdownoverlap == "up") {
      df2 <- (df == "up")
    }

    if (input$txtupdownoverlap == "down") {
      df2 <- (df == "down")
    }

    a<-as.matrix(df2)
    a[is.na(a)] = F
    v = gsub("sig_","",colnames(a)) #Specific nomenclature
    fit = eulerr::euler(a)
    #brewer.pal(ncol(a),"Set3")
    return(plot(fit,
         quantities = list(cex = 0.5),
         fills = list(fill = brewer.pal(ncol(a),"Set3"), alpha = 0.8, cex=0.7),
         # lty = 1:3,
         adjust_labels=T,
         labels = list(labels=v, font = 12, cex=1),
         main = paste("Venn Diagram:",input$txtupdownoverlap)))


  })
  
  output$printploteuler <- downloadHandler(
    filename = function() { paste("VennDiagram_",input$txtupdownoverlap, '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = ploteulerprint(), width = 8, height = 8)
    }
  )
  
  
  output$venndiagramtable <- renderReactable({
    
    if (input$customall == "All") {
      df<-tableresults()%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation,contains("sig_"))
    }
    
    if (input$customall == "Custom") {
      df<-tableresults()%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation,dplyr::one_of(inputvennsigcols()$condlist))
    }
    
    
    colnames(df) <- gsub("sig_","",colnames(df))
    
    if (input$txtupdownoverlap == "both") {
      df3 <- df[apply(df == "up",1,any),]
      df4 <- df[apply(df == "down",1,any),]
      df2 <- dplyr::bind_rows(df3,df4)
      # df2 <- (df == "up") | (df == "down")
    }
    
    if (input$txtupdownoverlap == "up") {
      df2 <- df[apply(df == "up",1,any),]
      # df2 <- (df == "up")
    }
    
    if (input$txtupdownoverlap == "down") {
      df2 <- df[apply(df == "down",1,any),]
      # df2 <- (df == "down")
    }
    
    reactable(df2,filterable = TRUE)
  })
  
  
  venndiagramtableexport <- reactive({
    if (input$customall == "All") {
      df<-tableresults()%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation,contains("sig_"))
    }
    
    if (input$customall == "Custom") {
      df<-tableresults()%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation,dplyr::one_of(inputvennsigcols()$condlist))
    }
    
    
    colnames(df) <- gsub("sig_","",colnames(df))
    
    if (input$txtupdownoverlap == "both") {
      df3 <- df[apply(df == "up",1,any),]
      df4 <- df[apply(df == "down",1,any),]
      df2 <- dplyr::bind_rows(df3,df4)
    }
    
    if (input$txtupdownoverlap == "up") {
      df2 <- df[apply(df == "up",1,any),]
    }
    
    if (input$txtupdownoverlap == "down") {
      df2 <- df[apply(df == "down",1,any),]
    }
    return(df2)
  })
  
  
  
  
  output$venndiagramtableexport2 <- downloadHandler(
    filename = function() {
      paste("venn_diagram_overlap_data_",input$txtupdownoverlap ,".csv", sep = "")
    },
    content = function(file) {
      write.csv(venndiagramtableexport(), file, row.names = FALSE)
    }
  )
  
  
  
  output$plotupset <-renderPlot({
    if (input$customall == "All") {
      df<-tableresults()%>%
        dplyr::ungroup()%>%
        dplyr::select(contains("sig_"))
    }
    
    if (input$customall == "Custom") {
      df<-tableresults()%>%
        dplyr::ungroup()%>%
        dplyr::select(dplyr::one_of(inputvennsigcols()$condlist))
    }

    colnames(df) <- gsub("sig_","",colnames(df))

    if (input$txtupdownoverlap == "both") {
      df2 <- (df == "up") | (df == "down")
    }

    if (input$txtupdownoverlap == "up") {
      df2 <- (df == "up")
    }

    if (input$txtupdownoverlap == "down") {
      df2 <- (df == "down")
    }

    a<-as.matrix(df2)
    a = a[as.logical(rowSums(a)),] #Remove proteins without group

    plot3.3 = ComplexUpset::upset(
      as.data.frame(a),colnames(a), name = paste("Upset Plot:",input$txtupdownoverlap),
       min_size = 1
    )
    #group_by= "sets",

    plot3.3
  })
  
  plotupsetprint <- reactive({
    
    
    if (input$customall == "All") {
      df<-tableresults()%>%
        dplyr::ungroup()%>%
        dplyr::select(contains("sig_"))
    }
    
    if (input$customall == "Custom") {
      df<-tableresults()%>%
        dplyr::ungroup()%>%
        dplyr::select(dplyr::one_of(inputvennsigcols()$condlist))
    }

    colnames(df) <- gsub("sig_","",colnames(df))

    if (input$txtupdownoverlap == "both") {
      df2 <- (df == "up") | (df == "down")
    }

    if (input$txtupdownoverlap == "up") {
      df2 <- (df == "up")
    }

    if (input$txtupdownoverlap == "down") {
      df2 <- (df == "down")
    }

    a<-as.matrix(df2)
    a = a[as.logical(rowSums(a)),] #Remove proteins without group
    plot3.3 = upset(
      as.data.frame(a),colnames(a), name = paste("Upset Plot:",input$txtupdownoverlap),
      min_size = 1
    )

    return(plot3.3)
  })
  
  output$printplotupset <- downloadHandler(
    filename = function() { paste("UpsetPlot_",input$txtupdownoverlap, '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotupsetprint(), width = 12, height = 6.5)
    }
  )
  
  
  
  
  
  
  
  output$plotVL <- renderPlot({
    req(input$maxNlabelVolcano)
    req(input$volcanoLabelOptions)
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse7vol])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse8vol])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    tableresults2<-tableresults()%>%
      dplyr::mutate(Reference = as.character(Reference))%>%
      tidyr::separate(Reference, into = c("type","Reference2","addit"),sep="\\|")%>%
      dplyr::mutate(ProtID = paste(Gene.Symbol,Reference2,sep="_"))

    if (input$volcanoLabelOptions == "All") {
      p1<-ggplot(tableresults2, aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]))+
        geom_point(alpha=input$alphanum2,size=0.5)+
        geom_text_repel(data=tableresults2%>%
                          dplyr::filter(.data[[sig_name]] != "n.s.")%>%
                          dplyr::mutate(label=ProtID),aes(.data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1)+
        theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste(foldchangevalue), y= "-Log10(q-value)")+
        theme(axis.text = element_text(size=16),
              axis.title = element_text(size = 18),
              legend.text = element_text(size=14),
              legend.title = element_text(size=16))
    }
    if (input$volcanoLabelOptions == "CustomN Log2FC") {
      p1<-ggplot(tableresults2, aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]))+
        geom_point(alpha=input$alphanum2,size=0.5)+
        geom_text_repel(data=tableresults2%>%
                          group_by(.data[[sig_name]])%>%
                          dplyr::filter(.data[[sig_name]] != "n.s.")%>%
                          dplyr::mutate(rankVal = dplyr::dense_rank(dplyr::desc(abs(.data[[foldchangevalue]]))))%>%
                          dplyr::filter(rankVal <= input$maxNlabelVolcano)%>%
                          dplyr::mutate(label=ProtID ),aes(.data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1)+
        theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste(foldchangevalue), y= "-Log10(q-value)")+
        theme(axis.text = element_text(size=16),
              axis.title = element_text(size = 18),
              legend.text = element_text(size=14),
              legend.title = element_text(size=16))
    }
    if (input$volcanoLabelOptions == "CustomN qvalue") {
      p1<-ggplot(tableresults2, aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]))+
        geom_point(alpha=input$alphanum2,size=0.5)+
        geom_text_repel(data=tableresults2%>%
                          group_by(.data[[sig_name]])%>%
                          dplyr::filter(.data[[sig_name]] != "n.s.")%>%
                          dplyr::mutate(rankVal = dplyr::dense_rank(abs(.data[[qvaluename]])))%>%
                          dplyr::filter(rankVal <= input$maxNlabelVolcano)%>%
                          dplyr::mutate(label=ProtID ),aes(.data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1)+
        theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste(foldchangevalue), y= "-Log10(q-value)")+
        theme(axis.text = element_text(size=16),
              axis.title = element_text(size = 18),
              legend.text = element_text(size=14),
              legend.title = element_text(size=16))
    }
    p1
  },width = 900,height=650)
  
  plotVLprint <- reactive({
    req(input$maxNlabelVolcano)
    req(input$volcanoLabelOptions)
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse7vol])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse8vol])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    tableresults2<-tableresults()%>%
      dplyr::mutate(Reference = as.character(Reference))%>%
      tidyr::separate(Reference, into = c("type","Reference2","addit"),sep="\\|")%>%
      dplyr::mutate(ProtID = paste(Gene.Symbol,Reference2,sep="_"))
    
    if (input$volcanoLabelOptions == "All") {
      p1<-ggplot(tableresults2, aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]))+
        geom_point(alpha=input$alphanum2,size=0.5)+
        geom_text_repel(data=tableresults2%>%
                          dplyr::filter(.data[[sig_name]] != "n.s.")%>%
                          dplyr::mutate(label=ProtID),aes(.data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1)+
        theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste(foldchangevalue), y= "-Log10(q-value)")+
        theme(axis.text = element_text(size=16),
              axis.title = element_text(size = 18),
              legend.text = element_text(size=14),
              legend.title = element_text(size=16))
    }
    if (input$volcanoLabelOptions == "CustomN Log2FC") {
      p1<-ggplot(tableresults2, aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]))+
        geom_point(alpha=input$alphanum2,size=0.5)+
        geom_text_repel(data=tableresults2%>%
                          group_by(.data[[sig_name]])%>%
                          dplyr::filter(.data[[sig_name]] != "n.s.")%>%
                          dplyr::mutate(rankVal = dplyr::dense_rank(dplyr::desc(abs(.data[[foldchangevalue]]))))%>%
                          dplyr::filter(rankVal <= input$maxNlabelVolcano)%>%
                          dplyr::mutate(label=ProtID ),aes(.data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1)+
        theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste(foldchangevalue), y= "-Log10(q-value)")+
        theme(axis.text = element_text(size=16),
              axis.title = element_text(size = 18),
              legend.text = element_text(size=14),
              legend.title = element_text(size=16))
    }
    if (input$volcanoLabelOptions == "CustomN qvalue") {
      p1<-ggplot(tableresults2, aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]))+
        geom_point(alpha=input$alphanum2,size=0.5)+
        geom_text_repel(data=tableresults2%>%
                          group_by(.data[[sig_name]])%>%
                          dplyr::filter(.data[[sig_name]] != "n.s.")%>%
                          dplyr::mutate(rankVal = dplyr::dense_rank(abs(.data[[qvaluename]])))%>%
                          dplyr::filter(rankVal <= input$maxNlabelVolcano)%>%
                          dplyr::mutate(label=ProtID ),aes(.data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1)+
        theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste(foldchangevalue), y= "-Log10(q-value)")+
        theme(axis.text = element_text(size=16),
              axis.title = element_text(size = 18),
              legend.text = element_text(size=14),
              legend.title = element_text(size=16))
    }
    return(p1)
    # con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse7vol])
    # con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse8vol])
    # 
    # foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    # qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    # sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    # 
    # 
    # 
    # p5<-ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]))+
    #   geom_point(alpha=input$alphanum2,size=0.5)+
    #   geom_text_repel(data=.%>%
    #                     mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(label=label),
    #                   max.overlaps = Inf,show.legend=FALSE,box.padding = 1)+
    #   theme_classic()+
    #   scale_color_viridis_d(end=0.8)+
    #   geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
    #   geom_vline(xintercept=log2(input$num78),linetype="dashed")+
    #   geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
    #   labs(x= paste(foldchangevalue), y= "-Log10(q-value)")
    # 
    # return(p5)

  })
  
  output$volcanoDownloadPlot <- downloadHandler(
    filename = function() { paste("Volcano_Label_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse8vol]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse7vol]), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotVLprint(), width = 12, height = 6.5)
    }
  )
  
  
  uniqueColProtein<-reactive({
    df_unique_val <-untidy_stats_all()%>%
      dplyr::ungroup()%>%
      dplyr::select(Gene.Symbol)%>%
      plotly::distinct()
  })
  

  
  
  
  output$plot5 <- renderPlot({
    req(input$text11)

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1vol])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2vol])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")

    gene_name_input <- input$text11
    plot_protein<-untidy_stats_all()%>%
      dplyr::filter(Gene.Symbol==gene_name_input)%>%
      dplyr::mutate(ProtID2 = paste(Gene.Symbol,Reference,sep="_"))

    plot_protein_v2<-plot_protein%>%dplyr::group_by(Condition,ProtID2)%>%dplyr::summarise(Abundance=median(Abundance))

    ggplot(plot_protein,aes(x=Condition,2^Abundance,fill = Condition,color=Condition))+
      geom_col(data=plot_protein_v2,alpha=0.5)+
      geom_jitter(size=3)+
      theme_classic()+scale_fill_viridis_d(end=0.8)+scale_color_viridis_d(end=0.8)+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=24,hjust=0.5))+
      labs(title = gene_name_input, x="",y="MS Intensity")+facet_wrap(~ProtID2,nrow=1,scales="free_y")

  })
  
  individProteinPlotInput <- reactive({
    req(input$text11)

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1vol])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2vol])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")

    gene_name_input <- input$text11
    plot_protein<-untidy_stats_all()%>%
      dplyr::filter(Gene.Symbol==gene_name_input)%>%
      dplyr::mutate(ProtID2 = paste(Gene.Symbol,Reference,sep="_"))

    plot_protein_v2<-plot_protein%>%dplyr::group_by(Condition,ProtID2)%>%dplyr::summarise(Abundance=median(Abundance))

    ggplot(plot_protein,aes(x=Condition,2^Abundance,fill = Condition,color=Condition))+
      geom_col(data=plot_protein_v2,alpha=0.5)+
      geom_jitter(size=3)+
      theme_classic()+scale_fill_viridis_d(end=0.8)+scale_color_viridis_d(end=0.8)+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=24,hjust=0.5))+
      labs(title = gene_name_input, x="",y="MS Intensity")+facet_wrap(~ProtID2,nrow=1,scales="free_y")
  })
  
  
  output$proteinBarplotDownloadPlot <- downloadHandler(
    filename = function() { paste(input$text11,"_SingleProtein_Barplot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = individProteinPlotInput(), width = 10, height = 7)
    }
  )
  
  volcanoPlot2Input <- reactive({
    x_axis <- paste(unique(norm_sum_final2()$Condition)[2], "-", unique(norm_sum_final2()$Condition)[1],sep=" ")

    p4000<-ggplot(data=results_msstats(), aes(log2_foldchange, -log2(q.val),color=Significant))+
      geom_point(alpha=input$alphanum1)+geom_text_repel(data=.%>%
                                     mutate(label=if_else(Significant != "n.s.",Gene.Symbol,"")),aes(label=label),
                                   max.overlaps = Inf,show.legend=FALSE,box.padding = 1)+theme_classic()+
      geom_vline(xintercept = log2(input$num4))+geom_vline(xintercept = -log2(input$num4))+
      geom_hline(yintercept = -log2(input$num3))+scale_color_viridis_d(end=0.8)+
      labs(x=paste("log2(",x_axis,")"), y="-log2(q-value)")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18))
  })
  

  
  output$volcanoDownloadPlot2 <- downloadHandler(
    filename = function() { paste("Volcano_Label_Plot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = volcanoPlot2Input(), width = 12, height = 6.5)
    }
  )
  
  
  output$text12 <- renderText({
    gene_name_input <- input$text11
    plot_protein<-untidy_stats_all()%>%
      dplyr::filter(Gene.Symbol==gene_name_input)
    annotation_plot<-unique(plot_protein$Annotation)
    return(annotation_plot)
  })
  output$text13 <- renderText({
    gene_name_input <- input$text11
    plot_protein<-untidy_stats_all()%>%
      dplyr::filter(Gene.Symbol==gene_name_input)
    ref_plot<-unique(plot_protein$Reference)
    return(ref_plot)
  })
  
  outputdata<-reactive({
    reps_df<-prot_biorep_final()%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")

    all_table<-left_join(results_msstats(),reps_df, by= c("Reference","Gene.Symbol","Annotation"))
  })
    
    
  output$table1 <- renderReactable({
    cutoff<-input$num1
    fc_cutoff<-input$num2

    reps_df<-prot_biorep_final()%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")

    all_table<-left_join(results_msstats(),reps_df, by= c("Reference","Gene.Symbol","Annotation"))

    slim_df_for_output<-all_table%>%
      dplyr::filter(q.val < cutoff, abs(log2_foldchange)>log2(fc_cutoff))
    reactable(slim_df_for_output, filterable = TRUE)
  })
  
  output$table16 <- renderReactable({
    
    clusterdesignation<-clusterTraceDFinitial()%>%
      dplyr::ungroup()%>%
      dplyr::select(ProtID,cluster)%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    
    all_table <- dplyr::left_join(results_msstats(),
                                  clusterdesignation,by=c("Reference","Gene.Symbol","Annotation"))%>%
      dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
      # the good stuff here
      dplyr::mutate_if(is.numeric,dplyr::coalesce,0)
    
    # all_table<-left_join(HotelClusterDF,reps_df, by= c("Reference","Gene.Symbol","Annotation"))

    reactable(all_table, filterable = TRUE)
  })
  
  tableresults <- reactive({

    clusterdesignation<-clusterTraceDFinitial()%>%
      dplyr::ungroup()%>%
      dplyr::select(ProtID,cluster)%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")

    all_table <- dplyr::left_join(results_msstats(),
                                  clusterdesignation,by=c("Reference","Gene.Symbol","Annotation"))%>%
      dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
      # the good stuff here
      dplyr::mutate_if(is.numeric,dplyr::coalesce,0)
    # print(head(all_table))

    return(all_table)
  })
  
  
  ### fix to table16 from df
  # tableresults <- reactive({
  #   
  #   all_table <- df()
  #   
  #   return(all_table)
  # })
  
  
  
  output$alldownload <- downloadHandler(
    filename = function() {
      paste("all_stats_results_",input$uniqueID,".csv", sep = "")
    },
    content = function(file) {
      write.csv(tableresults(), file, row.names = FALSE)
    }
  )
  ### localization start
  annotation_ttest<-reactive({

    if (input$txthumanmouse == "mouse") {
      ann_combine<-tableresults()%>%
        # dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale(),by="Gene.Symbol", relationship = "many-to-many")
    }

    if (input$txthumanmouse == "human") {
      ann_combine<-tableresults()%>%
        # dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale(),by="Gene.Symbol", relationship = "many-to-many")
    }
    if (input$txthumanmouse == "yeast") {
      ann_combine<-tableresults()%>%
        # dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale(),by="Gene.Symbol", relationship = "many-to-many")
    }
    return(ann_combine)
  })
  
  annotation_ttest2<-reactive({
    if (input$txthumanmouse == "mouse") {
      ann_combine<-tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale2,by="Gene.Symbol", relationship = "many-to-many")
    }

    if (input$txthumanmouse == "human") {
      ann_combine<-tableresults()%>%
        # dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale2,by="Gene.Symbol", relationship = "many-to-many")
    }
    if (input$txthumanmouse == "yeast") {
      ann_combine<-tableresults()%>%
        # dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale2,by="Gene.Symbol", relationship = "many-to-many")
    }
    return(ann_combine)
  })
  
  ###localization
  output$plot7 <- renderPlot({


    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))

    if (input$txtinversionloc == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con2,"-",con1,sep=""), y= "-Log10(q-value)",title=input$localization)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }

    if (input$txtinversionloc == "Yes") {
      # foldchangevalue2<-as.character(paste("log2FC_",con1,"-",con2,sep=""))
      p1<-ggplot()+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con1,"-",con2,sep=""), y= "-Log10(q-value)",title=input$localization)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    p1
  })
  
  plot7print <- reactive({


    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))

    if (input$txtinversionloc == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con2,"-",con1,sep=""), y= "-Log10(q-value)",title=input$localization)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }

    if (input$txtinversionloc == "Yes") {
      # foldchangevalue2<-as.character(paste("log2FC_",con1,"-",con2,sep=""))
      p1<-ggplot()+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con1,"-",con2,sep=""), y= "-Log10(q-value)",title=input$localization)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    return(p1)
  })
  
  
  output$printplot7 <- downloadHandler(
    filename = function() { paste("Volcano_Localizationv1_",input$localization,"_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc]),"-", as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plot7print(), width = 12, height = 6.5)
    }
  )
  
  
  
  output$plot7loc2 <- renderPlot({


    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc5])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc5])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))

    if (input$txtinversionloc5 == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con2,"-",con1,sep=""), y= "-Log10(q-value)",title=input$localization5)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }

    if (input$txtinversionloc5 == "Yes") {
      # foldchangevalue2<-as.character(paste("log2FC_",con1,"-",con2,sep=""))
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con1,"-",con2,sep=""), y= "-Log10(q-value)",title=input$localization5)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    p1
  })
  
  plot7loc2print <- reactive({


    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc5])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc5])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))

    if (input$txtinversionloc5 == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con2,"-",con1,sep=""), y= "-Log10(q-value)",title=input$localization5)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }

    if (input$txtinversionloc5 == "Yes") {
      # foldchangevalue2<-as.character(paste("log2FC_",con1,"-",con2,sep=""))
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con1,"-",con2,sep=""), y= "-Log10(q-value)",title=input$localization5)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    return(p1)
  })
  
  output$printplot7loc2 <- downloadHandler(
    filename = function() { paste("Volcano_Localizationv2_",input$localization5,"_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc5]),"-", as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc5]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plot7loc2print(), width = 12, height = 6.5)
    }
  )
  

  
  output$corrcolorcompartment<-renderPlot({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc])

    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))

    ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con1]], .data[[con2]]),color="black")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con1]], .data[[con2]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=con1, y=con2,title=input$localization)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
  })

  
  
  output$downloadlocalizationdata <- downloadHandler(
    filename = function() {
      paste("All_results_stats_localizationAdded", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(annotation_ttest(), file, row.names = FALSE)
    }
  )
  
  corrcolorcompartmentprint<-reactive({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc])

    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))


    return(ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con1]], .data[[con2]]),color="black")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con1]], .data[[con2]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=con1, y=con2,title=input$localization)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5)))
  })
  
  
  output$printcorrcolorcompartment <- downloadHandler(
    filename = function() { paste("Correlation_Localizationv1_",input$localization,"_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc]),"vs", as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corrcolorcompartmentprint(), width = 8, height = 8)
    }
  )
  
  
  output$corrcolorcompartmentloc5<-renderPlot({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc5])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc5])

    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))


    ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con1]], .data[[con2]]),color="black")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con1]], .data[[con2]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=con1, y=con2,title=input$localization5)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
  })
  
  corrcolorcompartmentloc5print<-reactive({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc5])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc5])

    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))


    return(ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con1]], .data[[con2]]),color="black")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con1]], .data[[con2]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=con1, y=con2,title=input$localization5)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5)))
  })
  
  output$printcorrcolorcompartmentloc5 <- downloadHandler(
    filename = function() { paste("Correlation_Localizationv2_",input$localization5,"_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc5]),"vs", as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc5]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corrcolorcompartmentloc5print(), width = 8, height = 8)
    }
  )
  
  output$corrcolorcompartmentmulti<-renderPlot({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc2])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc2])
    con3<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse3loc2])
    con4<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse4loc2])



    ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="black")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=paste(con2, "-",con1,sep=""), y= paste(con4,"-",con3,sep=""),title=input$localization)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
  })
  
  corrcolorcompartmentmultiprint<-reactive({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc2])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc2])
    con3<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse3loc2])
    con4<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse4loc2])



    return(ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="black")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=paste(con2, "-",con1,sep=""), y= paste(con4,"-",con3,sep=""),title=input$localization)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5)))
  })
  
  output$printcorrcolorcompartmentmulti <- downloadHandler(
    filename = function() { paste("Correlation_Localizationv1_Ratio_",input$localization,"_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc2]),"_DIV_", as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc2]),"_vs_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse4loc2]),"_DIV_", as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse3loc2]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corrcolorcompartmentmultiprint(), width = 8, height = 8)
    }
  )
  
  output$corrcolorcompartmentmultiloc5<-renderPlot({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc6])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc6])
    con3<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse3loc6])
    con4<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse4loc6])



    ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="black")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=paste(con2, "-",con1,sep=""), y= paste(con4,"-",con3,sep=""),title=input$localization5)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
  })
  
  corrcolorcompartmentmultiloc5print<-reactive({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc6])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc6])
    con3<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse3loc6])
    con4<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse4loc6])



    return(ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="black")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=paste(con2, "-",con1,sep=""), y= paste(con4,"-",con3,sep=""),title=input$localization5)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5)))
  })
  
  output$printcorrcolorcompartmentmultiloc5 <- downloadHandler(
    filename = function() { paste("Correlation_Localizationv2_Ratio_",input$localization5,"_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc6]),"_DIV_", as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc6]),"_vs_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse4loc6]),"_DIV_", as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse3loc6]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corrcolorcompartmentmultiloc5print(), width = 8, height = 8)
    }
  )
  
  
  
  
  output$plot8 <- renderPlot({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))

    x_axis <- paste(con2, "-", con1,sep=" ")
    ggplot(data=annotation_ttest()%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)), aes(reorder(Compartment, .data[[foldchangevalue]], FUN = median),.data[[foldchangevalue]],color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      # geom_jitter(alpha=0.1)+
      theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="Localization", y=paste("log2(",x_axis,")"),title="Fold change by localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+
      coord_flip()+
      ylim(c(input$numlocale1,input$numlocale2))
  })
  
  
  output$pointPlot<-renderPlot({
    
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc2])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc2])
    con3<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse3loc2])
    con4<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse4loc2])
    
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    fc2<-as.character(paste("log2FC_",con4,"-",con3,sep=""))
    
    
    
    if (input$ratioratioindicator == "Single Ratio") {

      
      pointplotinput<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_x = 5,
                      delta_y = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::group_by(Compartment)%>%
        dplyr::summarise(med_x=median(delta_x),
                         med_y = median(delta_y),
                         q1_x= quantile(delta_x,0.25),
                         q3_x= quantile(delta_x,0.75),
                         q1_y= quantile(delta_y,0.25),
                         q3_y= quantile(delta_y,0.75),
                         num_proteins = n())
      
      plot_output<-ggplot(pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment))+
        geom_point(size=4,alpha=0.6)+
        geom_segment(aes(x=q1_x,xend=q3_x,y=med_y,yend=med_y,color=Compartment),alpha=0.5,size=0.5)+
        geom_segment(aes(x=med_x,xend=med_x,y=q1_y,yend=q3_y,color=Compartment),alpha=0.5,size=0.5)+
        theme_classic()+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_vline(xintercept = 0,linetype="dashed")+
        geom_abline(slope = 1,linetype="dashed")+
        geom_text_repel(data=pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment, label=Compartment),
                        max.overlaps = Inf,show.legend=FALSE,size=4)+
        theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size=18))+
        labs(y=paste(fc2),x=paste(fc1))
    }
    
    if (input$ratioratioindicator == "X axis RoR") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse5loc2])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse6loc2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))
      
      join1_df<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc3) )%>%
        dplyr::rename(delta_w = 5)%>%
        dplyr::ungroup()
      
      
      pointplotinput<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_num = 5,
                      delta_y = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "Compartment"))%>%
        dplyr::mutate(delta_x = delta_num - delta_w)%>%
        dplyr::ungroup()%>%
        dplyr::group_by(Compartment)%>%
        dplyr::summarise(med_x=median(delta_x),
                         med_y = median(delta_y),
                         q1_x= quantile(delta_x,0.25),
                         q3_x= quantile(delta_x,0.75),
                         q1_y= quantile(delta_y,0.25),
                         q3_y= quantile(delta_y,0.75),
                         num_proteins = n())
      
      plot_output<-ggplot(pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment))+
        geom_point(size=4,alpha=0.6)+
        geom_segment(aes(x=q1_x,xend=q3_x,y=med_y,yend=med_y,color=Compartment),alpha=0.5,size=0.5)+
        geom_segment(aes(x=med_x,xend=med_x,y=q1_y,yend=q3_y,color=Compartment),alpha=0.5,size=0.5)+
        theme_classic()+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_vline(xintercept = 0,linetype="dashed")+
        geom_abline(slope = 1,linetype="dashed")+
        geom_text_repel(data=pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment, label=Compartment),
                        max.overlaps = Inf,show.legend=FALSE,size=4)+
        theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size=18))+
        labs(y=paste(fc2),x=paste(fc1,"-",fc3,sep=""))
    }
    
    if (input$ratioratioindicator == "Y axis RoR") {
      con7<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse7loc2])
      con8<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse8loc2])
      fc4<-as.character(paste("log2FC_",con8,"-",con7,sep=""))
      
      join1_df<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc4) )%>%
        dplyr::rename(delta_w = 5)%>%
        dplyr::ungroup()
      
      
      pointplotinput<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_x = 5,
                      delta_num = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "Compartment"))%>%
        dplyr::mutate(delta_y = delta_num - delta_w)%>%
        dplyr::ungroup()%>%
        dplyr::group_by(Compartment)%>%
        dplyr::summarise(med_x=median(delta_x),
                         med_y = median(delta_y),
                         q1_x= quantile(delta_x,0.25),
                         q3_x= quantile(delta_x,0.75),
                         q1_y= quantile(delta_y,0.25),
                         q3_y= quantile(delta_y,0.75),
                         num_proteins = n())
      
      plot_output<-ggplot(pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment))+
        geom_point(size=4,alpha=0.6)+
        geom_segment(aes(x=q1_x,xend=q3_x,y=med_y,yend=med_y,color=Compartment),alpha=0.5,size=0.5)+
        geom_segment(aes(x=med_x,xend=med_x,y=q1_y,yend=q3_y,color=Compartment),alpha=0.5,size=0.5)+
        theme_classic()+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_vline(xintercept = 0,linetype="dashed")+
        geom_abline(slope = 1,linetype="dashed")+
        geom_text_repel(data=pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment, label=Compartment),
                        max.overlaps = Inf,show.legend=FALSE,size=4)+
        theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size=18))+
        labs(y=paste(fc2,"-",fc4,sep=""),x=paste(fc1))
    }
    if (input$ratioratioindicator == "RoR same Denominator") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse5loc2])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse6loc2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))
      
      join1_df<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc3) )%>%
        dplyr::rename(delta_w = 5)%>%
        dplyr::ungroup()
      
      
      pointplotinput<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_numx = 5,
                      delta_numy = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "Compartment"))%>%
        dplyr::mutate(delta_x = delta_numx - delta_w,
                      delta_y = delta_numy - delta_w)%>%
        dplyr::ungroup()%>%
        dplyr::group_by(Compartment)%>%
        dplyr::summarise(med_x=median(delta_x),
                         med_y = median(delta_y),
                         q1_x= quantile(delta_x,0.25),
                         q3_x= quantile(delta_x,0.75),
                         q1_y= quantile(delta_y,0.25),
                         q3_y= quantile(delta_y,0.75),
                         num_proteins = n())
      
      plot_output<-ggplot(pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment))+
        geom_point(size=4,alpha=0.6)+
        geom_segment(aes(x=q1_x,xend=q3_x,y=med_y,yend=med_y,color=Compartment),alpha=0.5,size=0.5)+
        geom_segment(aes(x=med_x,xend=med_x,y=q1_y,yend=q3_y,color=Compartment),alpha=0.5,size=0.5)+
        theme_classic()+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_vline(xintercept = 0,linetype="dashed")+
        geom_abline(slope = 1,linetype="dashed")+
        geom_text_repel(data=pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment, label=Compartment),
                        max.overlaps = Inf,show.legend=FALSE,size=4)+
        theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size=18))+
        labs(y=paste(fc2,"-",fc3,sep=""),x=paste(fc1,"-",fc3,sep=""))
    }
    if (input$ratioratioindicator == "RoR different Denominator") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse5loc2])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse6loc2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))
      
      con7<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse7loc2])
      con8<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse8loc2])
      fc4<-as.character(paste("log2FC_",con8,"-",con7,sep=""))
      
      join1_df<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc3),all_of(fc4) )%>%
        dplyr::rename(delta_w = 5,
                      delta_z = 6)%>%
        dplyr::ungroup()
      
      
      pointplotinput<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_numx = 5,
                      delta_numy = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "Compartment"))%>%
        dplyr::mutate(delta_x = delta_numx - delta_w,
                      delta_y = delta_numy - delta_z)%>%
        dplyr::ungroup()%>%
        dplyr::group_by(Compartment)%>%
        dplyr::summarise(med_x=median(delta_x),
                         med_y = median(delta_y),
                         q1_x= quantile(delta_x,0.25),
                         q3_x= quantile(delta_x,0.75),
                         q1_y= quantile(delta_y,0.25),
                         q3_y= quantile(delta_y,0.75),
                         num_proteins = n())
      
      plot_output<-ggplot(pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment))+
        geom_point(size=4,alpha=0.6)+
        geom_segment(aes(x=q1_x,xend=q3_x,y=med_y,yend=med_y,color=Compartment),alpha=0.5,size=0.5)+
        geom_segment(aes(x=med_x,xend=med_x,y=q1_y,yend=q3_y,color=Compartment),alpha=0.5,size=0.5)+
        theme_classic()+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_vline(xintercept = 0,linetype="dashed")+
        geom_abline(slope = 1,linetype="dashed")+
        geom_text_repel(data=pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment, label=Compartment),
                        max.overlaps = Inf,show.legend=FALSE,size=4)+
        theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size=18))+
        labs(y=paste(fc2,"-",fc4,sep=""),x=paste(fc1,"-",fc3,sep=""))
    }
    plot_output
    
    

  })
  
  pointPlotForPrint<-reactive({
    
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc2])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc2])
    con3<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse3loc2])
    con4<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse4loc2])
    
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    fc2<-as.character(paste("log2FC_",con4,"-",con3,sep=""))
    
    
    
    if (input$ratioratioindicator == "Single Ratio") {
      
      
      pointplotinput<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_x = 5,
                      delta_y = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::group_by(Compartment)%>%
        dplyr::summarise(med_x=median(delta_x),
                         med_y = median(delta_y),
                         q1_x= quantile(delta_x,0.25),
                         q3_x= quantile(delta_x,0.75),
                         q1_y= quantile(delta_y,0.25),
                         q3_y= quantile(delta_y,0.75),
                         num_proteins = n())
      
      plot_output<-ggplot(pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment))+
        geom_point(size=4,alpha=0.6)+
        geom_segment(aes(x=q1_x,xend=q3_x,y=med_y,yend=med_y,color=Compartment),alpha=0.5,size=0.5)+
        geom_segment(aes(x=med_x,xend=med_x,y=q1_y,yend=q3_y,color=Compartment),alpha=0.5,size=0.5)+
        theme_classic()+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_vline(xintercept = 0,linetype="dashed")+
        geom_abline(slope = 1,linetype="dashed")+
        geom_text_repel(data=pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment, label=Compartment),
                        max.overlaps = Inf,show.legend=FALSE,size=4)+
        theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size=18))+
        labs(y=paste(fc2),x=paste(fc1))
    }
    
    if (input$ratioratioindicator == "X axis RoR") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse5loc2])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse6loc2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))
      
      join1_df<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc3) )%>%
        dplyr::rename(delta_w = 5)%>%
        dplyr::ungroup()
      
      
      pointplotinput<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_num = 5,
                      delta_y = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "Compartment"))%>%
        dplyr::mutate(delta_x = delta_num - delta_w)%>%
        dplyr::ungroup()%>%
        dplyr::group_by(Compartment)%>%
        dplyr::summarise(med_x=median(delta_x),
                         med_y = median(delta_y),
                         q1_x= quantile(delta_x,0.25),
                         q3_x= quantile(delta_x,0.75),
                         q1_y= quantile(delta_y,0.25),
                         q3_y= quantile(delta_y,0.75),
                         num_proteins = n())
      
      plot_output<-ggplot(pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment))+
        geom_point(size=4,alpha=0.6)+
        geom_segment(aes(x=q1_x,xend=q3_x,y=med_y,yend=med_y,color=Compartment),alpha=0.5,size=0.5)+
        geom_segment(aes(x=med_x,xend=med_x,y=q1_y,yend=q3_y,color=Compartment),alpha=0.5,size=0.5)+
        theme_classic()+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_vline(xintercept = 0,linetype="dashed")+
        geom_abline(slope = 1,linetype="dashed")+
        geom_text_repel(data=pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment, label=Compartment),
                        max.overlaps = Inf,show.legend=FALSE,size=4)+
        theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size=18))+
        labs(y=paste(fc2),x=paste(fc1,"-",fc3,sep=""))
    }
    
    if (input$ratioratioindicator == "Y axis RoR") {
      con7<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse7loc2])
      con8<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse8loc2])
      fc4<-as.character(paste("log2FC_",con8,"-",con7,sep=""))
      
      join1_df<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc4) )%>%
        dplyr::rename(delta_w = 5)%>%
        dplyr::ungroup()
      
      
      pointplotinput<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_x = 5,
                      delta_num = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "Compartment"))%>%
        dplyr::mutate(delta_y = delta_num - delta_w)%>%
        dplyr::ungroup()%>%
        dplyr::group_by(Compartment)%>%
        dplyr::summarise(med_x=median(delta_x),
                         med_y = median(delta_y),
                         q1_x= quantile(delta_x,0.25),
                         q3_x= quantile(delta_x,0.75),
                         q1_y= quantile(delta_y,0.25),
                         q3_y= quantile(delta_y,0.75),
                         num_proteins = n())
      
      plot_output<-ggplot(pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment))+
        geom_point(size=4,alpha=0.6)+
        geom_segment(aes(x=q1_x,xend=q3_x,y=med_y,yend=med_y,color=Compartment),alpha=0.5,size=0.5)+
        geom_segment(aes(x=med_x,xend=med_x,y=q1_y,yend=q3_y,color=Compartment),alpha=0.5,size=0.5)+
        theme_classic()+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_vline(xintercept = 0,linetype="dashed")+
        geom_abline(slope = 1,linetype="dashed")+
        geom_text_repel(data=pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment, label=Compartment),
                        max.overlaps = Inf,show.legend=FALSE,size=4)+
        theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size=18))+
        labs(y=paste(fc2,"-",fc4,sep=""),x=paste(fc1))
    }
    if (input$ratioratioindicator == "RoR same Denominator") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse5loc2])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse6loc2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))
      
      join1_df<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc3) )%>%
        dplyr::rename(delta_w = 5)%>%
        dplyr::ungroup()
      
      
      pointplotinput<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_numx = 5,
                      delta_numy = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "Compartment"))%>%
        dplyr::mutate(delta_x = delta_numx - delta_w,
                      delta_y = delta_numy - delta_w)%>%
        dplyr::ungroup()%>%
        dplyr::group_by(Compartment)%>%
        dplyr::summarise(med_x=median(delta_x),
                         med_y = median(delta_y),
                         q1_x= quantile(delta_x,0.25),
                         q3_x= quantile(delta_x,0.75),
                         q1_y= quantile(delta_y,0.25),
                         q3_y= quantile(delta_y,0.75),
                         num_proteins = n())
      
      plot_output<-ggplot(pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment))+
        geom_point(size=4,alpha=0.6)+
        geom_segment(aes(x=q1_x,xend=q3_x,y=med_y,yend=med_y,color=Compartment),alpha=0.5,size=0.5)+
        geom_segment(aes(x=med_x,xend=med_x,y=q1_y,yend=q3_y,color=Compartment),alpha=0.5,size=0.5)+
        theme_classic()+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_vline(xintercept = 0,linetype="dashed")+
        geom_abline(slope = 1,linetype="dashed")+
        geom_text_repel(data=pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment, label=Compartment),
                        max.overlaps = Inf,show.legend=FALSE,size=4)+
        theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size=18))+
        labs(y=paste(fc2,"-",fc3,sep=""),x=paste(fc1,"-",fc3,sep=""))
    }
    if (input$ratioratioindicator == "RoR different Denominator") {
      con5<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse5loc2])
      con6<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse6loc2])
      fc3<-as.character(paste("log2FC_",con6,"-",con5,sep=""))
      
      con7<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse7loc2])
      con8<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse8loc2])
      fc4<-as.character(paste("log2FC_",con8,"-",con7,sep=""))
      
      join1_df<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc3),all_of(fc4) )%>%
        dplyr::rename(delta_w = 5,
                      delta_z = 6)%>%
        dplyr::ungroup()
      
      
      pointplotinput<-annotation_ttest()%>%
        dplyr::select(Reference,Gene.Symbol,Annotation, Compartment,all_of(fc1), all_of(fc2) )%>%
        dplyr::rename(delta_numx = 5,
                      delta_numy = 6 )%>%
        dplyr::ungroup()%>%
        dplyr::left_join(.,join1_df , by=c("Reference","Gene.Symbol","Annotation", "Compartment"))%>%
        dplyr::mutate(delta_x = delta_numx - delta_w,
                      delta_y = delta_numy - delta_z)%>%
        dplyr::ungroup()%>%
        dplyr::group_by(Compartment)%>%
        dplyr::summarise(med_x=median(delta_x),
                         med_y = median(delta_y),
                         q1_x= quantile(delta_x,0.25),
                         q3_x= quantile(delta_x,0.75),
                         q1_y= quantile(delta_y,0.25),
                         q3_y= quantile(delta_y,0.75),
                         num_proteins = n())
      
      plot_output<-ggplot(pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment))+
        geom_point(size=4,alpha=0.6)+
        geom_segment(aes(x=q1_x,xend=q3_x,y=med_y,yend=med_y,color=Compartment),alpha=0.5,size=0.5)+
        geom_segment(aes(x=med_x,xend=med_x,y=q1_y,yend=q3_y,color=Compartment),alpha=0.5,size=0.5)+
        theme_classic()+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_vline(xintercept = 0,linetype="dashed")+
        geom_abline(slope = 1,linetype="dashed")+
        geom_text_repel(data=pointplotinput%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)),aes(med_x,med_y,color=Compartment, label=Compartment),
                        max.overlaps = Inf,show.legend=FALSE,size=4)+
        theme(legend.position = "none",axis.title = element_text(size=18),axis.text = element_text(size=18))+
        labs(y=paste(fc2,"-",fc4,sep=""),x=paste(fc1,"-",fc3,sep=""))
    }
    return(plot_output)
    
    
    
  })
  
  output$printpointPlot <- downloadHandler(
    filename = function() { paste("Correlation_Plot_Localization",'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = pointPlotForPrint(), width = 8, height = 8)
    }
  )
  
  plot8print <- reactive({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))

    x_axis <- paste(con2, "-", con1,sep=" ")
    return(ggplot(data=annotation_ttest()%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)), aes(reorder(Compartment, .data[[foldchangevalue]], FUN = median),.data[[foldchangevalue]],color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      # geom_jitter(alpha=0.1)+
      theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="Localization", y=paste("log2(",x_axis,")"),title="Fold change by localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+
      coord_flip()+
      ylim(c(input$numlocale1,input$numlocale2)))
  })
  
  
  output$printplot8 <- downloadHandler(
    filename = function() { paste("Localev1_FC_MultipleLocale_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plot8print(), width = 16, height = 12)
    }
  )
  
  output$plot8loc2 <- renderPlot({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc5])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc5])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))

    x_axis <- paste(con2, "-", con1,sep=" ")
    ggplot(data=annotation_ttest2()%>%dplyr::filter(Compartment %in% input$locales5,!is.na(Compartment)), aes(reorder(Compartment, .data[[foldchangevalue]], FUN = median),.data[[foldchangevalue]],color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      # geom_jitter(alpha=0.1)+
      theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="Localization", y=paste("log2(",x_axis,")"),title="Fold change by localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+
      coord_flip()+
      ylim(c(input$numlocale1loc2,input$numlocale2loc2))
  })
  
  plot8loc2print <- reactive({
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc5])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc5])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))

    x_axis <- paste(con2, "-", con1,sep=" ")
    return(ggplot(data=annotation_ttest2()%>%dplyr::filter(Compartment %in% input$locales5,!is.na(Compartment)), aes(reorder(Compartment, .data[[foldchangevalue]], FUN = median),.data[[foldchangevalue]],color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      # geom_jitter(alpha=0.1)+
      theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="Localization", y=paste("log2(",x_axis,")"),title="Fold change by localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+
      coord_flip()+
      ylim(c(input$numlocale1loc2,input$numlocale2loc2)))
  })
  
  output$printplot8loc2 <- downloadHandler(
    filename = function() { paste("Localev2_FC_MultipleLocale_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc5]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc5]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plot8loc2print(), width = 16, height = 12)
    }
  )
  
  
  
  output$plotmultilocale <- renderPlot({
    numcon<-length(unique(norm_sum_final2()$Condition))

    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])

      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))

      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]])
    }

    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))

      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    }

    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))

      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    }

    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))


      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    }


    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))

      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    }

    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      con8<-as.character(unique(norm_sum_final2()$Condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))

      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    }

    size_df <- dim(df)[2]

    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,Compartment)%>%
      tidyr::gather("Condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()


    dftidy<-dftidy%>%
      dplyr::mutate(Condition = str_remove(Condition,pattern=paste("-",con1,sep="")))


    if (numcon  == 3) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3))
    }

    if (numcon  == 4) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4))
    }

    if (numcon  == 5) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5))
    }

    if (numcon  == 6) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6))
    }

    if (numcon  == 7) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7))
    }

    if (numcon  == 8) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }


    ggplot(dftidy%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)), aes(Condition,log2_FC,color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),title="Foldchange by Localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))

  })
  
  plotmultilocaleprint <- reactive({
    numcon<-length(unique(norm_sum_final2()$Condition))

    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])

      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))

      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]])
    }

    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))

      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    }

    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))

      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    }

    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))


      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    }


    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))

      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    }

    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      con8<-as.character(unique(norm_sum_final2()$Condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))

      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    }

    size_df <- dim(df)[2]

    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,Compartment)%>%
      tidyr::gather("Condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()


    dftidy<-dftidy%>%
      dplyr::mutate(Condition = str_remove(Condition,pattern=paste("-",con1,sep="")))


    if (numcon  == 3) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3))
    }

    if (numcon  == 4) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4))
    }

    if (numcon  == 5) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5))
    }

    if (numcon  == 6) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6))
    }

    if (numcon  == 7) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7))
    }

    if (numcon  == 8) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }


    return(ggplot(dftidy%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)), aes(Condition,log2_FC,color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),title="Foldchange by Localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5)))

  })
  
  output$printplotmultilocale <- downloadHandler(
    filename = function() { paste("Localev1_FC_MultipleLocale_AlltoPriority1",'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotmultilocaleprint(), width = 16, height = 8)
    }
  )
  
  output$plotmultilocaleloc2 <- renderPlot({
    numcon<-length(unique(norm_sum_final2()$Condition))

    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])

      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))

      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]])
    }

    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))

      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    }

    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))

      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    }

    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))


      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    }

    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))

      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    }

    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      con8<-as.character(unique(norm_sum_final2()$Condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))

      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    }



    size_df <- dim(df)[2]

    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,Compartment)%>%
      tidyr::gather("Condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()


    dftidy<-dftidy%>%
      dplyr::mutate(Condition = str_remove(Condition,pattern=paste("-",con1,sep="")))


    if (numcon  == 3) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3))
    }

    if (numcon  == 4) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4))
    }

    if (numcon  == 5) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5))
    }

    if (numcon  == 6) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6))
    }

    if (numcon  == 7) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7))
    }

    if (numcon  == 8) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }


    ggplot(dftidy%>%dplyr::filter(Compartment %in% input$locales5,!is.na(Compartment)), aes(Condition,log2_FC,color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),title="Foldchange by Localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))

  })
  
  
  plotmultilocaleloc2print <- reactive({
    numcon<-length(unique(norm_sum_final2()$Condition))

    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])

      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))

      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]])
    }

    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))

      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    }

    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))

      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    }

    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))


      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    }

    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))

      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    }

    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final2()$Condition)[1])
      con2<-as.character(unique(norm_sum_final2()$Condition)[2])
      con3<-as.character(unique(norm_sum_final2()$Condition)[3])
      con4<-as.character(unique(norm_sum_final2()$Condition)[4])
      con5<-as.character(unique(norm_sum_final2()$Condition)[5])
      con6<-as.character(unique(norm_sum_final2()$Condition)[6])
      con7<-as.character(unique(norm_sum_final2()$Condition)[7])
      con8<-as.character(unique(norm_sum_final2()$Condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))

      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    }



    size_df <- dim(df)[2]

    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,Compartment)%>%
      tidyr::gather("Condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()


    dftidy<-dftidy%>%
      dplyr::mutate(Condition = str_remove(Condition,pattern=paste("-",con1,sep="")))


    if (numcon  == 3) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3))
    }

    if (numcon  == 4) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4))
    }

    if (numcon  == 5) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5))
    }

    if (numcon  == 6) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6))
    }

    if (numcon  == 7) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7))
    }

    if (numcon  == 8) {
      dftidy$Condition <- factor(dftidy$Condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }


    return(ggplot(dftidy%>%dplyr::filter(Compartment %in% input$locales5,!is.na(Compartment)), aes(Condition,log2_FC,color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),title="Foldchange by Localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5)))

  })
  
  
  output$printplotmultilocaleloc2 <- downloadHandler(
    filename = function() { paste("Localev2_FC_MultipleLocale_AlltoPriority1",'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotmultilocaleloc2print(), width = 16, height = 8)
    }
  )

  # 
  
  
  
  output$plot44 <- renderPlot({
    req(input$maxNlabelVolcano2)
    req(input$volcanoLabelOptions2)
    
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc3])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc3])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    
    tableresults2<-annotation_ttest()%>%
      dplyr::mutate(Reference = as.character(Reference))%>%
      tidyr::separate(Reference, into = c("type","Reference2","addit"),sep="\\|")%>%
      dplyr::mutate(ProtID = paste(Gene.Symbol,Reference2,sep="_"))
    
    if (input$volcanoLabelOptions2 == "All") {
      p1<-ggplot()+
        geom_point(data=tableresults2%>%dplyr::filter(Compartment != input$localization2), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]),color="black")+
        geom_point(data=tableresults2%>%dplyr::filter(Compartment == input$localization2), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]),color="red")+
        geom_text_repel(data=tableresults2%>%
                          dplyr::filter(Compartment == input$localization2)%>%
                          # dplyr::filter(.data[[sig_name]] != "n.s.")%>%
                          dplyr::mutate(label=ProtID),aes(.data[[foldchangevalue]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
      theme_classic()+
      geom_vline(xintercept = log2(input$num78))+
      geom_vline(xintercept = -log2(input$num78))+
      geom_hline(yintercept = -log10(input$num77))+
      labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title=input$localization2)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    if (input$volcanoLabelOptions2 == "CustomN Log2FC") {
      p1<-ggplot()+
        geom_point(data=tableresults2%>%dplyr::filter(Compartment != input$localization2), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]),color="black")+
        geom_point(data=tableresults2%>%dplyr::filter(Compartment == input$localization2), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]),color="red")+
        geom_text_repel(data=tableresults2%>%
                          dplyr::filter(Compartment == input$localization2)%>%
                          dplyr::mutate(upDown = dplyr::if_else(.data[[foldchangevalue]] >0,"up","down"))%>%
                          dplyr::group_by(upDown)%>%
                          dplyr::mutate(rankVal = dplyr::dense_rank(dplyr::desc(abs(.data[[foldchangevalue]]))))%>%
                          dplyr::filter(rankVal <= input$maxNlabelVolcano2)%>%
                          dplyr::mutate(label=ProtID ),aes(.data[[foldchangevalue]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+
        geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title=input$localization2)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    if (input$volcanoLabelOptions2 == "CustomN qvalue") {
      p1<-ggplot()+
        geom_point(data=tableresults2%>%dplyr::filter(Compartment != input$localization2), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]),color="black")+
        geom_point(data=tableresults2%>%dplyr::filter(Compartment == input$localization2), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]),color="red")+
        geom_text_repel(data=tableresults2%>%
                          dplyr::filter(Compartment == input$localization2)%>%
                          dplyr::mutate(upDown = dplyr::if_else(.data[[foldchangevalue]] >0,"up","down"))%>%
                          dplyr::group_by(upDown)%>%
                          dplyr::mutate(rankVal = dplyr::dense_rank(dplyr::desc(abs(.data[[foldchangevalue]]))))%>%
                          dplyr::filter(rankVal <= input$maxNlabelVolcano2)%>%
                          dplyr::mutate(label=ProtID ),aes(.data[[foldchangevalue]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+
        geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title=input$localization2)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    p1

  },width = 900,height=650)
  
  
  plot44print <- reactive({

    req(input$maxNlabelVolcano2)
    req(input$volcanoLabelOptions2)
    
    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc3])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc3])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    
    tableresults2<-annotation_ttest()%>%
      dplyr::mutate(Reference = as.character(Reference))%>%
      tidyr::separate(Reference, into = c("type","Reference2","addit"),sep="\\|")%>%
      dplyr::mutate(ProtID = paste(Gene.Symbol,Reference2,sep="_"))
    
    if (input$volcanoLabelOptions2 == "All") {
      p1<-ggplot()+
        geom_point(data=tableresults2%>%dplyr::filter(Compartment != input$localization2), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]),color="black")+
        geom_point(data=tableresults2%>%dplyr::filter(Compartment == input$localization2), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]),color="red")+
        geom_text_repel(data=tableresults2%>%
                          dplyr::filter(Compartment == input$localization2)%>%
                          # dplyr::filter(.data[[sig_name]] != "n.s.")%>%
                          dplyr::mutate(label=ProtID),aes(.data[[foldchangevalue]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+
        geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title=input$localization2)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    if (input$volcanoLabelOptions2 == "CustomN Log2FC") {
      p1<-ggplot()+
        geom_point(data=tableresults2%>%dplyr::filter(Compartment != input$localization2), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]),color="black")+
        geom_point(data=tableresults2%>%dplyr::filter(Compartment == input$localization2), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]),color="red")+
        geom_text_repel(data=tableresults2%>%
                          dplyr::filter(Compartment == input$localization2)%>%
                          dplyr::mutate(upDown = dplyr::if_else(.data[[foldchangevalue]] >0,"up","down"))%>%
                          dplyr::group_by(upDown)%>%
                          dplyr::mutate(rankVal = dplyr::dense_rank(dplyr::desc(abs(.data[[foldchangevalue]]))))%>%
                          dplyr::filter(rankVal <= input$maxNlabelVolcano2)%>%
                          dplyr::mutate(label=ProtID ),aes(.data[[foldchangevalue]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+
        geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title=input$localization2)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    if (input$volcanoLabelOptions2 == "CustomN qvalue") {
      p1<-ggplot()+
        geom_point(data=tableresults2%>%dplyr::filter(Compartment != input$localization2), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]),color="black")+
        geom_point(data=tableresults2%>%dplyr::filter(Compartment == input$localization2), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]),color="red")+
        geom_text_repel(data=tableresults2%>%
                          dplyr::filter(Compartment == input$localization2)%>%
                          dplyr::mutate(upDown = dplyr::if_else(.data[[foldchangevalue]] >0,"up","down"))%>%
                          dplyr::group_by(upDown)%>%
                          dplyr::mutate(rankVal = dplyr::dense_rank(dplyr::desc(abs(.data[[foldchangevalue]]))))%>%
                          dplyr::filter(rankVal <= input$maxNlabelVolcano2)%>%
                          dplyr::mutate(label=ProtID ),aes(.data[[foldchangevalue]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+
        geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title=input$localization2)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    p1

    return(p1)

  })
  
  
  output$printplot44 <- downloadHandler(
    filename = function() { paste("Localizationv1_Volcano_Labels_",input$localization2,"_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc3]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc3]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plot44print(), width = 12, height = 6.5)
    }
  )
  
  output$plot88 <- renderPlot({

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc8])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc8])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))


    if (input$txtinversionloc88 == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization8), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="red")+
        geom_text_repel(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8)%>%
                          mutate(label=if_else(Compartment == input$localization8,Gene.Symbol,"")),aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title=input$localization8)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }

    if (input$txtinversionloc88 == "Yes") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization8), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=Significant),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=Significant),color="red")+
        geom_text_repel(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8)%>%
                          mutate(label=if_else(Compartment == input$localization8,Gene.Symbol,"")),aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title=input$localization8)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }

    p1

  },width = 1200,height=650)
  
  
  plot88print <- reactive({

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc8])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc8])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))


    if (input$txtinversionloc88 == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization8), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="red")+
        geom_text_repel(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8)%>%
                          mutate(label=if_else(Compartment == input$localization8,Gene.Symbol,"")),aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title=input$localization8)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }

    if (input$txtinversionloc88 == "Yes") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization8), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=Significant),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=Significant),color="red")+
        geom_text_repel(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8)%>%
                          mutate(label=if_else(Compartment == input$localization8,Gene.Symbol,"")),aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title=input$localization8)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }

    return(p1)

  })
  
  output$printplot88 <- downloadHandler(
    filename = function() { paste("Localizationv2_Volcano_Labels_",input$localization8,"_",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse2loc8]),"-",as.character(unique(norm_sum_final2()$Condition)[input$numtimecourse1loc8]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plot88print(), width = 12, height = 6.5)
    }
  )

  
  dfuser <- reactive({
    req(input$file4)
    inFile <- input$file4
    if (is.null(inFile))
      return(NULL)
    tbl <- read.csv(inFile$datapath)

    return(tbl)
  })
  
  userinput<-reactive({
    req(input$file4)
    ann_combine<-tableresults()%>%
      dplyr::left_join(.,dfuser(),by="Gene.Symbol")%>%
      dplyr::mutate(UserAnnotation = if_else(!is.na(UserAnnotation),UserAnnotation,"__"))
    return(ann_combine)
  })
  
  heatmapUserdata <- reactive({
    req(input$file4)
    heatmapinputdataUser<-dplyr::inner_join(userheatmapreference(),dfuser(),by="Gene.Symbol", relationship="many-to-many")%>%
      # dplyr::select(-uservalue)%>%
      dplyr::arrange(UserAnnotation)
    return(heatmapinputdataUser)
  })
  
  output$heatmapUser <- renderPlot({
    heatmapuser_step1 <- heatmapUserdata()%>%
      tidyr::separate(Reference, into = c("type","Reference","gene_info"),sep="\\|")%>%
      dplyr::select(-type,-gene_info,-Annotation)%>%
      tidyr::unite("ProtID",c(Gene.Symbol,Reference),sep="_")
    
    heatmapuser_step2<-heatmapuser_step1%>%
      dplyr::select(-UserAnnotation)%>%
      tibble::column_to_rownames(var="ProtID")
      
    AnnoB<-data.frame(row.names=heatmapuser_step1$ProtID, Annotation=heatmapuser_step1$UserAnnotation)

    rg <- max(abs(heatmapuser_step2))
    
    names_vals<-heatmapuser_step2%>%
      tibble::rownames_to_column("ProtID")
    names_cols<-data.frame("sample" = colnames(names_vals)[2:dim(names_vals)[[2]]])%>%
      dplyr::mutate(Condition = str_sub(sample,start = 1,end = -3))%>%
      dplyr::left_join(.,norm_sum_final2()%>%dplyr::distinct(priority,Condition),by="Condition")%>%
      dplyr::arrange(priority,Condition,sample)%>%
      dplyr::mutate(sample = factor(sample, unique(sample)),
                    values = factor(Condition,unique(Condition)))
    
    annoD<-data.frame(row.names=names_cols$sample, Condition=names_cols$values)
    
      

    #, cluster_cols = FALSE,color=cols,
    if(input$txtclustercolsUser == "Both" && input$heatmapLabelRow == "Yes" ){
      p1<-pheatmap::pheatmap(heatmapuser_step2,   main = "Heatmap",annotation_row = AnnoB,annotation = annoD, breaks = seq(-rg, rg, length.out = 100))
    }

    if(input$txtclustercolsUser == "Rows Only" && input$heatmapLabelRow == "Yes"){
      p1<-pheatmap::pheatmap(heatmapuser_step2, main = "Heatmap",annotation_row = AnnoB,annotation = annoD, cluster_cols = FALSE, breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercolsUser == "Columns Only" && input$heatmapLabelRow == "Yes" ){
      p1<-pheatmap::pheatmap(heatmapuser_step2,   main = "Heatmap",annotation_row = AnnoB,annotation = annoD, cluster_rows = FALSE,breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercolsUser == "Neither" && input$heatmapLabelRow == "Yes"){
      p1<-pheatmap::pheatmap(heatmapuser_step2, main = "Heatmap",annotation_row = AnnoB,annotation = annoD, cluster_cols = FALSE,cluster_rows = FALSE, breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercolsUser == "Both" && input$heatmapLabelRow == "No" ){
      p1<-pheatmap::pheatmap(heatmapuser_step2,   main = "Heatmap",annotation_row = AnnoB,annotation = annoD, show_rownames=F,breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercolsUser == "Rows Only" && input$heatmapLabelRow == "No"){
      p1<-pheatmap::pheatmap(heatmapuser_step2,  main = "Heatmap",annotation_row = AnnoB,annotation = annoD, cluster_cols = FALSE, show_rownames=F,breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercolsUser == "Columns Only" && input$heatmapLabelRow == "No" ){
      p1<-pheatmap::pheatmap(heatmapuser_step2,   main = "Heatmap",annotation_row = AnnoB,annotation = annoD, cluster_rows = FALSE,show_rownames=F,breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercolsUser == "Neither" && input$heatmapLabelRow == "No"){
      p1<-pheatmap::pheatmap(heatmapuser_step2,  main = "Heatmap",annotation_row = AnnoB,annotation = annoD, cluster_rows = FALSE,cluster_cols = FALSE, show_rownames=F,breaks = seq(-rg, rg, length.out = 100))
    }
    
    
    p1


  })
  
  
  heatmapUserPrintInput <- reactive({
    heatmapuser_step1 <- heatmapUserdata()%>%
      tidyr::separate(Reference, into = c("type","Reference","gene_info"),sep="\\|")%>%
      dplyr::select(-type,-gene_info,-Annotation)%>%
      tidyr::unite("ProtID",c(Gene.Symbol,Reference),sep="_")
    
    heatmapuser_step2<-heatmapuser_step1%>%
      dplyr::select(-UserAnnotation)%>%
      tibble::column_to_rownames(var="ProtID")
    
    AnnoB<-data.frame(row.names=heatmapuser_step1$ProtID, Annotation=heatmapuser_step1$UserAnnotation)
    
    rg <- max(abs(heatmapuser_step2))
    
    names_vals<-heatmapuser_step2%>%
      tibble::rownames_to_column("ProtID")
    # names_cols<-data.frame("sample" = colnames(names_vals)[2:dim(names_vals)[[2]]])%>%
    #   dplyr::mutate(values = str_sub(sample,start = 1,end = -3))
    names_cols<-data.frame("sample" = colnames(names_vals)[2:dim(names_vals)[[2]]])%>%
      dplyr::mutate(Condition = str_sub(sample,start = 1,end = -3))%>%
      dplyr::left_join(.,norm_sum_final2()%>%dplyr::distinct(priority,Condition),by="Condition")%>%
      dplyr::arrange(priority,Condition,sample)%>%
      dplyr::mutate(sample = factor(sample, unique(sample)),
                    values = factor(Condition,unique(Condition)))
    
    annoD<-data.frame(row.names=names_cols$sample, Condition=names_cols$values)
    
    #, cluster_cols = FALSE,color=cols,
    if(input$txtclustercolsUser == "Both" && input$heatmapLabelRow == "Yes" ){
      p1<-pheatmap::pheatmap(heatmapuser_step2,   main = "Heatmap",annotation_row = AnnoB,annotation = annoD, breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercolsUser == "Rows Only" && input$heatmapLabelRow == "Yes"){
      p1<-pheatmap::pheatmap(heatmapuser_step2, main = "Heatmap",annotation_row = AnnoB,annotation = annoD, cluster_cols = FALSE, breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercolsUser == "Columns Only" && input$heatmapLabelRow == "Yes" ){
      p1<-pheatmap::pheatmap(heatmapuser_step2,   main = "Heatmap",annotation_row = AnnoB,annotation = annoD, cluster_rows = FALSE,breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercolsUser == "Neither" && input$heatmapLabelRow == "Yes"){
      p1<-pheatmap::pheatmap(heatmapuser_step2, main = "Heatmap",annotation_row = AnnoB,annotation = annoD, cluster_cols = FALSE,cluster_rows = FALSE, breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercolsUser == "Both" && input$heatmapLabelRow == "No" ){
      p1<-pheatmap::pheatmap(heatmapuser_step2,   main = "Heatmap",annotation_row = AnnoB,annotation = annoD, show_rownames=F,breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercolsUser == "Rows Only" && input$heatmapLabelRow == "No"){
      p1<-pheatmap::pheatmap(heatmapuser_step2,  main = "Heatmap",annotation_row = AnnoB,annotation = annoD, cluster_cols = FALSE, show_rownames=F,breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercolsUser == "Columns Only" && input$heatmapLabelRow == "No" ){
      p1<-pheatmap::pheatmap(heatmapuser_step2,   main = "Heatmap",annotation_row = AnnoB,annotation = annoD, cluster_rows = FALSE,show_rownames=F,breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercolsUser == "Neither" && input$heatmapLabelRow == "No"){
      p1<-pheatmap::pheatmap(heatmapuser_step2,  main = "Heatmap",annotation_row = AnnoB,annotation = annoD, cluster_rows = FALSE,cluster_cols = FALSE, show_rownames=F,breaks = seq(-rg, rg, length.out = 100))
    }
    
    
    return(p1)
    
    
  })
  
  
  output$printheatmapUser <- downloadHandler(
    filename = function() { paste("UserDefined_Heatmap", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = heatmapUserPrintInput(), width = 10, height = 10)
    }
  )
  
  output$plotuser <- renderPlot({
    req(input$file4)
    req(input$size1)
    req(input$size2)
    req(input$alpha5)
    req(input$alpha6)
    # x_axis <- paste(unique(norm_sum_final2()$Condition)[2], "-", unique(norm_sum_final2()$Condition)[1],sep=" ")

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconduser1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconduser2])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))

    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "none")%>%
                          mutate(label=if_else(UserAnnotation != "__",Gene.Symbol,"")),aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }

    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "none")%>%
                          mutate(label=if_else(UserAnnotation != "__",Gene.Symbol,"")),aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }

    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "No") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }

    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "No") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }
    p1


  })
  
  plotuserprint <- reactive({
    req(input$file4)
    req(input$size1)
    req(input$size2)
    req(input$alpha5)
    req(input$alpha6)
    # x_axis <- paste(unique(norm_sum_final2()$Condition)[2], "-", unique(norm_sum_final2()$Condition)[1],sep=" ")

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconduser1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconduser2])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))

    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "none")%>%
                          mutate(label=if_else(UserAnnotation != "__",Gene.Symbol,"")),aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }

    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "none")%>%
                          mutate(label=if_else(UserAnnotation != "__",Gene.Symbol,"")),aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }

    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "No") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }

    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "No") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }
    return(p1)


  })
  
  output$printplotuser <- downloadHandler(
    filename = function() { paste("UserDefined_VolcanoPlot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotuserprint(), width = 12, height = 8)
    }
  )
  
  
  
  
  output$categviolin <- renderPlot({

    req(input$file4)
    req(input$numconduser1)
    req(input$numconduser2)

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconduser1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconduser2])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))


    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con2]] - .data[[con1]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
                          mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con2]] - .data[[con1]],label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con2,"-",con1,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }

    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con1]] - .data[[con2]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
                          mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con1]] - .data[[con2]],label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con1,"-",con2,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }

    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "No") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con2]] - .data[[con1]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        # geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
        #                   mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con2]] - .data[[con1]],label=label,color=UserAnnotation),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con2,"-",con1,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }

    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "No") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con1]] - .data[[con2]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        # geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
        #                   mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con1]] - .data[[con2]],label=label,color=UserAnnotation),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con1,"-",con2,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }

    p1


  })
  
  categviolinprint <- reactive({

    req(input$file4)
    req(input$numconduser1)
    req(input$numconduser2)

    con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconduser1])
    con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconduser2])

    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))


    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con2]] - .data[[con1]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
                          mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con2]] - .data[[con1]],label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con2,"-",con1,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }

    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con1]] - .data[[con2]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
                          mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con1]] - .data[[con2]],label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con1,"-",con2,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }

    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "No") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con2]] - .data[[con1]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        # geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
        #                   mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con2]] - .data[[con1]],label=label,color=UserAnnotation),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con2,"-",con1,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }

    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "No") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con1]] - .data[[con2]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        # geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
        #                   mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con1]] - .data[[con2]],label=label,color=UserAnnotation),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con1,"-",con2,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }

    return(p1)


  })
  
  output$printcategviolin <- downloadHandler(
    filename = function() { paste("UserDefined_ViolinPlot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = categviolinprint(), width = 16, height = 8)
    }
  )
  
  
  # output$corruserplot <- renderPlot({
  # 
  #   req(input$file4)
  #   req(input$numconduser1)
  #   req(input$numconduser2)
  # 
  #   con1<-as.character(unique(norm_sum_final2()$Condition)[input$numconduser1])
  #   con2<-as.character(unique(norm_sum_final2()$Condition)[input$numconduser2])
  # 
  #   foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
  #   qvaluename<-paste("q.val_",con2,"-",con1,sep="")
  #   sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
  # 
  # 
  #   ggplot()+
  #     # geom_hline(yintercept = 0,linetype="dashed")+
  #     #geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(x=.data[[con2]] - .data[[con1]],y=uservalue),size=0.25,alpha=0.2,color="grey")+
  #     geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(x=.data[[con2]] - .data[[con1]],y=uservalue,color=UserAnnotation),size=1,alpha=0.5)+
  #     theme_classic()+scale_color_viridis_d(end=0.8)+
  #     labs(y="User-defined Values", x=paste("Log2FC(",con2,"-",con1,")",sep=""),title="User Correlation Plot")+
  #     theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
  # 
  # 
  # })
  
  dfuser2 <- reactive({
    req(input$file6)
    inFile <- input$file6
    if (is.null(inFile))
      return(NULL)
    tbl <- read.csv(inFile$datapath)

    return(tbl)
  })
  
  
  localizationvolcanomaker <-reactive({

    if (input$txthumanmouse == "mouse") {
      ann_combine<-dfuser2()%>%
        # dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale(),by="Gene.Symbol")
    }

    if (input$txthumanmouse == "human") {
      ann_combine<-dfuser2()%>%
        # dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale(),by="Gene.Symbol")%>%
        dplyr::mutate(Compartment = dplyr::if_else(is.na(Compartment), "__",Compartment))
    }
    return(ann_combine)
  })
  
  
  output$plotvolcanomaker<- renderPlot({
    req(input$file6)
    req(input$size10)
    req(input$size11)
    req(input$alpha10)
    req(input$alpha11)


    userinput3<-dfuser2()%>%
      dplyr::mutate(sig = dplyr::if_else(abs(Log2FC)>log2(input$num11),dplyr::if_else(p.value<input$num10,dplyr::if_else(Log2FC>0,"up","down"),"n.s."),"n.s."))
  ### need to finish the plotting aes

    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }

    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=labelplot),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }

    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }

    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }


    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha11,size=input$size11,color="red")+
        geom_text_repel(data=localizationvolcanomaker()%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color="red"),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }

    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha11,size=input$size11,color="red")+
        # geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
        #                   mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color="red"),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }



    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "no"), aes(Log2FC, p.value,color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "yes"), aes(Log2FC, p.value,color=sig),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, p.value,label=label,color=sig),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }

    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, p.value,label=label,color=labelplot),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }

    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, p.value,color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, p.value,color=sig),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }

    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }


    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha11,size=input$size11,color="red")+
        geom_text_repel(data=localizationvolcanomaker()%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, p.value,label=label,color="red"),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }

    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha11,size=input$size11,color="red")+
        # geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
        #                   mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color="red"),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }

    p1
  },width=1000,height=650)
  
  
  plotvolcanomakerprint<- reactive({
    req(input$file6)
    req(input$size10)
    req(input$size11)
    req(input$alpha10)
    req(input$alpha11)



    userinput3<-dfuser2()%>%
      dplyr::mutate(sig = dplyr::if_else(abs(Log2FC)>log2(input$num11),dplyr::if_else(p.value<input$num10,dplyr::if_else(Log2FC>0,"up","down"),"n.s."),"n.s."))

    ### need to finish the plotting aes

    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }

    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=labelplot),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }

    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }

    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }


    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha11,size=input$size11,color="red")+
        geom_text_repel(data=localizationvolcanomaker()%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color="red"),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }

    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha11,size=input$size11,color="red")+
        # geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
        #                   mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color="red"),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }



    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "no"), aes(Log2FC, p.value,color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "yes"), aes(Log2FC, p.value,color=sig),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, p.value,label=label,color=sig),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }

    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, p.value,label=label,color=labelplot),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }

    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, p.value,color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, p.value,color=sig),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }

    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }


    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha11,size=input$size11,color="red")+
        geom_text_repel(data=localizationvolcanomaker()%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, p.value,label=label,color="red"),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }

    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha11,size=input$size11,color="red")+
        # geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
        #                   mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color="red"),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }

    return(p1)
  })
  
  output$printplotvolcanomaker <- downloadHandler(
    filename = function() { paste("UserDefined_Volcano", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotvolcanomakerprint(), width = 12, height = 8)
    }
  )
  
  
  output$exampleUserInput <- downloadHandler(
    filename = function() {
      paste("Example_UserInput_file",".csv", sep = "")
    },
    content = function(file) {
      write.csv(userinput_example, file, row.names = FALSE)
    }
  )
  
  
  
  
  
}


# Run the application 
shinyApp(ui = ui, server = server)