library(shiny)#shiny library
library(DT)data tables library
chrlens<-c(249212500,243031394,197894711,190949272,180688758,170894141 ,158938491,
	146279898,141213431,135534747,135006516,133851895,115169878,107349540,
	102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,
	155270560,59373566)#declare chrmosome lenghs
shinyUI(navbarPage(theme = "bootstrap.css",windowTitle = "MethSurv-A web tool to perform multivariable survival analysis using DNA methylation data", footer = includeHTML("www/footer.html"),
                   headerPanel(""),
                   
                   tabPanel(h2("MethSurv : A web tool to perform multivariable survival analysis using DNA methylation data
                               "),title="MethSurv",
                            sidebarLayout(       
                              sidebarPanel(
                                bsTooltip("file", "Choose a cancer type from the list for survival analysis", placement = "bottom"),
                                bsTooltip("gene", "Type a gene symbol. Available genes will be automatically displayed", placement = "bottom"),
                                bsTooltip("Gene", "Type a gene symbol. Available genes will be automatically displayed for all the cancers", placement = "bottom"),
                                bsTooltip("island","Choose CpG island region",placement = "bottom"),
                                bsTooltip("region","Choose gene sub-region",placement = "bottom"),
                                bsTooltip("probe","Choose a CpG site (probe)",placement = "bottom"),
                                bsTooltip("split", "Choose a method for splitting the patients", placement = "bottom"),
                                bsTooltip("covariate","Choose the co-variate(s) for multivariate analysis",placement="bottom"),
                                
                                bsTooltip("cancer3", "Choose a cancer type from the list for gene visualization", placement="bottom"),
                                bsTooltip("Gene3","Type a gene symbol. Available genes will be automatically displayed", placement="bottom"),
                                bsTooltip("cancer4", "Choose a cancer type from the list for the Top biomakers view", placement="bottom"),
				bsTooltip("cancer5", "Choose a cancer type from the list for the region based analysis", placement="bottom"),
				bsTooltip("chr", "Specify a chromosome", placement="bottom"),
                                bsTooltip("adjustment", "Choose an adjustment type from the list", placement="bottom"),
                                conditionalPanel(condition = "input.tabs1 == 'Single CpG'",
                                                 selectInput('file', 'TCGA cancer datasets',list("Acute Myeloid Leukemia [LAML] TCGA March 2017"="LAML",
                                                                                                 "Adrenocortical carcinoma [ACC] TCGA March 2017"="ACC", 
                                                                                                 "Bladder Urothelial Carcinoma [BLCA] TCGA March 2017"="BLCA" ,
                                                                                                 "Brain Lower Grade Glioma [LGG] TCGA March 2017"= "LGG",
                                                                                                 "Breast invasive carcinoma [BRCA] TCGA March 2017"= "BRCA",
                                                                                                 "Cervical squamous cell carcinoma and endocervical adenocarcinoma [CESC] TCGA March 2017"= "CESC",
                                                                                                 "Colon adenocarcinoma [COAD] TCGA March 2017"= "COAD",
                                                                                                 "Esophageal carcinoma [ESCA] TCGA March 2017"= "ESCA",
                                                                                                 "Glioblastoma multiforme [GBM] TCGA March 2017"= "GBM",
                                                                                                 "Head and Neck squamous cell carcinoma  [HNSC] TCGA March 2017"= "HNSC",
                                                                                                 "Kidney Chromophobe  [KICH] TCGA March 2017"="KICH",
                                                                                                 "Kidney renal clear cell carcinoma  [KIRC] TCGA March 2017"="KIRC",
                                                                                                 "Kidney renal papillary cell carcinoma  [KIRP] TCGA March 2017"="KIRP",
                                                                                                 "Liver hepatocellular carcinoma [LIHC] TCGA March 2017"="LIHC",
                                                                                                 "Lung adenocarcinoma [LUAD] TCGA March 2017"="LUAD",
                                                                                                 "Lung squamous cell carcinoma  [LUSC] TCGA March 2017"="LUSC",
                                                                                                 "Mesothelioma [MESO] TCGA March 2017"= "MESO", 
                                                                                                 "Pancreatic adenocarcinoma [PAAD] TCGA March 2017"= "PAAD",
                                                                                                 "Rectum adenocarcinoma [READ] TCGA March 2017"= "READ",
                                                                                                 "Sarcoma [SARC] TCGA March 2017"= "SARC",
                                                                                                 "Skin Cutaneous Melanoma [SKCM] TCGA March 2017"= "SKCM",
												"Stomach Adenocarcinoma [STAD] TCGA March 2017"= "STAD",	
                                                                                                 "Uterine carcinosarcoma [UCS] TCGA March 2017"="UCS",
                                                                                                 "Uterine Corpus Endometrial Carcinoma [UCEC] TCGA March 2017"="UCEC",
                                                                                                 "Uveal Melanoma [UVM] TCGA March 2017"="UVM"),selected="UCS"
                                                 ),
                                                 helpText("Select a cancer type to explore"),
                                                 #Quality control exploration:</span></p>"),
                                                 selectizeInput("gene", "Gene", choices = "", selected = "",multiple=FALSE,options = list(placeholder = "Search for a gene")),              
                                                 
                                                 helpText("Query a gene"),
                                                 #selectInput('gene','Gene',choices=""),
                                                 
                                                 #selectizeInput("gene", "Gene", choices = "NULL", selected = NULL,multiple=FALSE,options = list(maxOptions = 0),),
                                                 #selectizeInput("gene", "Gene", choices = NULL, selected = NULL,multiple=FALSE,options = list(placeholder = "Search for a gene")),
                                                 selectizeInput('island','Relation to island',choices="",multiple=FALSE), 
                                                 #selectizeInput('genomic_region','Genomic Region',choices="",multiple=FALSE),
                                                 
                                                 #textInput("gene", label = "Enter gene", value = "MGMT"),
                                                 
                                                 
                                                 helpText("Choose island region from the list"),
                                                 #selectInput('split','Split by',choices=""),
                                                 #helpText("Split the patients based on"),
                                                 #selectizeInput('island','Relation to island',choices="",multiple=FALSE),
                                                 
                                                 #textInput("gene", label = "Enter gene", value = "MGMT"),
                                                 
                                                 selectizeInput('region','Genomic Region',choices="",multiple=FALSE),
                                                 helpText("Choose gene region from the list"),
                                                 selectizeInput('probe','CpG site',choices="",multiple=FALSE),
                                                 
                                                 #textInput("gene", label = "Enter gene", value = "MGMT"),
                                                 
                                                 
                                                 helpText("Choose CpG site from the list"),
                                                 selectizeInput('split','Split by',choices=""),
                                                 helpText("Choose splitting option to dichotomize methylation profiles of patients: mean, median, upper and lower quantiles (q25 and q75), maxstat and best among them"),
                                                 checkboxInput("includeCovariate", "include covariate", FALSE),
                                                 helpText("Would you like to adjust for co-variate(s)?"),
                                                 
                                                 conditionalPanel(condition = "input.includeCovariate",
                                                                  selectInput('covariate','Covariate adjustment',choices="",multiple = TRUE),
                                                                  #checkboxGroupInput('covariate','Covariate Adjustment',choices=NULL)					
                                                                  
                                                                  helpText("Add one or more covariates from the list")
                                                                  
                                                 )
                                                 
                                                 
                                                 
                                                 
                                ),


conditionalPanel(condition = "input.tabs1 == 'Region based analysis'",
                                                 
                                                selectInput('cancer5','Cancer' ,list(
                                                                                                 "Acute Myeloid Leukemia [LAML] TCGA March 2017"="LAML",
                                                                                                 "Adrenocortical carcinoma [ACC] TCGA March 2017"="ACC", 
                                                                                                 "Bladder Urothelial Carcinoma [BLCA] TCGA March 2017"="BLCA" ,
                                                                                                 "Brain Lower Grade Glioma [LGG] TCGA March 2017"= "LGG",
                                                                                                 "Breast invasive carcinoma [BRCA] TCGA March 2017"= "BRCA",
                                                                                                 "Cervical squamous cell carcinoma and endocervical adenocarcinoma [CESC] TCGA March 2017"= "CESC",
                                                                                                 "Colon adenocarcinoma [COAD] TCGA March 2017"= "COAD",
                                                                                                 "Esophageal carcinoma [ESCA] TCGA March 2017"= "ESCA",
                                                                                                 "Glioblastoma multiforme [GBM] TCGA March 2017"= "GBM",
                                                                                                 "Head and Neck squamous cell carcinoma  [HNSC] TCGA March 2017"= "HNSC",
                                                                                                 "Kidney Chromophobe  [KICH] TCGA March 2017"="KICH",
                                                                                                 "Kidney renal clear cell carcinoma  [KIRC] TCGA March 2017"="KIRC",
                                                                                                 "Kidney renal papillary cell carcinoma  [KIRP] TCGA March 2017"="KIRP",
                                                                                                 "Liver hepatocellular carcinoma [LIHC] TCGA March 2017"="LIHC",
                                                                                                 "Lung adenocarcinoma [LUAD] TCGA March 2017"="LUAD",
                                                                                                 "Lung squamous cell carcinoma  [LUSC] TCGA March 2017"="LUSC",
                                                                                                 "Mesothelioma [MESO] TCGA March 2017"= "MESO", 
                                                                                                 "Pancreatic adenocarcinoma [PAAD] TCGA March 2017"= "PAAD",
                                                                                                 "Rectum adenocarcinoma [READ] TCGA March 2017"= "READ",
                                                                                                 "Sarcoma [SARC] TCGA March 2017"= "SARC",
                                                                                                 "Skin Cutaneous Melanoma [SKCM] TCGA March 2017"= "SKCM",
												"Stomach Adenocarcinoma [STAD] TCGA March 2017"= "STAD",	
                                                                                                 "Uterine carcinosarcoma [UCS] TCGA March 2017"="UCS",
                                                                                                 "Uterine Corpus Endometrial Carcinoma [UCEC] TCGA March 2017"="UCEC",
                                                                                                 "Uveal Melanoma [UVM] TCGA March 2017"="UVM"),selected="UCS",multiple=FALSE) ,
                                helpText("Select a cancer type"),
selectInput("chr", "Chromosome", list( "chr1" = "1",
				"chr2" = "2",  "chr3" = "3",  "chr4" = "4", "chr5" = "5",
				"chr6" = "6",  "chr7" = "7",  "chr8" = "8",  "chr9" = "9",
				"chr10" = "10",  "chr11" = "11",  "chr12" = "12",  "chr13" = "13",
				"chr14" = "14",  "chr15" = "15",  "chr16" = "16",  "chr17" = "17",	
				"chr18" = "18" ,  "chr19" = "19",  "chr20" = "20",  "chr21" = "21",
				 "chr22" = "22",  "chrX" = "X",  "chrY" = "Y"),selected="21"),

		conditionalPanel(condition = "input.chr == '1'",
			sliderInput("range1", "Choose basepair position [hg19]",
                	min = 15865, max = 249212500, value = c(1,chrlens[1]),width="1000")),
		conditionalPanel(condition = "input.chr == '2'",
			sliderInput("range2", "Choose basepair position [hg19]",
                	min = 45264, max = 24303139, value = c(1,chrlens[2]),width="1000")),
		conditionalPanel(condition = "input.chr == '3'",
			sliderInput("range3", "Choose basepair position [hg19]",
                	min = 238048 , max = 197894711, value = c(1,chrlens[3]),width="1000")),
		conditionalPanel(condition = "input.chr == '4'",
			sliderInput("range4", "Choose basepair position [hg19]",
                	min = 53010, max = 190949272 , value = c(1,chrlens[4]),width="1000")),
		conditionalPanel(condition = "input.chr == '5'",
			sliderInput("range5", "Choose basepair position [hg19]",
                	min = 139839, max = 180688758 , value = c(1,chrlens[5]),width="1000")),
		conditionalPanel(condition = "input.chr == '6'",
			sliderInput("range6", "Choose basepair position [hg19]",
                	min = 290800, max = 170894141, value = c(1,chrlens[6]),width="1000")),
		conditionalPanel(condition = "input.chr == '7'",
			sliderInput("range7", "Choose basepair position [hg19]",
                	min = 191704, max = 158938491, value = c(1,chrlens[7]),width="1000")),
		conditionalPanel(condition = "input.chr == '8'",
			sliderInput("range8", "Choose basepair position [hg19]",
                	min = 176533, max = 146279898, value = c(1,chrlens[8]),width="1000")),
		conditionalPanel(condition = "input.chr == '9'",
			sliderInput("range9", "Choose basepair position [hg19]",
                	min = 117044, max = 141111396, value = c(1,chrlens[9]),width="1000")),
		conditionalPanel(condition = "input.chr == '10'",
			sliderInput("range10", "Choose basepair position [hg19]",
                	min = 93193, max = 135441372, value = c(1,chrlens[10]),width="1000")),
		conditionalPanel(condition = "input.chr == '11'",
			sliderInput("range11", "Choose basepair position [hg19]",
                	min = 127495, max = 134283252, value = c(1,chrlens[11]),width="1000")),
		conditionalPanel(condition = "input.chr == '12'",
			sliderInput("range12", "Choose basepair position [hg19]",
                	min = 149764, max = 133779169, value = c(1,chrlens[12]),width="1000")),
		conditionalPanel(condition = "input.chr == '13'",
			sliderInput("range13", "Choose basepair position [hg19]",
                	min = 19412184, max = 115092495, value = c(1,chrlens[13]),width="1000")),
		conditionalPanel(condition = "input.chr == '14'",
			sliderInput("range14", "Choose basepair position [hg19]",
                	min = 19377109 , max = 106941205 , value = c(1,chrlens[14]),width="1000")),
		conditionalPanel(condition = "input.chr == '15'",
			sliderInput("range15", "Choose basepair position [hg19]",
                	min = 20747200 , max = 102513445 , value = c(1,chrlens[15]),width="1000")),
		conditionalPanel(condition = "input.chr == '16'",
			sliderInput("range16", "Choose basepair position [hg19]",
                	min = 101260, max = 90143815, value = c(1,chrlens[16]),width="1000")),
		conditionalPanel(condition = "input.chr == '17'",
			sliderInput("range17", "Choose basepair position [hg19]",
                	min = 6018 , max = 81052243 , value = c(1,chrlens[17]),width="1000")),
		conditionalPanel(condition = "input.chr == '18'",
			sliderInput("range18", "Choose basepair position [hg19]",
                	min=158229,max=78005665, value = c(1,chrlens[18]),width="1000")),
		conditionalPanel(condition = "input.chr == '19'",
			sliderInput("range19", "Choose basepair position [hg19]",
                	min=200092,max= 59093449, value = c(1,chrlens[19]),width="1000")),
		conditionalPanel(condition = "input.chr == '20'",
			sliderInput("range20", "Choose basepair position [hg19]",
                	min = 68170, max = 62903185, value = c(1,chrlens[20]),width="1000")),
		conditionalPanel(condition = "input.chr == '21'",
			sliderInput("range21", "Choose basepair position [hg19]",
                	min = 10988071, max = 48084674, value = c(1,chrlens[21]),width="1000")),
		conditionalPanel(condition = "input.chr == '22'",
			sliderInput("range22", "Choose basepair position [hg19]",
                	min = 16256627, max = 51225561, value = c(1,chrlens[22]),width="1000")),
		conditionalPanel(condition = "input.chr == 'X'",
			sliderInput("range23", "Choose basepair position [hg19]",
                	min = 2715017, max = 154843228, value = c(1,chrlens[23]),width="1000")),
		conditionalPanel(condition = "input.chr == 'Y'",
			sliderInput("range24", "Choose basepair position [hg19]",
                	min = 2709598, max = 21238472, value = c(1,chrlens[24]),width="1000"))
),





                                conditionalPanel(condition = "input.tabs1 == 'Gene visualization'",selectInput('cancer3', 'Cancer',list("Acute Myeloid Leukemia [LAML] TCGA March 2017"="LAML",
                                                                                                                                        "Adrenocortical carcinoma [ACC] TCGA March 2017"="ACC", 
                                                                                                                                        "Bladder Urothelial Carcinoma [BLCA] TCGA March 2017"="BLCA" ,
                                                                                                                                        "Brain Lower Grade Glioma [LGG] TCGA March 2017"= "LGG",
                                                                                                                                        "Breast invasive carcinoma [BRCA] TCGA March 2017"= "BRCA",
                                                                                                                                        "Cervical squamous cell carcinoma and endocervical adenocarcinoma [CESC] TCGA March 2017"= "CESC",
                                                                                                                                        "Colon adenocarcinoma [COAD] TCGA March 2017"= "COAD",
                                                                                                                                        "Esophageal carcinoma [ESCA] TCGA March 2017"= "ESCA",
                                                                                                                                        "Glioblastoma multiforme [GBM] TCGA March 2017"= "GBM",
                                                                                                                                        "Head and Neck squamous cell carcinoma  [HNSC] TCGA March 2017"= "HNSC",
                                                                                                                                        "Kidney Chromophobe  [KICH] TCGA March 2017"="KICH",
                                                                                                                                        "Kidney renal clear cell carcinoma  [KIRC] TCGA March 2017"="KIRC",
                                                                                                                                        "Kidney renal papillary cell carcinoma  [KIRP] TCGA March 2017"="KIRP",
                                                                                                                                        "Liver hepatocellular carcinoma [LIHC] TCGA March 2017"="LIHC",
                                                                                                                                        "Lung adenocarcinoma [LUAD] TCGA March 2017"="LUAD",
                                                                                                                                        "Lung squamous cell carcinoma  [LUSC] TCGA March 2017"="LUSC",
                                                                                                                                        "Mesothelioma [MESO] TCGA March 2017"= "MESO", 
                                                                                                                                        "Pancreatic adenocarcinoma [PAAD] TCGA March 2017"= "PAAD",
                                                                                                                                        "Rectum adenocarcinoma [READ] TCGA March 2017"= "READ",
                                                                                                                                        "Sarcoma [SARC] TCGA March 2017"= "SARC",
                                                                                                                                        "Skin Cutaneous Melanoma [SKCM] TCGA March 2017"= "SKCM",
                                                                                                                                        "Stomach Adenocarcinoma [STAD] TCGA March 2017"="STAD",
                                                                                                                                        "Uterine carcinosarcoma [UCS] TCGA March 2017"="UCS",
                                                                                                                                        "Uterine Corpus Endometrial Carcinoma [UCEC] TCGA March 2017"="UCEC",
                                                                                                                                        "Uveal Melanoma [UVM] TCGA March 2017"="UVM"),selected="UCS"
                                ),
                                helpText("Select a cancer type to draw heatmap"),
                                #Quality control exploration:</span></p>"),
                                
                                selectizeInput('Gene3','Gene', choices = "", selected = "",multiple=FALSE,options = list(placeholder = "Search for a gene")) 
                                ),
                                
                                
                                
                                
                                #actionButton("update", "Update View"),
                                
                                conditionalPanel(condition = "input.tabs1 == 'All cancers'",
                                                 
                                                 selectizeInput('Gene','Gene', choices = "", selected = "",multiple=FALSE,options = list(placeholder = "Search for a gene")) 
                                ),
                                
                                conditionalPanel(condition = "input.tabs1 == 'Top biomarkers'",
                                                 
                                                 selectInput('cancer4','Cancer' ,list("Acute Myeloid Leukemia [LAML] TCGA March 2017"="LAML",
                                                                                      "Adrenocortical carcinoma [ACC] TCGA March 2017"="ACC", 
                                                                                      "Bladder Urothelial Carcinoma [BLCA] TCGA March 2017"="BLCA" ,
                                                                                      "Brain Lower Grade Glioma [LGG] TCGA March 2017"= "LGG",
                                                                                      "Breast invasive carcinoma [BRCA] TCGA March 2017"= "BRCA",
                                                                                      "Cervical squamous cell carcinoma and endocervical adenocarcinoma [CESC] TCGA March 2017"= "CESC",
                                                                                      "Colon adenocarcinoma [COAD] TCGA March 2017"= "COAD",
                                                                                      "Esophageal carcinoma [ESCA] TCGA March 2017"= "ESCA",
                                                                                      "Glioblastoma multiforme [GBM] TCGA March 2017"= "GBM",
                                                                                      "Head and Neck squamous cell carcinoma  [HNSC] TCGA March 2017"= "HNSC",
                                                                                      "Kidney Chromophobe  [KICH] TCGA March 2017"="KICH",
                                                                                      "Kidney renal clear cell carcinoma  [KIRC] TCGA March 2017"="KIRC",
                                                                                      "Kidney renal papillary cell carcinoma  [KIRP] TCGA March 2017"="KIRP",
                                                                                      "Liver hepatocellular carcinoma [LIHC] TCGA March 2017"="LIHC",
                                                                                      "Lung adenocarcinoma [LUAD] TCGA March 2017"="LUAD",
                                                                                      "Lung squamous cell carcinoma  [LUSC] TCGA March 2017"="LUSC",
                                                                                      "Mesothelioma [MESO] TCGA March 2017"= "MESO", 
                                                                                      "Pancreatic adenocarcinoma [PAAD] TCGA March 2017"= "PAAD",
                                                                                      "Rectum adenocarcinoma [READ] TCGA March 2017"= "READ",
                                                                                      "Sarcoma [SARC] TCGA March 2017"= "SARC",
                                                                                      "Skin Cutaneous Melanoma456789 [SKCM] TCGA March 2017"= "SKCM",
                                                                                      "Stomach Adenocarcinoma 
 [STAD] TCGA March 2017"="STAD",
                                                                                      "Uterine carcinosarcoma [UCS] TCGA March 2017"="UCS",
                                                                                      "Uterine Corpus Endometrial Carcinoma [UCEC] TCGA March 2017"="UCEC",
                                                                                      "Uveal Melanoma [UVM] TCGA March 2017"="UVM"), selected = "UCS",multiple=FALSE) ,
                                                 helpText("Select a cancer type to explore top CpG biomarkers")
                                )
                                
                              ),
                              #######
                              mainPanel(
                                #progressInit(),
                                
                                ##plotOutput("plot1"),
                                ##br(),br(),
                                ##tableOutput("estimates")
                                
                                ###here try tab###
                                
                                tabsetPanel(
                                  tabPanel("Single CpG",id = "Survival_plots",br(), h4(textOutput("caption3")),br(), 
                                           #plotOutput("plot1"),br(),br(),
                                           tags$style(type="text/css"
                                                      #".shiny-output-error { visibility: hidden; }",#hide errors from the front end
                                                      #".shiny-output-error:before { visibility: hidden; }"
                                           ),                plotOutput("plot1")  , br(),plotOutput("plot3"),br(),
                                           selectizeInput("Feature", "Display methylation profiles based on", choices = "", selected = "",multiple=FALSE,options = list(placeholder = "Search for a feature")), 
                                           plotOutput("plot4"),br(),br(),br(),br(),
                                           h4(textOutput("caption1")),br(),
                                           tableOutput("estimates"),br(),
                                           downloadButton("downloadkmpng", label = "Save KM plot as png"),
                                           downloadButton("downloadkmpdf", label = "Save KM plot as pdf"),
                                           downloadButton("downloaddenpng", label = "Save density plot as png"),
                                           downloadButton("downloaddenpdf", label = "Save density plot as pdf"),
                                           downloadButton("downloadviopng", label = "Save violin plot as png"),
                                           downloadButton("downloadviopdf", label = "Save violin plot as pdf"),
                                           downloadButton("downloadUserExpressionFile", label = "Save as tab delimited .txt"),
                                           uiOutput("link"),
                                           uiOutput("link2"),
                                           uiOutput("link3")
                                           
                                  ),
                                   tabPanel("Region based analysis",id="region_based_analysis",br(),br(),h4(textOutput("caption4")),br(),br(),
                                           dataTableOutput("region_res"),br(),br(), selectizeInput("probe_sel", "Draw survival curve for the CpG", choices = "", selected = "",multiple=FALSE,options = list(placeholder = "Choose a CpG")),plotOutput("plot5"), 
                                           #tableOutput("allcan"),br(),br(),
                                           downloadButton("regiontable", label = "Save summary from the chosen region as tab delimited .txt"),downloadButton("downloadkm2png",label = "Save KM plot as png"),downloadButton("downloadkm2pdf",label = "Save KM plot as pdf")
                                           #plotOutput("plot2")
                                  ),
                                  tabPanel("All cancers",id = "All_cancers",br(),br(),
                                           dataTableOutput("allcan"),br(),br(),uiOutput("popup2"),
                                           #tableOutput("allcan"),br(),br(),
                                           downloadButton("allcanTable", label = "Save as tab delimited .txt")
                                           #plotOutput("plot2")
                                  ),
                                  tabPanel("Top biomarkers",id = "Top_biomarkers",br(),br(),
                                           dataTableOutput("topbm"),br(),br(),uiOutput("popup1"),
                                           #tableOutput("allcan"),br(),br(),
                                           downloadButton("topbmtable", label = "Save as tab delimited .txt")
                                  ),
                                  tabPanel("Gene visualization",id = "Gene_visualization",br(),br(), h4(textOutput("caption2")),br(),br(),
                                           #dataTableOutput("table1")
                                           #verbatimTextOutput("test"),
                                           plotOutput("plot2",height="100%"),br(),br(),
                                           downloadButton("downloadHMpdf",label = "Save heatmap as pdf"),
                                           downloadButton("downloadHMpng",label = "Save heatmap as png"),
                                           actionButton("clustVisLink", label = "Browse in ClustVis"),
                                           conditionalPanel(condition = "input.clustVisLink",
                                                            uiOutput("clustVisLinkText")
                                           )
                                  ),
                                  
                                  
                                  id = "tabs1" 
                                  
                                )
                                
                                
                              )
                              
                            )),
                   tabPanel("About",
                            h1("About", align = "center"),
                            includeHTML("www/home.html"),
                            p(a(img(src = "UT_logo.png", height="100",width="100"), href = "http://www.ut.ee/en"),
                              a(img(src = "BIIT_logo.png", width="150"), href = "http://biit.cs.ut.ee/"), align = "center")
                   ), 
                   tabPanel("Quick start",
                            h1("Quick start", align = "center"),
                            includeHTML("www/quickstart.html"),
                            p(a(img(src = "UT_logo.png", height="100",width="100"), href = "http://www.ut.ee/en"),
                              a(img(src = "BIIT_logo.png", width="150"), href = "http://biit.cs.ut.ee/"), align = "center")
                   ),
                   tabPanel("Tutorial",
                            h1("Tutorial", align = "center"),
                            includeHTML("www/tutorial.html"),
                            p(a(img(src = "UT_logo.png", height="100",width="100"), href = "http://www.ut.ee/en"),
                              a(img(src = "BIIT_logo.png", width="150"), href = "http://biit.cs.ut.ee/"), align = "center")
                   ),
                   tabPanel("Download",
                            h1("Data download page", align = "center"),
                            includeHTML("www/download.html"),
                            p(a(img(src = "UT_logo.png", height="100",width="100"), href = "http://www.ut.ee/en"),
                              a(img(src = "BIIT_logo.png", width="150"), href = "http://biit.cs.ut.ee/"), align = "center")
                   ),
                   tabPanel("FAQ",
                            h1("FAQ", align = "center"),
                            includeHTML("www/FAQ.html"),
                            p(a(img(src = "UT_logo.png", height="100",width="100"), href = "http://www.ut.ee/en"),
                              a(img(src = "BIIT_logo.png",width="150"), href = "http://biit.cs.ut.ee/"),align = "center")
                   )
))
