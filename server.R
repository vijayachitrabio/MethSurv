#load the requied libaries
library(survminer)
library(shiny)
library(survival)
library(stringr)
library(shinyBS)
library(DT)
library("survMisc")
library(RColorBrewer)
options(shiny.sanitize.errors = FALSE)
shinyServer(function(input,output,session){

  
  #getting cancer genes from the data
  cang <- reactive({
    can_genes=fread("/srv/shiny-server/cancer_mat_genes",header=T,sep="\t")
    can_genes=as.character(unique(can_genes$gene))
    can_genes
  })
  #extract caner data for user selected cancer type
  filedata1 <- reactive({
    
    meth=get(paste0(input$file,"_meth")) 
    meth3=get(paste0(input$cancer3,"_meth"))
    meth3=subset(meth3,UCSC_RefGene_Name!="Unknown")
    meth4=get(paste0(input$cancer4,"_meth"))
    meth5=get(paste0(input$cancer5,"_meth"))
    cn=as.character(unique(meth$UCSC_RefGene_Name))
    cn_new=as.character(unique(meth3$UCSC_RefGene_Name))
    can_genes = cang()
    clinic=get(paste0(input$file,"_clinic_new"))
   colnames(clinic)=gsub("gender","sex",colnames(clinic)) 
    clinic4=get(paste0(input$cancer4,"_clinic_new"))
colnames(clinic4)=gsub("gender","sex",colnames(clinic4)) 
clinic5=get(paste0(input$cancer5,"_clinic_new")) 
    vio=get(paste0(input$file,"_vio"))
colnames(vio)=gsub("gender","sex",colnames(vio)) 
    vio=vio[,!colnames(vio) %in%"event"]				
    cov_vector<- get(paste0(input$file,"_covariate"))
    cov_vector=sort(cov_vector,decreasing=TRUE)
    annot=get(paste0(input$cancer3,"_annot"))	
    path="/srv/shiny-server/genes/"
    if(input$Gene == ""){
      canmat = NULL
    } else {
      canmat=fread(paste0(path,input$Gene,".txt"))
    }
    
 #get top biomarker list
    
    
    top=get(paste0(input$cancer4,"_top_2017"))
    #round the p.values with significant data points
    top$adj.P.value= signif(top$adj.P.value,2)
    top$P.value= signif(top$P.value,2)
    top$LR_test_pvalue= signif(top$LR_test_pvalue,2)
    top=subset(top,select=c("Name","HR","CI","P.value","adj.P.value","LR_test_pvalue","Best_split","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island"))
    meth_split_threshold=c("mean","median","q75","q25","maxstat","best")
gw=get(paste0(input$cancer5,"_null_2017"))
chrlens<-c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,
	146364022,141213431,135534747,135006516,133851895,115169878,107349540,
	102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,
	155270560,59373566)
    info=list(meth=meth,clinic=clinic,clinic4=clinic4,cn=cn,cov_vector=cov_vector,annot=annot,canmat=canmat,can_genes=can_genes,meth_split_threshold=meth_split_threshold,
              meth3=meth3,meth4=meth4,cn_new=cn_new,top=top,vio=vio,chrlens=chrlens,gw=gw,meth5=meth5,clinic5=clinic5)
    return(info)
  })
  observe({
    #can_genes<-filedata1()$can_genes
    can_genes = cang()
    if (is.null(can_genes))
      return(NULL)
    can_genes=as.character(can_genes)
    updateSelectizeInput(session, "Gene",selected=can_genes[1], choices = can_genes,server=TRUE)
  })
  
  
  
  
  
  
  ### for top biomarkers
  # observe({
  #   adj<-filedata1()$adj 
  #   if (is.null(adj))
  #     return(NULL)
  #   adj=as.character(adj)
  #   updateSelectizeInput(session, "adjustment",selected=NULL, choices = adj,server=TRUE)
  # })
  ####
  
  ### for violin plots
  observe({
    vio<-filedata1()$vio 
    if (is.null(vio))
      return(NULL)
    vio=as.character(names(vio))
    updateSelectizeInput(session, "Feature",selected="age", choices = vio,server=TRUE)
  })


### for kmplot in region-tab
  observe({
    region_meth<-filedata1()$meth5
    if (is.null(region_meth))
      return(NULL)
  
    if(input$chr == 1){range_start<-input$range1[1];range_end<-input$range1[2]}
		if(input$chr == 2){range_start<-input$range2[1];range_end<-input$range2[2]}
		if(input$chr == 3){range_start<-input$range3[1];range_end<-input$range3[2]}
		if(input$chr == 4){range_start<-input$range4[1];range_end<-input$range4[2]}
		if(input$chr == 5){range_start<-input$range5[1];range_end<-input$range5[2]}
		if(input$chr == 6){range_start<-input$range6[1];range_end<-input$range6[2]}
		if(input$chr == 7){range_start<-input$range7[1];range_end<-input$range7[2]}
		if(input$chr == 8){range_start<-input$range8[1];range_end<-input$range8[2]}
		if(input$chr == 9){range_start<-input$range9[1];range_end<-input$range9[2]}
		if(input$chr == 10){range_start<-input$range10[1];range_end<-input$range10[2]}
		if(input$chr == 11){range_start<-input$range11[1];range_end<-input$range11[2]}
		if(input$chr == 12){range_start<-input$range12[1];range_end<-input$range12[2]}
		if(input$chr == 13){range_start<-input$range13[1];range_end<-input$range13[2]}
		if(input$chr == 14){range_start<-input$range14[1];range_end<-input$range14[2]}
		if(input$chr == 15){range_start<-input$range15[1];range_end<-input$range15[2]}
		if(input$chr == 16){range_start<-input$range16[1];range_end<-input$range16[2]}
		if(input$chr == 17){range_start<-input$range17[1];range_end<-input$range17[2]}
		if(input$chr == 18){range_start<-input$range18[1];range_end<-input$range18[2]}
		if(input$chr == 19){range_start<-input$range19[1];range_end<-input$range19[2]}
		if(input$chr == 20){range_start<-input$range20[1];range_end<-input$range20[2]}
		if(input$chr == 21){range_start<-input$range21[1];range_end<-input$range21[2]}
		if(input$chr == 22){range_start<-input$range22[1];range_end<-input$range22[2]}
		if(input$chr == "X"){range_start<-input$range23[1];range_end<-input$range23[2]}
		if(input$chr == "Y"){range_start<-input$range24[1];range_end<-input$range24[2]}
    #region_meth<-subset(region_meth,CHR==input$chr & MAPINFO >= range_start & MAPINFO <= range_end )
region_meth=data.table(region_meth)
region_meth=region_meth[CHR==as.character(input$chr) & MAPINFO >= range_start & MAPINFO<=range_end,]
cpgs=as.character(region_meth$Name)
 
  
    updateSelectizeInput(session, "probe_sel", choices=cpgs, selected=cpgs[1],server=TRUE)
  })
  #####
  observe({
    cn<-filedata1()$cn  
    if (is.null(cn))
      return(NULL)
    if (is.null(filedata1))
      return(NULL)
    cn=as.character(cn)
    #updateSelectizeInput(session, "gene", choices = cn, selected=cn[1],server=TRUE)
    updateSelectizeInput(session, "gene", choices = cn, selected="NANOG",server=TRUE)
  })
  
  
  observe({
    meth_split_threshold<-filedata1()$meth_split_threshold
    if (is.null(meth_split_threshold))
      return(NULL)
    if (is.null(filedata1))
      return(NULL)
    #meth_split_threshold=as.character(cn)
    updateSelectizeInput(session, "split", choices = meth_split_threshold,selected="best")
  })
  ###  for heat map
  observe({
    if (is.null(filedata1))
      return(NULL)
    cn_new<-filedata1()$cn_new  
    if (is.null(cn_new))
      return(NULL)
    
    cn_new=as.character(cn_new)
    updateSelectizeInput(session, "Gene3",selected=cn_new[1], choices = cn_new,server=TRUE)
  })
 
  observe({
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)
    mat<-filedata1()$meth
    if (is.null(mat) | input$gene == "")
      return(NULL)
    
    #print(summary(mat))
    #mat3 = mat[mat$UCSC_RefGene_Name == input$gene, ]
    #mat3=mat[ which(mat$UCSC_RefGene_Name==input$gene),]
    #mat3=subset(mat,UCSC_RefGene_Name==as.character(input$gene))
    mat3=mat[mat$UCSC_RefGene_Name==input$gene,]
    #print(head(mat3))
    cn2=unique(mat3$Relation_to_UCSC_CpG_Island)
    cn2=as.character(cn2)
    if(!("row" %in% names(parseQueryString(session$clientData$url_search)))){
      #updateSelectizeInput(session, "island", selected =cn2[1],choices = cn2,  server = TRUE)
      updateSelectizeInput(session, "island", choices = cn2,  server = TRUE)
    }
    
    
  })
  
  
  observe({
    mat<-filedata1()$meth
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)
    if (is.null(mat) | input$gene == "") 
      return(NULL)
    
    mat2 = mat[mat$UCSC_RefGene_Name == input$gene, ]
    mat2 = mat2[mat2$Relation_to_UCSC_CpG_Island == input$island, ]
    cn1=unique(mat2$UCSC_RefGene_Group)
    cn1=as.character(cn1)
    if(!("row" %in% names(parseQueryString(session$clientData$url_search)))){
      updateSelectizeInput(session, "region", choices = cn1,  server = TRUE)
      #updateSelectizeInput(session, "region",selected =cn1[1], choices = cn1,  server = TRUE)
    }
  })
  
  
  observe({
    mat<-filedata1()$meth
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)
    if (is.null(mat)|input$gene == "")
      return(NULL)
    
    mat4 = mat[mat$UCSC_RefGene_Name == input$gene, ]
    mat4 = mat4[mat4$Relation_to_UCSC_CpG_Island == input$island, ]
    mat4 = mat4[mat4$UCSC_RefGene_Group == input$region, ]
    cn4=unique(mat4$Name)
    cn4=as.character(cn4)
    if(!("row" %in% names(parseQueryString(session$clientData$url_search)))){withProgress(message = 'Updating, please wait...', value = NULL, {
      updateSelectizeInput(session, "probe", choices = cn4, server = TRUE)
      #updateSelectizeInput(session, "probe",selected=cn4[1], choices = cn4, server = TRUE)
    })
      
      
    }
   
  })
  
  
  observe({
    if(!is.null(input$gene) && input$gene != ""){
      clinic<-filedata1()$clinic
      if (is.null(clinic))
        return(NULL)
      inFile<-filedata1()
      if (is.null(inFile))
        return(NULL)
      cov_vector=filedata1()$cov_vector
      cov_vector=sort(cov_vector,decreasing=TRUE)
      #updateCheckboxGroupInput(session, "covariate", choices = cov_vector)
      #print(cov_vector)
      updateSelectInput(session, "covariate", choices = cov_vector, selected = NULL)
      
    }
  })
  
  getPlot = reactive( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)  
    dat=filedata1()$meth
    clinic=filedata1()$clinic
    if (is.null(dat))
      return(NULL)
    if(!input$includeCovariate){
      covariate = NULL
    } else {
      covariate = input$covariate
    }
    cov_vector=filedata1()$cov_vector
    validate(
      need(input$gene != "", "Gene is empty, select a gene please")
    )
    validate(
      need(input$island != "", "Relation to island is empty, select an Island location please")
    )
    validate(
      need(input$region != "", "Genomic Region is empty, select a region please")
    )
    validate(
      need(input$probe != "", "CpG site is empty, select a CpG site please")
    )
validate(
      need(input$split != "", "Split is empty, choose a split site please")
    )
    cox_model_best_split_plot(cancer_data1=dat,cancer_data2=clinic,gene=as.character(input$gene),region=as.character(input$region),
                              group=as.character(input$island),probe=as.character(input$probe),cov_vector=input$covariate,threshold=0.25,split=input$split)
  })
  ##for gene region tab

 getPlot5 = reactive( {
    #region_tab<-getTable5()
    #if (is.null(region_tab))
      #return(NULL)
region_tab=filedata1()$meth5
#region_tab$Row.names=region_tab$Name
#colnames(region_tab)=gsub("Name","Row.names",colnames(region_tab))
  #region_tab=region_tab[,!colnames(region_tab) %in% c("CHR","MAPINFO","Name")]

    clinic=filedata1()$clinic5
    if (is.null(clinic))
      return(NULL)
    #updateSelectizeInput(session, "probe_sel",selected=cpgs[1], choices = cpgs,server=TRUE)
    cox_model_best_split_plot_v2(cancer_data1=region_tab,cancer_data2=clinic,cov_vector=NULL,threshold=0.25,probe=as.character(input$probe_sel),split="best")
  })

###fordownload
 km2_reac = reactive( {
    #region_tab<-getTable5()
    #if (is.null(region_tab))
      #return(NULL)
region_tab=filedata1()$meth5
#region_tab$Row.names=region_tab$Name
#colnames(region_tab)=gsub("Name","Row.names",colnames(region_tab))
  #region_tab=region_tab[,!colnames(region_tab) %in% c("CHR","MAPINFO","Name")]

    clinic=filedata1()$clinic5
    if (is.null(clinic))
      return(NULL)
    #updateSelectizeInput(session, "probe_sel",selected=cpgs[1], choices = cpgs,server=TRUE)
    cox_model_best_split_plot_v2(cancer_data1=region_tab,cancer_data2=clinic,cov_vector=NULL,threshold=0.25,probe=as.character(input$probe_sel),split="best")
  }) 
###
  

 
  output$downloadkm2png <- downloadHandler("survplot.png",
                                           content = function(file) {
                                             png(file)
                                             km2_reac()
                                             dev.off()
                                           },
                                           contentType = "image/png")

 
  output$downloadkm2pdf <- downloadHandler(filename ="survplot.pdf",
                                          content = function(file) {
                                            pdf(file)
                                            km2_reac()
                                            dev.off()
                                          },
                                          contentType = "application/pdf")
  
  output$downloadkmpdf <- downloadHandler(filename ="survplot.pdf",
                                          content = function(file) {
                                            pdf(file)
                                            km_reac()
                                            dev.off()
                                          },
                                          contentType = "application/pdf")
  ####reactive one for download
  km_reac = reactive( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)  
    dat=filedata1()$meth
    clinic=filedata1()$clinic
    if (is.null(dat)|input$gene == "")
      return(NULL)
    if(!input$includeCovariate){
      covariate = NULL
    } else {
      covariate = input$covariate
    }
    cov_vector=filedata1()$cov_vector
    cox_model_best_split_plot(cancer_data1=dat,cancer_data2=clinic,gene=as.character(input$gene),region=as.character(input$region),group=as.character(input$island),probe=as.character(input$probe),cov_vector=input$covariate,threshold=0.25,split=input$split)
  })
  ######
  output$plot1 <- renderPlot({
    withProgress(message = 'Processing, please wait...', value = NULL, {
      getPlot()
    })
    #dev.off()
  })	     
  ###densityplot
  getPlot3 =  reactive( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)  
    dat=filedata1()$meth
    clinic=filedata1()$clinic
    if (is.null(dat)|input$gene == "")
      return(NULL)
    
    if(!input$includeCovariate){
      covariate = NULL
    } else {
      covariate = input$covariate
    }
    cov_vector=filedata1()$cov_vector
    validate(
      need(input$gene != "", "Gene is empty, select a gene please")
    )
    validate(
      need(input$island != "", "Relation to island is empty, select an Island location please")
    )
    validate(
      need(input$region != "", "Genomic Region is empty, select a region please")
    )
    validate(
      need(input$probe != "", "CpG site is empty, select a CpG site please")
    )
validate(
      need(input$split != "", "Split is empty, choose a split site please")
    )
    distr_plot(cancer_data1=dat,cancer_data2=clinic,gene=as.character(input$gene),region=as.character(input$region),group=as.character(input$island),probe=as.character(input$probe),cov_vector=input$covariate,threshold=0.25,split=input$split) 
    
  })
  
  
  
 
  den_reac = reactive( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)  
    dat=filedata1()$meth
    clinic=filedata1()$clinic
    if (is.null(dat)|input$gene == "")
      return(NULL)
    if(!input$includeCovariate){
      covariate = NULL
    } else {
      covariate = input$covariate
    }
    cov_vector=filedata1()$cov_vector
    distr_plot(cancer_data1=dat,cancer_data2=clinic,gene=as.character(input$gene),region=as.character(input$region),group=as.character(input$island),probe=as.character(input$probe),cov_vector=input$covariate,threshold=0.25,split=input$split)
  })
  ####
  
  ###densityplot download
  
  output$downloaddenpng <- downloadHandler("density.png",
                                           content = function(file) {
                                             png(file)
                                             den_reac()
                                             dev.off()
                                           },
                                           contentType = "image/png")
  output$downloaddenpdf <- downloadHandler("density.pdf",
                                           content = function(file) {
                                             pdf(file)
                                             den_reac()
                                             #hm2 = den_reac()
                                             #grid.newpage()
                                             #grid.draw(hm2$gtable)
                                             dev.off()
                                           },
                                           contentType = "application/pdf")
  
  
  output$plot3 <- renderPlot({
    withProgress(message = 'Processing, please wait...', value = NULL, {
      getPlot3()
    })
    #dev.off()
  })
  
  ###
  ###here for violinplot####
  getPlot4 = reactive( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)  
    dat=filedata1()$meth
    if (is.null(dat)|input$gene == "")
      return(NULL)
    if(!input$includeCovariate){
      covariate = NULL
    } else {
      covariate = input$covariate
    }
    validate(
      need(input$gene != "", "Gene is empty, select a gene please")
    )
    validate(
      need(input$island != "", "Relation to island is empty, select an Island location please")
    )
    validate(
      need(input$region != "", "Genomic Region is empty, select a region please")
    )
    validate(
      need(input$probe != "", "CpG site is empty, select a CpG site please")
    )

    drawvp(meth_mat=dat, gene=as.character(input$gene),var=as.character(input$Feature), cpg= as.character(input$probe),annot=filedata1()$vio)  
    
  })
  
  output$plot4 <- renderPlot({
    withProgress(message = 'Processing, please wait...', value = NULL, {
      getPlot4()
    })
    #dev.off()
  })	
output$plot5 <- renderPlot({
    withProgress(message = 'Processing, please wait...', value = NULL, {
validate(
        need(input$probe_sel != "", "There are no CpGs, Choose a CpG please")
      )
      getPlot5()
    })
    #dev.off()
  })	
     
  #####ends here
  ###here for heatmap####
  getPlot2 = reactive( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)  
    mat=filedata1()$meth3
    #mat=inFile[[meth3]]
    #print(head(mat))
    if (is.null(mat)|input$Gene3 == "")
      return(NULL)
    
    drawhmap(mat=mat, gene=as.character(input$Gene3),annot=filedata1()$annot)  
    
  })
  
  output$plot2 <- renderPlot({
    withProgress(message = 'Processing, please wait...', value = NULL, {
      getPlot2()
    })
    #dev.off()
  },height=950,width=950)	     
  #####ends here
  ####reactive one for hm download
  hm_reac = reactive( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)  
    mat=filedata1()$meth3
    if (is.null(mat)|input$Gene3 == "")
      return(NULL)
    
    drawhmap(mat=mat, gene=as.character(input$Gene3),annot=filedata1()$annot)   
    
  })
 
  vp_reac = reactive( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)  
    dat=filedata1()$meth
    if (is.null(dat)|input$gene == "")
      return(NULL)
    if(!input$includeCovariate){
      covariate = NULL
    } else {
      covariate = input$covariate
    }
    
    drawvp(meth_mat=dat, gene=as.character(input$gene),var=as.character(input$Feature), cpg= as.character(input$probe),annot=filedata1()$vio)  
    
  })
  
  ####
  ###
  output$downloadviopng <- downloadHandler("violinplot.png",
                                           content = function(file) {
                                             png(file, width = 950, height = 950)
                                             print(vp_reac())
                                             #grid.newpage()
                                             #grid.draw(vp$gtable)
                                             dev.off()
                                           },
                                           contentType = "image/png")
  output$downloadviopdf <- downloadHandler("violinplot.pdf",
                                           content = function(file) {
                                             pdf(file,width=950/72,height=950/72)
                                            print(vp_reac())
                                             #vp = vp_reac()
                                             #grid.newpage()
                                             #grid.draw(vp$gtable)
                                             dev.off()
                                           },
                                           contentType = "application/pdf")
  ###
  
  # #download data here
  # output$downloadData <- downloadHandler(
  #     filename = function() {
  #       paste('data-', Sys.Date(), '.csv', sep='')
  #     },
  #     content = function(con) {
  #       write.csv(data, con)
  #     }
  #   )
  # ##
  
  ####render table here####
  
  getTable = reactive( {
    dat<-filedata1()$meth
    clinic=filedata1()$clinic
    cov_vector=filedata1()$cov_vector
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)
    mat=filedata1()$meth
    if (is.null(mat)|input$gene == "")
      return(NULL)
    if(!input$includeCovariate){
      covariate = NULL
    } else {
      covariate = input$covariate
    }
    validate(
      need(input$gene != "", "Gene is empty, select a gene please")
    )
    validate(
      need(input$island != "", "Relation to island is empty, select an Island location please")
    )
    validate(
      need(input$region != "", "Genomic Region is empty, select a region please")
    )
    validate(
      need(input$probe != "", "CpG site is empty, select a CpG site please")
    )
validate(
      need(input$split != "", "Split is empty, choose a split site please")
    )
    cox_model_best_split_table(cancer_data1=dat,cancer_data2=clinic,gene=as.character(input$gene),region=as.character(input$region),group=as.character(input$island),probe=as.character(input$probe),cov_vector=input$covariate,threshold=0.25,split=input$split)     
    #cox_model_best_split(cancer_data1=filedata1()$meth,cancer_data2=filedata1()$clinic,gene=input$gene,cov_vector=covariate)
    
  })
  
  
  output$estimates <- renderTable({ rownames = TRUE
  #print(row.names(getTable()))
  getTable()
  
  })
  
  
  ####table5 for chromosome based selection of survival plots
  getTable5= reactive( {
    region_res<-filedata1()$gw
    region_res=data.table(region_res)
    if (is.null(region_res))
      return(NULL)
    if(input$chr == 1){range_start<-input$range1[1];range_end<-input$range1[2]}
		if(input$chr == 2){range_start<-input$range2[1];range_end<-input$range2[2]}
		if(input$chr == 3){range_start<-input$range3[1];range_end<-input$range3[2]}
		if(input$chr == 4){range_start<-input$range4[1];range_end<-input$range4[2]}
		if(input$chr == 5){range_start<-input$range5[1];range_end<-input$range5[2]}
		if(input$chr == 6){range_start<-input$range6[1];range_end<-input$range6[2]}
		if(input$chr == 7){range_start<-input$range7[1];range_end<-input$range7[2]}
		if(input$chr == 8){range_start<-input$range8[1];range_end<-input$range8[2]}
		if(input$chr == 9){range_start<-input$range9[1];range_end<-input$range9[2]}
		if(input$chr == 10){range_start<-input$range10[1];range_end<-input$range10[2]}
		if(input$chr == 11){range_start<-input$range11[1];range_end<-input$range11[2]}
		if(input$chr == 12){range_start<-input$range12[1];range_end<-input$range12[2]}
		if(input$chr == 13){range_start<-input$range13[1];range_end<-input$range13[2]}
		if(input$chr == 14){range_start<-input$range14[1];range_end<-input$range14[2]}
		if(input$chr == 15){range_start<-input$range15[1];range_end<-input$range15[2]}
		if(input$chr == 16){range_start<-input$range16[1];range_end<-input$range16[2]}
		if(input$chr == 17){range_start<-input$range17[1];range_end<-input$range17[2]}
		if(input$chr == 18){range_start<-input$range18[1];range_end<-input$range18[2]}
		if(input$chr == 19){range_start<-input$range19[1];range_end<-input$range19[2]}
		if(input$chr == 20){range_start<-input$range20[1];range_end<-input$range20[2]}
		if(input$chr == 21){range_start<-input$range21[1];range_end<-input$range21[2]}
		if(input$chr == 22){range_start<-input$range22[1];range_end<-input$range22[2]}
		if(input$chr == "X"){range_start<-input$range23[1];range_end<-input$range23[2]}
		if(input$chr == "Y"){range_start<-input$range24[1];range_end<-input$range24[2]}

		# Subset of genome that you want to plot
region_res=region_res[region_res$CHR==input$chr,]
    		region_res<-subset(region_res,MAPINFO >= range_start & MAPINFO <= range_end )
 
region_res$P.value= signif(region_res$P.value,2)
    region_res$LR_test_pvalue= signif(region_res$LR_test_pvalue,2)
region_res=subset(region_res,select=-adj.P.value)
region_res=subset(region_res,select=c("Name","HR","CI","P.value","LR_test_pvalue","Best_split","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island"))
#region_res=as.data.frame(cbind(View = shinyInput(actionButton, nrow(region_res),'button_', label = " Click for KM Plot", onclick = 'Shiny.onInputChange(\"select_button3\",  this.id)' ),region_res))

    #canmat=as.data.frame(cbind(View = shinyInput(actionButton, nrow(canmat),'button_', label = "Click for KM Plot", onclick = 'Shiny.onInputChange(\"select_button2\",  this.id)' ),canmat))
    return(region_res)
  })
  
  
  
  output$region_res <- DT::renderDataTable({ 
    #output$allcan <-renderTable({
    withProgress(message = 'Fetching table, please wait...', value = NULL, {
      getTable5()
    })
  }, escape=FALSE,rownames= FALSE)
  #####
  ####new
  getTable2 = reactive( {
    canmat<-filedata1()$canmat
    
    if (is.null(canmat))
      return(NULL)
    
    #name = as.vector(canmat$Name)
    #hostname = "127.0.0.1:6696"
    #hostname = paste0(session$clientData$url_hostname, ":", session$clientData$url_port)
    #url = paste0("http://", hostname, "/?row=", 1:nrow(canmat), "&gene=", input$Gene)
    #canmat$Name = paste0("<a href=\"", url, "\" target=\"_blank\">", name, "</a>")
    canmat$P.value= signif(canmat$P.value,2)
    canmat$LR_test_pvalue= signif(canmat$LR_test_pvalue,2)
    canmat$PH_test_Pvalue=signif(canmat$PH_test_Pvalue,2)
    canmat=as.data.frame(cbind(View = shinyInput(actionButton, nrow(canmat),'button_', label = "Click for KM Plot", onclick = 'Shiny.onInputChange(\"select_button2\",  this.id)' ),canmat))
    return(canmat)
  })
  
  #for all cancers table
  
  output$allcan <- DT::renderDataTable({ 
    #output$allcan <-renderTable({
    withProgress(message = 'Fetching table, please wait...', value = NULL, {
      getTable2()
    })
  }, escape=FALSE,rownames= FALSE,selection='single',
  options = list(scrollX = TRUE))
  #####
  getlink = reactive( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)
    #     mat=filedata1()$meth
    #     if (is.null(mat)|input$gene == "")
    #       return(NULL)
    
    
    actionButton("link",
                 label=a( "Browse at genecards", target="_blank",href = paste0("https://genecards.weizmann.ac.il/v3/cgi-bin/carddisp.pl?gene=",input$gene )))  
  })
  
  
  output$link<-renderUI({
    getlink()
    
  })
  
  
  ###new link
  getlink2 = reactive( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)
    #     mat=filedata1()$meth
    #     if (is.null(mat)|input$gene == "")
    #       return(NULL)
    
    
    actionButton("link2",
                 label=a( "Browse at COSMIC", target="_blank",href = paste0("http://grch37-cancer-legacy.sanger.ac.uk/cosmic/gene/analysis?ln=",input$gene )))  
  })
  
  
  output$link2<-renderUI({
    getlink2()
    
  })
  ####
  
  
  ###new link
  getlink3 = reactive( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)
    #     mat=filedata1()$meth
    #     if (is.null(mat)|input$gene == "")
    #       return(NULL)
    
    
    actionButton("link3",
                 label=a( "Browse at GO", target="_blank",href = paste0("http://amigo.geneontology.org/amigo/search/annotation?q=",input$gene )))  
  })
  
  
  output$link3<-renderUI({
    getlink3()
    
  })
  ####
  cpg_reac = reactive( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)  
    mat=filedata1()$meth
    if (is.null(mat)|input$gene == "")
      return(NULL)
    
    mat = mat[mat$UCSC_RefGene_Name == input$gene, ]
    
  })
  ####for clustvis
  
  output$clustVisLinkText = renderPrint({
    d = filedata1()
    clicked = input$clustVisLink
    isolate({
      set.seed(NULL) #force random start
      sessId = paste(sample(c(letters, LETTERS), 15, replace = TRUE), collapse = "")
      if(clicked){
        f = paste0("/srv/clustvis-data/settings_", sessId, ".RData")
        mat = d$meth3[which(d$meth3$UCSC_RefGene_Name %in% input$Gene3), ]
	mat=as.data.frame(mat)
        rownames(mat)=mat$Name
        annotation_col =d$annot
        #annotation_col = subset(d$annot, select=-c(bcr_patient_barcode ,PID ,time_to_event))
        #annotation_col$age = as.numeric(as.vector(annotation_col$age))
        #quant =quantile(annotation_col$age)
        #annotation_col$age=cut(annotation_col$age,breaks=quant)
        
        annotation_row = subset(mat, select = c("UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island"))
        
        mat=mat[,!colnames(mat) %in% c("Name","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island","CHR","MAPINFO")]
        #rownames(annotation_col) = colnames(mat)
        rownames(annotation_row) = rownames(mat)
        data = list(annoCol = annotation_col, annoRow = annotation_row, mat = mat,inputSaved = list(procCentering = NULL, procScaling = "none"))
        save(data, file = f)
        link = paste0("http://biit.cs.ut.ee/clustvis_large/?e=", sessId)
        text = p("Please click here to go to ClustVis: ", a(link, href = link, target = "_blank"), style = "color:green;")
      }
    })
    text
  })
  
  
  
  output$caption1 <- renderText( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)
    
    "Survival analysis summary"
  })
  output$caption2 <- renderText( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)
    
    paste0("Heat map for ",input$Gene3 , " in ",input$cancer3)
  })
  output$caption3 <- renderText( {
    inFile<-filedata1()
    if (is.null(inFile))
      return(NULL)
    
    "Kaplan meier, density and violin plots"
  })
  
  output$table1<-renderDataTable({
    inFile<-filedata1()$meth
    #   if (is.null(inFile))
    #      return(NULL)
    head(inFile)
    #print(input$adjustment)
    #inFile[[input$adjustment]]
  })
  
  ###for clickable plots
  shinyInput <- function(FUN, len, id, ...) {inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))}
  inputs
  }
  # Here I created a reactive to save which row was clicked which can be stored for further analysis
  SelectedRow <- eventReactive(input$select_button,{
    as.numeric(strsplit(input$select_button, "_")[[1]][2])
  })
  
  
  SelectedRow2 <- eventReactive(input$select_button2,{
    as.numeric(strsplit(input$select_button2, "_")[[1]][2])
  })
 
  observeEvent(input$select_button, {
    toggleModal(session, "modalExample", "open")
  }) 
  observeEvent(input$select_button2, {
    toggleModal(session, "modalExample2", "open")
  })
  ###
  
  getTable3 = reactive( {
    top<-filedata1()$top
    if (is.null(top))
      return(NULL)
    
    #top=subset(top,Adjustment=input$adjustment)
    
    #name = as.vector(top$Name)
    #hostname = paste0(session$clientData$url_hostname, ":", session$clientData$url_port)
    #url = paste0("http://", hostname, "/?row=", 1:nrow(top), "&gene=", top$UCSC_RefGene_Name, "&island=", top$Relation_to_UCSC_CpG_Island,"&region=", top$UCSC_RefGene_Group,"&probe=", top$Name,"&cancer=", input$cancer4, "&top=1")
    #url = paste0("http://", hostname, "/?row=", 1:nrow(top), "&gene=", top$UCSC_RefGene_Name, "&top=1","&adjustment=",top$Adjustment,"&cancer=", input$cancer4)
    #top$Name = paste0("<a href=\"", url, "\" >", name, "</a>")
    #top=top[top$Adjustment=="NULL",]
    #top= subset(top, select=-Adjustment)
    top=as.data.frame(cbind(View = shinyInput(actionButton, nrow(top),'button_', label = " Click for KM Plot", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' ),top))
    return(top)
  })
  
 ####table5_reac
  getTable5_reac= reactive( {
    region_res<-filedata1()$gw
    
    if (is.null(region_res))
      return(NULL)
    if(input$chr == 1){range_start<-input$range1[1];range_end<-input$range1[2]}
		if(input$chr == 2){range_start<-input$range2[1];range_end<-input$range2[2]}
		if(input$chr == 3){range_start<-input$range3[1];range_end<-input$range3[2]}
		if(input$chr == 4){range_start<-input$range4[1];range_end<-input$range4[2]}
		if(input$chr == 5){range_start<-input$range5[1];range_end<-input$range5[2]}
		if(input$chr == 6){range_start<-input$range6[1];range_end<-input$range6[2]}
		if(input$chr == 7){range_start<-input$range7[1];range_end<-input$range7[2]}
		if(input$chr == 8){range_start<-input$range8[1];range_end<-input$range8[2]}
		if(input$chr == 9){range_start<-input$range9[1];range_end<-input$range9[2]}
		if(input$chr == 10){range_start<-input$range10[1];range_end<-input$range10[2]}
		if(input$chr == 11){range_start<-input$range11[1];range_end<-input$range11[2]}
		if(input$chr == 12){range_start<-input$range12[1];range_end<-input$range12[2]}
		if(input$chr == 13){range_start<-input$range13[1];range_end<-input$range13[2]}
		if(input$chr == 14){range_start<-input$range14[1];range_end<-input$range14[2]}
		if(input$chr == 15){range_start<-input$range15[1];range_end<-input$range15[2]}
		if(input$chr == 16){range_start<-input$range16[1];range_end<-input$range16[2]}
		if(input$chr == 17){range_start<-input$range17[1];range_end<-input$range17[2]}
		if(input$chr == 18){range_start<-input$range18[1];range_end<-input$range18[2]}
		if(input$chr == 19){range_start<-input$range19[1];range_end<-input$range19[2]}
		if(input$chr == 20){range_start<-input$range20[1];range_end<-input$range20[2]}
		if(input$chr == 21){range_start<-input$range21[1];range_end<-input$range21[2]}
		if(input$chr == 22){range_start<-input$range22[1];range_end<-input$range22[2]}
		if(input$chr == "X"){range_start<-input$range23[1];range_end<-input$range23[2]}
		if(input$chr == "Y"){range_start<-input$range24[1];range_end<-input$range24[2]}

		# Subset of genome that you want to plot
    		region_res<-subset(region_res,CHR==input$chr & MAPINFO >= range_start & MAPINFO <= range_end )
 
region_res$P.value= signif(region_res$P.value,2)
    region_res$LR_test_pvalue= signif(region_res$LR_test_pvalue,2)
region_res=subset(region_res,select=-adj.P.value)
region_res=subset(region_res,select=c("Name","HR","CI","P.value","LR_test_pvalue","Best_split","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island"))
#region_res=as.data.frame(cbind(View = shinyInput(actionButton, nrow(region_res),'button_', label = " Click for KM Plot", onclick = 'Shiny.onInputChange(\"select_button3\",  this.id)' ),region_res))

    #canmat=as.data.frame(cbind(View = shinyInput(actionButton, nrow(canmat),'button_', label = "Click for KM Plot", onclick = 'Shiny.onInputChange(\"select_button2\",  this.id)' ),canmat))
    return(region_res)
  })
  ###
  
  output$popup1 <- renderUI({
    
    sel=SelectedRow()
    top=getTable3()
    a=top[sel,]
    cancer_data1=filedata1()$meth4
    cancer_data2=filedata1()$clinic4
    gene=as.character(a$UCSC_RefGene_Name)
    group=as.character(a$Relation_to_UCSC_CpG_Island)
    probe=as.character(a$Name)
    region=as.character(a$UCSC_RefGene_Group)
    split=as.character(a$Best_split)
    bsModal("modalExample", paste0("Survival Plot for the CpG ",a$Name), "", size = "large",
            column(12,renderPlot(
              cox_model_best_split_plot(cancer_data1=cancer_data1,cancer_data2=cancer_data2,
                                        gene=gene,cov_vector=NULL,threshold=0.25
                                        ,region=region,group=group,probe=probe ,split=split))
            )
    )
  })
  #render topbiomarkers
  output$topbm <- renderDataTable({ 
    #output$allcan <-renderTable({
    withProgress(message = 'Fetching table, please wait...', value = NULL, {
      tab=getTable3()
    })
  }, escape=FALSE,rownames= FALSE)
  
  ###
  ####
  
  output$popup2 <- renderUI({
    
    sel=SelectedRow2()
    canmat=getTable2()
    a=canmat[sel,]
    cancer_data1=get(paste0(a$Cancer,"_meth"))
    cancer_data2=get(paste0(a$Cancer,"_clinic_new"))
    gene=as.character(a$UCSC_RefGene_Name)
    group=as.character(a$Relation_to_UCSC_CpG_Island)
    probe=as.character(a$Name)
    region=as.character(a$UCSC_RefGene_Group)
    split=as.character(a$Best_split)
    bsModal("modalExample2", paste0("Survival Plot for the CpG ",a$Name, " in ",a$Cancer), "", size = "large",
            column(12,renderPlot(
              cox_model_best_split_plot(cancer_data1=cancer_data1,cancer_data2=cancer_data2,
                                        gene=gene,cov_vector=NULL,threshold=0.25
                                        ,region=region,group=group,probe=probe ,split=split))
            )
    )
  })
  
  
  
  
  ####
  ####
  
  ###for downloads
  output$downloadHMpng <- downloadHandler("heatmap.png",
                                          content = function(file) {
                                            png(file, width = 950, height = 950)
                                            hm = hm_reac()
                                            grid.newpage()
                                            grid.draw(hm$gtable)
                                            dev.off()
                                          },
                                          contentType = "image/png")
  output$downloadHMpdf <- downloadHandler("heatmap.pdf",
                                          content = function(file) {
                                            pdf(file,onefile=FALSE,width=950/72,height=950/72)
                                            hm = hm_reac()
                                            grid.newpage()
                                            grid.draw(hm$gtable)
                                            dev.off()
                                          },
                                          contentType = "application/pdf")
  
  output$downloadkmpng <- downloadHandler(filename ="survplot.png",
                                          content = function(file) {
                                            png(file)
                                            km_reac()
                                            dev.off()
                                          },
                                          contentType = "image/png")
  
  
  output$downloadkmpdf <- downloadHandler(filename ="survplot.pdf",
                                          content = function(file) {
                                            pdf(file)
                                            km_reac()
                                            dev.off()
                                          },
                                          contentType = "application/pdf")
  
  output$allcanTable <- downloadHandler("all_cancers.txt",
                                        content = function(file) {
                                          write.table(filedata1()$canmat, file = file, row.names = FALSE, quote = FALSE, sep = "\t")
                                        },
                                        contentType = "text/tsv")
  
  output$downloadUserExpressionFile <- downloadHandler("Methylation data for entire gene",
                                                       content = function(file) {
                                                         
                                                         write.table(cpg_reac(),
                                                                     file = file, row.names = FALSE, quote = FALSE, sep = "\t"
                                                         )
                                                       },
                                                       contentType = "text/tsv")
  output$topbmtable <- downloadHandler("top_biomarkers.txt",
                                       content = function(file) {
                                         write.table(filedata1()$top, file = file, row.names = FALSE, quote = FALSE, sep = "\t")
                                       },
                                       contentType = "text/tsv")

output$regiontable<-downloadHandler("table_selected_region.txt",content = function(file) {write.table(getTable5_reac(), file = file, row.names = FALSE, quote = FALSE, sep = "\t")
                                       },
                                       contentType = "text/tsv")
  
  session$allowReconnect(TRUE)
  ###
  
  
    
  })
})   
