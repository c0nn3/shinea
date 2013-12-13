library(shiny);library(shinyIncubator);library(xtable);library(gplots);library(hwriter);library(reshape2);library(pheatmap);library(plyr)


shinyServer(function(input, output, session) {
   output$pheatmap<-renderUI({plotOutput("heatmap", height=paste0(input$heat_hight,"px"))})
      output$heatsliderfdrcut<-renderUI({sliderInput("heat_fdrcut", "FDR q-value cutoff:",min=0,max=1,step=0.001, value=1)})
   output$heatselectdb<-renderUI({selectInput("heat_selectdb", "Select database:", multiple=F, list('Show all'='', 'SWISS'='SWISS', 'PSP'='PSP', 'HPRD'='HPRD', 'PHOSPHO.ELM'='PHOSPHO.ELM'))})
   

   #conditional select human mouse
   output$orgName <- renderUI({
      if(input$orgtype=="Human"){
         checkboxGroupInput("database","",selected=c("PhosphoSitePlus (human)", "SWISS-Prot (human)", "HPRD (human)", "Phospho.ELM (human)") , list(
                        'PhosphoSitePlus (human)'='databases/PSP_human_phosphosite_set_edit_NAclean.gmx',
                        'SWISS-Prot (human)'='databases/SWISS_human_phosphosite_set_edit_NAclean.gmx',
                        'HPRD (human)'='databases/HPRD_human_phosphosite_set_clean_edit_NAclean.gmx',
                        'Phospho.ELM (human)'='databases/Phospho_ELM_human_phosphosite_set_clean_edit_NAclean.gmx'))} 
      else if(input$orgtype=="Mouse"){
         checkboxGroupInput("database","", selected=c("PhosphoSitePlus (mouse)", "SWISS-Prot (mouse)", "Phospho.ELM (mouse)"), list(
                        'PhosphoSitePlus (mouse)'='databases/PSP_mouse_phosphosite_set_edit_NAclean.gmx',
                        'SWISS-Prot (mouse)'='databases/SWISS_mouse_phosphosite_set_edit_NAclean.gmx',
                        'Phospho.ELM (mouse)'='databases/Phospho_ELM_mouse_phosphosite_set_clean_edit_NAclean.gmx'))}
   })
   
   #counter
   views <- as.numeric(read.table("viewFile.txt", header = FALSE)[1, 1]) + 1
   write(views, file = "viewFile.txt")
   output$count <- renderText({
      paste("Views:", views, sep = " ")
   })
   
   www<-"www/"
   
   data <- reactive({
      withProgress(session, min=1, max=15, expr={
         for (i in 1:15) {setProgress(message = 'Analysing data...', detail = 'This may take a while!', value = i) / Sys.sleep(0.01)} 
         setProgress(message = 'Analysing data...', detail = 'This may take a while!', value = 15)
         inFile <- input$files # bit unclear what about this
         if(input$exampledata==T){x3 <<- "52288"
         imput<-read.table(paste0(x3,".rnk"), header=input$header, sep=input$sep, dec=input$dec, skip=input$obsskip)
         imput<-imput[c(input$obspep,input$obsratio)]
         }
         else {x3 <<- sample(1:100000, 1)
               imput<-read.table(inFile$datapath, header=input$header, sep=input$sep, dec=input$dec, skip=input$obsskip)
               imput<-imput[c(input$obspep,input$obsratio)]}
         write.table(imput, file=paste(x3, "rnk", sep="."), quote=F, row.names=F, col.names=F, sep="\t")

         #Performing javaGSEA
         if(input$exampledata==F){  system(paste0("java -cp gsea2-2.0.13.jar -Xmx512m xtools.gsea.GseaPreranked -gmx ",  paste(input$database, collapse = ","), " -collapse false -mode Max_probe -norm meandiv -nperm ", input$obsperm, " -rnk ", paste0(x3, ".rnk"), " -scoring_scheme ", input$obsscoring, " -rpt_label my_analysis -include_only_symbols true -make_sets true -plot_top_x 9999 -rnd_seed timestamp -set_max ", input$obsmax, " -set_min ", input$obsmin, "  -zip_report false -out ", paste0(www,x3, "out"), " -gui false"))}
      
         for (i in 15:1) {setProgress(message = 'Loading data', detail = 'This may take a while!', value = i) / Sys.sleep(0.01)}
         
         system(paste("cp ",www,x3, "out/*/* ",www,x3,"out/", sep=""))
         system(paste("cp ",www, x3, "out/", "gsea_report_for_na_neg_*.xls ",www, x3, "out/gsea_neg.xls", sep=""))
         system(paste("cp ",www, x3, "out/", "gsea_report_for_na_pos_*.xls ",www, x3, "out/gsea_pos.xls", sep=""))    
         
         dodo_neg<-read.table(paste(www,x3,"out/gsea_neg.xls", sep=""), sep="\t", header=T)      
         dodo_neg<-dodo_neg[c(1,4:8)]
         dodo_neg[7]<-(1-dodo_neg[6])*dodo_neg[4]    
         names(dodo_neg)<-c("Kinase", "Size", "ES", "NES", "NOM p-val", "FDR q-val", "PHOSSEA Score")
         dodo_pos<-read.table(paste(www,x3,"out/gsea_pos.xls", sep=""), sep="\t", header=T)      
         dodo_pos<-dodo_pos[c(1,4:8)]
         dodo_pos[7]<-(1-dodo_pos[6])*dodo_pos[4]    
         names(dodo_pos)<-c("Kinase", "Size", "ES", "NES", "NOM p-val", "FDR q-val", "PHOSSEA Score")
         dodo<-"Bad file."
         dodo<-rbind(dodo_pos, dodo_neg)
         resh<-dodo
         resh<-data.frame(resh)

         # backmapping
         dio_map<-resh
         dio<-resh
         dio<-data.frame(dio, colsplit(dio$Kinase, "\\(", c("Kinase", "Db_Org")))
         dio$Kinase.1<-substr(dio$Kinase.1, 1, nchar(dio$Kinase.1)-1)
         dio$Db_Org<-(gsub("PHOSPHO.ELM", "PELM", dio$Db_Org))
         dio$Db_Org<-(gsub("SWISS-PROT", "Swiss", dio$Db_Org))
         dio$Db_Org<-(gsub("\\)", "", dio$Db_Org))
         dio<-data.frame(dio, colsplit(dio$Db_Org, "_", c("Db","Org")))
         dio$ES <-round(dio$ES,2); dio$NES <-round(dio$NES,2); dio$PHOSSEA.Score <-round(dio$PHOSSEA.Score,2)
         dio$NOM.p.val<-round(dio$NOM.p.val,4); dio$FDR.q.val<-round(dio$FDR.q.val,4);
         dio<-data.frame(dio$Kinase.1, dio$Db, dio$Size, dio$ES, dio$NES, dio$NOM.p.val, dio$FDR.q.val, dio$PHOSSEA.Score)
         # backmapping off        
         
         nu<-xtable(dio)
         nu_map<-xtable(dio_map)
         names(nu)<-c("Kinase", "Database","Size","ES","NES","NOM p-val","FDR q-val","PHOSSEA score")
         print.xtable(nu, type="html", file=paste0("www/",x3,"out/","Summary.html"), include.rownames=T, sanitize.text.function = force)
         write.table(x=nu, file=paste0("www/",x3,"out/","Summary.xls"), col.names=NA, quote=F, sep="\t", )
         write.table(x=nu_map, file=paste0("www/",x3,"out/","map_Summary.xls"), col.names=NA, quote=F, sep="\t", )
      dodo 
      })
   })

#### Heatmap
   
   heatmapfun<-function(){
      head(data(),5) #dirty start of javaGSEA # todo: differente aquire, start button
      dodo_neg<-read.table(paste(www,x3,"out/gsea_neg.xls", sep=""), sep="\t", header=T)      
      dodo_neg<-dodo_neg[c(1,4:8)]
      dodo_neg[7]<-(1-dodo_neg[6])*dodo_neg[4]    
      names(dodo_neg)<-c("Kinase", "Size", "ES", "NES", "NOM p-val", "FDR q-val", "PHOSSEA Score")
      dodo_pos<-read.table(paste(www,x3,"out/gsea_pos.xls", sep=""), sep="\t", header=T)      
      dodo_pos<-dodo_pos[c(1,4:8)]
      dodo_pos[7]<-(1-dodo_pos[6])*dodo_pos[4]    
      names(dodo_pos)<-c("Kinase", "Size", "ES", "NES", "NOM p-val", "FDR q-val", "PHOSSEA Score")
      dodo<-rbind(dodo_pos, dodo_neg)
      dodo <- dodo[grep(input$heat_filter, dodo$Kinase), ]
      dodo <- dodo[grep(input$heat_selectdb, dodo$Kinase), ]
      dodo<-subset(dodo, dodo[6] <= input$heat_fdrcut)
      dodo<-dodo[c(1,7)]
      names(dodo)<-c("name", "value")
      row.names(dodo)<-dodo$name
      dodo[1]<-dodo[2]
      dodo[2]<-dodo[2]
      dodo<-dodo[ order(-dodo[,1]), ]
      #     upper<-max(input$heat_bigger)
      #     lower<-min(input$heat_smaller)
      #     dodo <- subset(dodo, value >= upper | value < lower)
      row.names(dodo)<-gsub(pattern="_", replacement=" ", x=row.names(dodo))
      row.names(dodo)<-gsub(' HUMAN', '', x=row.names(dodo)); row.names(dodo)<-gsub(' MOUSE', '', x=row.names(dodo))
      
      output$heatfontsize<-renderUI({sliderInput("heat_rowsize","Lable font size:", min=1,max=100, value=((input$heat_hight/nrow(dodo[1]))/2))}) 
      
      #      if(input$heat_uniqnames==T){dodo$names<-row.names(dodo)
      #                           dodo$names<-gsub(pattern="_", replacement=" ", x=dodo$names)
      #                           dodo$names<-gsub('(*) .*', '\\1', dodo$names)
      #                           dodo<-dodo[!duplicated(dodo$names), ]
      #                           dodo$names<-NULL
      #                           }
      if(input$heat_showdb==F){row.names(dodo)<-gsub(".PSP.","",row.names(dodo));row.names(dodo)<-gsub(".SWISS-PROT."," ",row.names(dodo));row.names(dodo)<-gsub(".HPRD.","  ",row.names(dodo));row.names(dodo)<-gsub(".PHOSPHO.ELM.","   ",row.names(dodo))} # adjust for other databases later! uniqueness is gained by adding different numbers of " " to each database
      
      pheatmap(data.matrix(dodo), border_color = NA, clegend=T,show_colnames=F, display_numbers=F, show_rownames=T, fontsize=15, fontsize_row=input$heat_rowsize,  cellwidth=input$heat_cellwidth, cluster_rows=F, color=bluered(255), breaks=unique(c(seq(min(dodo[1]), -0.5, length.out=(255/2)),seq(0.5, max(dodo[1]), length.out=255/2))), main=input$heat_main)
      
   }
   output$heatmap <- renderPlot(function() {
      if (is.null(input$exampledata==F)) {return(NULL)}
      heatmapfun()
      })
   
   output$downloadData_heat <- downloadHandler(
      filename = function() { paste(x3, '.csv', sep=' ') },
      content = function(file) {
         png(file)
         heatmapfun()
         dev.off()}, contentType = 'image/png'
      )
   
   
#### HTML report
  #reactiveUI for html
  output$previewhtml <- reactiveUI(function() {
     if (is.null(input$exampledata==F)) {return(NULL)}
     data() # dirty starting analysis # todo: start with Button
      someout<-paste(www,x3,"out/AKT1.html", sep="")
      HTML(readLines(someout))
  })
  
  #### Input Details
  output$filetable2 <- renderTable(function() {
     if (is.null(input$exampledata==F)) {return(NULL)}
    input$files
  })
    
####  File preview detail on path etc
  output$filetable3 <- renderTable(function() {
     if (is.null(input$exampledata==F)) {return(NULL)}
    inFile <- input$files
    dio<-read.table(inFile$datapath, header=input$header, sep=input$sep, dec=input$dec, skip=input$obsskip)
    head(dio[c(input$obspep,input$obsratio)], 5)
  })
   
   ####  File raw preview
   output$RAWpreview <- renderPrint(function() {
      if (is.null(input$exampledata==F)) {return(NULL)}
      inFile <- input$files
      if (input$exampledata==T) {inFile<-paste0(x3,".rnk")} 
      diom<-read.table(inFile, sep="\t")
      names(diom)<-gsub("V","", names(diom))
      head(diom,5)
 #     diom
   })   
   
   ####  File preview
   output$preview <- renderPrint(function() {
      if (is.null(input$exampledata==F)) {return(NULL)}
      inFile <- input$files
      if (input$exampledata==T) {inFile<-paste0(x3,".rnk")} 
      dio<-read.table(inFile, header=input$header, sep=input$sep, dec=input$dec, skip=input$obsskip)
      names(dio)<-gsub("V","", names(dio))
      head(dio[c(input$obspep,input$obsratio)],5)
#      dio[c(input$obspep,input$obsratio)]                                      
   })   
   
  output$downloadData <- downloadHandler(
      filename = function() { paste(x3, '.csv', sep=' ') },
       content = function(file) {
        write.table(data(), file, quote=F, sep="\t", col.names=NA) 
  })
   
   output$downloadData_detail <- downloadHandler(
      filename = function() { paste(x3, '.csv', sep='') },
      content = function(file) {
         if (is.null(input$exampledata==F)) {return(NULL)}
         data()

         filenames<-list.files(paste0(www,x3,"out/"), pattern="*.xls")
         filenames<-dirselect[grep(c("Summary|index|gsea|snapshot|gene"), filenames, invert=T)]
         filenames<-na.omit(filenames)
         
         mydata = ldply(filenames, function(filename) {
            detkin = read.table(paste0(www,x3,"out/",filename), sep="\t", header=T, check.names=T)[,2:8]
            dbfrom<-gsub("", "", filename)
            dbfrom<-colsplit(dbfrom, "\\(", c("col1","col2"))
            dbfrom<-dbfrom$col2
            dbfrom<-substr(dbfrom, 1, nchar(dbfrom)-1)
            dbfrom<-colsplit(dbfrom, "_", c("db","org"))
            dbfrom<-dbfrom$db
            
            querydb<-read.table(paste0("backmapping/",dbfrom,"_human.txt"), sep="\t", header=T, check.names=F, fill=T)   
            names(querydb)<-c("PROBE", "Protname", "ResPos", "AccNo", "symbol")
            total<-merge(querydb,detkin, by="PROBE", all.y=T, all.x=T, )
            total<-total[order(total$RANK.IN.GENE.LIST),]
            finalout<-data.frame(total$PROBE, total$Protname, total$ResPos, total$AccNo, total$symbol, round(total$RANK.METRIC.SCORE,4), total$RANK.IN.GENE.LIST, round(total$RUNNING.ES,4), total$CORE.ENRICHMENT)
            finalout<-subset(finalout, finalout[6] != "")
            names(finalout)<-c("Phosphosite Sequence", "Protein Name", "Residue + Position", "Accession No.", "Gene Symbol", "Experimental value", "Rank in Phosphoprofile", "Running ES", "Core Enrichment")
            finalout$filename = gsub(".xls","",filename)
            return(finalout)
         })
         write.table(mydata, file, quote=F, sep="\t", col.names=NA)
      })
   
  
   
   # Check boxes for resultselect
   output$choose_columns <- reactiveUI(function(){
      # If missing input, return to avoid error later in function
      if (is.null(input$exampledata==F)) {return(NULL)}
      
     # dirselect<-NULL
      data()
      dirselect<-list.files(paste0(www,x3,"out/"), pattern="*.xls")
      dirselect<-dirselect[grep(c("Summary|index|gsea|snapshot|gene"), dirselect, invert=T)]
      dirselect=gsub(".xls","",dirselect)
      
      dada<-read.table(paste0("www/",x3,"out/","map_Summary.xls"), sep="\t", header=T, check.names=F, row.names=1)
      fileident<-dada$Kinase
      dada$Kinase<-gsub("\\(", "halligalli ", dada$Kinase)
      dada$Kinase<-gsub('HUMAN|\\)', "", dada$Kinase) 
      dada$Kinase<-gsub("_"," ", dada$Kinase)
      dada$Kinase<-gsub("halligalli"," | ", dada$Kinase)
      dada$Kinase<-(gsub("PHOSPHO.ELM", "PELM", dada$Kinase))
      dada$Kinase<-(gsub("SWISS-PROT", "Swiss", dada$Kinase))
      dtailedselect<-as.character(paste(dada$Kinase, " | ", dada$Size,  " | ",round(dada$ES,2), " | ",round(dada$NES,2), " | ",round(dada$NOM.p.val,4), " | ",round(dada$FDR.q.val,4), " | ", round(dada$PHOSSEA.Score,2), sep=""))
      dirselect<-as.list(setNames(as.character(fileident), dtailedselect))      
      select2Input("tool1s", "", choices = dirselect, multiple=F)
      })   
   
   
   
   # Output the data
   output$dynpreview <- reactiveUI(function(){
      # If missing input, return to avoid error later in function
      if (is.null(input$exampledata==F)) {return(NULL)}
      
      data()
      HTML(readLines(paste0(www,x3,"out/",input$tool1s,".html")))
      
   })
   
   output$dynpreviewxls_summary <- renderDataTable(function(){ # change to renderTable for std
      if (is.null(input$exampledata==F)) {return(NULL)}
      detkin<-read.table(paste0("www/",x3,"out/Summary.xls"), sep="\t", header=T, check.names=F, row.names=1)
      detkin
   }, options =list(aLengthMenu = c(10, 25, 50), iDisplayLength=10)
   )
   
   
   dynoutputfun<-function(){
      data()
      # get propper database name from the kinase in detailed results
      dbfrom<-input$tool1s
      dbfrom<-colsplit(dbfrom, "\\(", c("col1","col2"))
      dbfrom<-dbfrom$col2
      dbfrom<-substr(dbfrom, 1, nchar(dbfrom)-1)
      dbfrom<-colsplit(dbfrom, "_", c("db","org"))
      dbfrom<-dbfrom$db      
      
      #backmap details
      detkin<-read.table(paste0("www/",x3,"out/",input$tool1s,".xls"), sep="\t", header=T, check.names=T)[,2:8]
      querydb<-read.table(paste0("backmapping/",dbfrom,"_human.txt"), sep="\t", header=T, check.names=F, fill=T)
      names(querydb)<-c("PROBE", "Protname", "ResPos", "AccNo", "symbol")
      total<-merge(querydb,detkin, by="PROBE", all.y=T, all.x=T, )
      total<-total[order(total$RANK.IN.GENE.LIST),]
      
      if(dbfrom=="PSP"){
         total$Protname<-gsub("(.*)", '<a href=\"http://www.phosphosite.org/simpleSearchSubmitAction.do?queryId=-1&from=0&searchStr=\\1\" target=\"_blank\">\\1</a>', total$Protname)
         total$AccNo<-gsub("(.*)", '<a href=\"http://www.uniprot.org/uniprot/\\1#section_features\" target=\"_blank\">\\1</a>', total$AccNo)   
      }
      if(dbfrom=="HPRD"){
         total$symbol<-total$Protname
         total$Protname<-gsub("(.*)", '<a href=\"http://www.hprd.org/resultsQuery?multiplefound=&prot_name=\\1&external=Ref_seq&accession_id=&hprd=&gene_symbol=&chromo_locus=&function=&ptm_type=&localization=&domain=&motif=&expression=&prot_start=&prot_end=&limit=0&mole_start=&mole_end=&disease=&query_submit=Search\" target=\"_blank\">\\1</a>', total$Protname)
         total$AccNo<-gsub("(.*)", '<a href=\"http://www.ncbi.nlm.nih.gov/protein/\\1\" target=\"_blank\">\\1</a>', total$AccNo)
      }
      if(dbfrom=="PHOSPHO.ELM"){
         total$Protname<-total$symbol
         total$Protname<-gsub("(.*)", '<a href=\"http://phospho.elm.eu.org/bySubstrate/\\1.html\" target=\"_blank\">\\1</a>', total$Protname)
         total$AccNo<-gsub("(.*)", '<a href=\"http://www.uniprot.org/uniprot/\\1#section_features\" target=\"_blank\">\\1</a>', total$AccNo)   
      }
      if(dbfrom=="SWISS-PROT"){
         total$Protname<-total$symbol
         total$Protname<-gsub("(.*)", '<a href=\"http://www.uniprot.org/uniprot/?query=\\1&sort=score\" target=\"_blank\">\\1</a>', total$Protname)
         total$AccNo<-gsub("(.*)", '<a href=\"http://www.uniprot.org/uniprot/\\1#section_features\" target=\"_blank\">\\1</a>', total$AccNo)   
      }
      
      finalout<-data.frame(total$PROBE, total$Protname, total$ResPos, total$AccNo, total$symbol, round(total$RANK.METRIC.SCORE,4), total$RANK.IN.GENE.LIST, round(total$RUNNING.ES,4), total$CORE.ENRICHMENT)
      finalout<-subset(finalout, finalout[6] != "")
      names(finalout)<-c("Phosphosite Sequence", "Protein Name", "Residue + Position", "Accession No.", "Gene Symbol", "Experimental value", "Rank in Phosphoprofile", "Running ES", "Core Enrichment")
      #backmapping end
      #finalout<-finalout[!duplicated(finalout[c("Rank in Phosphoprofile")]),] #duplica filtering gives hint on redundancy problem??
      finalout
   }
   
   # Output the data only bottom table of detailed report  # initial try to inlcude xls
   output$dynpreviewxls <- renderDataTable(function(){ # change to renderTable for std
      if (is.null(input$exampledata==F)) {return(NULL)}
      
      dynoutputfun()
      
   }, options =list(aLengthMenu = c(10, 25, 50))
   )
   
   output$downloadData_single <- downloadHandler(
      filename = function() { paste(x3, '.csv', sep=' ') },
      content = function(file) {
         write.table(dynoutputfun(), file, quote=F, sep="\t", col.names=NA) 
      })
   
   
   
   ##### Results Summary page well1
   output$summaryresults1 <- renderText(function() {
      if (is.null(input$exampledata==F)) {return(NULL)}
      head(data(),5) #dirty start of javaGSEA # todo: differente aquire, start button
      htmlToText <- source("htmlToText.R")
      doc<-htmlToText(paste(www,x3,"out/index.html", sep=""))
      doc.df<-as.data.frame(doc)
      names(doc.df)<-"duuu"
      #substr(doc.df[64,1], start=1, stop=13)
      #nrofuniqgenes<-gsub(' features .genes.|The dataset has ', "", unique(grep("The dataset has", doc.df$duuu, value=T)))
      indexex1<-gsub("There were duplicate row identifiers in the specified ranked list. One id was arbitarilly choosen. Details are below. Generally, this is OK, but if you want to avoid this automagic, edit your ranked list so that all row ids are unique# of row ids in original dataset: |# of row UNIQUE ids in original dataset:|# The duplicates were", "", unique(grep("There were duplicate row identifiers", doc.df$duuu, value=T)))
      agid<-unlist(strsplit(indexex1," "))[1] #all genes in dataset
      ugid<-unlist(strsplit(indexex1," "))[2] #uniq genes in dataset
      nodg<-as.numeric(unlist(strsplit(indexex1," ")))[1]-as.numeric(unlist(strsplit(indexex1," "))[2]) # number of duplicated genes
      indexex2<-unlist(strsplit(unique(grep("Gene set size filters", doc.df$duuu, value=T)), " "))
      nops<-indexex2[13] # No. of phosphosite sets:
      nopswms<-as.numeric(indexex2[13])-as.numeric(indexex2[11]) # No. of phosphosite sets with min. size:
      
      dbshow<-gsub("databases/","", gsub(", ","",toString(input$database)))
      dbshow<-gsub("_human_phosphosite_set_edit_NAclean.gmx", " ", dbshow)
      dbshow<-gsub("_human_phosphosite_set_clean_edit_NAclean.gmx", " ", dbshow)
      dbshow<-gsub("_mouse_phosphosite_set_edit_NAclean.gmx", " ", dbshow)
      dbshow<-gsub("_mouse_phosphosite_set_clean_edit_NAclean.gmx", " ", dbshow)
      dbshow<-substr(dbshow, 1, nchar(dbshow)-1)
      dbshow<-gsub(" ", ", ", dbshow)
      
      if(input$exampledata==T){example.filename<-"Example data set"}
      else {example.filename <- input$files$name}   
        
         
     
      filename<-data.frame(
         c("Used Dataset:",
           "Name of uploaded file:", 
           "No. of phosphopeptides:", 
           "No. of unique phosphopeptides:", 
           "No. of duplicate phosphopeptides:", 
           "Used phosphosite set database:", 
           "Name of used database:",
           "Organism:",
           "No. of phosphosite sets:", 
           "No. of phosphosite sets with min. size:",
           "Used parameter settings:",
           "No of permutations:",
           "Min. size of phosphosite sets:"
         ),
         c(" ", 
           example.filename, 
           agid, 
           ugid, 
           nodg, 
           " ",
           dbshow,
           input$orgtype,
           nops, 
           nopswms,
           " ",
           input$obsperm,
           input$obsmin
         ))
   
      names(filename)<-c("ID", "Value")
      #filename
      hwrite(filename,
             style=matrix(c('font-weight:bold', NA, NA, NA,NA, 'font-weight:bold', NA, NA, NA, NA,  'font-weight:bold', NA, NA),nr=13,nc=2),
             border=0, col.names=FALSE)
   })
   
   ##### Results Summary page well2
   output$summaryresults2 <- renderText(function() {
      if (is.null(input$exampledata==F)) {return(NULL)}
      head(data(),5) #dirty start of javaGSEA # todo: differente aquire, start button
      htmlToText <- source("htmlToText.R")
      doc<-htmlToText(paste(www,x3,"out/index.html", sep=""))
      doc.df<-as.data.frame(doc)
      names(doc.df)<-"duuu"
      
      indexex3<-unique(grep("gene sets are upregulated in phenotype", doc.df$duuu, value=T))
      indexex3pos<-indexex3[1]
      indexex3pos<-unlist(strsplit(indexex3pos, " "))
      indexex3pos<-paste0(indexex3pos[1], " (of ", indexex3pos[3], ")")
      indexex3neg<-indexex3[2]
      indexex3neg<-unlist(strsplit(indexex3neg, " "))
      indexex3neg<-paste0(indexex3neg[1], " (of ", indexex3neg[3], ")")
      indexex4<-unique(grep("at FDR", doc.df$duuu, value=T))
      indexex4pos<-indexex4[1]
      indexex4pos<-unlist(strsplit(indexex4pos, " "))
      indexex4pos1<-indexex4pos[1]
      indexex4pos2<-indexex4pos[8]
      indexex4pos3<-indexex4pos[17]
      indexex4neg<-indexex4[2]
      indexex4neg<-unlist(strsplit(indexex4neg, " "))
      indexex4neg1<-indexex4neg[1]
      indexex4neg2<-indexex4neg[9]
      indexex4neg3<-indexex4neg[18]
      
      dbshow<-gsub("databases/","", gsub(", ","",toString(input$database)))
      dbshow<-gsub("_human_phosphosite_set_edit_NAclean.gmx", " ", dbshow)
      dbshow<-gsub("_human_phosphosite_set_clean_edit_NAclean.gmx", " ", dbshow)
      dbshow<-gsub("_mouse_phosphosite_set_edit_NAclean.gmx", " ", dbshow)
      dbshow<-gsub("_mouse_phosphosite_set_clean_edit_NAclean.gmx", " ", dbshow)
      dbshow<-substr(dbshow, 1, nchar(dbshow)-1)
      dbshow<-gsub(" ", ", ", dbshow)
      
      filename<-data.frame(
         c("PHOSSEA results:",
           "Positively enriched phosphosite sets:",
           "Total no.:",
           "Significant at FDR ≤ 25%:",
           "Significant at nom. P value ≤ 5%:",
           "Significant at nom. P value ≤ 1%:",
           "Negatively enriched phosphosite sets:",
           "Total no.:",
           "Significant at FDR ≤ 25%:",
           "Significant at nom. P value ≤ 5%:",
           "Significant at nom. P value ≤ 1%:"
         ),
         c(" ",
           " ",
           indexex3pos,
           indexex4pos1,
           indexex4pos2,
           indexex4pos3,
           " ",
           indexex3neg,
           indexex4neg1,
           indexex4neg2,
           indexex4neg3
         ))
      names(filename)<-c("ID", "Value")
      #filename
      hwrite(filename,
             style=matrix(c('font-weight:bold', 'font-weight:bold', NA,NA,NA,NA, 'font-weight:bold', NA,NA,NA,NA),nr=11,nc=2),
             border=0, col.names=FALSE)
   })
   
   
   
   
   
})
