#setwd("/Users/cfischer/Groundwork/PHOSSEA")
#runApp("/Users/cfischer/Groundwork/PHOSSEA") # to start locally
library(shiny)
library(shinyIncubator)
library(hwriter)
library(reshape2)
library(plyr)

#helptext
uploadhelp<-'Currently, choices are restricted to categorical variables for classification until I find time to generalize the app to include regression.'
settingshelp<-'keep calm its only settings'
databasehelp<-'keep calm this is different tabasases no more'

wellPanel2 <- function(...) {div(class="well", style ="background-color: #F2F5FF;min-height:200px;", ...)}
wellPanelheat <- function(...) {div(class="well", style ="background-color: #FFFFFF;min-height:200px;", ...)}
wellPanelsum <- function(...) {div(class="well", style ="background-color: #F2F5FF;min-height:300px;", ...)}

helpPopup <- function(title, content,
                      placement=c('left', 'right', 'top', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual')) {
   tagList(
      singleton(
         tags$head(
            tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })")
         )
      ),
      tags$a(
         href = "#", `data-toggle` = "popover", #class = "btn btn-mini",
         title = title, `data-content` = content, `data-animation` = T,
         `data-placement` = match.arg(placement, several.ok=T)[1],
         `data-trigger` = match.arg(trigger, several.ok=F)[1], # for click on/off
         tags$i(class="icon-question-sign")
      )
   )
}
shinyUI(basicPage(
   headerPanel(windowTitle="PHOSSEA  – Phosphosites Set Enrichment Analysis",
      HTML(
         '<div id="stats_header">
         PHOSSEA  <font size="5"> - <font color="#336666">P</font>hospho<font color="#336666">s</font>ite <font color="#336666">S</font>et <font color="#336666">E</font>nrichment <font color="#336666">A</font>nalysis</font>
         <a href="http://www.molgen.mpg.de" target="_blank">
         <img id="stats_logo" align="left" alt="Logo" src="site_header.jpg" hspace="7" vspace="2" width="45px"/></a>
          </div>'
         )
   ),
   mainPanel(
      progressInit(),
      div(class="row", div(class="span12",
         conditionalPanel("updateBusy() || $('html').hasClass('shiny-busy')",
            id='progressIndicator',
            "Processing data...",
            div(id='progress',includeHTML("timer.js"))
         ),
         tags$head(tags$style(type="text/css",
            '#progressIndicator {',
            '  position: fixed; top: 8px; right: 8px; width: 200px; height: 50px;',
            '  padding: 8px; border: 1px solid #CCC; border-radius: 8px;',
            '}'
         )),
         tabsetPanel(
            tabPanel("Upload and Settings",   
               div(class="row",
                  div(class="span4",
                     wellPanel(
                        div(class="row-fluid",
                           div(class="span11",           
                              h4("Data Upload"),
                              tags$div(title="Upload a dataset.", #Tooltip
                                 fileInput("files", "", multiple=FALSE)
                              ),
                              checkboxInput("exampledata", "Load example data.", T)
                           ),
                           div(class="span1",helpPopup("Uploadhelp",uploadhelp))
                        )
                     )
                  ),
                  div(class="span4", 
                     wellPanel(
                        tags$style(type='text/css', ".well {background-color: #F2F5FF; min-height:190px;}") ,
                        div(class="row-fluid",
                           div(class="span11",
                              h4("Settings"),
                               tags$style(type="text/css", '#obsperm {width: 60px;}'),
                              numericInput("obsperm", "No of permutations:", 1000, min=1,max=10000,step=1),
                               tags$style(type="text/css", '#obsmin {width: 60px;}'),
                              numericInput("obsmin", "Min. sites - exclude smaller sets:", 3, min=1,max=1000,step=1)
                           ),    
                           div(class="span1",helpPopup("Settingshelp",settingshelp))
                        )
                     )
                  ),
                  div(class="span4",    
                     wellPanel(
                        div(class="row-fluid",
                           div(class="span11", 
                              h4("Database selection"),
                               tags$style(type="text/css", '#orgtype {width: 100px;}'),
                               selectInput("orgtype","",list("Human","Mouse"),selected="Human"),
                               uiOutput("orgName")
                           ),
                           div(class="span1",helpPopup("databasehelp",databasehelp))
                        )      
                     )
                  )
               )
            ),
   
            tabPanel("Input data inspection",
                     tabPanel("Shape your data for PHOSSEA analysis",
                        div(class="row",
                            div(class="span4",
                                wellPanel2(
                                   h4("Data preview"),      
                                   verbatimTextOutput("preview")
                                )
                            ),     
                           div(class="span4",
                               wellPanel2(
                                  h4("Columns selection"),
                                  tags$style(type="text/css", '#obspep {width: 30px;}'),
                                  numericInput("obspep", "Peptide column:", 1), # to select columns from big data table
                                  tags$style(type="text/css", '#obsratio {width: 30px;}'),
                                  numericInput("obsratio", "Ratio column:", 2)
                               )
                           ),
                           div(class="span4",
                               wellPanel2(
                                  h4("Additional options"),
                                  checkboxInput('header', '1st row header', FALSE),
                                  tags$style(type="text/css", '#sep {width: 80px;}'),
                                  selectInput('sep', 'Separator:', c(Tab='\t',Comma=',', Semicolon=';'), 'Tab'),
                                  tags$style(type="text/css", '#obsskip {width: 80px;}'),
                                  numericInput("obsskip", "Skip first lines:", 0)
                               )
                           )
                        )      
                     ),
                  tabsetPanel(
                     tabPanel("The file you uploaded",
                        div(class="row",      
                           div(class="span12",
                              verbatimTextOutput("RAWpreview")
                           )
                        )
                     )
                  )   
            ),
            tabPanel("Start analysis",
               tabsetPanel(
                  
                  tabPanel("Summary",
                     div(class="row",
                        div(class="span5",      
                           wellPanelsum(
                              htmlOutput("summaryresults2")
                           )
                        ),
                        div(class="span7",
                           wellPanelsum(
                              htmlOutput("summaryresults1")
                           )
                        )
                     )
                  ),
                  tabPanel("Result tables",
                     div(class="row",
                        div(class="span12",
                           wellPanel(
                             h4("Overview"),
                             dataTableOutput("dynpreviewxls_summary"),
                             downloadButton('downloadData', 'Download')
                             
                           ),
                           wellPanel(
                             h4("Details"),
                    #         tags$style(type="text/css", '#choose_columns {width: 80px;}'), #size must be set in globalR-->select2:R
                             uiOutput("choose_columns"),
                             br(),
                             dataTableOutput("dynpreviewxls"),
                             downloadButton('downloadData_single', 'Download selected data'),
                             downloadButton('downloadData_detail', 'Download complete data set')                             
                           )
                        )
                     )
                  ),
                  tabPanel("Heatmap visualisation",
                        tags$style(type="text/css", # TO HIDE ERROR MESSAGES
                       ".shiny-output-error { visibility: hidden; }",
                       ".shiny-output-error:before { visibility: hidden; }"
                        ),       
                     div(class="row",
                        div(class="span4",
                           wellPanel(
                              h4("Visualisation"),
                              textInput("heat_main","Heatmap title:", value=""),
                              sliderInput("heat_cellwidth", "Heatmap width:",min = 1, max = 100, value =50),
                              sliderInput("heat_hight", "Heatmap hight:",min = 100, max = 1200, value =520, step=5),
                              uiOutput("heatfontsize"),
                              checkboxInput("heat_showdb", "Show database in row names?", value=T),
                              h4("Filtering"),
                              uiOutput("heatselectdb"),
                              uiOutput("heatsliderfdrcut"),
                              textInput("heat_filter", "Select kinases based on text string:", value="")
                         #    checkboxInput("heat_uniqnames", "Show only unique kinase names?", value=F)
                           ) 
                        ),
                        div(class="span8",
                           wellPanelheat(
                              uiOutput("pheatmap"),
                              downloadButton('downloadData_heat', 'Download')
                              
                           )
                        )
                     )
                  )
               )
            ),
            tabPanel("About",
               wellPanel(
                           
               h4("PHOSSEA tutorial"),
               p(style="text-align:justify",'For detailed description of PHOSSEA, please see the ',a("tuturial", href="http://www.rstudio.com/shiny/", target="_blank"),'.'),
               h4("Contact"),      
               p(style="text-align:justify",'If you have any comments, questions or suggestions, or if you are interested in collaborative projects, do not hesitate to contact',a("sauer@molgen.mpg.de", href="mailto:sauer@molgen.mpg.de"),'and to visit our research website at', a("www.molgen.mpg.de/nutrigenomics", href="http://www.molgen.mpg.de/nutrigenomics", target="_blank"),'.'),
               h4("Citation"),
               p(style="text-align:justify",'To cite the use of PHOSSEA, please reference Weidner et al. XXX.'),
               h4("Licensing"),      
               p(style="text-align:justify",'The use of this software is free to academic users. Commercial entities should contact ',a("sauer@molgen.mpg.de", href="mailto:sauer@molgen.mpg.de"),'for further information.'),
               h4("Disclaimer"),
               p(style="text-align:justify",'This software is supplied without any warranty. Neither the authors nor the Max Planck Society are responsible for its use, misuse, or functionality.'),
               br()
               )
#               div(class="row-fluid",        
#                  div(class="span4",
#                     wellPanel(
#                        HTML('<div style="clear: left;"><img src="http://www.molgen.mpg.de/68974/profile_image.jpg" width="60" alt="" style="float: left; margin-right:5px" /></div>'),
#                        strong('Dr. Sascha Sauer'),
#                        p('Group leader',br(),
#                           a('more', href="http://www.molgen.mpg.de/25022/Head_Contact", target="_blank"),'|',
#                           a("mail", href="mailto:sauer@molgen.mpg.de")
#                        )
#                     )
#                  ),
#                  div(class="span4",   
#                     wellPanel(
#                        HTML('<div style="clear: left;"><img src="http://www.molgen.mpg.de/280912/profile_image.jpg" width="60" alt="" style="float: left; margin-right:5px" /></div>'),
#                        strong('Dr. Christopher Weidner'),
#                        p('Biologist',br(),
#                           a('more', href="http://www.molgen.mpg.de/nutrigenomics/weidner", target="_blank"),'|',
#                           a("mail", href="mailto:weidner@molgen.mpg.de")
#                        )
#                     )      
#                  ),
#                  div(class="span4",   
#                     wellPanel(
#                        HTML('<div style="clear: left;"><img src="http://www.molgen.mpg.de/259767/Cornelius-Fischer.jpg" width="60" alt="" style="float: left; margin-right:5px" /></div>'),
#                        strong('Cornelius Fischer M.Sc.'),
#                        p('Computational Biologist', br(),
#                           a('more', href="http://www.molgen.mpg.de/228527/more_Fischer", target="_blank"),'|',
#                           a("mail", href="mailto:cfischer@molgen.mpg.de")
#                        )
#                     )      
#                  )                 
#              )
            )   
      ))
   ))
   ))


#I used the following to put two numericInput elements side-by-side inside a sidebar panel
#div(class="row-fluid",
#    div(class="span6",numericInput("viz_plot_height", label = "Plot height:", min = 100, step = 50, value = 650)),  
#    div(class="span6", numericInput("viz_plot_width", label = "Plot width:", min = 100, step = 50, value = 650))
#)
