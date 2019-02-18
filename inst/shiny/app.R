
library("shiny")
library("ggplot2")
library("scales")
library("shinythemes")
library("visNetwork")
###################################################
###############network summary#####################
load("data/Pavirnet_basicdata.Rda")
targets_summary <- netfordownload
colnames(targets_summary) <- c("TF ID","TF name","Target ID","Target name")
network_P5 <- list(network_summary_P5,targets_summary,vis.nodes,vis.links)
names(network_P5) <- c("network_summary","targets_summary","vis.nodes","vis.links")

GO_hit_set <- read.delim("data/mra/Hit_sets.txt",header=F,sep="\t")
whole_net <- read.delim("data/mra/PAVIRnet.txt",header=T,sep="\t")
  ##############ui#########################
  ui <- fixedPage(
  h3(strong("PAVIRnet:"),span(strong("P"),style = "word-spacing:-8px"),"seudomonas",
     span(strong("a"),style = "word-spacing:-8px"),"eruginosa",
     span(strong("V"),style = "word-spacing:-8px"),"irulence-related",
     span(strong("I"),style = "word-spacing:-8px"),"ntegrated",
     span(strong("R"),style = "word-spacing:-8px"),"egulatory",
     span(strong("net"),style = "word-spacing:-8px"), "work",
     style = "font-family: 'Source Sans Pro';
     color: #123456; text-align: center;
     background-image: url('background.png');
     padding: 20px"
  ),
  navbarPage("",
             tabPanel(h4("Introduction"),
                      tags$img(src = "homepage.png"),
                      h5("Copyright (c) 2018, Hao Huang of Dr. Xin Wang's Lab all rights reserved.",style = "text-align: center")
             ),
             ####################Network page###########################
             tabPanel(h4("PAVIRnet"),
                      fluidRow(
                        sidebarLayout(
                          sidebarPanel(
                            tags$h3("PAVIRnet"),
                            #network paramaters#
                            # checkboxInput("network_showall", "Display all Transcription Factors", TRUE),
                            radioButtons("network_show", "Subnetwork filtering",
                                         c("Full network" = "netshowall",
                                           "Sub-network" = "netshowpart"
                                         ),
                                         selected = "netshowall"
                            ),
                            selectizeInput(
                              "network_showpart", h5("Select transcription factors of interest to obtain a filtered subnetwork"), choices = tf[,1], multiple = TRUE
                            ),

                            tags$h4(""),

                            br(),
                            actionButton("update", "Update"),width = 2

                          ),
                          mainPanel(
                            tabsetPanel(
                              type = "tabs",
                              tabPanel("Visulazation",
                                       wellPanel(
                                         visNetworkOutput("network_plot", width = "100%", height = "800px")

                                       ),
                                       tags$img(src = "network_legend.png")

                              ),
                              tabPanel("Summary of TFs",
                                       DT::dataTableOutput("network_summary")
                              ),
                              tabPanel("Summary of network",
                                       DT::dataTableOutput("target_summary")
                              )
                            ),width = 10,height = "400px"
                          )
                        ))),
             tabPanel(h4("Master Regulator Analysis"),
                      sidebarLayout(
                        div(style = "wigth:100px;",
                            # Sidebar panel for inputs ----
                            sidebarPanel(

                              # Input: Select a file ----
                              radioButtons("Interest_set", "1. Choose PAVIRnet or upload your own network",
                                           c("PAVIRnet" = "online_interest_set",
                                             "Upload your own network" = "local_interest_set"
                                           ),
                                           selected = "online_interest_set"
                              ),

                              fileInput("file1", "Upload network file (in tab-deliminated format)",
                                        multiple = TRUE,
                                        accept = c("text/txt",
                                                   "text/comma-separated-values,text/plain",
                                                   ".txt")),
                              h6(HTML("<p>Check <a href = \"https://compbio-cityuhk.shinyapps.io/format/\" target = \"_blank\" >here</a> for format.")),
                              h6("Click",span(downloadLink("downloadexample1", "here")),"to download an example."),



                              ###################################
                              # Horizontal line ----
                              tags$hr(),
                              radioButtons("Hits_set", "2. Select or upload a gene signature ",
                                           c("Online gene signature" = "online_hit_set",
                                             "Upload a gene signature" = "local_hit_set"
                                           ),
                                           selected = "online_hit_set"
                              ),
                              selectizeInput(
                                "onlinehits", h5("Select from Gene Ontology or KEGG Pathways"), choices = GO_hit_set[,4], multiple = FALSE
                              ),
                              fileInput("file2", "Upload a gene signature",
                                        multiple = TRUE,
                                        accept = c("text/txt",
                                                   "text/comma-separated-values,text/plain",
                                                   ".txt")),
                              h6(HTML("<p> Check <a href = \"https://compbio-cityuhk.shinyapps.io/format/\" target = \"_blank\" >here</a> for format.")),
                              h6("Click",span(downloadLink("downloadexample2", "here")),"to download an example."),

                              br(),
                              actionButton("update_mra", "Analyze")

                            )),
                        mainPanel(
                          helpText("Master Regulator Analysis result"),
                          wellPanel(
                            #tabPanel("MRA Result",
                            DT::dataTableOutput("mra_result")
                            #)
                          )
                        )
                      )
             ),
             #           tabPanel("Search",
             #                   wellPanel(textInput("login", "Search", "")),
             #                   br(),
             #                  actionButton("update", "Search")
             #         ),
             tabPanel(h4("Download"),
                      titlePanel("PAVIRnet V1.0 (updated on 20180924)"),

                      # Button
                      downloadButton("downloadData", "Download")
             ),
             tabPanel(h4("Help"),
                      tags$img(src = "Help_MRA.png")
             ),
             tabPanel(h4("Contact us"),
                      h2("Contact Us"),
                      wellPanel(
                        h4("If you have any question or comments, please feel free to contact us."),
                        br(),
                        h4(strong("Dr. Xin Deng"),"(xindeng@cityu.edu.hk)"),
                        br(),
                        h4(strong("Address:"), "1B-106, 1/F, Block 1,To Yuen Building, City University of Hong Kong, 31 To Yuen Street ,
                           Kowloon Tong , Hong Kong SAR"),
                        br(),
                        h4(strong("Phone:"),"(852) 3442 5693"),
                        #tags$img(src = "cityu.png")
                        br(),
                        br(),
                        h4(strong("Dr. Xin Wang"),"(xin.wang@cityu.edu.hk)"),
                        br(),
                        h4(strong("Address:"), "1B-102, 1/F, Block 1,To Yuen Building, City University of Hong Kong, 31 To Yuen Street ,
                           Kowloon Tong , Hong Kong SAR"),
                        br(),
                        h4(strong("Phone:"),"(852) 3442 2367")
                        )
                      )

  )
  )

##############service####################
server <- function(input, output) {
  output$login_history <- renderText({
    input$Login
    isolate(input$login)
  })
  observeEvent(input$Login, {
    saveData(isolate(input$login))
  })
  saveData <- function(data) {
    write.table(
      x = data,
      file = "login_history.txt",
      row.names=F,col.names=F,quote=F,sep="\n",append=T)
  }

  net_to_visualization <- eventReactive(input$update, {
    if(input$network_show == "netshowall") {
      return(network_P5)
    }
    else {
      RN <- c()
      targets_count <- c()

      for(j in 1 : length(input$network_showpart)){
        iRN <- read.delim(paste("data/Regulon/",input$network_showpart[j],"_regulatory_network.txt",sep=""),header=F,sep="\t")
        targets_count_temp <- data.frame(as.matrix(iRN[1,1]),as.matrix(iRN[1,2]),dim(iRN)[1])
        RN <- rbind(RN,iRN)
        targets_count <- rbind(targets_count,targets_count_temp)

      }

      colnames(targets_count) <- c("TF ID", "TF name", "No. of targets")
      rownames(targets_count) <- c(1: dim(targets_count)[1])
      targets_summary <- RN[,c(1,2,3,4)]
      colnames(targets_summary) <- c("TF ID","TF name","Target ID","Target name")
      rownames(targets_summary) <- c(1: dim(targets_summary)[1])

      links_temp1 <- as.data.frame(RN[,c(1,3,5)])
      colnames(links_temp1) <- c("from","to","type")

      nodes_temp1_id <- unique(c(as.matrix(RN[,1]),as.matrix(RN[,3])))
      nodes_temp1 <- data.frame(nodes_temp[nodes_temp1_id,])
      rownames(nodes_temp1) <- c(1:dim(nodes_temp1)[1])
      colnames(nodes_temp1) <- c("id","name","expression","weight")

      vis.nodes <- nodes_temp1
      vis.links <- links_temp1

      vis.nodes$title  <- vis.nodes$name
      vis.nodes$shape <- c("dot")
      vis.nodes$shape[which(vis.nodes$weight == 10)] <- c("square")
      vis.nodes$size   <- vis.nodes$weight*2
      vis.nodes$borderWidth <- 2

      vis.nodes$label  <- vis.nodes$name
      vis.nodes$color.background <- c("slategrey")

      vis.nodes$color.background[which(vis.nodes$expression == -1)] <- c("#4682B4")
      vis.nodes$color.background[which(vis.nodes$expression == 1)] <- c("#F9CF9A")
      vis.nodes$color.background[which(vis.nodes$weight == 10)] <- c("darkred")
      vis.nodes$color.border <- "white"

      vis.links$color<- "gray"
      vis.links$color[which(links_temp1$type == -1)] <- c("#CD5C5C")
      vis.links$color[which(links_temp1$type == 1)] <- c("#2E8B57")
      # line color

      vis.links$arrows <- "to" # arrows: 'from', 'to', or 'middle'
      vis.links$smooth <- FALSE    # should the edges be curved?
      vis.links$shadow <- FALSE    # edge shadow





      sep_net <- list(targets_count,targets_summary,vis.nodes,vis.links)
      names(sep_net) <- c("network_summary","targets_summary","vis.nodes","vis.links")




      return(sep_net)

    }

  }, ignoreNULL = FALSE)


  output$network_plot <- renderVisNetwork({visNetwork(net_to_visualization()$vis.nodes, net_to_visualization()$vis.links)})


  #bounded = TRUE
  output$network_summary <- DT::renderDataTable(DT::datatable(
    net_to_visualization()$network_summary
  ))
  output$target_summary <- DT::renderDataTable(DT::datatable(
    net_to_visualization()$targets_summary
  ))

  ##############upload for MRA#####################
  mra_data <- eventReactive(input$update_mra, {
    if(input$Interest_set == "online_interest_set") {
      mra_tf_lists = tf
      mra_tf_list = mra_tf_lists[,2]
      content_interestset = whole_net[,c(1,3)]
    }
    else if(input$Interest_set == "local_interest_set") {
      req(input$file1)
      content_interestsets <-read.table(input$file1$datapath,header = T, sep = "\t")
      content_interestset = content_interestsets[,c(1,3)]
      mra_tf_lists <- unique(content_interestsets[,c(2,1)])
      mra_tf_list <- mra_tf_lists[,2]
    }

    if(input$Hits_set == "online_hit_set"){
      GO_hit_set1 <- as.matrix(GO_hit_set)
      content_hits <- unlist(strsplit(as.matrix(GO_hit_set[which(GO_hit_set[,4] == as.character(input$onlinehits)),2])," "))

    }
    else{
      #content_hits <- renderTable({
      req(input$file2)

      content_hits  <- read.table(input$file2$datapath,
                                  header = T,
                                  sep = "\n"
      )

    }

    universe <- unique(union(content_interestset[,1],content_interestset[,2]))
    universe.number <- length(universe)
    total.hits <- unique(intersect(as.matrix(universe),as.matrix(content_hits)))
    total.hits.number <- length(total.hits)
    mra_results<-c()
    # if(total.hits.number > 0){
    for(itf in 1 : length(mra_tf_list)){
      subnet <- content_interestset[which(content_interestset[,1] == as.character(mra_tf_list[itf])),]
      target.genes <- unique(subnet[,2])
      target.genes.number <- length(target.genes)
      observed.Hits <- length(intersect(target.genes,total.hits))
      pval <- round(phyper(observed.Hits - 1, m = total.hits.number, n = universe.number - total.hits.number, k = target.genes.number, lower.tail=F ),digits = 4)
      mra_results <- rbind(mra_results,data.frame(mra_tf_lists[itf,],universe.number,target.genes.number,total.hits.number,observed.Hits,pval))
    }

    #pvals.adj <- p.adjust(mra_results[,6], method="BH")
    #mra_results <- cbind(mra_results,pvals.adj)
    mra_results <- data.frame(mra_results[order(mra_results[,7]),])
    mra_results <- mra_results[,-3]
    colnames(mra_results) <- c("TF name","TF ID","No. of targets", "Total No. of hits", "Obseved hits", "P-value")
    mra_results[which(mra_results[,6] == 0),6] <- c("< 1e-4")
    mra_results <- mra_results[which(mra_results[,5]  != 0),]
    return(mra_results)

  })

  output$mra_result <- DT::renderDataTable(DT::datatable(
    mra_data(),rownames=FALSE
  ))



  ###################do mra#######################

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("PAVIRnet", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(netfordownload,file, row.names = FALSE,col.names =FALSE)
    }
  )
  #############################
  output$downloadexample1 <- downloadHandler(
    filename = function() {
      "Network.example.txt"
    },
    content = function(file) {
      write.table(network_example,file, row.names = FALSE,col.names =FALSE,quote=F,sep="\t")
    }
  )
  output$downloadexample2 <- downloadHandler(
    filename = function() {
      "Gene_signature.example.txt"
    },
    content = function(file) {
      write.table(signature_example,file, row.names = FALSE,col.names =FALSE,quote=F,sep="\t")
    }
  )
}
  ###################APP##################

  shinyApp(ui = ui, server = server)

