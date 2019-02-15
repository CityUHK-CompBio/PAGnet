
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

pagnet.mra.interface <- function(){
  ##############ui#########################
  ui <- fixedPage(
    h3(strong("PAGnet:"),"Pseudomonas aeruginosa genomic integrated regulatory network",
       style = "font-family: 'Source Sans Pro';
       color: #123456; text-align: center;
       background-image: url('background-2.png');
       padding: 20px"
    ),
    navbarPage("",
               tabPanel(h4("Master Regulator Analysis"),
                        sidebarLayout(
                          div(style = "wigth:100px;",
                              # Sidebar panel for inputs ----
                              sidebarPanel(

                                # Input: Select a file ----
                                radioButtons("Interest_set", "1. Choose PAGnet or upload your own network",
                                             c("PAGnet" = "online_interest_set",
                                               "Upload your own network" = "local_interest_set"
                                             ),
                                             selected = "online_interest_set"
                                ),

                                fileInput("file1", "Upload network file (in tab-deliminated format)",
                                          multiple = TRUE,
                                          accept = c("text/txt",
                                                     "text/tab-separated-values,text/plain",
                                                     ".txt")),


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
               tabPanel(h4("Help"),
                        tags$img(src = "data/Help_MRA.png")
               )

    )
    )

  ##############service####################
  server <- function(input, output) {

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

}

