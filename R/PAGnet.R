# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
library("shiny")
library("ggplot2")
library("scales")
library("networkD3")
library("shinythemes")
library("visNetwork")
###################################################
###############network summary#####################
#network_summary <- read.csv("data/network/network summary.csv",header=T)

#network_summary_P5 <- network_summary[,c(1,2,3)]
#colnames(network_summary_P5) <- c("TF ID", "TF name", "No. of targets")


#tf <- read.delim("data/network/TF_list.txt",sep="\t",header=F)

#network_inf <- read.delim("data/network/PAVIRnet_with_info.txt",header=F,sep="\t")

#nodes <- read.csv("data/network/PAVIRnet_nodes2.csv",header=T)
#links <- read.csv("data/network/PAVIRnet_edge_v5.csv",header=T)
#nodes_temp <- nodes
#rownames(nodes_temp) <- nodes_temp[,1]

#vis.nodes <- nodes
#vis.links <- links
#vis.nodes$title  <- vis.nodes$name
#vis.nodes$shape <- c("dot")
#vis.nodes$shape[which(vis.nodes$weight == 10)] <- c("square")
#vis.nodes$size   <- vis.nodes$weight*2
#vis.nodes$borderWidth <- 2

#vis.nodes$label  <- vis.nodes$name
#vis.nodes$color.background <- c("slategrey")

#vis.nodes$color.background[which(vis.nodes$expression == -1)] <- c("#4682B4")
#vis.nodes$color.background[which(vis.nodes$expression == 1)] <- c("#F9CF9A")
#vis.nodes$color.background[which(vis.nodes$weight == 10)] <- c("darkred")
#vis.nodes$color.border <- "white"
#vis.links$color <- "gray"    # line color
#vis.links$color[which(vis.links$type == "1")] <- c("#2E8B57")
#vis.links$color[which(vis.links$type == "-1")] <- c("#CD5C5C")
#vis.links$arrows <- "to" # arrows: 'from', 'to', or 'middle'
#vis.links$smooth <- FALSE    # should the edges be curved?
#vis.links$shadow <- FALSE    # edge shadow

#signature_example <- read.table("data/Gene_signature.example.txt",header=F,sep="\t")
#network_example <- read.table("data/Network.example.txt",header=F,sep="\t")

#netfordownload <- read.csv("data/network/PAVIRnet.csv",header=T,check.names = F)

#save(network_summary,network_summary_P5,tf,network_inf,nodes,links,nodes_temp,
#     vis.nodes,vis.links,signature_example,network_example,netfordownload,file="Pavirnet_basicdata.Rda")
load("data/Pavirnet_basicdata.Rda")
################q0.05#############################################

#node_Q5 <- read.csv("data/5net/net_v3_node_attr.csv",header=F)
#colnames(node_Q5) <- c("name","group","size")
#edge_code_Q5 <- read.csv("data/5net/net_v3_edge_attr_onlycode.csv",header=F)
# #colnames(edge_code_Q5) <- c("source","target","value")

##################################################################
targets_summary <- netfordownload
colnames(targets_summary) <- c("TF ID","TF name","Target ID","Target name")
network_P5 <- list(network_summary_P5,targets_summary,vis.nodes,vis.links)
names(network_P5) <- c("network_summary","targets_summary","vis.nodes","vis.links")

#network_P3 <- list(network_summary_P3,node_P3,edge_code_P3)
#names(network_P3) <- c("network_summary","node","edge_code")

#network_Q5 <- list(network_summary_Q5,node_Q5,edge_code_Q5)
#names(network_Q5) <- c("network_summary","node","edge_code")
##################################################################
#login_history <- c("login history")
##########################Go Hit sets############################
GO_hit_set <- read.delim("data/mra/Hit_sets.txt",header=F,sep="\t")
whole_net <- read.delim("data/mra/PAVIRnet.txt",header=T,sep="\t")
##############ui#########################

PAGnet_mra <- function() {

  if(input$Interest_set == "online_interest_set") {
    mra_tf_lists = tf
    mra_tf_list = mra_tf_lists[,2]
    content_interestset = whole_net[,c(1,3)]
  }else if(input$Interest_set == "local_interest_set") {
    req(input$file1)
    content_interestsets <-read.table(input$file1$datapath,header = T, sep = "\t")
    content_interestset = content_interestsets[,c(1,3)]
    mra_tf_lists <- unique(content_interestsets[,c(2,1)])
    mra_tf_list <- mra_tf_lists[,2]
  }

  if(input$Hits_set == "online_hit_set"){
    GO_hit_set1 <- as.matrix(GO_hit_set)
    content_hits <- unlist(strsplit(as.matrix(GO_hit_set[which(GO_hit_set[,4] == as.character(input$onlinehits)),2])," "))

  }else{
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




}
