library("ggplot2")
library("jsonlite")

rm(list=ls())
setwd("~/Desktop/biomodels_s2d")

source("0_Functions.R")

# model_ids <- c("BIOMD0000000049EGF","BIOMD0000000049NGF")
model_id <- "BIOMD0000000049"

reference <- read.table("sc_reference.tsv", sep = "\t" , header = F , stringsAsFactors = F)
test <- read.table("sc_test.tsv" , sep = "\t" , header = F , stringsAsFactors = F )

# for(model_id in model_ids){
  
  model <- SBMLDocument_getModel(readSBML(paste0("models/",model_id,"_url.xml")))
  
  new_concentrations <- new_initialConc(model,reference,test)
  # new_concentrations[c(96,97),] <- 0 # no stimulation
  simulate_model(model_id, new_concentrations, c(0,0,4800,1000))

  # post simulation ---------------------------------------------------------
  
  merged_data <- merge_reports(model_id)
  merged_data <- merged_data[ , c(1:2,grep("ERK",colnames(merged_data))) ]

  percent_pERK <- 100*rowSums(merged_data[,grepl("ppERK", colnames(merged_data))]) /rowSums(merged_data[,grepl("ERK", colnames(merged_data))])
  percent_pERK[is.nan(percent_pERK)] <- 0
  merged_data <- cbind(merged_data,percent_pERK)
  
  # generate_plots(merged_data)
  
  clustering <- hclust(dist(t(sapply(unique(merged_data$Id), function(i) percent_pERK[merged_data$Id == i]))))
  # plot(clustering)

# }

# section -----------------------------------------------------------------

is_sustained <- sapply(unique(merged_data$Id), function(id){
  temp <- percent_pERK[merged_data$Id == id]
  if(max(temp) == 0) FALSE else
  (temp[length(temp)]/max(temp) > 0.5) & (which(temp == max(temp)) < 1000)
})

data <- lapply(clustering$order , function(i) {
  data.frame(Time = merged_data$Time[merged_data$Id == i] , 
             Id = which(clustering$order %in% i) ,
             clust = cutree(clustering , k = 2)[i] ,
             sustained = is_sustained[i],
             pERK = percent_pERK[merged_data$Id == i])
})
data <- Reduce(rbind,data)

p1 <- ggplot(data,aes(x= Time ,y = pERK , color=clust, group = as.factor(Id)))+
  geom_line() +theme_bw()+guides(color=FALSE)

p2 <- ggplot(data,aes(x= Time ,y = pERK , color=sustained, group = as.factor(Id)))+
  geom_line() +theme_bw()+guides(color=FALSE)

# ggplot(data,aes(x= Time ,y = as.factor(Id) , fill = pERK,color=pERK))+ geom_tile(size=0.6) +
#   scale_color_viridis_c()+scale_fill_viridis_c()

# ggsave(paste0(gsub("0","",model_id),"_precent_phosphorylated_clustered.png"), p1, height = 10 , width = 10)
# ggsave(paste0(gsub("0","",model_id),"_precent_phosphorylated_sustain.png"), p2, height = 10 , width = 10)

# ppERK -------------------------------------------------------------------

meta <- fromJSON("extras/metadata.cart.2018-10-29.json")
barcode <- sapply(meta$associated_entities, function(i) i$entity_submitter_id)
barcode <- substring(barcode,1,12)

clinical <- read.table("extras/clinical.tsv" , sep = "\t" , header = TRUE , na.strings = "--")
clinical <- clinical[sapply(barcode, function(i) which(clinical$submitter_id == i)),]

length(unique(merged_data$Id)) == nrow(clinical)

temp <- sapply(sort(unique(clinical$days_to_death)) , function(day){
  g1 <- clinical[is_sustained,]
  g2 <- clinical[!is_sustained,]
  g1 <-100- sum(g1$vital_status[g1$days_to_death <= day] == "dead",na.rm = T)*100/nrow(g1)
  g2 <-100- sum(g2$vital_status[g2$days_to_death <= day] == "dead",na.rm = T)*100/nrow(g2)
  c(g1,g2)
})

data2 <- data.frame( days = sort(unique(clinical$days_to_death)) ,
                     sustained = temp[1,],
                     not_sustained = temp[2,])

ggplot(data2, aes(x = days))+ 
  geom_line(aes(y=sustained),color="red")+
  geom_line(aes(y=not_sustained) , color = "blue")+ 
  # ylim(0,100)+
  ylab("survival %")+
  ggtitle("KM plot sustained(r) transient(b)")+
  theme_bw()

is_sustained <- sapply(unique(merged_data$Id), function(id){
  temp <- merged_data$ppERK[merged_data$Id == id]
  temp[length(temp)] > 0.017
})

data <- data.frame(merged_data,
                   sustain = rep(is_sustained, each = sum(merged_data$Id ==1)))
ggplot(data, aes(x=Time, y = ppERK, color = sustain , group = as.factor(Id)))+ geom_line()


stage <- sub("stage ", "", clinical$tumor_stage)
stage <- sub("not reported" , NA , stage)
stage <- gsub("a|b|c" , "" , stage)
stage_34 <- grepl("iii|iv", stage)
stage_12 <- stage %in% c("ii","i","i/ii nos")

is_alive <- clinical$vital_status == "alive"

chisq.test(matrix(c(sum(stage_12 & is_sustained),
                    sum(stage_12 & !is_sustained),
                    sum(stage_34 & is_sustained),
                    sum(stage_34 & !is_sustained)), nrow = 2, byrow = TRUE),
           correct = FALSE)

chisq.test(matrix(c(sum(stage_12 & is_alive),
                    sum(stage_12 & !is_alive),
                    sum(stage_34 & is_alive),
                    sum(stage_34 & !is_alive)), nrow = 2, byrow = TRUE),
           correct = FALSE)

# sapply(unique(data$Id), function(i){
#   sum(data$ppERK[data$Id == i])
# })
plot(hist(sapply(unique(data$Id), function(i){
  sum(data$ppERK[data$Id == i])}),
  breaks = 50))

identical(data_egf$Id,data_ngf$Id)

data <- data.frame(id = data_egf$Id,
                   Time = data_egf$Time,
                   egf_erk = data_egf$ppERK,
                   ngf_erk = data_ngf$ppERK)
data <- melt(data,id.vars = c("Time","id"))
head(data)
data$id <- paste0(data$id,data$variable)
ggplot(data, aes(x=Time, y = id , color = value , fill = value))+
  geom_tile(size = 0.6)+scale_fill_viridis_c()+ scale_color_viridis_c()
ggsave("plot.png",height = 30 , width = 30)



# ngf vs egf --------------------------------------------------------------

egf_data <- read.csv("random/BIOMD49EGF.csv")
ngf_data <- read.csv("random/BIOMD49NGF.csv")

plots <- lapply(3:ncol(egf_data) , function(i){
  data <- data.frame(id = egf_data$Time,
                     egf = egf_data[,i],
                     ngf = ngf_data[,i])
  p <- ggplot(data,aes(x=egf,y=ngf, color=as.factor(id) ))+geom_point()+guides(color=FALSE)
  ggsave(paste0("plots/egf_ngs_",colnames(egf_data)[i],".png") , p ,height = 10 , width = 10)
})



# plotting- 26-11 ---------------------------------------------------------

is_sustained <- sapply(unique(merged_data$Id), function(id){
  temp <- percent_pERK[merged_data$Id == id]
  (temp[length(temp)]/max(temp) > 0.4) & all((which(temp == max(temp)) < 1500))
})

data <- lapply(clustering$order , function(i) {
  data.frame(Time = merged_data$Time[merged_data$Id == i] , 
             Id = which(clustering$order %in% i) ,
             clust = cutree(clustering , k = 2)[i] ,
             sustained = is_sustained[i],
             pERK = percent_pERK[merged_data$Id == i])
})
data <- Reduce(rbind,data)

p1 <- ggplot(data,aes(x= Time ,y = pERK , color=clust, group = as.factor(Id)))+
  geom_line() +theme_bw()+guides(color=FALSE)

p2 <- ggplot(data,aes(x= Time ,y = pERK , color=as.factor(sustained), group = as.factor(Id)))+
  geom_line() +theme_bw()+guides(color=FALSE)+ scale_color_hue( l = 55)+ 
  scale_x_continuous(limits = c(0,4500), expand=c(0,0)) +
  scale_y_continuous(limits = c(0,100), expand=c(0,0))+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size= 13),
        axis.title = element_text(size= 16,face = "bold"),
        plot.title = element_text(hjust = 0.5))+
  labs(x="Time(s)",y="% ERK phosphorylated",title = "NGF")
 ggsave("EGF_BIOMD49.png", p2)
 