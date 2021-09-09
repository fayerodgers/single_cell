m <- modules::use("~/git_repos/Trichuris_transwells/rnaseq.R")
library(ggplot2)
library(scales)

# paths
d1.go<-read.table('/Users/fr7/git_repos/single_cell/experiment_4/bulk_libraries/time.status.24h_control.24h_infected.bulk.sort/GO.BP_results.tsv',header=T, stringsAsFactors = FALSE)
d3.go<-read.table('/Users/fr7/git_repos/single_cell/experiment_4/bulk_libraries/time.status.72h_control.72h_infected.bulk.sort/GO.BP_results.tsv',header=T, stringsAsFactors = FALSE)
whole.go<-read.table('/Users/fr7/timecourse/kallisto_analysis/as_time_course/GO.BP_results.tsv',header=T, stringsAsFactors = FALSE)

figures_dir<-'/Users/fr7/git_repos/single_cell/figuresNov20'
files_dir<-'/Users/fr7/git_repos/single_cell/SupplementaryFiles'

# define GO terms to plot

go_to_plot<-c(
  "GO:0002221",
  "GO:0042742",
  "GO:0002218",
  "GO:0009615",
  "GO:0051707",
  "GO:0045087",
  "GO:0035456",
  "GO:0044419",
  "GO:0034340",
  "GO:0035455"
)

lookup_ngenes<-function(go,go.table){
  ngenes<-go.table[which(go.table$category == go),'numDEInCat']
  return(ngenes)
}

lookup_pval<-function(go,go.table){
  pval<-go.table[which(go.table$category == go),'over_represented_pvalue']
  return(pval)
}

lookup_description<-function(go,go.table){
  desc<-go.table[which(go.table$category == go),'term']
  return(desc)
}

# d1
go.df.d1 <- data.frame("pathway_id" = go_to_plot, stringsAsFactors = FALSE )
go.df.d1$time <- 'Epithelium, D1'
go.df.d1$ngenes <- sapply(go_to_plot,lookup_ngenes,d1.go)
go.df.d1$pval <- sapply(go_to_plot,lookup_pval,d1.go)
go.df.d1$description<-sapply(go_to_plot,lookup_description,d1.go)

# d3
go.df.d3 <- data.frame("pathway_id" = go_to_plot, stringsAsFactors = FALSE )
go.df.d3$time <- 'Epithelium, D3'
go.df.d3$ngenes <- sapply(go_to_plot,lookup_ngenes,d3.go)
go.df.d3$pval <- sapply(go_to_plot,lookup_pval,d3.go)
go.df.d3$description<-sapply(go_to_plot,lookup_description,d3.go)

# whole
go.df.whole <- data.frame("pathway_id" = go_to_plot, stringsAsFactors = FALSE )
go.df.whole$time <- 'Complete caecum'
go.df.whole$ngenes <- sapply(go_to_plot,lookup_ngenes,whole.go)
go.df.whole$pval <- sapply(go_to_plot,lookup_pval,whole.go)
go.df.whole$description<-sapply(go_to_plot,lookup_description,whole.go)

go.df <- rbind(go.df.d1,go.df.d3,go.df.whole)

# plot

F1<-ggplot(go.df,x=time,y=description)+
  geom_point(aes(x=time,y=description, size=ngenes, col=pval))+
  scale_color_gradient(limits = c(0.0001,0.01),
                       oob=squish, 
                       name = "Pathway p-value\n(adjusted)", 
                       breaks = c(0.0001,0.01),
                       labels = c("< 0.0001","> 0.01"),
                       guide = guide_colourbar())+
  scale_size_continuous(name="Number of\ngenes", range = c(0.1,2))+
  theme(text = element_text(size=5),
        axis.text.x=element_text(size=4, angle=45, hjust = 1),
        legend.title = element_text(size = 4), 
        legend.text = element_text(size = 4),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(-10,-10,-10,-10) )+
  labs(x=NULL,y=NULL)+
  guides(color = guide_colourbar(barwidth = 0.5, barheight = 2), 
         size = guide_legend(keywidth = 0.5, keyheight = 0.5))

pdf(file.path(figures_dir,'F1.pdf'),2.5,1.8)
print(F1)
dev.off()
