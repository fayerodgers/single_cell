m <- modules::use("~/git_repos/Trichuris_transwells/rnaseq.R")
library(ggplot2)
library(scales)

# paths
e1.go<-read.table('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/control/enterocytes_gprofiler/gProfiler_entero.1.csv',header=T, stringsAsFactors = FALSE, sep =",")
e2.go<-read.table('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/control/enterocytes_gprofiler/gProfiler_entero.2.csv',header=T, stringsAsFactors = FALSE, sep =",")
e3.go<-read.table('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/control/enterocytes_gprofiler/gProfiler_entero.3.csv',header=T, stringsAsFactors = FALSE, sep =",")
isg15.go<-read.table('/Users/fr7/git_repos/single_cell/experiment_4/FINAL/integrated/control/enterocytes_gprofiler/gProfiler_entero.Isg15.csv',header=T, stringsAsFactors = FALSE, sep =",")

  
figures_dir<-'/Users/fr7/git_repos/single_cell/figuresNov20'
files_dir<-'/Users/fr7/git_repos/single_cell/SupplementaryFiles'

# big table of all terms and descriptions
all.gos<-data.frame("term"=c(e1.go$term_id, e2.go$term_id, e3.go$term_id, isg15.go$term_id), "desc"=c(e1.go$term_name, e2.go$term_name, e3.go$term_name, isg15.go$term_name))
all.gos<-unique(all.gos)

# define GO terms to plot

go_to_plot<-c(
  "GO:0000103",
  "GO:1901615",
  "GO:0032528",
  "GO:0048002",
  "GO:0006811",
  "GO:0006820",
  "GO:0050892",
  "GO:0007588",
  "GO:0009617",
  "GO:0034340"
)

lookup_ngenes<-function(go,go.table){
  ngenes<-go.table[which(go.table$term_id == go),'intersection_size']
  if (length(ngenes) == 0 ){
    ngenes = 0
  }
  return(ngenes)
}

lookup_pval<-function(go,go.table){
  pval<-go.table[which(go.table$term_id  == go),'adjusted_p_value']
  if (length(pval) == 0 ){
    pval = 1
  }
  return(pval)
}

lookup_description<-function(go,go.table){
  desc<-go.table[which(go.table$term == go),'desc']
  return(desc)
}

# e1
go.df.e1 <- data.frame("pathway_id" = go_to_plot, stringsAsFactors = FALSE )
go.df.e1$cluster<- 'Entero.1'
go.df.e1$ngenes <- sapply(go_to_plot,lookup_ngenes,e1.go)
go.df.e1$pval <- sapply(go_to_plot,lookup_pval,e1.go)
go.df.e1$description<-sapply(go_to_plot,lookup_description,all.gos)

# e2
go.df.e2 <- data.frame("pathway_id" = go_to_plot, stringsAsFactors = FALSE )
go.df.e2$cluster<- 'Entero.2'
go.df.e2$ngenes <- sapply(go_to_plot,lookup_ngenes,e2.go)
go.df.e2$pval <- sapply(go_to_plot,lookup_pval,e2.go)
go.df.e2$description<-sapply(go_to_plot,lookup_description,all.gos)

# e3
go.df.e3 <- data.frame("pathway_id" = go_to_plot, stringsAsFactors = FALSE )
go.df.e3$cluster<- 'Entero.3'
go.df.e3$ngenes <- sapply(go_to_plot,lookup_ngenes,e3.go)
go.df.e3$pval <- sapply(go_to_plot,lookup_pval,e3.go)
go.df.e3$description<-sapply(go_to_plot,lookup_description,all.gos)

# isg15
go.df.isg15 <- data.frame("pathway_id" = go_to_plot, stringsAsFactors = FALSE )
go.df.isg15$cluster<- 'Entero.Isg15'
go.df.isg15$ngenes <- sapply(go_to_plot,lookup_ngenes,isg15.go)
go.df.isg15$pval <- sapply(go_to_plot,lookup_pval,isg15.go)
go.df.isg15$description<-sapply(go_to_plot,lookup_description,all.gos)

go.df <- rbind(go.df.e1,go.df.e2,go.df.e3,go.df.isg15)

# plot

F18<-ggplot(go.df,x=cluster,y=description)+
  geom_point(aes(x=cluster,y=description, size=ngenes, col=pval))+
  scale_color_gradient(limits = c(0.0001,0.05),
                       oob=squish, 
                       name = "Pathway p-value\n(adjusted)", 
                       breaks = c(0.0001,0.05),
                       labels = c("< 0.0001","> 0.05"),
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

pdf(file.path(figures_dir,'F18.pdf'),2.5,1.8)
print(F18)
dev.off()
