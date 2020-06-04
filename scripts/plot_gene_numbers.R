library(ggplot2)
m <- modules::use("~/git_repos/Trichuris_transwells/rnaseq.R")

genes<-read.table('~/git_repos/single_cell/experiment_4/regulated_genes.txt',header=T)
cbPalette <- c("#7851a9", "#FFFF99")
p<-ggplot(genes,x=time,y=n.genes)+
  geom_col(aes(time,n.genes,fill=regulation))+
  scale_fill_manual(values=cbPalette)+
#  facet_grid(~eval(substitute(facet), df),scales="free", space="free_x", labeller=labeller(labels))+
# scale_fill_manual(values=cbPalette)+
  theme(text = element_text(size=20))+
  ylab('Number of genes')+
  xlab('Timepoint')+
  labs(fill='Regulation (log2fc >1.5)')

pdf('~/git_repos/single_cell/experiment_4/ngenes.pdf')
print(p)
dev.off()
