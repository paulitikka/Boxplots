# Making a package for Steroid evaluations. Tikka, 1.11.24.

#https://ourcodingclub.github.io/tutorials/writing-r-package/
# https://docs.posit.co/ide/user/ide/guide/pkg-devel/writing-packages.html

# install.packages(c("devtools", "roxygen2", "testthat", "knitr"))

#' Boxplottings
#'
#' Converts data to boxplots.
#' Written by Pauli Tikka. University of Turku. Updated 1.11.24.
#' @param tvt The dataset
#' @return The plot
#' @examples
#' temp1 <- boxplots;
#' @export
boxplots=function(tvt,Group,Outcome,Out,oute,other) {

  if (Group=='Male') {tvt=tvt[tvt[,'Gender']==1,]} else if (Group=='Female')
  {tvt=tvt[tvt[,'Gender']==0,]} else if (Group=='All') {tvt=tvt}
  Steroid=rep(colnames(tvt[,9:28]), each=dim(tvt)[1])
  data2=rep('Control',dim(tvt)[1])
  num=min(tvt[,Outcome])
  if (Outcome=='HOMA-IR') {num=1.5}
  data2[tvt[,Outcome]>num]='Case' #'Steatosis.Grade.0.To.3' #
  Treatment=data2
  note=unlist(tvt[,9:28])
  Concentration=as.vector(note)
  data=data.frame(Steroid, Treatment ,  Concentration)
  data[,'Group'] = 0 # data$Steroid #check the Xs etc out..
  data$Steroid [data$Steroid  == '17aOH-P4']='17a-OHP4'
  # groups$Abbreviation[groups$Abbreviation == '17a-OHP4']='17aOH-P4'
  for (i in 1:21) {data[data$Steroid %in% groups$Abbreviation[i],'Group']=groups$Group[i]}
  title = paste(Out,"'s Effect to Concentrations of Steroids", ' in ',Group, ' Subjects',sep="")
  if (Group=='Male') {lep=theme(legend.position = "none")} else if (Group=='Female')
  {lep=theme(legend.position = "none")} else if (Group=='All') {lep=theme_classic2()+theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2,size = 14))}
  # lep=theme(legend.position = "none")
  if (num==1.5) {e1=paste('Case (>=',num,')',sep="");e2=paste('Control (<',num,')',sep="");
  e3=paste('****<0.001', '***<0.01', '**<0.05', '*<0.1')} else {e1=paste('Case (>',0,')',sep="");e2=paste('Control (=',0,')',sep="");e3=paste('****<0.001', '***<0.01', '**<0.05', '*<0.1')}

  e3=paste('****<0.001', '***<0.01', '**<0.05', '*<0.1')

  data=data[!is.na(data$Concentration),]
  # grouped boxplot: https://stackoverflow.com/questions/32539222/group-boxplot-data-while-keeping-their-individual-x-axis-labels-in-ggplot2-in-r
  p=ggplot(data, aes(x=Steroid, y=Concentration, fill=Treatment))+
    geom_boxplot(notch=F, notchwidth=0.5,outlier.shape=1,outlier.size=2, coef=1.5)+
    # coord_cartesian(ylim = c(0, 60000))+
    theme(axis.text=element_text(color="black"))+
    theme_classic2()+#https://stackoverflow.com/questions/34522732/changing-fonts-in-ggplot2
    # scale_x_discrete(guide = guide_axis(angle = 90))+ https://stackoverflow.com/questions/37488075/align-axis-label-on-the-right-with-ggplot2
    theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2,size = 14))+#annotate(geom="text",family="Broadway",size=20)+#annotate(geom="text",family="Calibri",size = 14)+
    theme(panel.grid.minor=element_blank())+ #http://www.sthda.com/english/wiki/ggplot2-rotate-a-graph-reverse-and-flip-the-plot
    labs(size= "Type",x = "Steroids",y = "Log2 of Picomolar Concentrations ", title=title,size = 14)+ #log2 Autoscaled
    scale_fill_manual(values=c("orange","blue"),name=c(oute),labels=c(e1,e2))+#abels=c("Case (>5)", "Control (=<5)"))
    facet_grid(~Group, scales = "free_x", space = "free")+lep+
    theme(text=element_text(size=10.5,family="Calibri"), #change font size of all text
          axis.text=element_text(size=14), #change font size of axis text
          axis.title=element_text(size=14), #change font size of axis titles
          plot.title=element_text(size=14), #change font size of plot title
          legend.text=element_text(size=14), #change font size of legend text
          legend.title=element_text(size=14))+theme(axis.text=element_text(color="black"))+stat_compare_means(hide.ns = TRUE,
                                                                                                              label= "p.signif",method = "wilcox.test",  symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                                                                                                                                                            symbols = c("****", "***", "**", "*", "ns")))#+annotation_custom( grid::textGrob(x =1, y = 0.1, label = "OK"))#+showSignificance( c(1,2), 10, -0.05, "*")
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  # https://stackoverflow.com/questions/76758153/is-there-a-way-to-change-the-asterisks-to-match-custom-p-values
  #+annotate(geom="text",family="Broadway",size=20) # add_pval(plot, pairs = list(c(1, 2)), test='wilcox.test');plot
  #https://datavizpyr.com/horizontal-boxplots-with-ggplot2-in-r/
  # https://stackoverflow.com/questions/72564551/a-custom-legend-unrelated-to-data-in-ggplot
  # p #+ annotation_custom( grid::textGrob(x =3, y = 5, label = "'****<0.001', '***<0.01', '**<0.05', '*<0.1'")).. adding a custom label (unrelated to data) to ggplot is again one of these not so easy but should be easy things... :)

  library(ragg)
  # path="C:/Users/patati/Documents/GitHub/new/" #oh, classical: https://forum.posit.co/t/r-markdown-html-document-doesnt-show-image/41629/2
  # pngfile <- fs::path(path,paste0(Group,Out,'box',".png"))#fs::path(knitr::fig_path(),  "theming2.png")
  # ragg::agg_png(pngfile, width = 60, height = 46, units = "cm", res = 300,scaling = 2)
  #https://stackoverflow.com/questions/66429500/how-to-make-raggagg-png-device-work-with-ggsave
  plot(p)
  # invisible(dev.off())
  # knitr::include_graphics(pngfile)
  }

#
# library(roxygen2); # Read in the roxygen2 R package
# gc()
# roxygenise();      # Builds the help files; read these separately not in a project I guess....
