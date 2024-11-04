# Making a package for Steroid evaluations. Tikka, 1.11.24.

#https://ourcodingclub.github.io/tutorials/writing-r-package/
# https://docs.posit.co/ide/user/ide/guide/pkg-devel/writing-packages.html

# install.packages(c("devtools", "roxygen2", "testthat", "knitr"))a

#' Forestplottings
#'
#' Converts data to forestplots.
#' Written by Pauli Tikka. University of Turku. Updated 1.11.24.
#' @param NAFLD The dataset with steroids and cases
#' @param Group All, Female, Male (typically a gender column)
#' @param Outcome The case (e.g. steatosis). Typically a column name in the dataset
#' @param name An header add (of the case) for the boxplot's header
#' @param oute The label title. Here, almost the same as case.
#' @return The boxplots. Yes.
#' @examples
#' ie=tv_half_log22;Outcome='Steatosis Grade';Out='Steatosis'; oute='Steatosis Grade';Group='All';
#' boxplots(ie,Group,Outcome,Out,oute)
#' @export
# Some further definitions are needed, yes

# This works with the autoscaled (raw if loge=1 and remove 1 in the means) data NAFLD as well...
forestplots=function(NAFLD,Outcome,Group,name,ordera,oute,first,e,xlim) { # Group='Female'

  if (Group=='Male') {NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==2,]} else if (Group=='Female')
  {NAFLDo=NAFLD[NAFLD[,'SEX.1F.2M']==1,]} else if (Group=='All') {NAFLDo=NAFLD}
  sample_data=c();n0=c();n1=c()

  for (i in 1:2) {
    if (i==1) {SG0=NAFLDo[NAFLDo[,Outcome] == 0,];n0=dim(SG0)[1]} else if (i==2) {SG0=NAFLDo[NAFLDo[,Outcome] > 0,];n1=dim(SG0)[1]}#Steatosis.Grade.0.To.3 Fibrosis.Stage.0.to.4
    #https://stats.stackexchange.com/questions/237256/is-the-shapiro-wilk-test-only-applicable-to-smaller-sample-sizes # https://statisticsbyjim.com/hypothesis-testing/nonparametric-parametric-tests/
    means=c();for (j in 9:28) {means=append(means,median(SG0[,j], na.rm=TRUE))}
    #https://www.statology.org/mean-standard-deviation-grouped-data/ # https://amsi.org.au/ESA_Senior_Years/SeniorTopic4/4h/4h_2content_11.html # https://www.themathdoctors.org/mean-and-standard-deviation-of-grouped-data/
    sds=c();for (j in 9:28) {sds=append(sds,sd(SG0[,j],na.rm=TRUE))} #here we are... :)
    error_lower=means-sds; error_upper=means+sds; error=sds
    sample_data <- append(sample_data,data.frame(study=colnames(NAFLD[,9:28]),index=colnames(NAFLD[,9:28]),result=means,error=error))} # cate<- (str_extract(colnames(SG0[,9:28]), "[aA-zZ]+"))    #https://stackoverflow.com/questions/29825537/group-categories-in-r-according-to-first-letters-of-a-string
  # https://datatofish.com/create-dataframe-in-r/#if you know: ga=groups[,'Abbreviation']; md=metad[,'name']; unique(c(ga,md)); ga[!(ga %in% md)]; ga[ga=="17a-OHP4"]="17aOH-P4"
  # https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot
  # http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
  # https://blogs.sas.com/content/iml/2011/04/27/log-transformations-how-to-handle-negative-data-values.html
  # https://stats.stackexchange.com/questions/155429/how-to-transform-negative-values-to-logarithms
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1534033/
  df=data.frame(sample_data) #
  ps=c();for (j in 9:28) {xnam <- colnames(NAFLDo)[j]; fmla <- as.formula(paste(xnam, "~",Outcome));
  ps=append(ps,wilcox.test(fmla, data = NAFLDo,exact = FALSE)$p.value)}#kruskal.test
  # http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
  # https://www.statisticshowto.com/probability-and-statistics/statistics-definitions/parametric-and-non-parametric-data/
  a=df[df[,1]==e,'result.1']/df[df[,1]==e,'result'];
  v2=data.frame(log(df$result.1/df$result)) #this should give the right order of the variables if not the absolute change
  v2[,'result']=v2[,1];v2[,'name']=df$study;v2=v2[,2:3]
  v2[,'name'] <- gsub("\\.", "-", v2[,'name']) #https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot
  v2[,'name'] <- gsub("X11", "11", v2[,'name'])
  v2[,'name'] <- gsub("X17", "17", v2[,'name'])
  v2[,'name'][v2[,'name']=="T-Epi-T"]="T/Epi-T"
  v2[,'pval']=ps# df2=df[rev(order(df[,'error'])),]
  # https://www.bmj.com/content/312/7038/1079.full# https://stats.stackexchange.com/questions/589920/how-can-i-back-transform-a-log-data-to-interpret-t-test-and-get-original-ci
  # https://www.biostars.org/p/16481/
  # https://whitlockschluter3e.zoology.ubc.ca/RLabs/R_tutorial_Contingency_analysis.html # https://sphweb.bumc.bu.edu/otlt/mph-modules/ep/ep713_randomerror/ep713_randomerror6.html
  # https://www.r-bloggers.com/2015/01/easy-error-propagation-in-r/# https://www.biostars.org/p/342756/ https://en.wikipedia.org/wiki/Tukey%27s_range_test
  #https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test
  v2[,'result_pure']=(df$result.1/df$result) #case and control
  v2[,'error']=(abs((1/df$result)*df$error.1)+abs((df$result.1/df$result^2)*df$error))/dim(NAFLDo)[1]*1.64
  # This 'error' definition above is ok... it is from the Jukka Vaari 1993 Fysiikan Laboratoriotyöt
  v2[,'error'][v2[,'error']>(median(v2[,'error'])+sd(v2[,'error']))]=median(v2[,'error'])*1.25
  v2[,'errord1a']=v2[,'result_pure']-v2[,'error']#(df$result.1-dy)/(df$result+dx)# #
  v2[,'errord2a']=v2[,'result_pure']+v2[,'error']#(df$result.1+dy)/(df$result-dx)##
  v2[,'errord1']=log(v2[,'errord1a'])
  v2[,'errord2']=log(v2[,'errord2a'])
  v2[,'result']=log(v2[,'result_pure'])
  v2[,'Control']=df$result
  v2[,'Case']=df$result.1 #
  #https://stackoverflow.com/questions/31518150/gsub-in-r-is-not-replacing-dot

  v2[,'pval0']=v2[,'pval']
  v2[,'pval1']=v2[,'pval']
  v2[,'Significance0']= v2[,'pval0']<0.1#
  v2[,'Significance0'][v2[,'Significance0']==TRUE]='Yes'
  v2[,'Significance0'][v2[,'Significance0']==FALSE]='No'
  v2[,'Color0']=v2[,'pval0'] < 0.1#(v2[,case1]  < 1 & v2[,case2] < 1) | (v2[,case1]  >= 1 & v2[,case2] >= 1)
  v2[,'Color0'][v2[,'Color0']==TRUE]='blue'
  v2[,'Color0'][v2[,'Color0']==FALSE]='grey'
  v2[,'Significance1']= v2[,'pval1']<0.1#(v2[,case1]  < 1 & v2[,case2] < 1) | (v2[,case1]  >= 1 & v2[,case2] >= 1)
  v2[,'Significance1'][v2[,'Significance1']==TRUE]='Yes'
  v2[,'Significance1'][v2[,'Significance1']==FALSE]='No'
  v2[,'Color1']=v2[,'pval1'] < 0.1#(v2[,case1]  < 1 & v2[,case2] < 1) | (v2[,case1]  >= 1 & v2[,case2] >= 1)
  v2[,'Color1'][v2[,'Color1']==TRUE]='blue'
  v2[,'Color1'][v2[,'Color1']==FALSE]='grey'

  gn=groups[,c('Group','Abbreviation')]
  gn=gn[gn[,'Abbreviation']!='F',]
  gn=gn[order(data.frame(gn[,'Abbreviation'])[,1]),]
  v2=v2[order(v2[,'name']),]
  v2=cbind(v2,gn[order(data.frame(gn[,'Abbreviation'])[,1]),])
  v2=v2[rev(order(v2[,'result'])),]

  xlab = "Autoscaled Concentrations (SE)" #xlab = "Raw Concentrations in Log10 Scale (SE)"}
  xlim=c(min(v2$errord1),max(v2$errord2)) #Occasionally: xlim=c(min(v2$result)*1.1,max(v2$result)*1.1) # if (xlim[2]>1) {xlim[2]=1};# if (xlim[1] < -0.75) {xlim[1]=-0.75};

  # Below perhaps wex, has been done already:
  plote2=forestplot(df = v2, #drive tba_example_v3_oh_tikka17823.R if not working via 'x' error #coef, #
                    estimate = result,
                    se=0,#abs(errord1-errord2)/4, #sterr,##this makes the significant value:
                    pvalue = pval1,psignif = 0.1,
                    xlim=xlim, xlab = 'Logged Ratio between Raw Concentrations of Case and Control with 90% CI',ylab='Steroid Groups',
                    title='',colour = Significance1 ) +#,colour = Significance
    ggforce::facet_col(facets = ~Group,scales = "free_y",space = "free", strip.position='left')+
  geom_errorbarh(aes(xmin = errord1, xmax = errord2,height = .0,colour=Significance1));#plote2

  if (sum(v2[,'Significance1']=='Yes')==20) {hp=c('blue','blue')} else {hp=c('#999999','blue')};#plote2
  if (Group=='All' & first==TRUE) {ordera=v2$name[order(v2$result)]; #
  plote2[["data"]][["name"]]=factor(plote2[["data"]][["name"]], levels = ordera)} else if
  (Group=='All' & first==FALSE) {plote2[["data"]][["name"]]=factor(plote2[["data"]][["name"]], levels = ordera)} else if
  (Group=='Female') {plote2[["data"]][["name"]]=factor(plote2[["data"]][["name"]], levels = ordera)} else if
  (Group=='Male') {plote2[["data"]][["name"]]=factor(plote2[["data"]][["name"]], levels = ordera)}
  #https://www.r-bloggers.com/2020/03/how-to-standardize-group-colors-in-data-visualizations-in-r/
  plote2$layers[[1]]$aes_params$odd <- "#00000000" #https://stackoverflow.com/questions/71745719/how-to-control-stripe-transparency-using-ggforestplot-geom-stripes
  v2$Group2=v2$Group
  v2 <- transform(v2,Group2 = as.numeric(as.factor(Group2)))
  v2$facet_fill_color <- c("red", "green", "blue", "yellow", "brown")[v2$Group2]
  jopon=plote2  +theme(axis.text.y=element_blank()) +theme_classic2();   #theme(axis.text.y = element_text(lineheight=.05));
  jopon2=jopon+geom_point(aes(colour = factor(Significance1)),colour = v2[,'Color1']) +
    scale_color_manual(values=hp)+theme(legend.position = "none")+theme(strip.text.y = element_text(size=-Inf)) #ggtext::element_markdown(size = 12)
  # https://stackoverflow.com/questions/10547487/remove-facet-wrap-labels-completely
  g <- ggplot_gtable(ggplot_build(jopon2))
  stripr <- which(grepl('strip-l', g$layout$name)); fills <- c("red","green","blue","yellow",'brown'); k <- 1;
  for (i in stripr) {j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]; k <- k+1}
  grid::grid.draw(g)
  # https://stackoverflow.com/questions/24169675/multiple-colors-in-a-facet-strip-background-in-ggplot
  # http://www.sthda.com/english/wiki/ggplot2-axis-scales-and-transformations
  # https://www.statology.org/geom_point-fill/
  # https://stackoverflow.com/questions/62093084/set-geom-vline-line-types-and-sizes-with-aes-mapping-in-ggplot2
  # https://en.wikipedia.org/wiki/Standard_error
  # This has been already done:
  jpeg(paste(name ,"divi.jpg"), width = 7500, height = 11000, quality = 100,pointsize = 16, res=1000); print(grid::grid.draw(g));dev.off();

  return(ordera) #If you do not want to have 'null' to the Rmarkdown/html take this away
}




#
# library(roxygen2); # Read in the roxygen2 R package
# gc()
# roxygenise();      # Builds the help files; read these separately not in a project I guess....
