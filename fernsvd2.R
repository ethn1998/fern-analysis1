
#Import fern data

env_data <- read.csv("Env-Data-Elaborated-SORTED.csv", header=TRUE, fileEncoding="UTF-8-BOM") #Recorded environmental data
colnames(env_data) <- make.names(colnames(env_data)) #Clean up column names
#Bug in read.csv making 4 extra columns. We only need first 18 columns
env_data <- env_data[,c(1:18)]

env_data <- subset(env_data, select=-c(N.tot.individuals,N.species)) 
#Number of individuals and number of species is already encoded in fern numbers anyways so these are somewhat redundant

bianca <- read.csv("Bianca-plot-June-cleaned-SORTED.csv",header=TRUE, fileEncoding="UTF-8-BOM") #Observed fern numbers
colnames(bianca) <- make.names(colnames(bianca)) #Clean up column names
class(bianca$FOREST) = "character" #Index different types of forests.

#get numeric columns only
num.env_data <- env_data[sapply(env_data,is.numeric)] 
num.bianca <- bianca[sapply(bianca,is.numeric)]

#get names for convenience
env_vars <- colnames(num.env_data)
fern_names <- colnames(num.bianca)

# Step 1: Normalize and calculate cross covariance matrix

normalize <- function(x) {
  if(is.numeric(x)){ 
    scale(x, center=TRUE, scale=TRUE) #make sure center=TRUE to calculate deviation values, set scale=TRUE for unit variance.
  } else x
}

mat.bianca <- t(apply(num.bianca,1,function(x)(x/sum(x)))) #Matrix of fern proportions
mat.env_data <- apply(num.env_data,2,normalize) #Matrix of normalized environmental data

cov.fern <- cov(x=mat.env_data,y=mat.bianca) #covariance matrix between normalized environmental data and fern proportions


# Step 2: Calculate SVD

svd.fern <- svd(cov.fern)

#Export SVD matrices into .csv for reference
sv_ids <- 1:ncol(svd.fern$v)
colnames(svd.fern$v) <- paste(rep_len("fernsv",ncol(svd.fern$v)),sv_ids,sep="") #Must set row and column names prior to write.csv
rownames(svd.fern$v) <- fern_names
write.csv(t(svd.fern$v), file="fern-sv.csv",row.names=TRUE)
colnames(svd.fern$u) <- paste(rep_len("envsv",ncol(svd.fern$u)),sv_ids,sep="") #Must set row and column names prior to write.csv
rownames(svd.fern$u) <- env_vars
write.csv(t(svd.fern$u), file="env-sv.csv",row.names=TRUE)

df.sv <- data.frame(c(1:length(svd.fern$d)),svd.fern$d)
colnames(df.sv) <- c("index","sv")
write.csv(df.sv, file="sv.csv",row.names=FALSE)

#percents <- c("0%","20%","40%","60%","80%","100%")

mar.default <- c(5,4,4,2)+0.1 #Default margin parameter
par(mar=mar.default)

plot(cumsum(svd.fern$d^2)/sum(svd.fern$d^2),main="Cumulative percentage of total covariance",xlab="Index",ylab="",ylim=c(0,1),xlim=c(0,length((svd.fern$d))),yaxt="n",las=1)
axis(side=2,at=seq(0.0,1.0,0.2),labels=c("0%","20%","40%","60%","80%","100%"),las=1)
abline(h=0.9,col="red")
abline(h=0.8,col="red")
#First 2 singular vectors describe almost 80% of the variation


#u: columns are left singular vectors (environment)
#v: columns are right singular vectors (fern composition)

#First 4 singular vectors describe 90% of variation.

#Horizontal barplot for environmental information

mar.default <- c(5,4,4,2)+0.1 #Default margin parameter
par(mar=mar.default)
bp.env <- barplot(height=svd.fern$u[,1],names.arg=env_vars,xlim=c(-1,1),las=1,main="1st environmental axis",horiz=TRUE,yaxt="n",xlab="Coefficient",ylab="Environmental factor")
axis(side=2,at=bp.env,labels=env_vars,pos=-0.5,las=1) #Axis on left
axis(side=4,at=bp.env,labels=env_vars,pos=0.5,las=1) #Axis on right

bp.env <- barplot(height=svd.fern$u[,2],names.arg=env_vars,xlim=c(-1,1),las=1,main="2nd environmental axis",horiz=TRUE,yaxt="n",xlab="Coefficient",ylab="Environmental factor")
axis(side=2,at=bp.env,labels=env_vars,pos=-0.5,las=1) #Axis on left
axis(side=4,at=bp.env,labels=env_vars,pos=0.5,las=1) #Axis on right


barplot(height=svd.fern$u[,3],names.arg=env_vars,xlim=c(-1,1),las=1,main="3rd environmental axis",horiz=TRUE)

#Pick out most important components of left (environment) singular vectors
sort.u1 <- sort(svd.fern$u[,1]^2,decreasing=TRUE,index.return=TRUE)
sort.u1$ivars <- env_vars[sort.u1$ix]


mar.default <- c(5,4,4,2)+0.1 #Default margin parameter


par(mar = c(5,0,0,0)+mar.default)
plot(sort.u1$x,ylab="sq coeff",main="1st environmental axis",xaxt="n",xlab="")
#Should I put percentages for cumsum of square coeffs? 
axis(1, at = c(1:length(sort.u1$x)),labels = sort.u1$ivars,las=2)
plot(cumsum(sort.u1$x),ylim=c(0,1),ylab="cumulative sum of squares",main="1st environmental axis",xaxt="n",xlab="")
axis(1, at = c(1:length(sort.u1$x)),labels = sort.u1$ivars,las=2)
abline(h=0.9,col="red")
abline(h=0.8,col="purple")
#1: GWC, OM, P, Silt, N


sort.u2 <- sort(svd.fern$u[,2]^2,decreasing=TRUE,index.return=TRUE)
sort.u2$ivars <- env_vars[sort.u2$ix]
plot(sort.u2$x,ylab="sq coeff",main="2nd environmental axis",xaxt="n",xlab="")
axis(1, at = c(1:length(sort.u2$x)),labels = sort.u2$ivars,las=2)
plot(cumsum(sort.u2$x),ylim=c(0,1),ylab="cumulative sum of squares",main="2nd environmental axis",xaxt="n",xlab="")
axis(1, at = c(1:length(sort.u2$x)),labels = sort.u2$ivars,las=2)
abline(h=0.9,col="red")
abline(h=0.8,col="purple")
#2: Clay, Sand, canopy, pH, Silt, Distance from water
###Silt important in both

#barplot(height=svd.fern$v[,1],names.arg=fern_names,ylim=c(-1,1),las=2,main="1st fern composition axis")
#abline(h=0)
###Fern species names are probably not a good idea, instead index because they are ordered alphabetically following bianca anyways

#Pick out important components of right (fern composition) singular vectors

sort.v1 <- sort(svd.fern$v[,1]^2,decreasing=TRUE,index.return=TRUE)
sort.v1$ivars <- fern_names[sort.v1$ix]

par(mar=mar.default)
plot(sort.v1$x[1:15],ylim=c(0,1),ylab="Sq coeff.",main="1st fern axis",xaxt="n",xlab="")
axis(1, at = c(1:15),labels = sort.v1$ivars[1:15],las=2)

#plot(cumsum(sort.v1$x[1:15]),ylim=c(0,1),ylab="cumulative sum of squares",main="1st fern axis",xaxt="n",xlab="") #Do top 15 only for simplicity
#axis(1, at = c(1:15),labels = sort.v1$ivars[1:15],las=2) #Don't do this unless we use abbreviated fern species names. Also need to clean up trailing numbers from species names.

plot(cumsum(sort.v1$x[1:15]),ylim=c(0,1),ylab="cumulative sum of squares",main="1st fern axis",xaxt="n",xlab="Index",las=2) #Do top 15 only for simplicity
axis(1,at=c(1:15),labels=c(1:15))
abline(h=0.9,col="red")
abline(h=0.8,col="purple")

topindex.v1 <- sort(sort.v1$ix[1:15]) #Top 15 most important ferns on 1st axis
top.v1 <- svd.fern$v[topindex.v1,1]

par(mar=c(0,0,0,10)+mar.default)
bp <- barplot(height=top.v1,xlim=c(-1,1),las=1,main="1st fern composition axis",horiz=TRUE,yaxt="n",xlab="Coefficient",ylab="Top 15 Fern Species")
axis(side=4,at=bp,labels=names(top.v1),pos=0.5,las=1)


#top.v1 <- data.frame(coeff=svd.fern$v[topindex.v1],)

sort.v2 <- sort(svd.fern$v[,2]^2,decreasing=TRUE,index.return=TRUE)
sort.v2$ivars <- fern_names[sort.v2$ix]
par(mar=mar.default)
plot(sort.v2$x[1:15],ylim=c(0,1),ylab="Sq coeff.",main="2nd fern axis",xaxt="n",xlab="")
axis(1, at = c(1:15),labels = sort.v2$ivars[1:15],las=2)
plot(cumsum(sort.v2$x[1:15]),ylim=c(0,1),ylab="cumulative sum of squares",main="2nd fern axis",xaxt="n",xlab="",las=2) #Do top 15 only for simplicity
#axis(1, at = c(1:15),labels = sort.v2$ivars[1:15],las=2)
axis(1, at = c(1:15),labels = c(1:15))

topindex.v2 <- sort(sort.v2$ix[1:15]) #Top 15 most important ferns on 1st axis
top.v2 <- svd.fern$v[topindex.v2,1]

par(mar=c(0,0,0,10)+mar.default)
bp <- barplot(height=top.v1,xlim=c(-1,1),las=1,main="2nd fern composition axis",horiz=TRUE,yaxt="n",xlab="Coefficient",ylab="Top 15 Fern Species")
axis(side=4,at=bp,labels=names(top.v2),pos=0.5,las=1)



abline(h=0.9,col="red")
abline(h=0.8,col="purple")

###Schizaea dicotoma, Selliguea heterocarpa, Asplenium longissimum, Asplenium phyllitidis, Drynaria Quercifolia appear in both 1st and 2nd axes.


#Maybe multiple plots per principal axis?! I can't squeeze all 89 species into one plot.

#Step 3: Projection plots
ftoc <- function(ftype){ #Color based on forest types
  switch(ftype, "K"="magenta", "PW"="navyblue", "MDF"="limegreen")
}
pcolors <- sapply(env_data$forest.type,ftoc)
ftos <- function(ftype){ #Symbol based on forest types  
  switch(ftype, "K"=1, "PW"=2, "MDF"=3)
}
psymbs <- sapply(env_data$forest.type,ftos)


projplot <- function(j, xmat, ymat) { #data projection plot along jth principal axis of svd, specialized for this dataset
  cov.mats <- cov(x=xmat,y=ymat)
  svd.mats <- svd(cov.mats)
  yproj <- c()
  xproj <- c()
  for(i in 1:length(mat.bianca[,1])){
    yproj <- append(yproj,sum(ymat[i,]*svd.mats$v[,j])) #Right projections on vertical axis
    xproj <- append(xproj,sum(xmat[i,]*svd.mats$u[,j])) #Left projections on horizontal axis
  }
  plottitle <- paste("Projection along axis:",j,sep="")
  plot(xproj,yproj,xlab="Env",ylab="Fern composition",main=plottitle,col=pcolors,pch=psymbs) 
  #text(xproj,yproj,labels=env_data$plot,col=pcolors)#LABEL PLOT IDS, THIS CAN MAKE PLOT A BIT MESSY
  legend("bottomright",legend=c("HF","PSF","MDF"),col=c("magenta","navyblue","limegreen"),pch=psymbs,title="Forest Types")
  abline(h=0)
  abline(v=0)
} #BUG: Plot IDs are not labelled FIXED

projplot(1, xmat=mat.env_data, ymat=mat.bianca)
projplot(2, xmat=mat.env_data, ymat=mat.bianca)

