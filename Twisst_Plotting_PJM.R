
library(ape)
library(parallel)
library(data.table)

#wrapper function for loess smoothing
simple.loess.predict <- function(x, y, span, new_x=NULL, weights = NULL, max = NULL, min = NULL, degree=NULL, family=NULL){
    y.loess <- loess(y ~ x, span = span, weights = weights, degree=degree, family=family)
    if (is.null(new_x)) {y.predict <- predict(y.loess,x)}
    else {y.predict <- predict(y.loess,new_x)}
    if (is.null(min) == FALSE) {y.predict = ifelse(y.predict > min, y.predict, min)}
    if (is.null(max) == FALSE) {y.predict = ifelse(y.predict < max, y.predict, max)}
    y.predict
    }

#function for smoothing multiple columns of a data frame at once
smooth_df <- function(x, df, span, new_x = NULL, col.names=NULL, weights=NULL, min=NULL, max=NULL, degree=NULL, family=NULL){
    if (is.null(new_x)) {smoothed <- df}
    else smoothed = df[1:length(new_x),]
    if (is.null(col.names)){col.names=colnames(df)}
    for (col.name in col.names){
        print(paste("smoothing",col.name))
        smoothed[,col.name] <- simple.loess.predict(x,df[,col.name],span = span, new_x = new_x, max = max, min = min, weights = weights, degree=degree, family=family)
        }
    smoothed
    }

#function to held stack polygons for plotting
stack <- function(mat){
    upper <- t(apply(mat, 1, cumsum))
    lower <- upper - mat
    list(upper=upper,lower=lower)
    }

#function for interleaving start and end positions (for plotting)
interleave <- function(x1,x2){
    output <- vector(length= length(x1) + length(x2))
    output[seq(1,length(output),2)] <- x1
    output[seq(2,length(output),2)] <- x2
    output
    }

#function for ploting weights as either stacked or overlayed polygons
plot_weights <- function(weights_dataframe,positions=NULL,line_cols=NULL,fill_cols=NULL,density=NULL,xlim=NULL,ylim=c(0,1),stacked=FALSE,
                                        ylab="Weights", xlab = "Position", main="",xaxt=NULL,yaxt=NULL,bty="n",lwd=1, add=FALSE){
    #get x axis
    x = positions
    #if a two-column matrix is given - plot step-like weights with start and end of each window    
    if (is.matrix(x)==TRUE) {
        x = interleave(positions[,1],positions[,2])
        yreps=2
        }
    else {
        if (is.null(x)==FALSE) x = positions
        else x = 1:nrow(weights_dataframe)
        yreps=1
        }
    
    #set x limits
    if(is.null(xlim)) xlim = c(min(x), max(x))
    
    #if not adding to an old plot, make a new plot
    if (add==FALSE) plot(0, pch = "", xlim = xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main,xaxt=xaxt,yaxt=yaxt,bty=bty)
    
    if (stacked == TRUE){
        y_stacked <- stack(weights_dataframe)
        for (n in 1:ncol(weights_dataframe)){
            y_upper = rep(y_stacked[["upper"]][,n],each=yreps)
            y_lower = rep(y_stacked[["lower"]][,n],each = yreps)
            polygon(c(x,rev(x)),c(y_upper, rev(y_lower)), col = fill_cols[n], density=density[n], border=NA)
            }
        }
    else{
        for (n in 1:ncol(weights_dataframe)){
            y = rep(weights_dataframe[,n],each=yreps)
            x_poly = c(x[1],x,x[length(x)])
            y_poly= c(0,y,0)
            polygon(x_poly, y_poly, col=fill_cols[n], border=NA,density=density[n])
            lines(x,y, type = "l", col = line_cols[n],lwd=lwd)
            }
        }
    
    }


options(scipen = 5)

########### input data ################


fsc2.best=tapply(weights_by_chrom_smooth$Lhood_est,fsc2$key,max)
fsc2.best=data.frame(key=names(fsc2.best),value=fsc2.best)
fsc2.1=merge(fsc2,fsc2.best[,c("key","value")],by=c("key"))
fsc3=cbind(fsc3,best=apply(fsc3[,c(2,3,4,5,6,7)],1,which.min))

#there are separate weights input files and window data files for each chromosome
chromNames = 1:8

weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst1.weights.tsv")
weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst2.weights.tsv")
weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst3.weights.tsv")
weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst3b.weights.tsv")
weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst3c.weights.tsv")
weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst3d.weights.tsv")
weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst3e.weights.tsv")
weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst4.weights.tsv")
weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst4b.weights.tsv")
weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst4c.weights.tsv")
weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst4d.weights.tsv")
weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst4e.weights.tsv")
weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst5.weights.tsv")

window_data_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.data.tsv")


#make list of weight and window data for each chromosme, throwing out windows with NA values
weights_by_chrom = list()
window_data_by_chrom = list()

for (i in 1:length(weights_files)){
    weights <- read.table(weights_files[i], header = T, as.is = T)
    weights <- weights/apply(weights,1,sum,na.rm=T)
    good_rows = which(is.na(weights[,1]) == FALSE)
    weights_by_chrom[[i]] <- weights[good_rows,]
    window_data_by_chrom[[i]] <- read.table(window_data_files[i], header = T, as.is = T)[good_rows,]
    }


#merge lists to get all weights and overall mean for each
weights_all=data.frame()
for(i in 1:8){
weights_all=rbind(weights_all,weights_by_chrom[[i]])}
mean_weights <- apply(weights_all, 2, mean)

##########################

#positions of centromeres. (APPROXIMATE, this was based on where SNPs seem to disappear).

centromeres = list(c(18443862, 21579352),
                   c(5043517,7259459),
                   c(15338332,17391288),
                   c(7158340,9655602),
                   c(7223722,9310246),
                   c(14632935,17138420),
                   c(15211666,17689663),
                   c(7520452,9782543)
                   )


###################################### colours ##############################################

#define colours and semi-transparent versions of each colour. And name colours according to topologies

#### cols for 3 topos

topoNames = names(weights)

cols = c(
"#0075DC", #Blue
"#FF0010", #Red
"#FFA405") #Orpiment (mustard)

trans_cols = paste0(cols, "25")

names(cols) = topoNames[c(1,2,3)]
names(trans_cols) = topoNames[c(1,2,3)]


#### cols for 15 topos

topoNames = names(weights)

cols = c(
"#2BCE48", #Green
"#005C31", #Forest
"#94FFB5", #Jade
"#9DCC00", #Lime
"#426600", #Quagmire (olive)
"#00998F", #Turquoise
"#5EF1F2", #Sky
"#0075DC", #Blue
"#003380", #Navy
"#740AFF", #Violet
"#FF5005", #Zinnia (orange)
"#F0A3FF", #Amethyst (pink)
"#C20088", #Mallow (red purple)
"#FF0010", #Red
"#FFA405") #Orpiment (mustard)

trans_cols = paste0(cols, "25")

names(cols) = topoNames[c(7,8,9,12,13,1,2,3,10,15,4,5,6,11,14)]
names(trans_cols) = topoNames[c(7,8,9,12,13,1,2,3,10,15,4,5,6,11,14)]

#### cols for 12 topos

topoNames = names(weights)

cols = c(
"#2BCE48", #Green
"#005C31", #Forest
"#94FFB5", #Jade
"#9DCC00", #Lime
"#426600", #Quagmire (olive)
"#00998F", #Turquoise
"#5EF1F2", #Sky
"#0075DC", #Blue
"#003380", #Navy
"#740AFF", #Violet
"#FF5005", #Zinnia (orange)
"#F0A3FF" #Amethyst (pink)
)

trans_cols = paste0(cols, "25")

names(cols) = topoNames[c(10,11,12,61,62,63,48,55,775,768,430,423)]
names(trans_cols) = topoNames[c(10,11,12,61,62,63,48,55,775,768,430,423)]


################################### smoothing #########################################

#define a set of positions for the smoothed data
#(reduces the number of data points without losing resolution, as long as spacing is not too big).
smoothed_positions <- lapply(chromNames, function(chromName) seq(min(window_data_by_chrom[[chromName]]$mid),max(window_data_by_chrom[[chromName]]$mid), 10000))

#get smoothed values for all weights for all chroms
weights_by_chrom_smooth <- mclapply(chromNames,function(chromName){smooth_df(x=window_data_by_chrom[[chromName]]$mid,
                                                                weights_by_chrom[[chromName]],
                                                                span=100000/max(window_data_by_chrom[[chromName]]$end),
                                                                new_x = smoothed_positions[[chromName]],
                                                                min = 0, max = 1, degree = 2, family = "gaus")}, mc.cores = 4)

names(weights_by_chrom_smooth) <- chromNames

#rescale values so they sum to 1
for (chromName in chromNames){
    weights_by_chrom_smooth[[chromName]][,topoNames] <- weights_by_chrom_smooth[[chromName]][,topoNames] / 
                                                    apply(weights_by_chrom_smooth[[chromName]][,topoNames], 1, sum)}


################################# plots across chromosomes ##############################
### SMOOTHED

png(file = "~/Documents/Research/300/Twisst/DM_HM_BP_BI.4dg.DP8MIN230MAC2.LDphase.phyml_bionj.Twisst4b.STACKED.weights.png", width = 3000, height = 3000, res=300)

plot_order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)

#plot_order = c(3,1,2)


par(mfrow=c(8,1), mar = c(3,3,1,1), cex.main = 1)

for (chromName in chromNames) {
    plot_weights(weights_dataframe=weights_by_chrom_smooth[[chromName]][,topoNames[plot_order]], positions=smoothed_positions[[chromName]]/1e6,
                 line_cols=NULL, fill_cols=cols[topoNames[plot_order]], stacked=TRUE, xlim = c(0, 32))
    rect(centromeres[[chromName]][1]/1e6,0,centromeres[[chromName]][2]/1e6,1, col="white", border="white")
    }

# for (chromName in chromNames) {
#     plot_weights(weights_dataframe=weights_by_chrom_smooth[[chromName]][,topoNames[plot_order]], positions=smoothed_positions[[chromName]]/1e6,
#                  line_cols=NULL, fill_cols=cols[topoNames[plot_order]], stacked=TRUE, xlim = c(0, 32),ylim = c(0,max(apply(weights_by_chrom_smooth[[chromName]][,topoNames[plot_order]],1,sum))))
#     rect(centromeres[[chromName]][1]/1e6,0,centromeres[[chromName]][2]/1e6,1, col="white", border="white")
#     }


axis(2,cex.axis = 1, line = -1.5, at = c(0,0.5,1))
mtext(side = 2, text = "Weighting", line = 1)

dev.off()

### SMOOTHED, overlayed

png(file = "~/Documents/Research/300/Twisst/DM_HM_BP_BI.4dg.DP8MIN230MAC2.LDphase.phyml_bionj.Twisst4b.OVLAY.weights.png", width = 3000, height = 3000, res=300)
# pdf(file = "twisst/arn_lyr_72.DP5MIN58MAC2.w100m40.LDphase.phyml_bionj.arn2_and4_lyr2_lyr4_halleri.weights2.pdf", width = 8, height = 10)
# plot_order = c(14,11)

# png(file = "twisst/arn_lyr_72.DP5MIN58MAC2.w100m40.LDphase.phyml_bionj.arn4_lyr4_OCH2_halleri.weights.png", width = 3000, height = 3000, res=300)
# png(file = "twisst/arn_lyr_72.DP5MIN58MAC2.w100m40.LDphase.phyml_bionj.arn4_lyr4_OCH6_halleri.weights.png", width = 3000, height = 3000, res=300)
# png(file = "twisst/arn_lyr_72.DP5MIN58MAC2.w100m40.LDphase.phyml_bionj.arn4_lyr4_ROK2_halleri.weights.png", width = 3000, height = 3000, res=300)
# png(file = "twisst/arn_lyr_72.DP5MIN58MAC2.w100m40.LDphase.phyml_bionj.arn4_lyr4_ROK4_halleri.weights.png", width = 3000, height = 3000, res=300)

plot_order <- order(mean_weights, decreasing=T)

par(mfrow=c(8,1), mar = c(4,4,1,1), cex.main = 1)

for (chromName in chromNames) {
    plot_weights(weights_dataframe=weights_by_chrom_smooth[[chromName]][,topoNames[plot_order]], positions=smoothed_positions[[chromName]]/1e6,
                 line_cols=cols[topoNames[plot_order]], fill_cols=trans_cols[topoNames[plot_order]], stacked=FALSE, xlim = c(0, 32), bty="l",
                 ylab = "", xlab = "", yaxt="n", lwd=0.5, ylim = c(0,max(apply(weights_by_chrom_smooth[[chromName]][,topoNames[plot_order]],1,max))))
    rect(centromeres[[chromName]][1]/1e6,0,centromeres[[chromName]][2]/1e6,1, col="white", border="white")
    axis(2, at = c(0,0.5,1), las = 2)
    mtext(2,text="Weighting", line = 2.5)
    mtext(3,text=chromName, line = -.1)
    }
mtext(1,text="Position (Mb)", line = 2.5)


dev.off()

png(file = "~/Documents/Research/300/Twisst/DM_HM_BP_BI.4dg.DP8MIN230MAC2.LDphase.phyml_bionj.Twisst3.OVLAY.weights.topo13n3.png", width = 3000, height = 3000, res=300)
plot_order <- order(mean_weights, decreasing=T)

par(mfrow=c(8,1), mar = c(4,4,1,1), cex.main = 1)

for (chromName in chromNames) {
    plot_weights(weights_dataframe=weights_by_chrom_smooth[[chromName]][,topoNames[c(3,13)]], positions=smoothed_positions[[chromName]]/1e6,
                 line_cols=cols[topoNames[c(3,13)]], fill_cols=trans_cols[topoNames[c(3,13)]], stacked=FALSE, xlim = c(0, 32), bty="l",
                 ylab = "", xlab = "", yaxt="n", lwd=0.5, ylim = c(0,max(apply(weights_by_chrom_smooth[[chromName]][,topoNames[c(3,13)]],1,max))))
    rect(centromeres[[chromName]][1]/1e6,0,centromeres[[chromName]][2]/1e6,1, col="white", border="white")
    axis(2, at = c(0,0.5,1), las = 2)
    mtext(2,text="Weighting", line = 2.5)
    mtext(3,text=chromName, line = -.1)
    }
mtext(1,text="Position (Mb)", line = 2.5)


dev.off()

############# RAW

png(file = "twisst/arn_lyr_72.DP5MIN58MAC2.w100m40.LDphase.phyml_bionj.arn2_and4_lyr2_lyr4_halleri.weights_raw.png", width = 3000, height = 3000, res=300)
plot_order <- order(mean_weights, decreasing=T)

png(file = "twisst/arn_lyr_72.DP5MIN58MAC2.w100m40.LDphase.phyml_bionj.arn2_and4_lyr2_lyr4_halleri.weights_raw3.png", width = 3000, height = 3000, res=300)
plot_order = c(11,14,6)

png(file = "twisst/arn_lyr_72.DP5MIN58MAC2.w100m40.LDphase.phyml_bionj.arn2_arn4p_lyr2_lyr4p_halleri.weights_raw3.png", width = 3000, height = 3000, res=300)
plot_order = c(11,14,6)

png(file = "twisst/arn_lyr_72.DP5MIN58MAC2.w100m40.LDphase.phyml_bionj.arn4_lyr4_OCH2_halleri.weights_raw2.png", width = 3000, height = 3000, res=300)
plot_order <- c(2,3)

png(file = "twisst/arn_lyr_72.DP5MIN58MAC2.w100m40.LDphase.phyml_bionj.arn4_lyr4_OCH6_halleri.weights_raw2.png", width = 3000, height = 3000, res=300)
plot_order <- c(2,3)


par(mfcol=c(8,1), mar = c(4,4,1,1), cex.main = 1)

for (chromName in chromNames) {
    plot(0, cex = 0, ylim = c(0,1), main = chromName, bty="l", xlim = c(0, 32), yaxt="n",xlab = "",ylab = "")
    for (topoName in topoNames[plot_order]){
        points(window_data_by_chrom[[chromName]]$mid/1000000, weights_by_chrom[[chromName]][,topoName],
                    col=cols[topoName], type = "h", lwd=1.2)
        }
    rect(centromeres[[chromName]][1]/1e6,0,centromeres[[chromName]][2]/1e6,1, col="white", border="white")
    axis(2, at = c(0,1), labels = c(0,1), las = 2)
    mtext(2,text="Weighting", line = 2.5)
    }
mtext(1,text="Position (Mb)", line = 2.5)

dev.off()



######### raw, selected regions

png(file = "~/Documents/Research/300/Twisst/DM_HM_BP_BI.4dg.DP8MIN230MAC2.LDphase.phyml_bionj.Twisst4e.weights_raw_selected.png", width = 3000, height = 3000, res=300)

plot_order = c(1,2,3)

regions <- matrix(c(1,0,0.5e6,220193,225616,
                    1,9.5e6,10.5e6,9779387,9791542,
                    2,17e6,18e6,17706052,17715900,
                    2,12e6,13e6,12427364,12431683,
                    4,11e6,12e6,11123922,11131646,
                    4,22.5e6,23.5e6,22845898,22850613,
                    6, 1.75e6,2.25e6,2001440,2005979
                    ), ncol=5, byrow=T)
reg_names=c("PRD3","ZYP1a/b","PDS5","ASY1","SMC3","ASY3","SYN1")



par(mfrow=c(7,2), mar = c(4,4,1,1), cex.main = 1)

for (x in 1:nrow(regions)) {
    plot_weights(weights_dataframe=weights_by_chrom[[regions[x,1]]][,topoNames[plot_order]],
                 positions=as.matrix(window_data_by_chrom[[regions[x,1]]][,c("start","end")]),
                 line_cols=cols[topoNames[plot_order]], fill_cols=trans_cols[topoNames[plot_order]], stacked=FALSE, xlim = regions[x,c(2,3)], bty="l",
                 ylab = "", xlab = "", yaxt="n", lwd=1)
    polygon(regions[x,c(4,5)],c(0,0))
    axis(2, at = c(0,0.5,1), las = 2)
    mtext(2,text="Weighting", line = 2.5)
#    mtext(3,text=paste(c("Chromosome", " ", regions[x,1], " ", regions[x,2]/1e6, "-", regions[x,3]/1e6), collapse=""), line = -.1)
    mtext(3,text=paste("Chromosome", regions[x,1],"--",reg_names[x]), line = -.1)
    }

dev.off()

plot_weights(weights_dataframe=weights_by_chrom[[5]][,topoNames[c(1,2,3)]],
             positions=as.matrix(window_data_by_chrom[[5]][,c("start","end")]),
             line_cols=cols[topoNames[c(1,2,3)]], fill_cols=trans_cols[topoNames[c(1,2,3)]], stacked=FALSE, xlim = c(15000000,16000000), bty="l",
             ylab = "", xlab = "", yaxt="n", lwd=1)

#Calculate the mean_weight of meiosis genes. Two genes span windows, PDS5 and SMC3.
#PDS5 is included because 94.8% of gene space is included in one window
#SMC3 is excluded because only 59.9% of gene is included in single window
m_winds=data.frame()
for (x in 1:nrow(regions)){
    df = window_data_by_chrom[[regions[x,1]]][,c("start","end")]
    print(head(df))
    print(regions[x,4])
    print(regions[x,5])
    row = rownames(df[df$start < regions[x,4] & df$end > regions[x,5],])
    print(row)
    m_winds = rbind(m_winds,weights_by_chrom[[regions[x,1]]][row,topoNames[plot_order]])
}
m_winds=rbind(m_winds,weights_by_chrom[[2]]["447",topoNames[plot_order]]) 

> colMeans(m_winds)
    topo1     topo2     topo3 
0.6500189 0.1846463 0.1653349 

genes=read.table("~/Documents/Research/300/AccessoryInfo/Gene_coords.txt",sep=",")
g_winds=data.frame()
for (x in 1:nrow(genes)){
    df = window_data_by_chrom[[genes[x,1]]][,c("start","end")]
    row = rownames(df[df$start < genes[x,2] & df$end > genes[x,3],])
    g_winds = rbind(g_winds,weights_by_chrom[[genes[x,1]]][row,topoNames[plot_order]])
}

############################################### boxplot ########################################

library(data.table)


weights_all <- rbindlist(weights_by_chrom)

mean_weights <- apply(weights_all, 2, mean)

plot_order <- order(mean_weights, decreasing=T)

png(file = "twisst/arn_lyr_72.DP5MIN58MAC2.w100m40.LDphase.phyml_bionj.arn2_and4_lyr2_lyr4_halleri.weights_box.png", width = 3000, height = 2000, res = 400)
par(mar = c(4,4,1,1), bty="n")
boxplot(as.data.frame(weights_all)[,plot_order], outline = F, col = cols[topoNames[plot_order]], ylab = "Weighting", las=2)
dev.off()





############################################# plot topologies ####################################


########## 15 topos

text = system(paste("cat", weights_files[[1]], "| head -n 15"), intern = T)
topos <- read.tree(text = text)
names(topos) = topoNames

#plot_order = c(3,10,15,1,2,7,8,9,12,13,11,14,6,4,5)
#plot_order = 1:15
plot_order <- order(mean_weights, decreasing=T)

png("~/Documents/Research/300/Twisst/Twisst3b.topos.png", width = 3000, height = 1500, res = 400)

par(mfrow = c(3,5), mar = c(1,1,2,1),xpd=NA)

for (topoName in topoNames[plot_order]){
  plot.phylo(topos[[topoName]], type = "cl", cex = 1,edge.color=cols[topoName], edge.width=5,label.offset=0.3, main = topoName)
  #text(1,1,paste("topo",n,sep=""))
  }

dev.off()


########## 3 topos

text = system(paste("cat", weights_files[[1]], "| head -n 3"), intern = T)
topos <- read.tree(text = text)
names(topos) = topoNames

plot_order = c(1,2,3)

png("~/Documents/Research/300/Twisst/Twisst3e.topos.png", width = 2000, height = 400, res = 400)

par(mfrow = c(1,3), mar = c(1,2,2,1),xpd=NA)

for (topoName in topoNames[plot_order]){
  plot.phylo(topos[[topoName]], type = "cl", cex = 1,edge.color=cols[topoName], edge.width=5,label.offset=0.5, adj=0.5,
             use.edge.length=F, rotate.tree=90)
  #text(1,1,paste("topo",n,sep=""))
  }

dev.off()


####### 12 topos from interspersed lines

lines=c(10,11,12,61,62,63,48,55,775,768,430,423)
text1=""
for(i in 1:(length(lines)-1)){
    text1=paste(text1,"NR==",lines[i],"||",sep="")
}
text1=paste(text1,"NR==",lines[length(lines)],"{print}NR>",max(lines),"{exit}",sep="")

text = system(paste("awk '", text1, "' ", weights_files[[1]],sep=""), intern = T)
topos <- read.tree(text = text)
names(topos) = names(cols)

png("~/Documents/Research/300/Twisst/Twisst1.topos.png", width = 200, height = 8000, res = 400)

par(mfrow = c(4,3),xpd=NA)

for (topoName in names(topos)){
  plot.phylo(topos[[topoName]], type = "cl", cex = 1,edge.color=cols[topoName], edge.width=5,label.offset=0.5, adj=0.5,
             use.edge.length=F, rotate.tree=90)
  #text(1,1,paste("topo",n,sep=""))
  }

dev.off()




################################## pull out outliers ###############################################


weights_all <- as.data.frame(rbindlist(weights_by_chrom))

window_data_all <- as.data.frame(rbindlist(window_data_by_chrom))


window_data_all[which(weights_all[,"topo6"] >= 0.5),][,1:3]

window_data_all[which(weights_all[,"topo11"] >= 0.5),][,1:3]

window_data_all[which(weights_all[,"topo14"] >= 0.5),][,1:3]


# 
# ################### plot some actual trees
# 
# 
# groups_table <- read.table("arn_lyr.big_groups.txt", as.is = T)
# groups <- groups_table[,2]
# samples <- groups_table[,1]
# 
# group_cols = c("red","orange","blue","skyblue","violet","gray80")
# names(group_cols) <- c("arenosa2", "arenosa4", "lyrata2", "lyrata4", "hybrid4", "halleri")
# 
# sample_cols <- group_cols[groups]
# names(sample_cols) <- samples
# 
# 
# 
# text = system("zcat twisst/arn_lyr_72.scf1.DP5MIN58MAC2.w100m40.LDphase.phyml_bionj.trees.gz | head -n 50", intern = T)
# trees <- read.tree(text = text)
# 
# 
# i = 10
# 
# plot.phylo(trees[[i]], type = "un", cex = 0.7, rotate.tree=90, no.margin=T, show.tip.label=T )
# tiplabels(pch = 19, col = sample_cols[trees[[i]]$tip.label], cex = 1, adj = 0.5)


chromNames = 1:8
window_data_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.data.tsv")
centromeres = list(c(18443862, 21579352),
                   c(5043517,7259459),
                   c(15338332,17391288),
                   c(7158340,9655602),
                   c(7223722,9310246),
                   c(14632935,17138420),
                   c(15211666,17689663),
                   c(7520452,9782543)
                   )

topoNames = names(weights)

cols = c(
"#0075DC", #Blue
"#FF0010", #Red
"#FFA405") #Orpiment (mustard)

trans_cols = paste0(cols, "25")

names(cols) = topoNames[c(1,2,3)]
names(trans_cols) = topoNames[c(1,2,3)]

regions <- matrix(c(1,0,0.5e6,220193,225616,
                    1,9.5e6,10.5e6,9779387,9791542,
                    2,17e6,18e6,17706052,17715900,
                    2,12e6,13e6,12427364,12431683,
                    4,11e6,12e6,11123922,11131646,
                    4,22.5e6,23.5e6,22845898,22850613,
                    6, 1.75e6,2.25e6,2001440,2005979
                    ), ncol=5, byrow=T)
reg_names=c("PRD3","ZYP1a/b","PDS5","ASY1","SMC3","ASY3","SYN1")
plot_order = c(1,2,3)

aa=c("VEL.DRA.VID","VEL.LAC.VID","VEL.TZI.VID","HOC.DRA.VID","HOC.LAC.VID","HOC.TZI.VID","SPI.DRA.VID","SPI.LAC.VID","SPI.TZI.VID","TKO.DRA.VID","TKO.LAC.VID","TKO.TZI.VID","TRE.DRA.VID","TRE.LAC.VID","TRE.TZI.VID")
aa=c("VEL.KOW.MIE","VEL.STE.MIE","VEL.TBG.MIE","HOC.KOW.MIE","HOC.STE.MIE","HOC.TBG.MIE","SPI.KOW.MIE","SPI.STE.MIE","SPI.TBG.MIE","TKO.KOW.MIE","TKO.STE.MIE","TKO.TBG.MIE","TRE.KOW.MIE","TRE.STE.MIE","TRE.TBG.MIE")

Mean_Weights_all=data.frame("Pops"=as.character(),"topo1"=as.numeric(),"topo2"=as.numeric(),"topo3"=as.numeric())
Mei_mean_weights=data.frame("Pops"=as.character(),"topo1"=as.numeric(),"topo2"=as.numeric(),"topo3"=as.numeric())
Weights_all=data.frame()
for (j in 1:length(aa)){
    weights_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.",aa[j],".weights.tsv")
    weights_by_chrom = list()
    window_data_by_chrom = list()

    for (i in 1:length(weights_files)){
        weights <- read.table(weights_files[i], header = T, as.is = T)
        weights <- weights/apply(weights,1,sum,na.rm=T)
        good_rows = which(is.na(weights[,1]) == FALSE)
        weights_by_chrom[[i]] <- weights[good_rows,]
        window_data_by_chrom[[i]] <- read.table(window_data_files[i], header = T, as.is = T)[good_rows,]
        }
    weights_all=data.frame()
    for(i in 1:8){

    weights_all=rbind(weights_all,weights_by_chrom[[i]])}

    mean_weights <- apply(weights_all, 2, mean)
    Mean_Weights_all=rbind(Mean_Weights_all,data.frame("Pops"=aa[j],"topo1"=mean_weights[1],"topo2"=mean_weights[2],"topo3"=mean_weights[3]))
    
    weights_all["Pops"]=aa[j]
    Weights_all=rbind(Weights_all,weights_all)

    m_winds=data.frame()
    for (x in 1:nrow(regions)){
        df = window_data_by_chrom[[regions[x,1]]][,c("start","end")]
        row = rownames(df[df$start < regions[x,4] & df$end > regions[x,5],])
        m_winds = rbind(m_winds,weights_by_chrom[[regions[x,1]]][row,topoNames[plot_order]])
    }
    m_winds=rbind(m_winds,weights_by_chrom[[2]]["447",topoNames[plot_order]]) 
    m_winds_mean=colMeans(m_winds)
    m_winds_mean["Pops"]=aa[j]
    Mei_mean_weights=rbind(Mei_mean_weights,data.frame("Pops"=aa[j],"topo1"=m_winds_mean[1],"topo2"=m_winds_mean[2],"topo3"=m_winds_mean[3]))
}

Weights_all.m=melt(Weights_all,id.var="Pops")
Mei_mean_weights.m=melt(Mei_mean_weights,id.var="Pops")
Mei_mean_weights.m$value=as.numeric(Mei_mean_weights.m$value)
ggplot(Weights_all.m,aes(x=variable, y=value,fill=Pops))+geom_boxplot(outlier.shape = NA)+geom_point(data=Mei_mean_weights.m,aes(x=variable,y=value),position = position_dodge(width = -.75))+scale_x_discrete(limits=c("topo1","topo3","topo2"),labels=c("Topology 1","Topology 3","Topology 2"))+ylab("Weight")+theme(axis.title.x=element_blank())


## get outliers ##
chromNames = 1:8
weights_files3 = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst4c.weights.tsv")
weights_files4 = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.Twisst3c.weights.tsv")
window_data_files = paste0("~/Documents/Research/300/Twisst/CurrentResults/DM_HM_BP_BI.4dg.scf", chromNames,".DP8MIN230MAC2.LDphase.phyml_bionj.data.tsv")


#make list of weight and window data for each chromosme, throwing out windows with NA values
weights_by_chrom3 = list()
window_data_by_chrom3 = list()

for (i in 1:length(weights_files3)){
    weights <- read.table(weights_files3[i], header = T, as.is = T)
    weights <- weights/apply(weights,1,sum,na.rm=T)
    good_rows = which(is.na(weights[,1]) == FALSE)
    weights_by_chrom3[[i]] <- weights[good_rows,]
    window_data_by_chrom3[[i]] <- read.table(window_data_files[i], header = T, as.is = T)[good_rows,]
    }

weights_all3=data.frame()
for(i in 1:8){
weights_all3=rbind(weights_all3,weights_by_chrom3[[i]])}
mean_weights3 <- apply(weights_all3, 2, mean)

weights_by_chrom4 = list()
window_data_by_chrom4 = list()

for (i in 1:length(weights_files4)){
    weights <- read.table(weights_files4[i], header = T, as.is = T)
    weights <- weights/apply(weights,1,sum,na.rm=T)
    good_rows = which(is.na(weights[,1]) == FALSE)
    weights_by_chrom4[[i]] <- weights[good_rows,]
    window_data_by_chrom4[[i]] <- read.table(window_data_files[i], header = T, as.is = T)[good_rows,]
    }

weights_all4=data.frame()
for(i in 1:8){
weights_all4=rbind(weights_all4,weights_by_chrom4[[i]])}
mean_weights4 <- apply(weights_all4, 2, mean)

twisst3=data.frame()
for (i in 1:8){
    df=cbind(window_data_by_chrom3[[i]],weights_by_chrom3[[i]])
    twisst3=rbind(twisst3,df)
}

twisst4=data.frame()
for (i in 1:8){
    df=cbind(window_data_by_chrom4[[i]],weights_by_chrom4[[i]])
    twisst4=rbind(twisst4,df)
}

twisst3e.topo1=twisst3[twisst3$topo1 > quantile(twisst3$topo1,0.99,na.rm=T),]
twisst3e.topo2=twisst3[twisst3$topo2 > quantile(twisst3$topo2,0.99,na.rm=T),]
twisst3e.topo3=twisst3[twisst3$topo3 > quantile(twisst3$topo3,0.99,na.rm=T),]

twisst4e.topo1=twisst4[twisst4$topo1 > quantile(twisst4$topo1,0.99,na.rm=T),]
twisst4e.topo2=twisst4[twisst4$topo2 > quantile(twisst4$topo2,0.99,na.rm=T),]
twisst4e.topo3=twisst4[twisst4$topo3 > quantile(twisst4$topo3,0.99,na.rm=T),]

c1=intersect(twisst3e.topo1$id,twisst4e.topo1$id)
c3=intersect(twisst3e.topo3$id,twisst4e.topo3$id)

wpAll.5["L"]=wpAll.5$start+1000
wpAll.5["U"]=wpAll.5$end-1000

aa=list(twisst3e.topo1,twisst3e.topo2,twisst3e.topo3)
outs.twisst3=data.frame()
for (i in 1:3){
    xx=data.frame("H"=as.numeric(),"E"=as.numeric(),"D"=as.numeric(),"topo"=as.numeric())
    for (j in 1:nrow(aa[[i]])){
        pp=wpAll.5[wpAll.5$scaff==aa[[i]][j,1] & ((wpAll.5$L <= aa[[i]][j,2] & wpAll.5$U >= aa[[i]][j,2]) | (wpAll.5$L <= aa[[i]][j,3] & wpAll.5$U >= aa[[i]][j,3]) | (wpAll.5$L >= aa[[i]][j,2] & wpAll.5$U <= aa[[i]][j,3])) & wpAll.5$pop == "DRA",]
        xx=rbind(xx,data.frame("H"=mean(pp[pp$num_snps>20,]$H),"E"=mean(pp[pp$num_snps>20,]$E),"D"=mean(pp[pp$num_snps>20,]$D),"topo"=i))
    }
    outs.twisst3=rbind(outs.twisst3,xx)
}

aa=list(twisst4e.topo1,twisst4e.topo2,twisst4e.topo3)
outs.twisst4=data.frame()
for (i in 1:3){
    xx=data.frame()
    for (j in 1:nrow(aa[[i]])){
        pp=wpAll.5[wpAll.5$scaff==aa[[i]][j,1] & ((wpAll.5$L <= aa[[i]][j,2] & wpAll.5$U >= aa[[i]][j,2]) | (wpAll.5$L <= aa[[i]][j,3] & wpAll.5$U >= aa[[i]][j,3]) | (wpAll.5$L >= aa[[i]][j,2] & wpAll.5$U <= aa[[i]][j,3])) & wpAll.5$pop == "KOW",]
        xx=rbind(xx,data.frame("H"=mean(pp[pp$num_snps>20,]$H),"E"=mean(pp[pp$num_snps>20,]$E),"D"=mean(pp[pp$num_snps>20,]$D),"topo"=i))
    }
    outs.twisst4=rbind(outs.twisst4,xx)
}


NEED TO USE QUANTILE OF H WITHIN DRA OR TSC NOT OVERALL.

HT.dra.overlap=data.frame()
    for (j in 1:nrow(twisst3e.topo3)){
        pp=wpAll.5[wpAll.5$scaff==twisst3e.topo3[j,1] & ((wpAll.5$L <= twisst3e.topo3[j,2] & wpAll.5$U >= twisst3e.topo3[j,2]) | (wpAll.5$L <= twisst3e.topo3[j,3] & wpAll.5$U >= twisst3e.topo3[j,3]) | (wpAll.5$L >= twisst3e.topo3[j,2] & wpAll.5$U <= twisst3e.topo3[j,3])) & wpAll.5$pop == "DRA" & wpAll.5$H < quantile(wpAll.5[wpAll.5$pop == "DRA",]$H,0.01,na.rm=T),]
        HT.dra.overlap=rbind(HToverlap,pp)
    }

HT.kow.overlap=data.frame()
    for (j in 1:nrow(twisst4e.topo3)){
        pp=wpAll.5[wpAll.5$scaff==twisst4e.topo3[j,1] & ((wpAll.5$L <= twisst3e.topo3[j,2] & wpAll.5$U >= twisst4e.topo3[j,2]) | (wpAll.5$L <= twisst4e.topo3[j,3] & wpAll.5$U >= twisst4e.topo3[j,3]) | (wpAll.5$L >= twisst4e.topo3[j,2] & wpAll.5$U <= twisst4e.topo3[j,3])) & wpAll.5$pop == "KOW" & wpAll.5$H < quantile(wpAll.5[wpAll.5$pop == "KOW",]$H,0.01,na.rm=T),]
        HT.kow.overlap=rbind(HToverlap,pp)
    }

H.dra.outs=wpAll.5[wpAll.5$pop == "DRA" & wpAll.5$H < quantile(wpAll.5$H,0.01,na.rm=T) & wpAll.5$num_snps>20,]
H.kow.outs=wpAll.5[wpAll.5$pop == "KOW" & wpAll.5$H < quantile(wpAll.5$H,0.01,na.rm=T) & wpAll.5$num_snps>20,]

WpAll["L"]=WpAll$end-9000
WpAll["U"]=WpAll$end-1000

HT.dra.overlap=data.frame()
    for (j in 1:nrow(twisst3e.topo3)){
        pp=WpAll[WpAll$scaff==twisst3e.topo3[j,1] & ((WpAll$L <= twisst3e.topo3[j,2] & WpAll$U >= twisst3e.topo3[j,2]) | (WpAll$L <= twisst3e.topo3[j,3] & WpAll$U >= twisst3e.topo3[j,3]) | (WpAll$L >= twisst3e.topo3[j,2] & WpAll$U <= twisst3e.topo3[j,3])) & WpAll$lin1 == "TSC" & WpAll$H < quantile(WpAll[WpAll$lin1 == "TSC",]$H,0.01,na.rm=T),]
        HT.dra.overlap=rbind(HT.dra.overlap,pp)
    }

HT.kow.overlap=data.frame()
    for (j in 1:nrow(twisst4e.topo3)){
        pp=WpAll[WpAll$scaff==twisst4e.topo3[j,1] & ((WpAll$L <= twisst4e.topo3[j,2] & WpAll$U >= twisst4e.topo3[j,2]) | (WpAll$L <= twisst4e.topo3[j,3] & WpAll$U >= twisst4e.topo3[j,3]) | (WpAll$L >= twisst4e.topo3[j,2] & WpAll$U <= twisst4e.topo3[j,3])) & WpAll$lin1 == "TRU" & WpAll$H < quantile(WpAll[WpAll$lin1 == "TRU",]$H,0.01,na.rm=T),]
        HT.kow.overlap=rbind(HT.kow.overlap,pp)
    }

#find average topo weights of H outliers
outs.H=data.frame()
for (j in 1:nrow(H.dra.outs)){
    pp=twisst3e.topo3[twisst3e.topo3$scaffold==H.dra.outs[j,1] & (!(twisst3e.topo3$start <= H.dra.outs[j,"L"] & twisst3e.topo3$end <= H.dra.outs[j,"L"]) | !(twisst3e.topo3$end >= H.dra.outs[j,"U"] & twisst3e.topo3$start >= H.dra.outs[j,"U"])),]
    outs.H=rbind(outs.H,pp)
}

wpAll.5['topoOut']=0
aa=list(twisst3e.topo1,twisst3e.topo2,twisst3e.topo3)
outs.twisst3=data.frame()
for (i in 1:3){
    for (j in 1:nrow(aa[[i]])){
        wpAll.5[wpAll.5$scaff==aa[[i]][j,1] & ((wpAll.5$L <= aa[[i]][j,2] & wpAll.5$U >= aa[[i]][j,2]) | (wpAll.5$L <= aa[[i]][j,3] & wpAll.5$U >= aa[[i]][j,3]) | (wpAll.5$L >= aa[[i]][j,2] & wpAll.5$U <= aa[[i]][j,3])) & wpAll.5$pop == "DRA",]$topoOut=i

        }}

write.table(twisst3e.topo1[,c("scaffold","start","end")],"~/Documents/Research/300/Twisst/Twisst3c_topo1_0.99outliers.txt",sep="\t",row.names = F,col.names = F,quote=F)

#need to find command used to get .annot files...must've used bedtools...should be in history after logging onto cluster.

bedtools intersect -a /usr/users/JIC_c1/monnahap/FileTransfer/Twisst4c_topo1_0.99outliers.txt -b /nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/Ly2.0/annotation/LyV2.gff -f 0.1 -wo | grep transcript | grep -v transcription | sort -u | awk '{print $1,$2,$3,$7,$8,$9,$10,$12}' | tr ' ' '\t' > /usr/users/JIC_c1/monnahap/FileTransfer/Twisst4c_topo1_0.99outliers_genes.gff
bedtools intersect -a /usr/users/JIC_c1/monnahap/FileTransfer/Twisst4c_topo3_0.99outliers.txt -b /nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/Ly2.0/annotation/LyV2.gff -f 0.1 -wo | grep transcript | grep -v transcription | sort -u | awk '{print $1,$2,$3,$7,$8,$9,$10,$12}' | tr ' ' '\t' > /usr/users/JIC_c1/monnahap/FileTransfer/Twisst4c_topo3_0.99outliers_genes.gff
bedtools intersect -a /usr/users/JIC_c1/monnahap/FileTransfer/Twisst3c_topo1_0.99outliers.txt -b /nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/Ly2.0/annotation/LyV2.gff -f 0.1 -wo | grep transcript | grep -v transcription | sort -u | awk '{print $1,$2,$3,$7,$8,$9,$10,$12}' | tr ' ' '\t' > /usr/users/JIC_c1/monnahap/FileTransfer/Twisst3c_topo1_0.99outliers_genes.gff
bedtools intersect -a /usr/users/JIC_c1/monnahap/FileTransfer/Twisst3c_topo3_0.99outliers.txt -b /nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/Ly2.0/annotation/LyV2.gff -f 0.1 -wo | grep transcript | grep -v transcription | sort -u | awk '{print $1,$2,$3,$7,$8,$9,$10,$12}' | tr ' ' '\t' > /usr/users/JIC_c1/monnahap/FileTransfer/Twisst3c_topo3_0.99outliers_genes.gff

awk '{print substr($8,4,9)}' ~/Documents/Research/300/Twisst/Twisst3c_topo3_0.99outliers_genes.gff | grep -f /dev/stdin ~/Documents/Research/300/AccessoryInfo/LyV2_TAIR10orth_des_20150927.txt | awk '{print $2}' > ~/Documents/Research/300/Twisst/Twisst3c_topo3_0.99outliers_ThalIDs.txt
awk '{print substr($8,4,9)}' ~/Documents/Research/300/Twisst/Twisst3c_topo1_0.99outliers_genes.gff | grep -f /dev/stdin ~/Documents/Research/300/AccessoryInfo/LyV2_TAIR10orth_des_20150927.txt | awk '{print $2}' > ~/Documents/Research/300/Twisst/Twisst3c_topo1_0.99outliers_ThalIDs.txt
awk '{print substr($8,4,9)}' ~/Documents/Research/300/Twisst/Twisst4c_topo3_0.99outliers_genes.gff | grep -f /dev/stdin ~/Documents/Research/300/AccessoryInfo/LyV2_TAIR10orth_des_20150927.txt | awk '{print $2}' > ~/Documents/Research/300/Twisst/Twisst4c_topo3_0.99outliers_ThalIDs.txt
awk '{print substr($8,4,9)}' ~/Documents/Research/300/Twisst/Twisst4c_topo1_0.99outliers_genes.gff | grep -f /dev/stdin ~/Documents/Research/300/AccessoryInfo/LyV2_TAIR10orth_des_20150927.txt | awk '{print $2}' > ~/Documents/Research/300/Twisst/Twisst4c_topo1_0.99outliers_ThalIDs.txt

for file in ./*annot.txt; do awk '{print $2}' "$file" > "$file"_ThalIDs.txt; done

twisst4['data']='4c'
twisst3['data']='3c'
twisst4['source']='All'
twisst3['source']='All'

topoNames = names(weights)
regions <- matrix(c("scaffold_1",0,0.5e6,220193,225616,
                    "scaffold_1",9.5e6,10.5e6,9779387,9791542,
                    "scaffold_2",17e6,18e6,17706052,17715900,
                    "scaffold_2",12e6,13e6,12427364,12431683,
                    "scaffold_4",11e6,12e6,11123922,11131646,
                    "scaffold_4",22.5e6,23.5e6,22845898,22850613,
                    "scaffold_6", 1.75e6,2.25e6,2001440,2005979
                    ), ncol=5, byrow=T)
reg_names=c("PRD3","ZYP1a/b","PDS5","ASY1","SMC3","ASY3","SYN1")
plot_order = c(1,2,3)

twisst3$start=as.numeric(twisst3$start)
twisst3$end=as.numeric(twisst3$end)
twisst4$start=as.numeric(twisst4$start)
twisst4$end=as.numeric(twisst4$end)

m_winds3=data.frame()
for (x in 1:nrow(regions)){
    df = twisst3[twisst3$scaffold==regions[x,1] & twisst3$start < as.numeric(regions[x,4]) & twisst3$end > as.numeric(regions[x,5]),]
    m_winds3 = rbind(m_winds3,df)
}
m_winds3=rbind(m_winds3,twisst3[twisst3$scaffold=="scaffold_2" & twisst3$topo1 == weights_by_chrom3[[2]]["447",1],]) 

m_winds4=data.frame()
for (x in 1:nrow(regions)){
    df = twisst4[twisst4$scaffold==regions[x,1] & twisst4$start < as.numeric(regions[x,4]) & twisst4$end > as.numeric(regions[x,5]),]
    m_winds4 = rbind(m_winds4,df)
}
m_winds4=rbind(m_winds4,twisst4[twisst4$scaffold=="scaffold_2" & twisst4$topo1 == weights_by_chrom4[[2]]["447",1],]) 

m_winds4['data']='4c'
m_winds3['data']='3c'
m_winds4['source']='Mei'
m_winds3['source']='Mei'

twisst3=rbind(twisst3,m_winds3)
twisst4=rbind(twisst4,m_winds4)
twisst=rbind(twisst3,twisst4)
twisst[twisst$data=="4c",]$data="Baltic"
twisst[twisst$data=="3c",]$data="S. Carpathians"

m.twisst=melt(twisst[,c("source","data","topo1","topo2","topo3")],id.vars=c("source","data"))
ggplot(m.twisst,aes(x=variable,y=value,fill=source))+geom_boxplot(outlier.size=0.5)+theme_bw()+ylab("Weights")+facet_grid(~data)+theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1,size=14),legend.title=element_blank(),legend.text=element_text(size=14),axis.title.y=element_text(size=14),strip.text.x=element_text(size=14),axis.text.y=element_text(size=12))


