################################################################################
# This script was used for the PC analyses of Ferreira et al. (2021)
# 

getwd()

# load raw data
#raw.data = read.csv("rawdata.csv", header = T, sep = ",")

# load dataset of log-transformed linear measurements
morpho.data = read.csv("morphospace-measurements.csv", header = T, row.names = 1, sep = ",")


# scatterplots of each variable
dev.new(width = 4, height = 6)

pdf("results/variable_plot.pdf", width = 6, height = 6, family = "Times")
tiff("results/variable_plot.tif", family = "serif", width = 2000, height = 1800, 
      res = 300)

plot(x = morpho.data$SH, y = rep(1, length(morpho.data$SH)), ylim = c(0.5, 7.5), 
     xlim = c(0, 1.3), type = "n", xlab = "log-transformed length", ylab = "",
     yaxt = "n")

axis(2, labels = c("MW","MiL", "TSL", "TSW", "TSML", "SH", "HC"), 
     at = c(1, 2, 3, 4, 5, 6, 7))

# add MW
# all points in gray
points(x = morpho.data$MW, y = rep(1, length(morpho.data$MW)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.3))
# Peltocephalus dumerilianus
points(x = morpho.data$MW[which(morpho.data$species == "pd")], 
       y = rep(1, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = "#bebebe", bg = "#CC6677", cex = 1.5)
# Stupendemys
points(x = morpho.data$MW[which(morpho.data$species == "sg")], 
       y = rep(1, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 25, col = "#bebebe", bg = "#fdae6b", cex = 2)
# Rondonia
points(x = morpho.data$MW[which(morpho.data$species == "it")], 
       y = rep(1, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 25, col = "#bebebe", bg = "#AA4499", cex = 2)

# add MiL
points(x = morpho.data$MiL, y = rep(2, length(morpho.data$MiL)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.3))
points(x = morpho.data$MiL[which(morpho.data$species == "pd")], 
       y = rep(2, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = "#bebebe", bg = "#CC6677", cex = 1.5)
# Stupendemys
points(x = morpho.data$MiL[which(morpho.data$species == "sg")], 
       y = rep(2, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 25, col = "#bebebe", bg = "#fdae6b", cex = 2)
# Rondonia
points(x = morpho.data$MiL[which(morpho.data$species == "it")], 
       y = rep(2, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 25, col = "#bebebe", bg = "#AA4499", cex = 2)

# add TSL
points(x = morpho.data$TSL, y = rep(3, length(morpho.data$TSL)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.3))
points(x = morpho.data$TSL[which(morpho.data$species == "pd")], 
       y = rep(3, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = "#bebebe", bg = "#CC6677", cex = 1.5)
# Stupendemys
points(x = morpho.data$TSL[which(morpho.data$species == "sg")], 
       y = rep(3, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 25, col = "#bebebe", bg = "#fdae6b", cex = 2)

# Rondonia
points(x = morpho.data$TSL[which(morpho.data$species == "it")], 
       y = rep(3, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 25, col = "#bebebe", bg = "#AA4499", cex = 2)

# add TSW
points(x = morpho.data$TSW, y = rep(4, length(morpho.data$TSW)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.3))
points(x = morpho.data$TSW[which(morpho.data$species == "pd")], 
       y = rep(4, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = "#bebebe", bg = "#CC6677", cex = 1.5)
# Stupendemys
points(x = morpho.data$TSW[which(morpho.data$species == "sg")], 
       y = rep(4, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 25, col = "#bebebe", bg = "#fdae6b", cex = 2)
# Rondonia
points(x = morpho.data$TSW[which(morpho.data$species == "it")], 
       y = rep(4, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 25, col = "#bebebe", bg = "#AA4499", cex = 2)

# add TSML
points(x = morpho.data$TSML, y = rep(5, length(morpho.data$TSML)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.3))
points(x = morpho.data$TSML[which(morpho.data$species == "pd")], 
       y = rep(5, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = "#bebebe", bg = "#CC6677", cex = 1.5)
# Stupendemys
points(x = morpho.data$TSML[which(morpho.data$species == "sg")], 
       y = rep(5, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 25, col = "#bebebe", bg = "#fdae6b", cex = 2)
# Rondonia
points(x = morpho.data$TSML[which(morpho.data$species == "it")], 
       y = rep(5, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 25, col = "#bebebe", bg = "#AA4499", cex = 2)

# add SH
points(x = morpho.data$SH, y = rep(6, length(morpho.data$SH)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.3))
points(x = morpho.data$SH[which(morpho.data$species == "pd")], 
       y = rep(6, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = "#bebebe", bg = "#CC6677", cex = 1.5)
# Stupendemys
points(x = morpho.data$SH[which(morpho.data$species == "sg")], 
       y = rep(6, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 25, col = "#bebebe", bg = "#fdae6b", cex = 2)
# Rondonia
points(x = morpho.data$SH[which(morpho.data$species == "it")], 
       y = rep(6, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 25, col = "#bebebe", bg = "#AA4499", cex = 2)

# add HC
points(x = morpho.data$HC, y = rep(7, length(morpho.data$HC)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.3))
points(x = morpho.data$HC[which(morpho.data$species == "pd")], 
       y = rep(7, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = "#bebebe", bg = "#CC6677", cex = 1.5)
# Stupendemys
points(x = morpho.data$HC[which(morpho.data$species == "sg")], 
       y = rep(7, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 25, col = "#bebebe", bg = "#fdae6b", cex = 2)
# Rondonia
points(x = morpho.data$HC[which(morpho.data$species == "it")], 
       y = rep(7, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 25, col = "#bebebe", bg = "#AA4499", cex = 2)

legend(x = "topright", legend = c("Peltocephalus maturin", "Peltocephalus dumerilianus", 
                                  "Stupendemys geographica"), 
       cex = 0.7, bty = "n", pch = c(25, 21, 25), pt.cex = 1, col = "#bebebe", 
       pt.bg = c("#AA4499","#CC6677", "#fdae6b"), text.font = 3)



dev.off()

################################################################################
# Principal Component Analysis
pca.res = prcomp(morpho.data[, 3:10], center = T, scale = T)
summary(pca.res)

var = pca.res$sdev^2/sum(pca.res$sdev^2)
barplot(var[1:8], main = "Explained Variance", axisnames = T, ylim = c(0, 0.5))
var = as.matrix(var)
rownames(var) = colnames(pca.res$rotation)

cumvar = cumsum(pca.res$sdev^2/sum(pca.res$sdev^2))
plot(cumvar[0:8], xlab = "PCs", ylab = "Explained variance", 
     main = "Cumulative variance plot")

# check loadings
biplot(pca.res)

# export to pdf
pdf("results/morphospace_biplot.pdf", height = 6, width = 8, family = "Times")
biplot(pca.res)
dev.off()

###############################
# Morphospace plot

# create hulls
tab = matrix(c(pca.res$x[, 1], pca.res$x[, 2]), ncol = 2)

tab.em = tab[(which(morpho.data$species == "em")), ]
tab.pd = tab[(which(morpho.data$species == "pd")), ]
tab.pe = tab[(which(morpho.data$species == "pe")), ]
tab.pl = tab[(which(morpho.data$species == "pl")), ]
tab.pu = tab[(which(morpho.data$species == "pu")), ]
tab.ps = tab[(which(morpho.data$species == "ps")), ]
tab.pv = tab[(which(morpho.data$species == "pv")), ]
tab.px = tab[(which(morpho.data$species == "px")), ]
tab.pl = tab[(which(morpho.data$genus == "pl")), ]

ch.em = chull(tab.em)
ch.pd = chull(tab.pd)
ch.pe = chull(tab.pe)
ch.pl = chull(tab.pl)
ch.pu = chull(tab.pu)
ch.ps = chull(tab.ps)
ch.pv = chull(tab.pv)
ch.px = chull(tab.px)
ch.pl = chull(tab.pl)

# plot PC1 x PC2 without points
pdf("results/morphospace_jaws.pdf", height = 6, width = 8, family = "Times")

plot(pca.res$x[, 1], pca.res$x[, 2], cex = 0, ylim = c(-3.7, 3.5), 
     xlim = c(-3.5, 4.2), xlab = "PC1 (39.88%)", ylab = "PC2 (22.26%)")

#plot hulls
polygon(tab.px[ch.px, ], border = NA, col = rgb(51/255, 34/255, 136/255, 0.3), 
        lwd = 1)
polygon(tab.em[ch.em, ], border = NA, col = rgb(136/255, 204/255, 238/255, 0.3), 
        lwd = 1)
polygon(tab.pd[ch.pd, ], border = NA, col = rgb(204/255, 102/255, 119/255, 0.3),
        lwd = 1)
polygon(tab.pu[ch.pu, ], border = NA, col = rgb(68/255, 170/255, 153/255, 0.3), 
        lwd = 1)
polygon(tab.ps[ch.ps, ], border = NA, col = rgb(68/255, 170/255, 153/255, 0.3), 
        lwd = 1)
polygon(tab.pv[ch.pv, ], border = NA, col = rgb(68/255, 170/255, 153/255, 0.3), 
        lwd = 1)


# plot points
points(pca.res$x[, 2][which(morpho.data$species == "pd")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "pd")], pch = 21, col = "#bebebe", 
       bg = "#CC6677", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "em")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "em")], pch = 21, col = "#bebebe", 
       bg = "#88CCEE", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "pe")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "pe")], pch = 21, col = "#bebebe", 
       bg = "#332288", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "pl")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "pl")], pch = 23, col = "#bebebe", 
       bg = "#332288", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "pu")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "pu")], pch = 21, col = "#bebebe", 
       bg = "#44AA99", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "ps")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "ps")], pch = 22, col = "#bebebe", 
       bg = "#44AA99", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "pv")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "pv")], pch = 23, col = "#bebebe", 
       bg = "#44AA99", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "px")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "px")], pch = 22, col = "#bebebe", 
       bg = "#332288", cex = 1.5)

points(pca.res$x[, 2][which(morpho.data$species == "it")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "it")], pch = 25, col = "#bebebe", 
       bg = "#AA4499", cex = 2.5)
points(pca.res$x[, 2][which(morpho.data$species == "sg")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "sg")], pch = 25, col = "#bebebe", 
       bg = "#fdae6b", cex = 2)

points(pca.res$x[, 2][which(morpho.data$species == "ca")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "ca")], pch = 22, col = "#bebebe", 
       bg = "#DDCC77", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "si")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "si")], pch = 22, col = "#bebebe", 
       bg = "#DDCC77", cex = 1.5)


legend(x = "topright", legend = c("Peltocephalus maturin", "Stupendemys geographica",
        "Peltocephalus dumerilianus", "Erymnochelys madagascariensis", "Pelusios spp.",
        "Podocnemis sextuberculata", "Podocnemis unifilis", "Podocnemis vogli",
        "Podocnemis erythrocephala", "Podocnemis expansa", "Podocnemis lewiana"), 
       cex = 0.7, bty = "n", pch = c(25, 25, 21, 21, 22, 22, 21, 23, 21, 22, 23), 
       pt.cex = 1, col = "#bebebe", pt.bg = c("#AA4499",
        "#fdae6b", "#CC6677", "#88CCEE", "#DDCC77", "#44AA99", "#44AA99", 
        "#44AA99", "#332288", "#332288", "#332288"), ncol = 2, text.font = 3)



dev.off()


################################################################################
# Regression analyses for size estimation

size.raw = read.table("size-meas.csv", header = T, sep = ",")
attach(size.raw)
head(size.raw)
dim(size.raw)

################################################################################
## predicting JL from MiL in Podocnemididae

size.raw$MiL
size.raw[,4]

# log-transform all metrics

for (i in 4:7) {
  size.raw[, i] = log(size.raw[, i])
}

size.raw  # check if everything is log-transformed

################################################################################
## predicting JL from MiL in Podocnemididae

temp = NULL
for (i in 1:length(size.raw$Clade)) {
  if(size.raw$Clade[i] == "Podocnemididae"){
    x = size.raw[i,]
    temp = rbind(temp, x)
    
  }
  else(
    next
  )
  
}

length(temp$Clade)

temp = temp[,-(6:7)]
den = na.omit(temp)
den
cor.test(den$MiL, den$JL)

# correlation between MiL & JL = 0.9949353

length(den$Specimen)
rownames(den) = 1:39

attach(den)
cor.test(MiL, JL)
den.model = lm(JL ~ MiL)
summary(den.model)  # R-squared = 0.9896

par(mfrow = c(2,2))
plot(den.model)


new.turtle = data.frame(MiL = log(266))
pred.JL = predict(object = den.model, newdata = new.turtle, interval = "confidence")
pred.JL
# prediction:       
#       fit      lwr      upr
#   5.749245  5.689218  5.809271 

lwr.JL = data.frame(pred.JL[2])
upr.JL = data.frame(pred.JL[3])
fit.JL = data.frame(pred.JL[1])
colnames(lwr.JL) = colnames(upr.JL) = colnames(fit.JL) = "JL"
lwr.JL
upr.JL
fit.JL
exp(lwr.JL)  # 295.6623
exp(upr.JL)  # 333.3761
exp(fit.JL)  # 313.9534

################################################################################

################################################################################
# predicting SCL from JL
attach(size.raw)
jaw = data.frame(Clade, Species, JL, SCL)
jaw = na.omit(jaw)
jaw
length(jaw$JL)
cor.test(jaw$JL, jaw$SCL)

# correlation JL ~ SCL = 0.8860217 

attach(jaw)
cor.test(JL, SCL)
jaw.model = lm(SCL ~ JL)
summary(jaw.model)  # R-squared = 0.784

plot(jaw.model)

jaw.model


# lower prediction
pred.SCL = predict(object = jaw.model, newdata = lwr.JL, interval = "confidence")
lwr.SCL1 = data.frame(pred.SCL[2])
colnames(lwr.SCL1) = "SCL"

# upper prediction
pred.SCL = predict(object = jaw.model, newdata = upr.JL, interval = "confidence")
upr.SCL1 = data.frame(pred.SCL[3])
colnames(upr.SCL1) = "SCL"

# mean prediction
pred.SCL = predict(object = jaw.model, newdata = fit.JL, interval = "confidence")
fit.SCL1 = data.frame(pred.SCL[1])
colnames(fit.SCL1) = "SCL"

lwr.SCL1
fit.SCL1
upr.SCL1
exp(lwr.SCL1)  # 1407.85
exp(upr.SCL1)  # 2078.15
exp(fit.SCL1)  # 1704.204


################################################################################

################################################################################
# predicting SCm from JL
jaw.skull = data.frame(size.raw$Clade, size.raw$Species, size.raw$JL, size.raw$SCm)
jaw.skull = na.omit(jaw.skull)
colnames(jaw.skull) = c("Clade", "Species", "JL", "SCm")
head(jaw.skull)
length(jaw.skull$Species)

attach(jaw.skull)
length(SCm)
cor.test(SCm, JL)
# correlation SCm ~ JL = 0.9916569

jaw.skull.model = lm(SCm ~ JL)
summary(jaw.skull.model)  # R-squared = 0.9833
jaw.skull.model

# lower prediction
pred.SCm = predict(object = jaw.skull.model, newdata = lwr.JL, interval = "confidence")
lwr.SCm = data.frame(pred.SCm[2])
colnames(lwr.SCm) = "SCm"

# upper prediction
pred.SCm = predict(object = jaw.skull.model, newdata = upr.JL, interval = "confidence")
upr.SCm = data.frame(pred.SCm[3])
colnames(upr.SCm) = "SCm"

# mean prediction
pred.SCm = predict(object = jaw.skull.model, newdata = fit.JL, interval = "confidence")
fit.SCm = data.frame(pred.SCm[1])
colnames(fit.SCm) = "SCm"

lwr.SCm
upr.SCm
fit.SCm
exp(lwr.SCm)  # 328.8021
exp(upr.SCm)  # 394.5754
exp(fit.SCm)  # 359.8795


################################################################################
# predicting SCL from SCm

attach(size.raw)
skull = data.frame(Clade, Species, SCm, SCL)
skull = na.omit(skull)

length(skull$SCm)
cor.test(skull$SCm, skull$SCL)
# correlation SCm ~ SCL = 0.88702014

length(skull$SCm)

attach(skull)

skull.model = lm(SCL ~ SCm)
summary(skull.model)  # R-squared = 0.786
skull.model


par(mfrow = c(2,2))
plot(skull.model)


# lower prediction
pred.SCL2 = predict(object = skull.model, newdata = lwr.SCm, interval = "confidence")
lwr.SCL2 = data.frame(pred.SCL2[2])
colnames(lwr.SCL2) = "SCL"

# upper prediction
pred.SCL2 = predict(object = skull.model, newdata = upr.SCm, interval = "confidence")
upr.SCL2 = data.frame(pred.SCL2[3])
colnames(upr.SCL2) = "SCL"

# mean prediction
pred.SCL2 = predict(object = skull.model, newdata = fit.SCm, interval = "confidence")
fit.SCL2 = data.frame(pred.SCL2[1])
colnames(fit.SCL2) = "SCL"

lwr.SCL2
upr.SCL2
fit.SCL2
exp(lwr.SCL2)  # 1407.40
exp(upr.SCL2)  # 2169.768
exp(fit.SCL2)  # 1736.961


################################################################################
# sensitivity analyses
SCm.expansa = data.frame(SCm = log(85.6))
pred.SCL2.expansa = predict(object = skull.model, newdata = SCm.expansa, interval = "confidence")
exp(pred.SCL2.expansa)  # SCL = 470.3765, upr = 496.4139

SCm.unifilis = data.frame(SCm = log(57.5))
pred.SCL2.unifilis = predict(object = skull.model, newdata = SCm.unifilis, interval = "confidence")
exp(pred.SCL2.unifilis)  # SCL = 327.5289, upr = 339.6576

SCm.geographicus = data.frame(SCm = log(190))
pred.SCL2.geographicus = predict(object = skull.model, newdata = SCm.geographicus, interval = "confidence")
exp(pred.SCL2.geographicus)  # SCL = 971.5051, upr = 1070.406

SCm.dumerilianus = data.frame(SCm = log(81.6))
pred.SCL2.dumerilianus = predict(object = skull.model, newdata = SCm.dumerilianus, interval = "confidence")
exp(pred.SCL2.dumerilianus)  # SCL = 450.3389, upr 474.1474

par(mfrow = c(1,1))
pdf("results/size_regression.pdf", height = 6, width = 8, family = "Times")
plot(size.raw$SCL ~ size.raw$SCm, col = rgb(255/255, 255/255, 255/255, 0), 
     bg = rgb(252/255, 205/255, 229/255, 0), ylim = c(4.2, 8), xlim = c(2.7, 6),
     xlab = "log(skull length)", ylab = "log(carapace length)")
#abline(skull.model, col = "#bebebe", lwd = 2)

# confidence interval
conf.SCm.SCL = seq(2, 7, length.out = 251)
preds = predict(skull.model, newdata = data.frame(SCm = conf.SCm.SCL), interval = "confidence", level = 0.95)

#matlines(conf.SCm.CML, preds[,2:3], col = "blue")

# fill in area between regression line and confidence interval
polygon(c(rev(conf.SCm.SCL), conf.SCm.SCL), c(rev(preds[,3]), preds[,2]), 
        col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)

# plot background points
# plot background points
points(y = size.raw$SCL[which(!size.raw$Clade == "Podocnemididae")], 
       x = size.raw$SCm[which(!size.raw$Clade == "Podocnemididae")], pch = 21,
       col = rgb(225/255, 225/255, 225/255, 0.5), bg = rgb(225/255, 225/255, 225/255, 0.5), cex = 1.5)


points(y = size.raw$SCL[which(size.raw$Clade == "Podocnemididae" & 
                                !size.raw$Species == "Peltocephalus_dumerilianus" & !size.raw$Species == "Stupendemys_geographicus")], 
       x = size.raw$SCm[which(size.raw$Clade == "Podocnemididae" & 
                                !size.raw$Species == "Peltocephalus_dumerilianus" & !size.raw$Species == "Stupendemys_geographicus")],
       pch = 21, col = "#bebebe", bg = "#44AA99", cex = 2)

points(y = size.raw$SCL[which(size.raw$Species == "Peltocephalus_dumerilianus")], 
       x = size.raw$SCm[which(size.raw$Species == "Peltocephalus_dumerilianus")], pch = 21,
       col = "#bebebe", bg = "#CC6677", cex = 2)

points(y = size.raw$SCL[which(size.raw$Species == "Stupendemys_geographicus")], 
       x = size.raw$SCm[which(size.raw$Species == "Stupendemys_geographicus")], pch = 25,
       col = "#bebebe", bg = "#AA4499", cex = 2.5)

points(y = fit.SCL2, x = fit.SCm, pch = 25,
       col = "#bebebe", bg = "#fdae6b", cex = 2.5)


dev.off()

################################################################################
# Building Figure 2

pdf("results/figure2.pdf", height = 4.5, width = 8, family = "Times")
tiff("results/figure2.tif", family = "serif", width = 3000, height = 1687, 
     res = 300)

layout(matrix(c(1,2), ncol = 2, byrow = TRUE), widths = c(7,5), heights = c(2,2))
par(mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1.5,0.4,0))
plot(pca.res$x[, 1], pca.res$x[, 2], cex = 0, cex.axis = 0.8, ylim = c(-3.7, 4), 
     xlim = c(-3.5, 4.2), xlab = "PC1 (39.88%)", ylab = "PC2 (22.26%)")
text(x = -3.4, y = 4, labels = "(a)", xpd = NA)

#plot hulls
polygon(tab.px[ch.px, ], border = NA, col = rgb(51/255, 34/255, 136/255, 0.3), 
        lwd = 1)
polygon(tab.em[ch.em, ], border = NA, col = rgb(136/255, 204/255, 238/255, 0.3), 
        lwd = 1)
polygon(tab.pd[ch.pd, ], border = NA, col = rgb(204/255, 102/255, 119/255, 0.3),
        lwd = 1)
polygon(tab.pu[ch.pu, ], border = NA, col = rgb(68/255, 170/255, 153/255, 0.3), 
        lwd = 1)
polygon(tab.ps[ch.ps, ], border = NA, col = rgb(68/255, 170/255, 153/255, 0.3), 
        lwd = 1)
polygon(tab.pv[ch.pv, ], border = NA, col = rgb(68/255, 170/255, 153/255, 0.3), 
        lwd = 1)


# plot points
points(pca.res$x[, 2][which(morpho.data$species == "pd")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "pd")], pch = 21, col = "#bebebe", 
       bg = "#CC6677", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "em")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "em")], pch = 21, col = "#bebebe", 
       bg = "#88CCEE", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "pe")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "pe")], pch = 21, col = "#bebebe", 
       bg = "#332288", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "pl")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "pl")], pch = 23, col = "#bebebe", 
       bg = "#332288", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "pu")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "pu")], pch = 21, col = "#bebebe", 
       bg = "#44AA99", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "ps")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "ps")], pch = 22, col = "#bebebe", 
       bg = "#44AA99", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "pv")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "pv")], pch = 23, col = "#bebebe", 
       bg = "#44AA99", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "px")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "px")], pch = 22, col = "#bebebe", 
       bg = "#332288", cex = 1.5)

points(pca.res$x[, 2][which(morpho.data$species == "it")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "it")], pch = 25, col = "#bebebe", 
       bg = "#AA4499", cex = 2)
points(pca.res$x[, 2][which(morpho.data$species == "sg")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "sg")], pch = 25, col = "#bebebe", 
       bg = "#fdae6b", cex = 1.5)

points(pca.res$x[, 2][which(morpho.data$species == "ca")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "ca")], pch = 22, col = "#bebebe", 
       bg = "#DDCC77", cex = 1.5)
points(pca.res$x[, 2][which(morpho.data$species == "si")] ~ pca.res$x[, 1]
       [which(morpho.data$species == "si")], pch = 22, col = "#bebebe", 
       bg = "#DDCC77", cex = 1.5)


legend(x = "topright", legend = c("Peltocephalus maturin", "Stupendemys geographica",
                                  "Peltocephalus dumerilianus", "Erymnochelys madagascariensis", "Pelusios spp.",
                                  "Podocnemis sextuberculata", "Podocnemis unifilis", "Podocnemis vogli",
                                  "Podocnemis erythrocephala", "Podocnemis expansa", "Podocnemis lewiana"), 
       cex = 0.6, bty = "n", pch = c(25, 25, 21, 21, 22, 22, 21, 23, 21, 22, 23), 
       pt.cex = 0.8, col = "#bebebe", pt.bg = c("#AA4499",
                                              "#fdae6b", "#CC6677", "#88CCEE", "#DDCC77", "#44AA99", "#44AA99", 
                                              "#44AA99", "#332288", "#332288", "#332288"), ncol = 2, text.font = 3)

#screen(2)
par(mai = c(0.5, 0.5, 0.2, 0.2), mgp = c(1.5,0.4,0))
plot(size.raw$SCL ~ size.raw$SCm, col = rgb(255/255, 255/255, 255/255, 0), 
     bg = rgb(252/255, 205/255, 229/255, 0), ylim = c(4.2, 8), xlim = c(2.7, 6),
     xlab = "log(skull length)", ylab = "log(carapace length)")
text(x = 2.8, y = 8, labels = "(b)", xpd = NA)

# confidence interval
conf.SCm.SCL = seq(2, 7, length.out = 251)
preds = predict(skull.model, newdata = data.frame(SCm = conf.SCm.SCL), interval = "confidence", level = 0.95)

#matlines(conf.SCm.CML, preds[,2:3], col = "blue")

# fill in area between regression line and confidence interval
polygon(c(rev(conf.SCm.SCL), conf.SCm.SCL), c(rev(preds[,3]), preds[,2]), 
        col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)

# plot background points
# plot background points
points(y = size.raw$SCL[which(!size.raw$Clade == "Podocnemididae")], 
       x = size.raw$SCm[which(!size.raw$Clade == "Podocnemididae")], pch = 21,
       col = rgb(225/255, 225/255, 225/255, 0.5), bg = rgb(225/255, 225/255, 225/255, 0.5), cex = 1)


points(y = size.raw$SCL[which(size.raw$Clade == "Podocnemididae" & 
                                !size.raw$Species == "Peltocephalus_dumerilianus" & !size.raw$Species == "Stupendemys_geographicus")], 
       x = size.raw$SCm[which(size.raw$Clade == "Podocnemididae" & 
                                !size.raw$Species == "Peltocephalus_dumerilianus" & !size.raw$Species == "Stupendemys_geographicus")],
       pch = 21, col = "#bebebe", bg = "#44AA99", cex = 1.5)

points(y = size.raw$SCL[which(size.raw$Species == "Peltocephalus_dumerilianus")], 
       x = size.raw$SCm[which(size.raw$Species == "Peltocephalus_dumerilianus")], pch = 21,
       col = "#bebebe", bg = "#CC6677", cex = 1.5)

points(y = size.raw$SCL[which(size.raw$Species == "Stupendemys_geographicus")], 
       x = size.raw$SCm[which(size.raw$Species == "Stupendemys_geographicus")], pch = 25,
       col = "#bebebe", bg = "#AA4499", cex = 1.5)

points(y = fit.SCL2, x = fit.SCm, pch = 25,
       col = "#bebebe", bg = "#fdae6b", cex = 2)



dev.off()



save.image("giantpodc.RData")
