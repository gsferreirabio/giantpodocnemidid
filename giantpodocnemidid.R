################################################################################
# This script was used for the PC analyses of Ferreira et al. (2021)
# 


getwd()

# load raw data
raw.data = read.csv("rawdata.csv", header = T, sep = ",")

# load unitless linear measurements dataset
morpho.data = read.csv("morphospace.csv", header = T, row.names = 1, sep = ",")

# log-transform data
log.morpho.data = morpho.data
for (i  in 3:10){
  log.morpho.data[, i] = log(log.morpho.data[, i])
}

# scatterplots of each variable
dev.new(width = 4, height = 6)

pdf("results/variable_plot.pdf", width = 6, height = 6)
plot(x = morpho.data$SH, y = rep(1, length(morpho.data$SH)), ylim = c(0.5, 7.5), 
     xlim = c(0, 1.3), type = "n", xlab = "log-transformed length", ylab = "",
     yaxt = "n")

axis(2, labels = c("MW","MiL", "TSL", "TSW", "TSML", "SH", "HC"), 
     at = c(1, 2, 3, 4, 5, 6, 7))

# add MW
# all points in gray
points(x = morpho.data$MW, y = rep(1, length(morpho.data$MW)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.5))
# Peltocephalus dumerilianus
points(x = morpho.data$MW[which(morpho.data$species == "pd")], 
       y = rep(1, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = rgb(255/255, 87/255, 51/255, 0.5), bg = rgb(255/255, 87/255, 
       51/255, 0.5))
# Stupendemys
points(x = morpho.data$MW[which(morpho.data$species == "sg")], 
       y = rep(1, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 23, col = rgb(184/255, 16/255, 81/255), bg = rgb(184/255, 16/255, 
       81/255), cex = 1.5)
# Rondonia
points(x = morpho.data$MW[which(morpho.data$species == "it")], 
       y = rep(1, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 22, col = rgb(88/255, 24/255, 69/255), bg = rgb(88/255, 24/255, 
       69/255), cex = 1.5)

# add MiL
points(x = morpho.data$MiL, y = rep(2, length(morpho.data$MiL)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.5))
points(x = morpho.data$MiL[which(morpho.data$species == "pd")], 
       y = rep(2, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = rgb(255/255, 87/255, 51/255, 0.5), bg = rgb(255/255, 87/255, 
                                                                   51/255, 0.5))
# Stupendemys
points(x = morpho.data$MiL[which(morpho.data$species == "sg")], 
       y = rep(2, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 23, col = rgb(184/255, 16/255, 81/255), bg = rgb(184/255, 16/255, 
                                                              81/255), cex = 1.5)
# Rondonia
points(x = morpho.data$MiL[which(morpho.data$species == "it")], 
       y = rep(2, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 22, col = rgb(88/255, 24/255, 69/255), bg = rgb(88/255, 24/255, 
                                                             69/255), cex = 1.5)

# add TSL
points(x = morpho.data$TSL, y = rep(3, length(morpho.data$TSL)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.5))
points(x = morpho.data$TSL[which(morpho.data$species == "pd")], 
       y = rep(3, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = rgb(255/255, 87/255, 51/255, 0.5), bg = rgb(255/255, 87/255, 
                                                                   51/255, 0.5))
# Stupendemys
points(x = morpho.data$TSL[which(morpho.data$species == "sg")], 
       y = rep(3, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 23, col = rgb(184/255, 16/255, 81/255), bg = rgb(184/255, 16/255, 
                                                              81/255), cex = 1.5)
# Rondonia
points(x = morpho.data$TSL[which(morpho.data$species == "it")], 
       y = rep(3, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 22, col = rgb(88/255, 24/255, 69/255), bg = rgb(88/255, 24/255, 
                                                             69/255), cex = 1.5)

# add TSW
points(x = morpho.data$TSW, y = rep(4, length(morpho.data$TSW)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.5))
points(x = morpho.data$TSW[which(morpho.data$species == "pd")], 
       y = rep(4, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = rgb(255/255, 87/255, 51/255, 0.5), bg = rgb(255/255, 87/255, 
                                                                   51/255, 0.5))
# Stupendemys
points(x = morpho.data$TSW[which(morpho.data$species == "sg")], 
       y = rep(4, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 23, col = rgb(184/255, 16/255, 81/255), bg = rgb(184/255, 16/255, 
                                                              81/255), cex = 1.5)
# Rondonia
points(x = morpho.data$TSW[which(morpho.data$species == "it")], 
       y = rep(4, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 22, col = rgb(88/255, 24/255, 69/255), bg = rgb(88/255, 24/255, 
                                                             69/255), cex = 1.5)

# add TSML
points(x = morpho.data$TSML, y = rep(5, length(morpho.data$TSML)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.5))
points(x = morpho.data$TSML[which(morpho.data$species == "pd")], 
       y = rep(5, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = rgb(255/255, 87/255, 51/255, 0.5), bg = rgb(255/255, 87/255, 
                                                                   51/255, 0.5))
# Stupendemys
points(x = morpho.data$TSML[which(morpho.data$species == "sg")], 
       y = rep(5, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 23, col = rgb(184/255, 16/255, 81/255), bg = rgb(184/255, 16/255, 
                                                              81/255), cex = 1.5)
# Rondonia
points(x = morpho.data$TSML[which(morpho.data$species == "it")], 
       y = rep(5, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 22, col = rgb(88/255, 24/255, 69/255), bg = rgb(88/255, 24/255, 
                                                             69/255), cex = 1.5)

# add SH
points(x = morpho.data$SH, y = rep(6, length(morpho.data$SH)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.5))
points(x = morpho.data$SH[which(morpho.data$species == "pd")], 
       y = rep(6, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = rgb(255/255, 87/255, 51/255, 0.5), bg = rgb(255/255, 87/255, 
                                                                   51/255, 0.5))
# Stupendemys
points(x = morpho.data$SH[which(morpho.data$species == "sg")], 
       y = rep(6, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 23, col = rgb(184/255, 16/255, 81/255), bg = rgb(184/255, 16/255, 
                                                              81/255), cex = 1.5)
# Rondonia
points(x = morpho.data$SH[which(morpho.data$species == "it")], 
       y = rep(6, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 22, col = rgb(88/255, 24/255, 69/255), bg = rgb(88/255, 24/255, 
                                                             69/255), cex = 1.5)

# add HC
points(x = morpho.data$HC, y = rep(7, length(morpho.data$HC)), pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.5))
points(x = morpho.data$HC[which(morpho.data$species == "pd")], 
       y = rep(7, length(morpho.data$species[which(morpho.data$species == "pd")])),
       pch = 21, col = rgb(255/255, 87/255, 51/255, 0.5), bg = rgb(255/255, 87/255, 
                                                                   51/255, 0.5))
# Stupendemys
points(x = morpho.data$HC[which(morpho.data$species == "sg")], 
       y = rep(7, length(morpho.data$species[which(morpho.data$species == "sg")])),
       pch = 23, col = rgb(184/255, 16/255, 81/255), bg = rgb(184/255, 16/255, 
                                                              81/255), cex = 1.5)
# Rondonia
points(x = morpho.data$HC[which(morpho.data$species == "it")], 
       y = rep(7, length(morpho.data$species[which(morpho.data$species == "it")])),
       pch = 22, col = rgb(88/255, 24/255, 69/255), bg = rgb(88/255, 24/255, 
                                                             69/255), cex = 1.5)

legend(x = "topright", legend = c("Rondonia Jaw", "Stupendemys geographicus",
                                  "Peltocephalus dumerilianus"), 
       cex = 0.7, bty = "n", pch = c(22, 23, 21), 
       pt.cex = 1, col = c("#581845", "#b81051", "#FF5733"), pt.bg = c("#581845", 
        "#b81051", "#FF5733"), text.font = c(1, rep(3, 2)))


dev.off()

# PCA
pca.res = prcomp(log.morpho.data[, 3:10], center = T, scale = T)
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
pdf("results/morphospace_biplot.pdf", height = 6, width = 8)
biplot(pca.res)
dev.off()

###############################
# Morphospace plot

# create hulls
tab = matrix(c(pca.res$x[, 1], pca.res$x[, 2]), ncol = 2)

tab.em = tab[(which(log.morpho.data$species == "em")), ]
tab.pd = tab[(which(log.morpho.data$species == "pd")), ]
tab.pe = tab[(which(log.morpho.data$species == "pe")), ]
tab.pl = tab[(which(log.morpho.data$species == "pl")), ]
tab.pu = tab[(which(log.morpho.data$species == "pu")), ]
tab.ps = tab[(which(log.morpho.data$species == "ps")), ]
tab.pv = tab[(which(log.morpho.data$species == "pv")), ]
tab.px = tab[(which(log.morpho.data$species == "px")), ]

ch.em = chull(tab.em)
ch.pd = chull(tab.pd)
ch.pe = chull(tab.pe)
ch.pl = chull(tab.pl)
ch.pu = chull(tab.pu)
ch.ps = chull(tab.ps)
ch.pv = chull(tab.pv)
ch.px = chull(tab.px)

# plot PC1 x PC2 without points
pdf("results/morphospace_jaws.pdf", height = 6, width = 8)

plot(pca.res$x[, 1], pca.res$x[, 2], cex = 0, ylim = c(-2.2, 4), xlab = "PC1 (43.5%)",
     ylab = "PC2 (21.4%)")

#plot hulls
polygon(tab.px[ch.px, ], border = NA, col = rgb(63/255, 95/255, 5/255, 0.3), 
        lwd = 1)
polygon(tab.em[ch.em, ], border = NA, col = rgb(178/255, 84/255, 64/255, 0.3), 
        lwd = 1)
polygon(tab.pd[ch.pd, ], border = NA, col = rgb(255/255, 87/255, 51/255, 0.3), 
        lwd = 1)
polygon(tab.pu[ch.pu, ], border = NA, col = rgb(166/255, 203/255, 99/255, 0.3), 
        lwd = 1)
polygon(tab.ps[ch.ps, ], border = NA, col = rgb(133/255, 172/255, 63/255, 0.3), 
        lwd = 1)
polygon(tab.pv[ch.pv, ], border = NA, col = rgb(102/255, 144/255, 26/255, 0.3), 
        lwd = 1)

# plot points
points(pca.res$x[, 2][which(log.morpho.data$species == "pd")] ~ pca.res$x[, 1]
       [which(log.morpho.data$species == "pd")], pch = 22, col = "#FF5733", 
       bg = "#FF5733")
points(pca.res$x[, 2][which(log.morpho.data$species == "em")] ~ pca.res$x[, 1]
       [which(log.morpho.data$species == "em")], pch = 23, col = "#B25440", 
       bg = "#B25440")
points(pca.res$x[, 2][which(log.morpho.data$species == "pe")] ~ pca.res$x[, 1]
       [which(log.morpho.data$species == "pe")], pch = 21, col = "#c4d8a0", 
       bg = "#c4d8a0")
points(pca.res$x[, 2][which(log.morpho.data$species == "pl")] ~ pca.res$x[, 1]
       [which(log.morpho.data$species == "pl")], pch = 21, col = "#c1e383", 
       bg = "#c1e383")
points(pca.res$x[, 2][which(log.morpho.data$species == "pu")] ~ pca.res$x[, 1]
       [which(log.morpho.data$species == "pu")], pch = 21, col = "#a6cb63", 
       bg = "#a6cb63")
points(pca.res$x[, 2][which(log.morpho.data$species == "ps")] ~ pca.res$x[, 1]
       [which(log.morpho.data$species == "ps")], pch = 21, col = "#85ac3f", 
       bg = "#85ac3f")
points(pca.res$x[, 2][which(log.morpho.data$species == "pv")] ~ pca.res$x[, 1]
       [which(log.morpho.data$species == "pv")], pch = 21, col = "#66901a", 
       bg = "#66901a")
points(pca.res$x[, 2][which(log.morpho.data$species == "px")] ~ pca.res$x[, 1]
       [which(log.morpho.data$species == "px")], pch = 21, col = "#3f5f05", 
       bg = "#3f5f05")

points(pca.res$x[, 2][which(log.morpho.data$species == "it")] ~ pca.res$x[, 1]
       [which(log.morpho.data$species == "it")], pch = 25, col = "#581845", 
       bg = "#581845", cex = 1.5)
points(pca.res$x[, 2][which(log.morpho.data$species == "sg")] ~ pca.res$x[, 1]
       [which(log.morpho.data$species == "sg")], pch = 25, col = "#b81051", 
       bg = "#b81051", cex = 1.5)

legend(x = "topright", legend = c("Rondonia Jaw", "Stupendemys geographicus",
        "Erymnochelys madagascariensis", "Peltocephalus dumerilianus", 
        "Podocnemis erythrocephala", "Podocnemis expansa", "Podocnemis lewiana",
        "Podocnemis sextuberculata", "Podocnemis unifilis", "Podocnemis vogli"), 
       cex = 0.7, bty = "n", pch = c(25, 25, 23, 22, 21, 21, 21, 21, 21, 21), 
       pt.cex = 1, col = c("#581845", "#b81051", "#B25440", "#FF5733", "#c4d8a0",
        "#3f5f05", "#c1e383", "#85ac3f", "#a6cb63", "#66901a"), pt.bg = c("#581845",
        "#b81051", "#B25440", "#FF5733", "#c4d8a0", "#3f5f05", "#c1e383", 
        "#85ac3f", "#a6cb63", "#66901a"), ncol = 2, text.font = c(1, rep(3, 9)))


dev.off()

# check data points labels
#text(pca.res$x[, 2] ~ pca.res$x[, 1], labels = rownames(pca.res$x))

# size estimation
size.data = read.csv("size_estimation.csv", header = T, row.names = 1, sep = ",")
reg.size = lm(size.data$CML~size.data$MiL)
summary(reg.size)

plot(size.data$CML~size.data$MiL)
abline(reg.size)
text(size.data$CML~size.data$MiL, labels = rownames(size.data))


# size estimation based on Peltocephalus alone
reg.size.pd = lm(size.data$CML[which(size.data$species == "pd")]
                 ~size.data$MiL[which(size.data$species == "pd")])

summary(reg.size.pd)
plot(size.data$CML[which(size.data$species == "pd")]
     ~size.data$MiL[which(size.data$species == "pd")])
abline(reg.size.pd)


# carapace length plots
carapace = read.csv("carapace.csv", header = T, sep = ",")

plot(y = carapace$CML, x = rep(1, length(carapace$Taxon)), xlim = c(0.5, 2), 
     ylim = c(250, 1000), type = "n", ylab = "carapace length", xlab = "",
     xaxt = "n")

axis(1, labels = c("Peltocephalus dumerilianus","Podocnemis expansa"), 
                   at = c(0.8, 1.7))
points(y = carapace$CML[which(carapace$species == "pd")], pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.5), 
       x = rep(0.8, length(carapace$species[which(carapace$species == "pd")])))

points(y = carapace$CML[which(carapace$species == "pe")], pch = 16, 
       col = rgb(192/255, 192/255, 192/255, alpha = 0.5), 
       x = rep(1.7, length(carapace$species[which(carapace$species == "pe")])))



save.image("giantpodc.RData")
