nmds <- read.table(file="CityofCalgary.braycurtis.0.03.lt.nmds.axes",sep="\t", head=TRUE, row.names="group")
#plot(nmds, xlim = xlim, ylim = ylim, asp = 1, type = "p")
points(nmds$axis1, nmds$axis2, pch = 1)
text(nmds, labels = rownames(nmds), col = 1)
