library("bio3d")
dcd <- read.dcd("pro_pca.dcd")
pdb <- read.pdb("pro_pca_end.pdb")
ca.inds <- atom.select(pdb, elety="CA")
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd, fixed.inds=ca.inds$xyz,mobile.inds=ca.inds$xyz)
print(dim(xyz) == dim(dcd))

cij <- dccm(xyz[,ca.inds$xyz])
png("pca_corr.png", width=5000,height=4000, res=600)
plot(cij)
dev.off()

pc <- pca.xyz(xyz[,ca.inds$xyz])
png("pca_pc.png", width=5000,height=4000, res=600)
plot(pc, col=bwr.colors(nrow(xyz)))
dev.off()

print("all done ! good day ! ")