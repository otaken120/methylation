source("https://bioconductor.org/biocLite.R")
biocLite("conumee")
biocLite("minfiData")
biocLite("missMethyl")
biocLite("Gviz")
biocLite("DMRcate")
biocLite("ChAMP")
biocLite("CopyNumber450kData")

# Tutorail: 
# The conumee vignette https://www.bioconductor.org/packages/devel/bioc/vignettes/conumee/inst/doc/conumee.html

library(minfi)
library(conumee)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICmanifest)

RGsetTNBC <- read.metharray.exp(base = "/users/kendeng/desktop/EPIC_Methylation/")  # adjust path
#RGsetTCGA <- read.450k.url()  # use default parameters for vignette examples
MsetTNBC <- preprocessIllumina(RGsetTNBC)
MsetTNBC

data(exclude_regions)
data(detail_regions)
head(detail_regions, n = 2)

#EPIC annotation
anno <- CNV.create_anno(array_type = "EPIC", exclude_regions = exclude_regions, detail_regions = detail_regions)
anno

#Get my TNBC EPIC methylation data
tnbc.data <- CNV.load(MsetTNBC)
tnbc.data
names(tnbc.data) # names of samples

# Get Control normal breast EPIC data 
# GSE81939; 850k of nine breast HMECs (Human Mammary Epithelial Cells) samples (184AA3, 184AA2, 184AA4, 184A1, 184D, 184BE1, 184B5, 184B5ME, 184FMY2). 
# Data downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81939
RGsetControl <- read.metharray.exp(base = "/users/kendeng/desktop/GSE81939_RAW/")  # adjust path
MsetControl <- preprocessIllumina(RGsetControl)
controls.data <- CNV.load(MsetControl)
controls.data

# To create a loop to go through each patient sample "names(tnbc.data)", and save the following plots with unique file names without overwriting, and save each sample's CNV in a specific file.
# Name of 24 samples:  
# [1] "200511490097_R01C01" "200511490097_R02C01" "200511490097_R03C01" "200511490097_R04C01" "200511490097_R05C01" "200511490097_R06C01" "200511490097_R07C01"
# [8] "200511490097_R08C01" "200511490136_R01C01" "200511490136_R02C01" "200511490136_R03C01" "200511490136_R04C01" "200511490136_R05C01" "200511490136_R06C01"
# [15] "200511490136_R07C01" "200511490136_R08C01" "200511490150_R01C01" "200511490150_R02C01" "200511490150_R03C01" "200511490150_R04C01" "200511490150_R05C01"
# [22] "200511490150_R06C01" "200511490150_R07C01" "200511490150_R08C01"

fileDir='/users/kendeng/desktop/'

pdf(sprintf("%s%s_genomeplot.pdf",fileDir,names(tnbc.data)[1]))
CNV.genomeplot(y)
CNV.detailplot(y, name = "CCNE1")
dev.off()

pdf("/users/kendeng/desktop/CNV_plot.pdf")    
for(i in 1:length(names(tnbc.data)))
  {
    names(tnbc.data)[i]
    y <- CNV.fit(tnbc.data[i], controls.data, anno)
    y <- CNV.bin(y)
    y <- CNV.detail(y)
    y <- CNV.segment(y)
    
    #pdf(sprintf("%s%s_genomeplot.pdf",fileDir,names(tnbc.data)[i]))
    CNV.genomeplot(y)
    CNV.genomeplot(y, chr = "chr19")
    CNV.detailplot(y, name = "CCNE1")
    CNV.detailplot_wrap(y)

#head(CNV.write(y, what = "segments"), n = 5)

# To add a sample-specific file name 
#CNV.write(y, what = "segments")
}
dev.off()