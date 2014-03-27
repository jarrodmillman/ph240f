kirp_normalized <- "../data/kirp/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/unc.edu.0057e13d-5489-4c48-9772-9fd45db4d84f.1041387.rsem.genes.normalized_results"
kirc_normalized <- "../data/kirc/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/unc.edu.005a0c90-a762-4df7-8b44-b4facd06e9b1.1230067.rsem.genes.normalized_results"

kirp <- read.table(kirp_normalized, sep="\t", header=TRUE)
kirc <- read.table(kirc_normalized, sep="\t", header=TRUE)

hist(log(kirp$normalized_count+1))
hist(log(kirc$normalized_count+1))
