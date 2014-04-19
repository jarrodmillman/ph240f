#https://tcga-data.nci.nih.gov/tcga/
#Kidney renal clear cell carcinoma [KIRC]
#Kidney renal papillary cell carcinoma [KIRP]
#https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/kirc/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/
#https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/kirp/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/

HTTPBASE=https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/
KIRCBASE=${HTTPBASE}/kirc/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2
KIRPBASE=${HTTPBASE}/kirp/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/

PATTERN="unc.edu.0*.rsem.genes.results"


cd kirc
wget -np -nd -r -A ${PATTERN}  ${KIRCBASE}/unc.edu_KIRC.IlluminaHiSeq_RNASeqV2.Level_3.1.4.0/
cd ../kirp
wget -np -nd -r -A ${PATTERN}  ${KIRPBASE}/unc.edu_KIRP.IlluminaHiSeq_RNASeqV2.Level_3.2.2.0/
cd ..
