metadata <- data.frame(
  sample = c("21-1A-AD","20-1T-AD","23-2A-AD","22-2T-AD",
             "26-3A-AD","24-3T-AD","27-5A-AD","25-5T-AD",
             "29-6T-AD","31-7T-AD","28-8T-AD","30-9T-AD",
             "18-10A-Old","11-10T-Old","19-11A-Old","13-11T-Old",
             "15-13T-Old","16-14T-Old","12-6A-Old","14-7A-Old",
             "10-8A-Old","17-9A-Old",
             "2-12A-Young","4-13A-Young","6-14A-Young","8-15A-Young",
             "9-16A-Young","3-17T-Young","5-18T-Young","7-19T-Young"),
  
  diagnosis = c(rep("AD", 12),
                rep("Old", 10),
                rep("Young", 8)),
  
  region = c("A","T","A","T","A","T","A","T","T","T","T","T",
             "A","T","A","T","T","T","A","A","A","A",
             "A","A","A","A","A","T","T","T"),
  
  stringsAsFactors = FALSE
)
write.csv(metadata,
          "\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\data\\metadata.csv",
          row.names = FALSE)
# Load the metadata
metadata <- read.csv("\\\\wsl.localhost\\Ubuntu\\home\\maulibhavsar\\AD_RNAseq_Project\\data\\metadata.csv")

# Check it
dim(metadata)        # Should say 30 rows, 3 columns
head(metadata)       # First 6 rows
table(metadata$diagnosis)  # Should show AD=12, Old=10, Young=8

