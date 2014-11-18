# resampling figure with correlation with shannon index instead of overall complexity

# loading df with factors
alldata <- read.csv("all_factor_data.csv", row.names=1)[1:41, ]


scent_diversity <- read.csv(paste("C:\\Users\\Martin\\Studium\\",
                                  "projects\\sealscent\\data_files\\",
                                  "Rdata\\csv_files\\",
                                  "scent diversity.csv", sep = ""),
                                  row.names=1)