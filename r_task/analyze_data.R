data <- read.csv("r_task/sample_data.csv")

cat("Mean Score:", mean(data$Score), "\n")

treatment <- data[data$Group == "Treatment", ]
cat("Max Score in Treatment:", max(treatment$Score), "\n")

png("r_task/score_boxplot.png")
boxplot(Score ~ Group, data = data, main = "Score Distribution by Group")
dev.off()

