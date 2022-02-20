library(tidyverse)
#Fig1A
myData <- read.csv("../data/First3BatchesHsMetrics.csv")
smallerTable <- myData[, c(1:2, 48:57)]
toBePlotted <- smallerTable %>% gather('PCT_TARGET_BASES_2X', 'PCT_TARGET_BASES_10X', 
                                       'PCT_TARGET_BASES_20X', 'PCT_TARGET_BASES_30X',
                                       'PCT_TARGET_BASES_40X', 'PCT_TARGET_BASES_50X',
                                       'PCT_TARGET_BASES_100X',  'PCT_TARGET_BASES_250X',
                                       'PCT_TARGET_BASES_500X', 'PCT_TARGET_BASES_1000X',
                                       key = Coverage, value = PCT)
#to get Coverage into digits
tempValue <- str_replace_all(toBePlotted$Coverage, "PCT_TARGET_BASES_", "")
tempValue2 <- str_replace_all(tempValue, "X", "")
toBePlotted$NewCoverage <- factor(tempValue2, levels = c("2", "10", "20", "30", "40", "50", "100",
                                                         "250", "500", "1000"))
toBePlotted$PCT <- toBePlotted$PCT * 100
ggplot(data = toBePlotted) +
    geom_boxplot(map = aes(x= NewCoverage, y = PCT)) +
    ylim(0, 100) +
    xlab("Depth of Coverage") +
    ylab("Percentage of Targets") +
    theme_bw()

#Fig2A
myDataRaw <- read.delim("../data/Fresh_FFPE_hsMetrics.txt", sep = "\t")
par(mfrow=c(1,3))
boxplot(myDataRaw$PF_READS/1000000 ~ myDataRaw$TissueType,
        ylim = c(10, 60), 
        xlab = "", 
        ylab = "Total Reads (Million)")
stripchart(myDataRaw$PF_READS/1000000 ~ myDataRaw$TissueType,
           data = myDataRaw,
           method = "jitter",
           pch = 19,
           vertical = TRUE,
           add = TRUE)
boxplot(myDataRaw$PF_UNIQUE_READS/1000000 ~ myDataRaw$TissueType,
        ylim = c(10, 60),
        xlab = "",
        ylab = "Unique Reads (Million)")
stripchart(myDataRaw$PF_UNIQUE_READS/1000000 ~ myDataRaw$TissueType,
           data = myDataRaw,
           method = "jitter",
           pch = 19,
           vertical = TRUE,
           add = TRUE)
boxplot(myDataRaw$MEAN_TARGET_COVERAGE ~ myDataRaw$TissueType, 
        ylim = c(500, 1000),
        xlab = "",
        ylab = "Depth of Coverge on Targets")
stripchart(myDataRaw$MEAN_TARGET_COVERAGE ~ myDataRaw$TissueType,
           data = myDataRaw,
           method = "jitter",
           pch = 19,
           vertical = TRUE,
           add = TRUE)
par(mfrow=c(1,1))

#p-values for tissue types
t.test(myDataRaw$PF_READS/1000000 ~ myDataRaw$TissueType)
t.test(myDataRaw$PF_UNIQUE_READS/1000000 ~ myDataRaw$TissueType)
t.test(myDataRaw$MEAN_BAIT_COVERAGE ~ myDataRaw$TissueType)

#Fig2B
BAT4samples <-read.csv("../data/DNAamountSamplesHsmetrics.csv")
par(mfrow=c(1,3))
boxplot(BAT4samples$PF_READS/1000000 ~ BAT4samples$StartingDNA,
        ylim = c(10, 60), 
        xlab = "", 
        ylab = "Reads per Sample (Million)")
stripchart(BAT4samples$PF_READS/1000000 ~ BAT4samples$StartingDNA,
           data = BAT4samples,
           method = "jitter",
           pch = 19,
           vertical = TRUE,
           add = TRUE)
boxplot(BAT4samples$PF_UNIQUE_READS/1000000 ~ BAT4samples$StartingDNA,
        ylim = c(10, 60),
        xlab = "",
        ylab = "Unique Reads (Million)")
stripchart(BAT4samples$PF_UNIQUE_READS/1000000 ~ BAT4samples$StartingDNA,
           data = BAT4samples,
           method = "jitter",
           pch = 19,
           vertical = TRUE,
           add = TRUE)
boxplot(BAT4samples$MEAN_TARGET_COVERAGE ~ BAT4samples$StartingDNA, 
        ylim = c(200, 1000),
        xlab = "",
        ylab = "Depth of Coverge on Targets")
stripchart(BAT4samples$MEAN_TARGET_COVERAGE ~ BAT4samples$StartingDNA,
           data = BAT4samples,
           method = "jitter",
           pch = 19,
           vertical = TRUE,
           add = TRUE)
par(mfrow=c(1,1))

#p-values for 50ng v 100ng
t.test(BAT4samples$PF_READS/1000000 ~ BAT4samples$StartingDNA)
t.test(BAT4samples$PF_UNIQUE_READS/1000000 ~ BAT4samples$StartingDNA)
t.test(BAT4samples$MEAN_TARGET_COVERAGE ~ BAT4samples$StartingDNA)

#Fig3
#for Fig3 panel A only
Clinicalcases <- read.delim("../data/Table_4_Clinical_Samples_OncoFocus_ActSeq.txt")
PanelA <- Clinicalcases
PA_pct <- round(table(PanelA$Source)/sum(table(PanelA$Source))*100, digits = 1)
PA_lbls <- paste(names(table(PanelA$Source)), " ", PA_pct, "%", sep = "")
pie(table(PanelA$Source), labels = PA_lbls)
#for Figure 3, panel B, C, and D
PanelBCD <- read.delim("../data/VariantsDetectedInPatients.txt")
#combine inframe indels
PanelBCD[PanelBCD == "inframe_deletion"] <- "Inframe_indel"
PanelBCD[PanelBCD == "inframe_insertion"] <- "Inframe_indel"
#shorten the variant category names
PanelBCD[PanelBCD == "frameshift_variant"] <- "Frameshift"
PanelBCD[PanelBCD == "missense_variant"] <- "Missense"
PanelBCD[PanelBCD == "stop_gained"] <- "Stop_gain"
#remove synonymous 
PanelBCD <- subset(PanelBCD, PanelBCD$Sequence.Ontology..Combined. != "synonymous_variant")
#Panel B: VAF distribution
hist(PanelBCD$AF, main ="", xlab = "VAF", ylab = "Count")
#Panel C: top 20 genes 
barplot(sort(table(PanelBCD$Gene.Names), decreasing = TRUE)[1:20], las = 2, ylab = "Count")
#Panel D: variant types
Variant_pct <- round((table(PanelBCD$Sequence.Ontology..Combined.))/sum(table(PanelBCD$Sequence.Ontology..Combined.))*100, digits = 1)
Variant_lbls <- paste(names(table(PanelBCD$Sequence.Ontology..Combined.)), " ", Variant_pct, "%", sep = "")
pie(table(PanelBCD$Sequence.Ontology..Combined.), labels = Variant_lbls)
#Some numbers cited in the last paragraph of Results
#total valid non-synonymous variants
dim(PanelBCD)[1]
#number of genes with variants
length(unique(PanelBCD$Gene.Names))
#list of cases and number of variants
sort(table(PanelBCD$SampleName))
#clinical cases
unique(PanelBCD$SampleName)
#variant types
paste(names(table(PanelBCD$Sequence.Ontology..Combined.)), " ", Variant_pct, "%", sep = "")

