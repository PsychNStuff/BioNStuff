
#make into percentage; change values as needed
dataframe$waterfallpercent <- (dataframe$rawnumber-0.5)*100 #numbers between 0-1
dataframe$waterfallchange <- (dataframe$rawnumber-dataframe$baseline) #if you want to show a change from baseline
dataframe$waterfalldiff <- 100-(dataframe$rawnumber) #maybe you have numbers between 0-50 that would be better shown subtracted from 100

#add case/control info using a matching variable
dataframe <- rownames_to_column(dataframe, var = "Var1")
dataframe <-inner_join(sampledata[,c("Var1", "diagnosis")], dataframe)

#sort by decreasing numbers so your plot is up on the left and down on the right. I will use waterfallpercent as an example.
dataframe <- dataframe[order(dataframe$waterfallpercent, decreasing = TRUE), ]

#assign colors
col <- ifelse(dataframe$diagnosis == "case", 
              "steelblue",
              "cadetblue")

#make the plot
p <- barplot(dataframe$waterfallpercent,
        space=0.5, ylim=c(-50,50),
        main = "Percent difference",
        ylab="Individual score (-0.5)*100",
        cex.axis=1.2, cex.lab=1.4,
        col = col,
        border = col,
        legend.text = c("Case", "Control"),
        args.legend = list(title = "Diagnosis", fill = c("steelblue", "cadetblue"),
                           border = NA, cex = 0.9))
png("yourdirectory/filename.png")
print(p)
dev.off()
