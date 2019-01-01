library(reshape)
library(ggplot2)


raw_counts <- readRDS('raw_counts.rds')
raw_counts <- melt(raw_counts)
raw_counts$value <- raw_counts$value/1e6



ggplot(raw_counts, aes(x = variable, y = value, group = variable, fill = variable)) +
       stat_boxplot(geom = 'errorbar', width = 0.25) +
       geom_boxplot(width = 0.5) +
       stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .5, col = 'orange', lwd = .75) +
       stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .5, col = 'black', lwd = 1) +
       ggtitle('Total Raw Read Counts by Weeks') +
       xlab('') +
       ylab('Raw Counts in Millions') +
       scale_y_continuous(breaks = seq(0,12, 1)) +
       theme(plot.title = element_text(hjust = 0.5, size = 25),
             panel.background = element_blank(),
             axis.line = element_line(size = 0.5),
             axis.text = element_text(size = 15),
             axis.title = element_text(size = 20),
             legend.position = 'none')







