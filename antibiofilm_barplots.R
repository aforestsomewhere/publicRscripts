#Script for EPS-F antibiofilm analysis
#Katie O'Mahony
#09/03/24

library(ggplot2)
library(dplyr)
library(multcompView)
library(multcomp)

#import data
df <- read.csv("data/epsf_lmono6179.csv")
df <- na.omit(df)
df$group <-factor(df$group)

# analysis of variance
anova.df <- aov(od ~ group, data = df)
summary(anova.df)
#check homogeneity of variance
plot(anova.df, 1)
#check normality
plot(anova.df, 2)
# Check normality with Shapiro-Wilk test
aov_residuals <- residuals(object = anova.df )
shapiro.test(x = aov_residuals ) #departure from normality, skew (kruskal wallis)

# table with factors, means and standard deviation
df_summary <- group_by(df, group) %>%
  summarise(mean=mean(od), sd=sd(od)) %>%
  arrange(desc(mean))
print(df_summary)

# Simple Tukey's test for compact letter display
tukey.df <- TukeyHSD(anova.df)
cld.df <- multcompLetters4(anova.df, tukey.df)
cld <- as.data.frame.list(cld.df$group)
df_summary$Tukey <- cld$Letters

#plot with compact letter display
df_summary$group <- factor(df_summary$group, levels=c("control","5 mg/mL", "10 mg/mL", "20 mg/mL"))

ggplot(df_summary, aes(group, mean, fill=group)) + 
  geom_bar(stat="identity", width=0.8, alpha=0.8, colour="black") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, colour="black") +
  scale_fill_brewer(palette = "Greys") +
  xlab("") + ylab(expression(OD[ 570])) + 
  geom_text(aes(label=Tukey), nudge_y = 0.2, size = 5) + theme_classic()

#More comprehensive Tukeys for p values
tukey_results <- glht(anova.df, linfct = mcp(group = "Tukey"))
summary(tukey_results) #use these for table of p values

       
