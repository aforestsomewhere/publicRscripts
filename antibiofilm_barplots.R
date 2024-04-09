#Antibiofilm analysis for EPS-F
#Katie O'Mahony
#29/03/24 - updated 09/04/24. X axis is now set at 0, labels added.

library(ggplot2)
library(dplyr)
library(multcompView)
library(multcomp)
library(cowplot)
library(broom)

#Multipanel, one panel for each species
#1. L. monocytogenes 6179
#import data
rm(df)
df <- read.csv("data/epsf_lmono6179.csv")
df <- na.omit(df)
df$group <-factor(df$group)

#divide od by mean od of control
mean_value <- df %>% filter(group=="control") %>% summarise(mean=mean(od))
mean_value <- mean_value[,1]
#df$relative <- (mean_value - df$od)/mean_value*100
df$relative <- (1-df$od/mean_value)*100
# table with factors, means and standard deviation
df_summary <- group_by(df, group) %>%
  summarise(mean=mean(relative), sd=sd(relative)) %>%
  arrange(desc(mean))

# analysis of variance
anova.df.rel <- aov(relative ~ group, data = df)
summary(anova.df.rel)
#check homogeneity of variance
plot(anova.df.rel, 1)
#check normality
plot(anova.df.rel, 2)
# Check normality with Shapiro-Wilk test
aov_residuals_rel <- residuals(object = anova.df.rel )
shapiro.test(x = aov_residuals_rel ) #departure from normality, skew (kruskal wallis)

# Simple Tukey's test for compact letter display
tukey.df.rel <- TukeyHSD(anova.df.rel)
cld.df.rel <- multcompLetters4(anova.df.rel, tukey.df.rel)
cld.rel <- as.data.frame.list(cld.df.rel$group)
df_summary$Tukey <- cld.rel$Letters

#reorder
df_summary$group <- factor(df_summary$group, levels=c("control","5 mg/mL", "10 mg/mL", "20 mg/mL"))

p1 <- ggplot(df_summary %>% filter(group!="control"), aes(group, mean, fill=group)) + 
  geom_bar(stat="identity", width=0.4, alpha=0.8, colour="black") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, colour="black") +
  scale_fill_brewer(palette = "Greys") +
  xlab("") + ylab("Biofilm Inhibition (%)") + 
  geom_text(aes(label=Tukey), nudge_y = 5, size = 5) + theme_classic() +
  ggtitle("L. monocytogenes 6179") + theme(legend.position="bottom") +
  guides(fill=guide_legend(title="EPS-F Concentration")) + 
  ylim(-100,100)

p1
#isolate legend and remove
#leg <- get_legend(p1)
leg = cowplot::get_plot_component(p1, 'guide-box-bottom', return_all = TRUE)

#replot without legend
p1_title <- expression(paste(italic("L. monocytogenes "), "6179"))

p1 <- ggplot(df_summary, aes(group, mean, fill=group)) + 
  geom_bar(stat="identity", width=0.5, alpha=0.8, colour="black") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, colour="black") +
  scale_fill_brewer(palette = "Greys") +
  xlab("") + ylab("Biofilm Inhibition (%)") + 
  geom_hline(yintercept = 0) +
  geom_text(aes(label=Tukey), nudge_y = 30, size = 5) + theme_classic() +
  ggtitle(p1_title) + theme(legend.position = "none",
                            axis.text.x=element_text(angle=45, vjust=-.01, size=12),
                            axis.title.y=element_text(size=16),
                            plot.title = element_text(hjust = 0.5),
                            axis.line.x = element_blank(),
                            axis.ticks.x = element_blank()) + 
  ylim(-38,100)
p1

#write stats to output
write.csv(df_summary, "data/df_summary_lm6179.csv")
write.csv(tidy(anova.df.rel), "data/anova_lm6179.csv")
write.csv(tidy(tukey.df.rel), "data/tukey_lm6179.csv")

#2. L. monocytogenes DSA_248
#import data
rm(df)
df <- read.csv("data/epsf_lmono_scott.csv")
df <- na.omit(df)
df$group <-factor(df$group)
#remove the blank readings
df <- df %>% filter(group!="blank")

#divide od by mean od of control
mean_value <- df %>% filter(group=="control") %>% summarise(mean=mean(od))
mean_value <- mean_value[,1]
df$relative <- (1-df$od/mean_value)*100

# table with factors, means and standard deviation
df_summary <- group_by(df, group) %>%
  summarise(mean=mean(relative), sd=sd(relative)) %>%
  arrange(desc(mean))

# analysis of variance
anova.df.rel <- aov(relative ~ group, data = df)
summary(anova.df.rel)
#check homogeneity of variance
plot(anova.df.rel, 1)
#check normality
plot(anova.df.rel, 2)
# Check normality with Shapiro-Wilk test
aov_residuals_rel <- residuals(object = anova.df.rel )
shapiro.test(x = aov_residuals_rel ) #departure from normality, skew (kruskal wallis)

# Simple Tukey's test for compact letter display
tukey.df.rel <- TukeyHSD(anova.df.rel)
cld.df.rel <- multcompLetters4(anova.df.rel, tukey.df.rel)
cld.rel <- as.data.frame.list(cld.df.rel$group)
df_summary$Tukey <- cld.rel$Letters

#reorder
df_summary$group <- factor(df_summary$group, levels=c("control","5 mg/mL", "10 mg/mL", "20 mg/mL"))

#replot without legend
p2_title <- expression(paste(italic("L. monocytogenes "), "DSA_248"))

p2 <- ggplot(df_summary, aes(group, mean, fill=group)) + 
  geom_bar(stat="identity", width=0.5, alpha=0.8, colour="black") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, colour="black") +
  scale_fill_brewer(palette = "Greys") +
  xlab("") + ylab("") + 
  geom_hline(yintercept = 0) +
  geom_text(aes(label=Tukey), nudge_y = 20, size = 5) + theme_classic() +
  ggtitle(p2_title) + theme(legend.position = "none",
                            axis.text.x=element_text(angle=45, vjust=-.01, size=12),
                            axis.title.y=element_text(size=16),
                            plot.title = element_text(hjust = 0.5),
                            axis.line.x = element_blank(),
                            axis.ticks.x = element_blank()) + 
  ylim(-38,100)
  ylim(-38,100)
p2
#write stats to output
write.csv(df_summary, "data/df_summary_lmd.csv")
write.csv(tidy(anova.df.rel), "data/anova_lmd.csv")
write.csv(tidy(tukey.df.rel), "data/tukey_lmd.csv")

#3. P. fluorescens
#import data
rm(df)
df <- read.csv("data/epsf_pfluorescens.csv")
df <- na.omit(df)
df$group <-factor(df$group)
#remove the blank readings
df <- df %>% filter(group!="blank")

#divide od by mean od of control
mean_value <- df %>% filter(group=="control") %>% summarise(mean=mean(od))
mean_value <- mean_value[,1]
df$relative <- (1-df$od/mean_value)*100

# table with factors, means and standard deviation
df_summary <- group_by(df, group) %>%
  summarise(mean=mean(relative), sd=sd(relative)) %>%
  arrange(desc(mean))

# analysis of variance
anova.df.rel <- aov(relative ~ group, data = df)
summary(anova.df.rel)
#check homogeneity of variance
plot(anova.df.rel, 1)
#check normality
plot(anova.df.rel, 2)
# Check normality with Shapiro-Wilk test
aov_residuals_rel <- residuals(object = anova.df.rel )
shapiro.test(x = aov_residuals_rel ) #departure from normality, skew (kruskal wallis)

# Simple Tukey's test for compact letter display
tukey.df.rel <- TukeyHSD(anova.df.rel)
cld.df.rel <- multcompLetters4(anova.df.rel, tukey.df.rel)
cld.rel <- as.data.frame.list(cld.df.rel$group)
df_summary$Tukey <- cld.rel$Letters

#reorder
df_summary$group <- factor(df_summary$group, levels=c("control","5 mg/mL", "10 mg/mL", "20 mg/mL"))

#replot without legend
p3_title <- expression(paste(italic("P. fluorescens "), "DSA_L22"))


p3 <- ggplot(df_summary, aes(group, mean, fill=group)) + 
  geom_bar(stat="identity", width=0.5, alpha=0.8, colour="black") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, colour="black") +
  scale_fill_brewer(palette = "Greys") +
  xlab("") + ylab("") + 
  geom_hline(yintercept = 0) +
  geom_text(aes(label=Tukey), nudge_y = 25, size = 5) + theme_classic() +
  ggtitle(p3_title) + theme(legend.position = "none",
                            axis.text.x=element_text(angle=45, vjust=-.01, size=12),
                            axis.title.y=element_text(size=16),
                            plot.title = element_text(hjust = 0.5),
                            axis.line.x = element_blank(),
                            axis.ticks.x = element_blank()) + 
  ylim(-38,100)
p3
#write stats to output
write.csv(df_summary, "data/df_summary_pf.csv")
write.csv(tidy(anova.df.rel), "data/anova_pf.csv")
write.csv(tidy(tukey.df.rel), "data/tukey_pf.csv")


p4 <- plot_grid(p1, p2, p3, nrow=1, labels="AUTO", label_size = 20)
plot_grid(p4, leg, nrow=2, rel_heights = c(1,0.4))
#ggsave("figures/EPSF_antibiofilm.png")
ggsave(filename = "figures/EPSF_antibiofilm2.png", width = 10, height = 6, dpi = 600)
