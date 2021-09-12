
# -- 
require(ssizeRNA)
require(exprso)
require(ggplot2)
require(plyr)
require(RColorBrewer)
require(data.table)

topfeatures = seq(from = 5, to = 100, by=5)
noise = c(0.1, 0.5, 1.0)
fc = c(1.05, 1.1, 1.2)
grdSearch = expand.grid(fc, noise)
colnames(grdSearch) = c("fc", "noise")

allStats = list()
for(y in 1:20){
  
  cat("\rRunning iteration:", y)
  cat("\n")
  
  saveData = list();
  for(x in 1:nrow(grdSearch)){
    
    d = sim.counts(nGenes=10e3, 
                   m = 50, 
                   mu = 10, 
                   disp = grdSearch$noise[[x]], 
                   up = 0.5, 
                   pi0 = 0.80, 
                   fc = grdSearch$fc[[x]])
    
    groups = d$group
    # datExpr = log(d$counts,2)
    saveData[[x]] = t(d$counts)
    saveData[[x]] = data.frame(Dx = groups, saveData[[x]])
    
  }
  grdSearch$var = unlist(lapply(saveData, function(x) sum(var(x[,-1]))))
  
  save_stats = list();
  for( x in 1:length(saveData)){
    
    exprs = saveData[[x]]
    
    training = exprs[sample(1:nrow(exprs), 0.66*nrow(exprs)), ]
    training = training[order(training$Dx, decreasing = T), ]
    trainGrp = data.frame(Dx = as.character(training$Dx))
    testing = exprs[!rownames(exprs) %in% rownames(training), ]
    testing = testing[order(testing$Dx, decreasing = T), ]
    testGrp = data.frame(Dx = as.character(testing$Dx))
    
    array.train = exprso(x = training[,-1],
                         y = trainGrp,
                         label ='Dx')
    
    array.test = exprso(x = testing[,-1], 
                        y = testGrp,
                        label = 'Dx')
    
    array.train <- fsStats(array.train, top = 0, how = "t.test")
    
    
    feature_count = length(array.train@preFilter[[1]])
    if(feature_count < 1){stop("Warning! Expected feature labels, but none provided.")}
    
    
    # build and test performance of classifiers over a grid search
    COST = seq(from=10^-2, to = 10^2, length.out = 5)
    
    svm_gs <- plGrid(array.train = array.train,
                     array.valid = array.test,
                     top = topfeatures,
                     how = "buildSVM",
                     fold = NULL, 
                     # gamma = c(0.5, 1.0, 2.0),
                     kernel = c("linear"),
                     cost = COST)
    
    svm_stats = data.frame(svm_gs@summary)
    
    rf_gs <- plGrid(array.train = array.train,
                    array.valid = array.test,
                    top = topfeatures,
                    how=c("buildRF"),
                    fold = NULL,
                    trees = 500)
    
    lr_gs = plGrid(array.train,
                   array.test,
                   top = topfeatures,
                   how=c("buildLR"), fold = NULL)
    
    lda_gs = plGrid(array.train, 
                    array.test,
                    top = topfeatures,
                    how='buildLDA', fold=NULL)
    
    ann_gs = plGrid(array.train,
                    array.test,
                    top=topfeatures,
                    how='buildANN', fold=NULL)
    
    nb_gs = plGrid(array.train, array.test,
                   top = topfeatures,
                   how='buildNB', fold=NULL)
    
    
    rf_stats = data.frame(rf_gs@summary)
    lr_stats = data.frame(lr_gs@summary)
    lda_stats = data.frame(lda_gs@summary)
    ann_stats = data.frame(ann_gs@summary)
    nb_stats = data.frame(nb_gs@summary)
    
    ttest_selected = ldply(list(rf_stats, svm_stats, lr_stats, lda_stats, ann_stats, nb_stats))
    ttest_selected$fs = "ttest"
    
    save_stats[[x]] = ldply(list(ttest_selected))
    save_stats[[x]]$noise = grdSearch$noise[[x]]
    save_stats[[x]]$effect_size = grdSearch$fc[[x]]
    
  }
  allStats[[y]] = ldply(save_stats)
  
}

names(allStats) = paste("Iteraction.",1:length(allStats), sep="")
all_stats_df = ldply(allStats, .id='iteration')

# --- export performance statistics to text file
fwrite(all_stats_df, "SimDE_MLA_perf.txt", quote = F, row.names=F, sep="\t")

# all_stats_df = all_stats_df[all_stats_df$noise > 0.1, ]

# Function to run correlation on each group
corByGroup = function(df, g, x, y) {
  split = split(df, df[,colnames(df) %in% g])
  stats = lapply(split, function(i) data.frame(broom::tidy(cor.test(i[,colnames(i) %in% x], i[,colnames(i) %in% y]))))
  return(ldply(stats, .id = g))
}
corByGroup(df = all_stats_df, g = "build", x = "valid.auc", y = "train.auc")


all_stats_df = all_stats_df[order(all_stats_df$valid.auc, decreasing = TRUE), ]
all_stats_df$set = paste(all_stats_df$build, "_", all_stats_df$top, "_",all_stats_df$noise, "_", all_stats_df$kernel,"_", all_stats_df$effect_size, "_", all_stats_df$iteration, sep="")
all_stats_df = all_stats_df[!duplicated(all_stats_df$set), ]

sumstats <- ddply(all_stats_df, .(build,top,noise,effect_size,kernel), 
                  summarize,
                  MeanAUC = mean(valid.auc), 
                  MeanAUC_train  = mean(train.auc),
                  SE = sd(valid.auc)/sqrt(length(valid.auc)))

sumstats$param = paste("Noise = ", sumstats$noise, ", Effect size = ", sumstats$effect_size, sep="")
sumstats$build = gsub("build", "", sumstats$build)
sumstats$type = paste(sumstats$build, " ",sumstats$kernel, sep="")
sumstats$type = gsub(" NA", "", sumstats$type)

split = split(sumstats, sumstats$param)
split = lapply(split, function(x) x[order(x$MeanAUC, decreasing = T), ])

percent_boost = list();
split_stats = list()
for(x in 1:length(split)){
  
  temp = split[[x]]
  
  overall = ddply(temp, .(type,param), summarize, grandMean = mean(MeanAUC), Variance = var(MeanAUC))
  overall = overall[order(overall$grandMean, decreasing = T), ]
  split_stats[[x]] = overall[order(overall$Variance, decreasing = F), ]
  best = overall[!duplicated(overall$param), ]
  best = best[order(best$param), ]
  
  percent_boost[[x]] = sapply(overall$grandMean, function(x) round(overall$grandMean[[1]] - x, 3))
  percent_boost[[x]] = data.frame(Boost = percent_boost[[x]], Model = overall$type)
  
  temp$alpha = NA
  temp$alpha = ifelse(temp$type %in% best$type[[1]], 1.0, 0.25)
  split[[x]] = temp
  
}
names(percent_boost) = names(split)
pdf = ldply(percent_boost, .id='param')
ddply(pdf, ~Model, summarize, MeanDiff = mean(Boost))

sumstats = ldply(split)
sumstats$LABEL = ifelse(sumstats$alpha == 1.0, paste("AUC = ", round(sumstats$MeanAUC,2), sep=""), NA)

label_df = sumstats[sumstats$alpha == 1.0, ]
label_df = label_df[order(label_df$MeanAUC, decreasing = T), ]
label_df = label_df[!duplicated(label_df$param), ]
label_df = label_df[order(label_df$param), ]

png("Plot_MLA_Accuracy_SimDE.png",res=300,units="in",height=7,width=8)
ggplot(sumstats, aes(x = factor(top), y = MeanAUC, col = type, alpha = factor(alpha))) +
  facet_wrap(~param, ncol = 3) +
  scale_color_manual('Model', values = c(rev(brewer.pal(n = 5, 'Set1')), 'black')) +
  # geom_errorbar(aes(ymin = MeanAUC - (1.96*SE), ymax = MeanAUC + (1.96*SE)), position=position_dodge(0.9), show.legend = F, lwd = 0.2, width = 0.5) +
  # geom_point(position=position_dodge(0.9), pch = 15) +
  scale_alpha_manual(values = c(0.3, 1.0)) +
  geom_line(aes(group = type), position=position_dodge(0.9), lwd=0.5) +
  geom_text(data = label_df, aes(x = 10, y = 1.1, label = paste("\nMax ",LABEL,sep="")), size=3, show.legend = FALSE) +
  geom_hline(aes(yintercept = 0.5), alpha=0.5, lty = 1, lwd = 0.3, col = 'black') +
  xlab("Number of genes in classifier") +
  ylab("Mean Area Under the Curve") +
  theme_grey() +
  guides(alpha=FALSE) +
  theme(panel.border=element_rect(size=1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(color='black'),
        axis.text=element_text(size = 10),
        axis.text.x=element_text(size = 8, angle = 90, vjust = 0.5))
dev.off()

split = lapply(split, function(x) x[1,])
print(split)

overall = ddply(sumstats, .(type,param), summarize, grandMean = mean(MeanAUC))
overall = overall[order(overall$grandMean, decreasing = T), ]
best = overall[!duplicated(overall$param), ]
best = best[order(best$param), ]
print(best)
worst = overall[order(overall$grandMean, decreasing = FALSE), ]
worst = worst[!duplicated(worst$param), ]
worst = worst[order(worst$param), ]
print(worst)

