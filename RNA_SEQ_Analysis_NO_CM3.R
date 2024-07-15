---
title: "Mar_RNASEQ_NO_CM3"
author: "Sebastian Vanin"
date: "2023-06-29"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Removing the CM3 cell line from analysis

```{r}
Counts_filt<- Counts[,-c(4,8,12,16)]
Sample_Info_df_filt <- Sample_Info_df[-c(4,8,12,16),]
CPM_filt <- cpm(Counts_filt)

#filtering ou tlowly expressed genes
thresh_filt <- CPM_filt > 0.4
filtered <- rowSums(thresh_filt) >= 3
summary(filtered)
Counts_keep_filt <- Counts_filt[filtered,]
dim(Counts_keep_filt)


Gene_info_filt <- gconvert(query = rownames(Counts_keep_filt),
                      organism = "hsapiens",
                      target = "ENSG",
                      filter_na = F,
                      mthreshold = 1)


Counts_Final_filt <- as_tibble(Counts_keep_filt)
Counts_Final_filt$ENSEMBL <- Gene_info_filt$input
Counts_Final_filt$Symbol <- Gene_info_filt$name
Counts_Final_filt$Description <- Gene_info_filt$description

Counts_Final_filt <- Counts_Final_filt %>%
  drop_na(Symbol)

dim(Counts_Final_filt)





```

```{r}
dgeFILT <- DGEList(counts = Counts_Final_filt[,1:28], 
                 group = Sample_Info_df_filt$Meta.Group, 
                 genes = Counts_Final_filt$Symbol)

dgeFILT$genes$ENSEMBL <- Counts_Final_filt$ENSEMBL
dgeFILT$genes$Description <- Counts_Final_filt$Description

dgeFILT <- calcNormFactors(dgeFILT)

design_filt <- model.matrix(~0 +Sample_Info_df_filt$Meta.Group,
                            data = dgeFILT$samples )
colnames(design_filt) <- gsub("Sample_Info_df_filt\\$Meta.Group", "", colnames(design_filt))
design_filt


dgeFILT<- estimateDisp(dgeFILT, design_filt)
plotBCV(dgeFILT)

fit_filt <- glmQLFit(dgeFILT, design = design_filt)


```

```{r}
color_MDS <- c("black", "red4", "blue4", "purple4", "darkgrey", "red1", "cyan4", "purple1")[Sample_Info_df_filt$Meta.Group]
par(mar = c(5,5,5,10), xpd = TRUE)
plotMDS(dgeFILT , main = "PCA Colored By Group (No CM3)- Dim 1&2", cex.main = 2, cex.lab = 1.2, col = color_MDS)
legend("topright", 
       legend=levels(Sample_Info_df_filt$Meta.Group), 
       col = c("black", "red4", "blue4", "purple4", "darkgrey", "red1", "cyan4", "purple1"), 
       pch =16,
       bty = "n",
       inset = c(-0.2, 0),
       title = "Treatment",
       title.cex = 1.5)

#Doesn't seem like there's much grouping by the meta groups
```

### Differential Expression Analysis - Within condition comparisons

First, the comparisons of treatment vs vehicle will be examined within the control and schizophrenia condition groups.

CONTROL: THC VS VEHICLE

```{r}
Cntrl_THC_Vs_VEH_filt_qlf <- glmQLFTest(fit_filt,
                                   contrast = c(-1, 1,0,0,0,0,0,0))
Cntrl_THC_Vs_VEH_filt_qlf

Cntrl_THC_Vs_VEH_filt_df <- as.data.frame(topTags(Cntrl_THC_Vs_VEH_filt_filt_qlf,
                                             n = Inf,
                                             sort.by = "PValue",
                                             adjust.method = "fdr"))
Cntrl_THC_Vs_VEH_filt_df <- Cntrl_THC_Vs_VEH_filt_df %>% 
  mutate(Significance = case_when(FDR < 0.05 & logFC > 0 ~ "Up",
                                  FDR < 0.05 & logFC < 0 ~ "Down",
                                  FDR > 0.05 ~ "No"))

table(Cntrl_THC_Vs_VEH_filt_df$Significance)
#no DE Genes
```

CONTROL: CBD VS VEHICLE

```{r}
Cntrl_CBD_Vs_VEH_qlf_filt <- glmQLFTest(fit_filt,
                                   contrast = c(-1, 0,1,0,0,0,0,0))
Cntrl_CBD_Vs_VEH_qlf_filt

Cntrl_CBD_Vs_VEH_df <- as.data.frame(topTags(Cntrl_CBD_Vs_VEH_qlf_filt,
                                             n = Inf,
                                             sort.by = "PValue",
                                             adjust.method = "fdr"))
Cntrl_CBD_Vs_VEH_df <- Cntrl_CBD_Vs_VEH_df %>% 
  mutate(Significance = case_when(FDR < 0.05 & logFC > 0 ~ "Up",
                                  FDR < 0.05 & logFC < 0 ~ "Down",
                                  FDR > 0.05 ~ "No"))

table(Cntrl_CBD_Vs_VEH_df$Significance)

#no DE Genes
```

CONTROL: COMBO VS VEHICLE

```{r}
Cntrl_COMBO_Vs_VEH_filt_qlf <- glmQLFTest(fit_filt,
                                   contrast = c(-1, 0,0,1,0,0,0,0))
Cntrl_COMBO_Vs_VEH_filt_qlf

Cntrl_COMBO_Vs_VEH_filt_df <- as.data.frame(topTags(Cntrl_COMBO_Vs_VEH_filt_qlf,
                                             n = Inf,
                                             sort.by = "PValue",
                                             adjust.method = "fdr"))
Cntrl_COMBO_Vs_VEH_filt_df <- Cntrl_COMBO_Vs_VEH_filt_df %>% 
  mutate(Significance = case_when(FDR < 0.05 & logFC > 0 ~ "Up",
                                  FDR < 0.05 & logFC < 0 ~ "Down",
                                  FDR > 0.05 ~ "No"))

table(Cntrl_COMBO_Vs_VEH_filt_df$Significance)
#no DE Genes
```

Now the comparisons between treatments and vehicle in the schizophrenic condition

SCZ: THC VS VEHICLE

```{r}
SCZ_THC_Vs_VEH_filt_qlf <- glmQLFTest(fit_filt,
                                   contrast = c(0,0,0,0,-1,1,0,0))
SCZ_THC_Vs_VEH_filt_qlf

SCZ_THC_Vs_VEH_filt_df <- as.data.frame(topTags(SCZ_THC_Vs_VEH_filt_qlf,
                                             n = Inf,
                                             sort.by = "PValue",
                                             adjust.method = "fdr"))
SCZ_THC_Vs_VEH_filt_df <- SCZ_THC_Vs_VEH_filt_df %>% 
  mutate(Significance = case_when(FDR < 0.05 & logFC > 0 ~ "Up",
                                  FDR < 0.05 & logFC < 0 ~ "Down",
                                  FDR > 0.05 ~ "No"))

table(SCZ_THC_Vs_VEH_filt_df$Significance)
#no DE Genes
```

SCZ: CBD VS VEHICLE

```{r}
SCZ_CBD_Vs_VEH_filt_qlf <- glmQLFTest(fit_filt,
                                   contrast = c(0,0,0,0,-1,0,1,0))
SCZ_CBD_Vs_VEH_filt_qlf

SCZ_CBD_Vs_VEH_filt_df <- as.data.frame(topTags(SCZ_CBD_Vs_VEH_filt_qlf,
                                             n = Inf,
                                             sort.by = "PValue",
                                             adjust.method = "fdr"))
SCZ_CBD_Vs_VEH_filt_df <- SCZ_CBD_Vs_VEH_filt_df %>% 
  mutate(Significance = case_when(FDR < 0.05 & logFC > 0 ~ "Up",
                                  FDR < 0.05 & logFC < 0 ~ "Down",
                                  FDR > 0.05 ~ "No"))

table(SCZ_CBD_Vs_VEH_filt_df$Significance)

#no DE Genes
```

SCZ: COMBO VS VEHICLE

```{r}
SCZ_COMBO_Vs_VEH_filt_qlf <- glmQLFTest(fit_filt,
                                   contrast = c(0,0,0,0,-1,0,0,1))
SCZ_COMBO_Vs_VEH_filt_qlf

SCZ_COMBO_Vs_VEH_filt_df <- as.data.frame(topTags(SCZ_COMBO_Vs_VEH_filt_qlf,
                                             n = Inf,
                                             sort.by = "PValue",
                                             adjust.method = "fdr"))
SCZ_COMBO_Vs_VEH_filt_df <- SCZ_COMBO_Vs_VEH_filt_df %>% 
  mutate(Significance = case_when(FDR < 0.05 & logFC > 0 ~ "Up",
                                  FDR < 0.05 & logFC < 0 ~ "Down",
                                  FDR > 0.05 ~ "No"))

table(SCZ_COMBO_Vs_VEH_filt_df$Significance)
#no DE Genes
```

### Differential Expression Analysis - Between condition comparisons

First, the contrasts of interest will be generated.

```{r}
myContrasts_filt<- makeContrasts(SZC.THCvsCNT.THC = (ST-SV)-(CT-CV),
                            SZC.CBDvsCNT.CBD = (SC-SV)-(CC-CV),
                            SZC.COMBOvsCNT.COMBO = (STC-SV)-(CTC-CV),
                            SZC.VEHvsCNT.VEH = SV-CV,
                            levels = design_filt)

```

SCZ Vehicle vs CNTRL Vehicle

```{r}
SCZ.VEH_Vs_CNTRL.VEH_filt_qlf <- glmQLFTest(fit_filt,
                                   contrast = myContrasts[,"SZC.VEHvsCNT.VEH"])
SCZ.VEH_Vs_CNTRL.VEH_filt_qlf

SCZ.VEH_Vs_CNTRL.VEH_filt_df <- as.data.frame(topTags(SCZ.VEH_Vs_CNTRL.VEH_filt_qlf,
                                             n = Inf,
                                             sort.by = "PValue",
                                             adjust.method = "fdr"))
SCZ.VEH_Vs_CNTRL.VEH_filt_df <- SCZ.VEH_Vs_CNTRL.VEH_filt_df %>% 
  mutate(Significance = case_when(FDR < 0.05 & logFC > 0 ~ "Up",
                                  FDR < 0.05 & logFC < 0 ~ "Down",
                                  FDR > 0.05 ~ "No"))

table(SCZ.VEH_Vs_CNTRL.VEH_filt_df$Significance)

```

SCZ THC vs CNTRL THC

```{r}
SZC.THCvsCNT.THC_filt_qlf <- glmQLFTest(fit_filt,
                                   contrast = myContrasts[,"SZC.THCvsCNT.THC"])
SZC.THCvsCNT.THC_filt_qlf

SZC.THCvsCNT.THC_filt_df <- as.data.frame(topTags(SZC.THCvsCNT.THC_filt_qlf,
                                             n = Inf,
                                             sort.by = "PValue",
                                             adjust.method = "fdr"))
SZC.THCvsCNT.THC_filt_df <- SZC.THCvsCNT.THC_filt_df %>% 
  mutate(Significance = case_when(FDR < 0.05 & logFC > 0 ~ "Up",
                                  FDR < 0.05 & logFC < 0 ~ "Down",
                                  FDR > 0.05 ~ "No"))

table(SZC.THCvsCNT.THC_filt_df$Significance)

```

SCZ CBD vs CNTRL CBD

```{r}
SZC.CBDvsCNT.CBD_filt_qlf <- glmQLFTest(fit_filt,
                                   contrast = myContrasts[,"SZC.CBDvsCNT.CBD"])
SZC.CBDvsCNT.CBD_filt_qlf

SZC.CBDvsCNT.CBD_filt_df <- as.data.frame(topTags(SZC.CBDvsCNT.CBD_filt_qlf,
                                             n = Inf,
                                             sort.by = "PValue",
                                             adjust.method = "fdr"))
SZC.CBDvsCNT.CBD_filt_df <- SZC.CBDvsCNT.CBD_filt_df %>% 
  mutate(Significance = case_when(FDR < 0.05 & logFC > 0 ~ "Up",
                                  FDR < 0.05 & logFC < 0 ~ "Down",
                                  FDR > 0.05 ~ "No"))

table(SZC.CBDvsCNT.CBD_filt_df$Significance)

```

SCZ COMBO vs CNTRL COMBO

```{r}
SZC.COMBOvsCNT.COMBO_filt_qlf <- glmQLFTest(fit_filt,
                                   contrast = myContrasts[,"SZC.COMBOvsCNT.COMBO"])
SZC.COMBOvsCNT.COMBO_filt_qlf

SZC.COMBOvsCNT.COMBO_filt_df <- as.data.frame(topTags(SZC.COMBOvsCNT.COMBO_filt_qlf,
                                             n = Inf,
                                             sort.by = "PValue",
                                             adjust.method = "fdr"))
SZC.COMBOvsCNT.COMBO_filt_df <- SZC.COMBOvsCNT.COMBO_filt_df %>% 
  mutate(Significance = case_when(FDR < 0.05 & logFC > 0 ~ "Up",
                                  FDR < 0.05 & logFC < 0 ~ "Down",
                                  FDR > 0.05 ~ "No"))

table(SZC.COMBOvsCNT.COMBO_filt_df$Significance)

```
