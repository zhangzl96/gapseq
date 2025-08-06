# 两种胃肠道细菌的交叉喂养

> <mark>:heavy_exclamation_mark: 请注意</mark>：本教程适用于gapseq版本≤v1.3.1。针对gapseq版本≥v1.4.0的更新版教程正在制作中。

### 背景

肠道细菌直肠真杆菌（*Eubacterium rectale*）已知能够在厌氧条件下利用乙酸盐作为能源，并生成丁酸盐作为代谢终产物（[Rivère *et al.* (2015) Appl Envrion Microbiol](https://pubmed.ncbi.nlm.nih.gov/26319874/)）。乙酸盐是多种其他肠道细菌（包括双歧杆菌属，如长双歧杆菌（*Bifidobacterium longum*））的常见发酵终产物。本教程中，将使用**gapseq**工具重建直肠真杆菌和长双歧杆菌的基因组规模代谢模型，随后模拟两者的共培养过程并探究其相互作用。

*注：以下所有由命令生成的中间文件均存储在[GitHub仓库](https://github.com/Waschina/gapseq.tutorial.data)中，若您希望从教程的后续步骤开始而非从头操作，可下载或克隆该仓库。*

### 输入

- 基因组文件：

  - *Eubacterium rectale* ATCC 33656

    RefSeq: `GCF_000020605.1`

  - *Bifidobacterium longum* NCC2705: 

    RefSeq: `GCF_000007525.1`

- 生长培养基文件：`gf_medium.csv` 

  这本质上是一种葡萄糖和乙酸盐的最低限度培养基。基于根据[GapMind](http://papers.genomics.lbl.gov/cgi-bin/gapView.cgi)的预测（[Price et al. (2019) mSystems](https://doi.org/10.1101/741918 )），这两种微生物很可能具备原养型特性（能够自主合成所有蛋白质合成所需的氨基酸），因此培养基中未额外添加氨基酸。

  E. rectale: [View Gapmind results](http://papers.genomics.lbl.gov/cgi-bin/gapView.cgi?orgs=NCBI__GCF_000020605.1&set=aa); B. longum: [View Gapmind results](http://papers.genomics.lbl.gov/cgi-bin/gapView.cgi?orgs=NCBI__GCF_000007525.1&set=aa)

### 准备工作

下载基因组文件和gapfill培养基，重命名文件。

```sh
#!/bin/bash

# Download genome assemblies 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/605/GCF_000020605.1_ASM2060v1/GCF_000020605.1_ASM2060v1_protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/525/GCF_000007525.1_ASM752v1/GCF_000007525.1_ASM752v1_protein.faa.gz

# Download gapfill-medium file
wget https://github.com/Waschina/gapseq.tutorial.data/raw/master/CF_eure_bilo/gf_medium.csv

# Rename genomes to "eure" (E. rectale) and "bilo" (B. longum) 
mv GCF_000020605.1_ASM2060v1_protein.faa.gz eure.faa.gz
mv GCF_000007525.1_ASM752v1_protein.faa.gz bilo.faa.gz
```

### gapseq重建 

Now we have the genome sequences and a gapfill medium. That is all we need. Lets reconstruct models:

```sh
#!/bin/bash

modelA="eure"
modelB="bilo"

# Reaction & Pathway prediction
gapseq find -p all -b 200 -m Bacteria $modelA.faa.gz
gapseq find -p all -b 200 -m Bacteria $modelB.faa.gz

# Transporter prediction
gapseq find-transport -b 200 $modelA.faa.gz 
gapseq find-transport -b 200 $modelB.faa.gz

# Building Draft Model - based on Reaction-, Pathway-, and Transporter prediction
gapseq draft -r $modelA-all-Reactions.tbl -t $modelA-Transporter.tbl -p $modelA-all-Pathways.tbl -c $modelA.faa.gz -u 200 -l 100
gapseq draft -r $modelB-all-Reactions.tbl -t $modelB-Transporter.tbl -p $modelB-all-Pathways.tbl -c $modelB.faa.gz -u 200 -l 100

# Gapfilling
gapseq fill -m $modelA-draft.RDS -n gf_medium.csv -c $modelA-rxnWeights.RDS -g $modelA-rxnXgenes.RDS -b 100
gapseq fill -m $modelB-draft.RDS -n gf_medium.csv -c $modelB-rxnWeights.RDS -g $modelB-rxnXgenes.RDS -b 100
```

最终模型存储为R-目标文件：`eure.RDS`和`bilo.RDS`，可以直接使用`readRDS()`命令加载到R语言中。


### 群落模拟

在此，将使用R包`BacArena`对*B. longum*和*E. rectale*的共生长进行基于基质的模拟。以下代码块展示了用于简单群落代谢模拟的R源代码。

```R
# Load R-packages
library(BacArena)
library(data.table)

# Load reconstructed models
er <- readRDS("eure.RDS") # E. rectale
bl <- readRDS("bilo.RDS") # B. longum

# Small fix to D/L-Lactate secretion (*) and model names
bl <- rmReact(bl, react = "EX_cpd00221_e0")
er@mod_desc <- "E. rectale"
bl@mod_desc <- "B. longum"

# Construct the organism objects for BacArena simulations
eure <- Bac(er)
bilo <- Bac(bl)

# Construct the arena size 10x10 grid cells
arena <- Arena(n = 10, m = 10)

# For each organism, populate randomly 2 grid cells in the Arena as 
# 'starter culture'
arena <- addOrg(arena, eure, amount = 2)
arena <- addOrg(arena, bilo, amount = 2)

# add substrates to arena
arena_subs <- fread("gf_medium.csv") # same as gapfill medium
arena_subs[, ex.rxn := paste0("EX_", compounds, "_e0")]

arena <- addSubs(arena, smax = arena_subs$maxFlux, 
                 mediac = arena_subs$ex.rxn, unit = "mM", addAnyway = T)
# Remove acetate from initial substrate list to see effect of Cross-Feeding
arena <- rmSubs(arena, mediac = "EX_cpd00029_e0") 

```

（*gapseq工具常预测：若生物体将乳酸作为发酵终产物，其最优代谢路径可能同时涉及两种对映异构体——D-乳酸和L-乳酸。为便于绘图和分析，我们强制禁止D-乳酸的产生，以确保观测到的乳酸产量仅以单一代谢物L-乳酸的形式呈现。此操作仅影响数据可视化方式，不会改变模拟结果的本质。*）

现在已经准备好进行实际的群落模拟并绘制结果：

```R
# Simulation for 13 time steps
CF_sim <- simEnv(arena, time = 13, sec_obj = "mtf")

# Plot levels of Acetate, Buyrate, and Lactate as well as growth
par(mfrow = c(1, 2))
plotCurves2(CF_sim, legendpos = "topleft",
            subs = c("cpd00211_e0", "cpd00029_e0", "cpd00159_e0"),
            dict = list(cpd00211_e0 = "Butyrate", 
                        cpd00029_e0 = "Acetate", 
                        cpd00159_e0 = "Lactate"))
```

![](https://github.com/Waschina/gapseq.tutorial.data/raw/master/CF_eure_bilo/CF_eure_bilo.png)

模拟预测显示，在直肠真杆菌（*E. rectale*）与长双歧杆菌（*B. longum*）共生长过程中，乙酸盐、丁酸盐和乳酸盐会被产生。

接下来，我们需验证部分发酵产物是否被其中一种微生物部分消耗。这需要几行数据处理代码的辅助：

```R
# Lets get the exchange fluxs for each grid cell at time step 11
dt.cf <- CF_sim@exchangeslist[[11]]

# Limit output to Acetate, Butyrate, and Lactate
dt.cf <- as.data.table(dt.cf[,c("species",
                                "EX_cpd00029_e0",
                                "EX_cpd00211_e0",
                                "EX_cpd00159_e0")])

# Rename column names (this is just aestetics)
dt.cf <- dt.cf[,.(species, 
                  Acetate = EX_cpd00029_e0, 
                  Butyrate = EX_cpd00211_e0,
                  Lactate = EX_cpd00159_e0)]

# Wide-To-Long table transformation
dt.cf <- melt(dt.cf, 
              id.vars = "species", 
              variable.name = "Metabolite", 
              value.name = "Flux")
dt.cf <- dt.cf[!is.na(Flux)] # rm NA entries (no exchange reaction in grid cell)

# Sum exchanges for each species and metabolite over all 100 grid cells
dt.cf[, .(summed.Flux = sum(Flux)), by = c("species", "Metabolite")]
```

输出信息：

```R
      species Metabolite summed.Flux
1: E. rectale    Acetate  -1066.7647
2:  B. longum    Acetate   1925.5036
3: E. rectale   Butyrate   2694.5406
4:  B. longum   Butyrate      0.0000
5:  B. longum    Lactate    806.4733
```

我们可以看到，长双歧杆菌（*B. longum*）将乙酸盐（Acetate）和乳酸盐（Lactate）作为主要代谢终产物分泌。其中约55%（1066.76/1925.50=0.554）的乙酸盐被直肠真杆菌（*E. rectale*）消耗。此外，长双歧杆菌还额外产生乳酸盐，而直肠真杆菌则分泌丁酸盐（Butyrate）。