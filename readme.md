# 标题：

作者：

肖瑶	    521010910058

赵苑羽	521141910004



## 摘要

本流程利用RNA-seq对肺炎克雷伯菌（Klebsiella pneumoniae）在庆大霉素处理前后的基因表达差异进行了分析。通过比较经过庆大霉素处理与未经过处理的六组样本，揭示了基因表达水平的显著变化，旨在理解抗生素对细菌基因表达的影响。工作流程主要包括利用anaconda构建计算环境，进行质量控制、序列剪切、比对等RNA序列上游分析，进行RNA序列下游分析筛选显著差异表达的基因和展示基因表达的整体变化趋势。

## 前言

肺炎克雷伯菌（Klebsiella pneumoniae）是一种重要的病原菌，常导致医院感染及社区获得性感染，具有较强的抗药性。庆大霉素是一种广泛使用的氨基糖苷类抗生素，对多种细菌具有强效的抑制作用。RNA-seq是一种高通量测序技术，用于研究转录组的复杂性和动态变化。作为一种强大的工具，RNA-seq能够同时测量数千个基因的表达水平，提供关于基因表达、可变剪接、新的转录本和转录起始位点的信息。本流程通过比较庆大霉素处理前后的基因表达水平，揭示基因表达的动态变化，从而帮助理解基因调控机制。

## 数据与方法

### 分析工具

1. SRA Toolkit：一系列用于检索、转换、处理和分析来自Sequence Read Archive（SRA）数据的工具。在此工作流程中，使用了`prefetch`下载NCBI中的SRR数据，并用`fasterq-dump`将SRA数据转换为fastq格式。
2. Trim Galore：包含`FastQC`，用于高通量序列数据的质量控制程序，使用fastqc对下载的六组SRR数据进行质量分析。且可用于对fastq格式的文件进行序列剪切，去除低质量末端和接头序列。
3. BWA：一个用于DNA序列比对的软件包，它主要用于将短读序列比对到参考基因组上，利用BWT算法和suffix array实现快速匹配。
4. Samtools：用于将`.sam`和`.bam`格式的文件相互转换，利用`.bam`格式的文件提高计算效率，且可用于构建`.bam`文件的索引。
5. Featurecounts：利用参考基因组注释文件，对感兴趣的基因组特征（通常为外显子、基因、启动子等区域）进行记数，用于下游分析。
6. R：利用DESeq2、EnhancedVolcano、pheatmat、ggplot2等R包，通过输入的基因表达量矩阵进行计算，得到差异基因表达计算结果，并进行绘图。

### 原始数据

1. RNA测序结果原始数据：共六组，分别为经过或未经过庆大霉素处理的肺炎克雷伯菌RNA测序结果。下载自美国国家生物技术信息中心（NCBI）数据库，样本编号：SRR15174659，SRR15174661，SRR15174663，SRR15174670，SRR15174671，SRR15174672，均来自于项目PRJNA747165（GEO: GSE180237）。数据信息如下：

   | 编号        | 物种                  | 处理情况   | 重复序号 |
   | ----------- | --------------------- | ---------- | -------- |
   | SRR15174659 | Klebsiella pneumoniae | untreated  | rep1     |
   | SRR15174661 | Klebsiella pneumoniae | untreated  | rep2     |
   | SRR15174663 | Klebsiella pneumoniae | untreated  | rep3     |
   | SRR15174670 | Klebsiella pneumoniae | gentamicin | rep1     |
   | SRR15174671 | Klebsiella pneumoniae | gentamicin | rep2     |
   | SRR15174672 | Klebsiella pneumoniae | gentamicin | rep3     |

2. 参考基因组数据：用于RNA测序数据的比对，包括参考基因组序列（`.fna`）和基因组注释（`.gtf`），下载自NCBI，序号为GCF_000694555.1。

### 计算环境

1. rna_seq环境：主要用于数据下载和上游分析。添加上海交通大学anaconda镜像后配置该环境。conda创建rna_seq虚拟环境，随后通过conda安装FastQC、BWA、Samtools、SRA Toolkit、Trim Galore、Subread。
2. r_env环境：用于运行R脚本，分析基因表达差异，并绘图。依次安装基础R包，以及DESeq2、EnhancedVolcano、pheatmat、ggplot2等差异分析需要的包。

### 分析流程

* RNA序列上游分析
  * 构建工作目录
  * 利用prefetch下载六组RNA测序数据，并利用fasterq-dump将格式转换为`.fastq`
  * 下载参考基因组序列与基因组注释文件，并解压
  * 将两个参考基因组文件利用fasterq-dump转换为`.fastq`
  * 利用FastQC检测六组RNA序列的质量，生成质量分析报告
  * 利用Trim Galore对原始序列进行剪切，去除低质量序列与接头
  * 利用BWA将六组测序数据对应到参考基因组序列上，生成`.sam`比对文件
  * 利用Samtools将`.sam`文件转化为`.bam`文件，进行排序，并构建索引
  * 利用FeatureCounts，对排序后的`.bam`文件进行基因表达记数，得到表达量矩阵
  * 去除中间文件
* RNA序列下游分析
  * 将表达量矩阵进行处理，去除不需要的列，作为后续R脚本的输入
  * 运行R脚本，加载所需要的包
  * 输入原始数据，设置样本信息，构建DataFrame
  * 使用表达量数据和样本信息，构建DESeq数据集
  * 数据预处理，包括去除总数过低的基因，利用DESeq函数将数据归一化，并进行差异分析
  * 得到差异分析结果，进行排序
  * 设置参数，绘制火山图
  * 筛选显著差异表达的基因，绘制热图
  * 绘制MA图
  * 输出差异分析结果的数据和图

## 分析结果



## 讨论



## 参考文献

[1] Peijun Ma, Haley M. Amemiya, Lorrie L. He, et al. Bacterial droplet-based single-cell RNA-seq reveals antibiotic-associated heterogeneous cellular states. *Cell*, 186(4), pp 877-891.

[2] Stark, R., Grzelak, M. & Hadfield, J. RNA sequencing: the teenage years. *Nat Rev Genet* 20, 631–656 (2019).

[3] Conesa, A., Madrigal, P., Tarazona, S. *et al.* A survey of best practices for RNA-seq data analysis. *Genome Biol* 17, 13 (2016).



## 附录

核心代码和脚本文件地址：[Cannizzaro-reaction/RNA-seq_workflow: course project for BIO2503 (github.com)](https://github.com/Cannizzaro-reaction/RNA-seq_workflow)