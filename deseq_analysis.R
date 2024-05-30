# 设置CRAN镜像
options(repos = c(CRAN = "https://cloud.r-project.org"))

# 安装和加载必要的R包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install("limma")

library("limma")
library("edgeR")

# 读取计数数据
countData <- read.table("result/expression/counts.txt", header = TRUE, sep = "\t", comment.char = "#", row.names = 1)

# 查看数据框的列名称
print(colnames(countData))

# 选择计数数据列
countData <- countData[, c("..bam.SRR15174659_sorted.bam", "..bam.SRR15174664_sorted.bam")]

# 创建样本信息数据框
group <- factor(c("untreated", "meropenem"))

# 创建DGEList对象
y <- DGEList(counts = countData, group = group)

# 计算规范化因子
y <- calcNormFactors(y)

# 创建设计矩阵
design <- model.matrix(~ group)

# 使用voom转换数据
v <- voom(y, design, plot = TRUE)

# 线性模型拟合
fit <- lmFit(v, design)

# 计算差异表达
fit <- eBayes(fit)

# 提取结果
res <- topTable(fit, coef = 2, number = nrow(countData))

# 查看结果
head(res)

# 保存结果
write.csv(res, file = "differential_expression_results_limma_voom.csv")
