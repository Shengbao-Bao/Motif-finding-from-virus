# 病毒 RNA 候选调控 motif 发现

从公开病毒 RNA 序列中识别并排序候选调控 motif 的最小可运行原型。流程是 k-mer 全枚举、
一阶 Markov 背景期望、Poisson 右尾检验、BH-FDR 校正，最后综合显著性、覆盖范围与支持簇数进行排序。

## 目录

```
.
├─ README.md
├─ task.md                        任务原文
├─ run_all.ps1                    一键复现（含自测）
├─ examples/sample.fasta          离线演示输入：16 条流感 A 参考段，29.8 kb
├─ src/
│  ├─ cli.py                      入口：fetch | qc | motifs | known | all
│  ├─ fetch.py                    E-utilities 抓取、GenBank 解析、来源清单
│  ├─ qc.py                       规范化、过滤、去重、同源聚类、簇权重
│  ├─ stats.py                    Markov 背景、Poisson 尾、BH-FDR
│  ├─ motifs.py                   k-mer 富集、家族归并、排名输出
│  ├─ known_motifs.py             已知元件在同一排名中的位置
│  └─ common.py                   FASTA/TSV/md5 等 IO
├─ tests/test_stats.py            统计核心自测
├─ data/
│  ├─ raw/manifest.tsv            36 条记录：accession、物种、taxid、类型、注释、版本、md5
│  ├─ raw/raw.fasta               原始序列（QC 输入，可离线）
│  └─ clean/                      sequences.fasta、seqinfo.tsv、qc_report.tsv、clusters.tsv
└─ results/
   ├─ motifs.tsv                  候选排名（Top 50）
   ├─ known_motifs.tsv            已知 RNA 元件的判定（5 个可用短序列代理的元件）
   └─ summary.md                  结果速览
```

## 快速开始

```powershell
python src\cli.py all                                                          # 全流程
python src\cli.py motifs --fasta examples\sample.fasta --results-dir results\demo
python tests\test_stats.py                                                     # 统计自测
```

`all` 依次执行 fetch → qc → motifs → known；缓存存在时自动离线，`--refresh` 重新抓取。
`motifs --fasta` 可直接对任意 FASTA 运行，`--seqinfo` 可选（缺省时权重为 1，类群取 header 第二段）。

## 任务 1：数据获取与质量控制

### 来源与规模

用 NCBI Nucleotide，RefSeq 与 GenBank 并用，通过 E-utilities 的 esearch + efetch
（`rettype=gbwithparts`）取 36 条记录、323 255 bp：

类群          | 记录 | 分子类型 / 架构  | 长度           
SARS-CoV-2    | 6   | +ssRNA，单片段   | 29.9 kb       
Hepacivirus C | 8   | +ssRNA，单片段   | 9.0–9.7 kb    
HIV-1         | 4   |  +ssRNA         | 9.0–9.3 kb    
Zika virus    | 3   | +ssRNA          | 10.3–10.8 kb  
Influenza A   | 16  | −ssRNA，8 段     | 0.86–2.34 kb  

每条记录写入 `data/raw/manifest.tsv`：`accession.version`、`organism`、`taxid`、`mol_type`、
`topology`、`segment`、`length`、`gc_pct`、`release_date`（数据库版本）、`fetch_date`、序列 `md5`，
以及解析出的 `5'UTR`、`3'UTR`、`CDS` 区间。GenBank 原文不入库，
`fetch` 按 `src/fetch.py:GROUPS` 里的检索词重抓。

### 解析的字段
字段        | 来源                                          | 说明                           
序列类型    | LOCUS 行与 `source` 的 `/mol_type`            | genomic RNA、viral cRNA        
物种与分类  | `source` 的 `/organism`、`/db_xref="taxon:…"` | 按类群分层统计                   
分段        | `source` 的 `/segment`                        | 流感 8 个节段                   
已有注释    | `5'UTR`、`3'UTR`、`CDS` 的 location            | 支持 complement、join、partial 

### 质量控制
规则                         | 处理                                     | 记录字段            
非核酸记录                    | 丢弃                                     | drop_reason 
小写、空白、数字、间隔符       | 移除                                     | n_nonstd 
U → T                        | 替换，内部字母表统一为 ACGT               | n_u_to_t 
IUPAC 模糊碱基                | 映射到代表碱基，`--iupac drop` 可改为丢弃 | n_iupac 
N 占比 > 5% 或连续 N ≥ 10     | 丢弃                                     | n_frac 
长度 < 100 bp 或 > 50 kb     | 丢弃                                     | len_flag 
长度 Tukey IQR               | 默认关闭，`--iqr` 开启                    | len_flag 
完全重复（md5 相同）          | 保留 1 条，RefSeq 优先                    | duplicate_of 
近同源（16-mer Jaccard ≥ 0.9）| 单链接聚类，保留簇代表                    | cluster_id、cluster_size 
元数据缺失                    | 从 GenBank 补齐                          | meta_missing 

本次 36 条中剔除 1 条：`MN908947.3` 与 `NC_045512.2` 序列完全相同，判为完全重复。
IQR 默认关闭的原因是本数据集同时含 30 kb 冠状病毒与 0.9–2.3 kb 流感分段，
全局 IQR 假设长度单峰，会把流感分段判为离群；更合理的做法是按类群分别计算。

### 同源重复采样的偏差
公共库中同一毒株常有上万条几乎相同的基因组，若按记录统计，毒株特异的 k-mer 会被放大成
极显著富集。本项目的处理分三层：
1. 序列层：以 16-mer 集合的 Jaccard ≥ 0.9 判定同一谱系，单链接聚类，每簇只留一个代表。
   先用 `min|A|/max|A|` 做上界剪枝（Jaccard 不可能超过该值），再对幸存配对算精确交集。
2. 统计层：簇代表权重取 1，观测与期望用同一权重。簇成员之间的 `1/cluster_size` 只写进
   `qc_report.tsv` 表示冗余程度，不进入统计，否则 z 会被缩小 sqrt(cluster_size) 倍。
3. 结果层：输出支持该 motif 的独立簇数 `n_clusters_with` 与覆盖率 `coverage_pct`，
   用于判断显著性是否来自个别序列。

结果：35 条合并为 31 个簇（4 条 SARS-CoV-2 合一簇，2 条 HIV-1 合一簇），
代表序列合计 194 086 bp。统计单位是簇，不是记录。

## 任务 2：motif 发现算法

### 方法与统计量

k-mer 全枚举（k = 4…8），背景为一阶 Markov 模型，显著性用 Poisson 右尾检验，
多重检验用 Benjamini-Hochberg FDR。选 k-mer 的理由是能穷举（k = 8 时 65 536 种），
零假设期望可以精确计算，且显著 k-mer 可直接作为后续 PWM 的种子。

```
obs_w = Σ_i w_i · count_i(motif)                   加权观测次数
L_w   = Σ_i w_i · (L_i − k + 1)                    加权可用窗口数
exp_w = L_w · P_bg(motif)                          零假设期望
fold  = (obs_w + 0.5) / (exp_w + 0.5)              富集倍数
z     = (obs_w − exp_w) / sqrt(max(exp_w, 1e-9))   Poisson 近似 z
p     = P(X ≥ obs_w), X ~ Poisson(exp_w)           右尾 p
q     = BH(p)，假设族为该 k 的全部 4^k 个 k-mer     FDR
```

### 频率与覆盖范围

除观测次数与富集倍数外，每个候选还输出四项：支持它的独立簇数与占比
（`n_clusters_with`、`coverage_pct`）、跨类群数（`n_taxa_with`）、平均每簇拷贝数
（`mean_copies_per_cluster`）。这四项用来区分单拷贝的类群特异信号与广泛存在的重复或组成信号。

### 背景模型

默认 `--bg markov1`，`P(motif) = f(b0) · Π f(b_i | b_{i−1})`，条件概率用 α = 1 的 Dirichlet
先验平滑，使未出现过的二核苷酸也有非零概率。保留二核苷酸结构可以吸收密码子偏好，
只留下超出二核苷酸水平的超额；`--bg composition` 是零阶对照。

### 多重检验与冗余家族

- 每个 k 单独做 BH 校正，假设族是该 k 的全部 4^k 个 k-mer（k = 4…8 依次为 256、1 024、4 096、16 384、65 536，合计 87 296）；`motifs_all.tsv` 只写出实际出现过的 73 889 个。
- 短 k-mer 必然被长 k-mer 包含，排名会被同一家族占满，因此按排名顺序做子串与超串归并，
  只留家族代表，成员数记在 `note=family_size`。

### 输出

`results/motifs.tsv` 共 17 列：`rank`、`motif`、`k`、`gc_pct`、`obs_weighted`、`exp_weighted`、
`fold`、`log2_fold`、`z`、`p_value`、`q_value`、`n_clusters_with`、`coverage_pct`、`n_taxa_with`、
`mean_copies_per_cluster`、`family_of`、`note`；排名键为 q、fold、支持簇数。

### 结果

在输出的 73 889 个 k-mer 中，1 928 个 q 小于 0.05（假设族见上）。前 5 名：

| rank  | motif    | k | fold | q       | 簇 | 类群 |
| 1     | TTTTTTTT | 8 | 67.69| 6.6e-296| 6  | 2   |
| 2     | AAGAA    | 5 | 1.75 | 1.6e-34 | 31 | 5    |
| 3     | GAAGA    | 5 | 1.68 | 1.5e-25 | 29 | 5    |
| 4     | AGAAG    | 5 | 1.55 | 9.5e-18 | 29 | 5    |
| 5     | AAAGA    | 5 | 1.52 | 1.3e-17 | 31 | 5    |

第 1 名覆盖 19% 的簇、只跨 2 个类群，对应流感基因组末端的 poly(U) 与 poly(A) 摆动区。
第 2 至 5 名属同一家族（family_size 84、49、23、25），fold 只有 1.5 至 1.8 倍，
但覆盖 94% 到 100% 的簇、跨全部 5 个类群，是密码子与核苷酸组成偏好，不是具体调控位点。

### 候选与已验证元件的区别

`results/known_motifs.tsv` 只收录能用短序列代理的 5 个元件，并记录 DNA/RNA 核心序列、适用范围、
检测方式与是否依赖结构上下文：

| 元件                | motif | fold | q        | 同类排名 | 判定          
| 流感 poly(U) 摆动区 | TTTTTT | 5.70 | 1.1e-105 | 1       | 富集       
| Kozak 核心          | GCCACC | 1.31 | 0.132   | 513      | 与背景不可区分 
| polyA 信号          | AATAAA | 1.27 | 0.154   | 551      | 与背景不可区分 
| SARS-CoV-2 TRS-B   | ACGAAC | 0.99 | 1.0      | 2002     | 与背景不可区分 
| ARE 核心            | ATTTA | 0.75 | z = −3.2 | 912      | 显著缺失 

5 个里只有 1 个被判为富集，原因是富集衡量的是相对背景的超额，不是存在与否：
polyA 信号、TRS、Kozak 上下文的出现次数本来就接近组成期望，ARE 在病毒中常被负选择压低而表现为缺失，
只有局部成簇的 poly(U) 摆动区才显著。

结构依赖的元件不进入这张表：HIV-1 TAR 的功能来自整段茎环及其凸起，HCV IRES 来自多结构域折叠，
用其中任意一个六聚体作代理既定位不到元件，也不承载活性；这类元件只能靠协方差模型或结构探测识别，
因此放在任务 3 讨论。

两者的差别是结构性的。本工具的候选由滑窗得到，没有位置先验，不建二级结构，
证据是富集倍数与 FDR；已验证元件的证据是突变、回补、报告基因与结构探测。
因此候选只能表述为"值得进一步验证"，不能称为调控元件。

## 任务 3：结构与功能推断规划
### 3.1 结构与功能信息的引入
当前算法主要依赖 k-mer 富集程度识别候选 RNA 调控 motif，尚未考虑 RNA 二级结构、跨类群保守性及已知功能元件等信息。为提高候选 motif 的生物学可信度，后续计划从以下四个方面引入结构与功能证据。
（1）RNA 二级结构
使用 RNAfold-p预测 RNA 二级结构，获取最小自由能（MFE）及碱基配对概率。根据 motif 命中位点的结构状态，将其划分为 loop（环区）和 stem（茎区），分别计算富集程度，并结合碱基可及性（accessibility）评估候选位点的结构开放程度。

该分析可以通过计算完成，但预测结构与RNA在体内的真实构象可能存在差异。
（2）跨类群保守性
在现有类群覆盖率统计（`n_taxa_with`）的基础上，对同源簇内序列进行多序列比对（MSA），比较 motif 区域与侧翼区域的保守程度。对于茎环结构，还可以进一步分析配对碱基的协同变化（covariation）。

该分析可以通过序列比对与统计分析完成，但需要足够的同源序列支持。

（3）已知 RNA motif 与 RNA 家族
使用 Rfam 数据库及 `cmsearch` 检测 IRES、TAR、RRE、CRE、s2m 等已知 RNA 家族，并分析其与候选 motif 命中区域的重叠情况。

此外，可以利用 MEME 根据显著 k-mer 构建位置权重矩阵（PWM），再通过 FIMO 对病毒 RNA 序列进行扫描，以识别具有一定序列变异的潜在 motif。

该分析可以通过现有软件和数据库完成，检测效果取决于数据库覆盖范围及候选 motif 的质量。

（4）已有注释与实验数据
利用 GenBank 中的 5′UTR、3′UTR 和 CDS 注释开展区域特异性富集分析；进一步引入 SHAPE/DMS 结构探测数据、CLIP-seq 蛋白结合位点及 Ribo-seq 翻译效率等外部数据，评估候选 motif 的潜在功能。

其中，GenBank 注释可以直接用于计算分析；公开实验数据需要从相关数据库获取，部分功能证据则需要通过额外实验产生。

3.2 优先引入 RNA 二级结构信息
在上述四个方向中，RNA 二级结构是优先考虑的改进方向。
真实的 RNA 顺式调控元件通常具有一定的结构依赖性，例如位于茎环结构的环区、凸起区或其他结构开放区域。相比之下，密码子偏好和核苷酸组成偏差等信号通常不会表现出明显的结构区域特异性。
因此，可以利用 RNAfold 预测结果，将候选 motif 的命中位点按照结构状态进行分层统计，分别计算 loop 区域与 stem 区域的富集程度，并结合碱基可及性进行进一步筛选。
如果某个候选 motif 在 loop 区域表现出显著富集，而在 stem 区域不显著，同时具有较高的结构可及性，则说明其可能具有一定的结构依赖性，值得优先开展后续研究。

需要注意的是，这种结构分层结果只能作为候选功能的计算证据，不能直接证明该 motif 具有调控活性。

3.3 计算分析与实验验证的边界
本项目需要区分能够通过生物信息学分析获得的证据，以及必须通过实验才能确认的功能结论。

可通过计算分析获得的证据包括：
1. RNA 二级结构预测及候选位点的结构可及性。
2. 候选 motif 在不同病毒类群及同源序列中的保守性。
3. 候选区域与 Rfam 已知 RNA 家族的匹配情况。
4. 候选 motif 在 5′UTR、3′UTR 和 CDS 等不同区域的分布差异。
5. 候选 motif 与公开 RNA 结构探测、蛋白结合及翻译相关数据的关联。

需要实验验证的功能结论包括：
1. 候选 motif 是否真正影响病毒 RNA 的翻译效率或稳定性。
2. 候选 motif 是否能够被特定宿主蛋白或其他调控因子识别和结合。
3. RNA 在细胞内的真实二级结构是否与计算预测一致。
4. 候选 motif 对相关生物学功能是否具有必要性和充分性。
5. 候选 motif 的功能是否依赖特定宿主因子或细胞环境。

因此，本工具的输出定位为具有统计支持的候选调控 motif，而不是已经得到功能验证的 RNA 调控元件。

### 3.4 候选 motif 的可靠性分级
为避免将统计显著性直接等同于生物学功能，本项目采用分级证据体系评估候选 motif 的可靠性。

L0：统计支持
判定条件：
- `q < 0.05`。
- 独立支持簇数不少于 3。

该级别说明候选 motif 具有统计显著性，且不是由单条序列的偶然重复造成的。

L1：结构支持
判定条件：
- 候选 motif 命中位点位于 RNA 预测结构的开放区域。
- 具有相对较高的碱基可及性。
该级别说明候选 motif 具有潜在的结构依赖性。

L2：进化保守性支持
判定条件：
- `n_taxa_with ≥ 3`，即候选至少在 3 个病毒类群中出现。
- 在同源簇的多序列比对中，motif 所在位点在类群之间相互对应。
- 该位点比其侧翼更保守，且配对位点之间存在协同变化。

`n_taxa_with` 只反映跨类群的分布范围，出现不等于保守；位点对应关系与保守程度必须由 MSA 给出，
本工具尚未实现，只提供前一项。该级别说明候选 motif 可能具有跨类群保守的生物学意义。

L3：功能验证支持
判定条件：
- 获得公开数据库中的可靠功能注释；或
- 通过突变、回补、报告基因实验、结构探测等实验获得正向证据。
该级别说明候选 motif 获得了功能层面的支持。

其中，当前工具主要提供 L0 级统计证据，并输出支持后续 L1 和 L2 级分析所需的相关字段。L3 级证据需要依赖外部实验数据或已有功能验证结果。

需要强调的是，证据等级代表支持信息的类型和强度，不意味着所有候选 motif 都必须依次通过前一级才能获得更高级别的证据。


## 任务 4：工具实现与展示

可运行的程序是 `src/cli.py`，只用标准库，Python 3.8 以上。两条命令即可跑通：

```powershell
python src\cli.py all
python src\cli.py motifs --fasta examples\sample.fasta --results-dir results\demo
```

交付项                                                                                                     
1.FASTA 数据输入 
`examples/sample.fasta`（离线，16 条流感 A 段）；`data/clean/sequences.fasta`（QC 后）；`motifs --fasta` 可接任意 FASTA 
2.候选 motif 排名文件（TSV） `results/motifs.tsv`（Top 50，17 列） 


### 结果说明
`results/summary.md` 记录输入序列数、k-mer 假设数、显著数与 Top 表。
主要结论：榜首对应流感基因组末端的 poly(U) 摆动区，其余高排名项是密码子与核苷酸组成偏好；
5 个已知元件中只有 1 个被判为富集，说明单靠自背景富集无法替代功能注释。

### 局限
1. 自背景法只能识别组成离群值，天然高频的真实元件会被漏掉，被负选择压低的元件会显示缺失。
2. k-mer 没有错配容忍，一个 SNP 就打散 8-mer。
3. 不使用结构信息，同一茎环在不同病毒中序列不同，无法合并。
4. Poisson 的独立性假设在重复区被违反，k-mer 成簇出现时 p 值偏小，只能靠覆盖度指标间接识别。
5. 只做右尾检验，缺失只能从 z 判读，没有对应的 q 值。
6. 样本量小（31 个簇），未做按类群的独立富集矩阵与系统发育校正。

### 下一步
1. 前景与背景对比：以 3′UTR、5′UTR 为前景，CDS 为背景，使位置特异的元件显形。
2. 用显著 k-mer 的命中位点对齐建 PWM 并迭代，替代固定长度的 k-mer。
3. 接入 RNAfold，按 loop 与 stem 分层重算富集。
4. 用 Rfam cmsearch 与 MEME fimo 做交叉验证。
5. 若要下系统发育层面的结论，引入 PhyloGLM 或独立对比。
