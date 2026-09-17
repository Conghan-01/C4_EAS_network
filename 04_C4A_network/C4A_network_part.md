# Construction and Comparative Analysis of the C4A Seed-Co-expression Network
---
## Core Objectives

- **C4A 为核心的正相关和负相关的基因表达网络各自执行什么样的功能/Developmental Dynamics: Determine whether the C4A-centered network changes across developmental stages and undergoes functional transitions.**
- **C4A 为核心的网络是否随着发育阶段变化，是否发生功能转换/Developmental Dynamics: Determine whether the C4A-centered network changes across developmental stages and undergoes functional transitions.**
- **C4A 为核心的网络与神经精神疾病的关联，比如SCZ和BP，介导网络风险的发育阶段/Disease Associations: Investigate the association between C4A-centered networks and neuropsychiatric disorders (e.g., SCZ, BP), and identify the critical developmental stages mediating genetic risk.**

## Main Method

1. *Gene Expression Preprocessing*
2. *Compute correlations for C4A-associated genes to define gene sets at different R2 thresholds.*
3. *GO Functional Enrichment Analysis*
4. *MAGMA GWAS Heritability Enrichment: Assess heritability enrichment for SCZ (Schizophrenia) and BP (Bipolar Disorder).*

# Comparative Co-expression Network Analysis Framework
### Preprocessing Principle: The goal is not to force expression matrices from different stages to be identical. Instead, ensure high internal matrix quality for each stage, use consistent gene definitions and statistical pipelines, and preserve true biological differences across developmental stages.

prenatal RNA-seq
        |
        |-- QC
        |-- gene filtering
        |-- normalization
        |-- covariate correction
        |
        ↓
C4A network (prenatal)


postnatal RNA-seq
        |
        |-- QC
        |-- gene filtering
        |-- normalization
        |-- covariate correction
        |
        ↓
C4A network (postnatal)

raw counts
 |
 |-- TMM normalization
 |
 |-- log2 CPM
 |
 |-- gene filtering
 |
 |-- remove chrM chrY chrX
 |
 |-- quantile normalization
 |
 |-- ComBat batch correction
 |
 ↓
expr.qn.comb

expr.qn.comb
 |
 regression:
 sex
 RIN
 GW/Age
 PCs
 |
 residual
 |
 C4A correlation

