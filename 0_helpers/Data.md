### **Understanding the Data Files:**
The dataset consists of **FASTQ** files which are raw sequencing data files.

From the paper:
1. **4f-SAMMY-seq** (4 fractions Sequential Analysis of Macromolecules accessibility sequencing) involves **S2, S3, and S4 chromatin fractions**, which are processed into **S2S (< 300bp) and S2L (> 300bp) fractions**.
2. **RNA-seq** involves **total RNA extraction** followed by sequencing.
3. **ChIP-seq** (Chromatin Immunoprecipitation sequencing) is done with specific antibodies for histone marks like **H3K9me3, H3K27ac, and H3K4me3**.

---

### **Grouping the Files into Clusters:**
#### **Cluster 1: 4f-SAMMY-seq Data**
These files correspond to chromatin fraction sequencing using the 4f-SAMMY-seq technique. The fractions (S2, S2S, S2L, S3, S4) reflect different chromatin states.

- **Files starting with "4F_" and "3F_" in `Sammy_Seq_fastq/`**
- **Files in `Sammy_Seq2_fastq/fastq` and `fastq2/` with S2S, S2L, S3, or S4 in the names**
- Example:  
  - `4F_Neu1_S2S_A1_S1_L001_R1_001.fastq.gz`
  - `3F_Neu1_M2_S3_H1_S8_L001_R1_001.fastq.gz`
  - `Neu1_M2_4F_S2S_G1_S43_L001_R1_001.fastq.gz`
  - `Neu1_CTRL_S2S_A1_S36_L001_R1_001.fastq.gz`

---

#### **Cluster 2: RNA-seq Data**
RNA-seq files contain transcriptome sequencing data.

- **Files with "RNA" or "GFP" in the name**
- **Files in `Sammy_Seq_fastq/` where no fractionation (S2, S3, S4) is specified**
- Example:
  - `Neu1_GFP_S2L_4F_D11_S20_L001_R1_001.fastq.gz`
  - `Neu1_GFP_S3_4F_E11_S21_L001_R1_001.fastq.gz`
  - `NSC1_GFP_S2L_4F_B4_S32_L001_R1_001.fastq.gz`

---

#### **Cluster 3: ChIP-seq Data**
ChIP-seq is used to analyze protein-DNA interactions and histone modifications.

- **Files corresponding to chromatin-associated proteins (e.g., H3K9me3, H3K27ac, H3K4me3)**
- **Potentially found in the `Sammy_Seq_fastq/` dataset, especially for NSC (Neural Stem Cell) samples**
- Example:
  - `NSC1_GFP_S3_B4_S40_L001_R1_001.fastq.gz`
  - `NSC1_M2_S3_A4_S41_L001_R1_001.fastq.gz`

---

#### **Cluster 4: Controls (Non-Treated or Wild-Type Samples)**
These files likely represent **wild-type (WT) or untreated control samples** used for comparison.

- **Files containing “CTRL” in the name**
- Example:
  - `Neu1_CTRL_S2S_A1_S36_L001_R1_001.fastq.gz`
  - `Neu2_CTRL_S2S_A6_S30_L001_R1_001.fastq.gz`
  - `Neu3_CTRL_S2S_A3_S33_L001_R1_001.fastq.gz`

---

### **Summary of Clusters**
| Cluster Name         | File Naming Pattern | Description |
|----------------------|--------------------|-------------|
| **4f-SAMMY-seq** | `4F_*`, `3F_*`, `S2S`, `S2L`, `S3`, `S4` | Chromatin fraction sequencing |
| **RNA-seq** | `GFP`, `RNA` | Transcriptome sequencing |
| **ChIP-seq** | `H3K9me3`, `H3K27ac`, `H3K4me3` | Chromatin immunoprecipitation |
| **Controls (WT/Untreated)** | `CTRL` | Baseline comparison |