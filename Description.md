# Summary of conditions for both NPCs and Neurons:
- Untreated control (at this point, I would exclude it from the analyses)
- Infected with GFP
- Infected with Mecp2

**The phases to be analyzed are S2S (euchromatin) and S3 (heterochromatin).**

# Tasks to be done:
1. Check how many chromatin regions change moving from phase S2 (euchromatin) to S3 (heterochromatin) and vice versa after Mecp2 overexpression compared to the GFP control (in both NPCs and Neurons).  
    From the preliminary analyses by Giosuè, very few regions changed their chromatin state, suggesting that Mecp2 overexpression does not radically alter the chromatin state.
2. Check the phase of the target genes of endogenous Mecp2 by comparing neurons and NPCs (only in GFP group).
3. Check the phase of genes where endogenous Mecp2 is enriched (FC>2, Neurons vs NPCs ) in NPCs and neurons (only in GFP-treated samples).
4. Check the phase of genes where exogenous Mecp2 is enriched (FC>2) in NPCs and neurons, and see if their chromatin state changes compared to the GFP control.


# Selection of Files for Analysis

Data directory: `./DATA/Sammy_Seq_fastq/`

Focus on **S2S (euchromatin) and S3 (heterochromatin) phases** in samples that are either **infected with GFP or Mecp2** (excluding untreated controls). 
Analyze **Neurons and NPCs** separately.
---

## **Files to Use for Analysis**

### **1. Exclude Untreated Controls**
- Exclude untreated control samples, **any file containing "CTRL" should be removed**.

### **2. Select Only S2S and S3 Phases**
- **Keep files containing `S2S` and `S3`**.
- **Exclude `S2L`, `S4`, or non-fractionated samples**.

### **3. Select Only GFP and Mecp2 Samples**
- **Keep files with `GFP` (control) and `M2` (Mecp2-infected)**.
- **Exclude files that do not contain either `GFP` or `M2`**.

### **4. Select NPCs and Neurons**
From the file names, **NSC refers to Neural Stem Cells (NPCs)**, and **Neu refers to Neurons**.

---

### **Final File Selection**
#### **Neurons (Neu)**
- **GFP-Infected (Control)**
  - `Neu1_GFP_S2S_*`
  - `Neu2_GFP_S2S_*`
  - `Neu3_GFP_S2S_*`
- **Mecp2-Infected**
  - `Neu1_M2_S2S_*`
  - `Neu2_M2_S2S_*`
  - `Neu3_M2_S2S_*`
- **S3 Phase (for both conditions)**
  - `Neu1_S3_*`
  - `Neu2_S3_*`
  - `Neu3_S3_*`

#### **NPCs (NSC)**
- **GFP-Infected (Control)**
  - `NSC1_GFP_S2S_*`
  - `NSC2_GFP_S2S_*`
  - `NSC3_GFP_S2S_*`
- **Mecp2-Infected**
  - `NSC1_M2_S2S_*`
  - `NSC2_M2_S2S_*`
  - `NSC3_M2_S2S_*`
- **S3 Phase (for both conditions)**
  - `NSC1_S3_*`
  - `NSC2_S3_*`
  - `NSC3_S3_*`

---

## **Breakdown by Task**
| Task | Required Files |
|------|--------------|
| **1. Chromatin state change (S2S ↔ S3) in Mecp2 vs GFP** | `Neu*_GFP_S2S_*`, `Neu*_M2_S2S_*`, `Neu*_S3_*`, `NSC*_GFP_S2S_*`, `NSC*_M2_S2S_*`, `NSC*_S3_*` |
| **2. Phase of endogenous Mecp2 target genes (GFP group only, Neurons vs NPCs)** | `Neu*_GFP_S2S_*`, `NSC*_GFP_S2S_*`, `Neu*_S3_*`, `NSC*_S3_*` |
| **3. Phase of genes enriched for endogenous Mecp2 (GFP group, FC>2 Neurons vs NPCs)** | `Neu*_GFP_S2S_*`, `NSC*_GFP_S2S_*`, `Neu*_S3_*`, `NSC*_S3_*` |
| **4. Phase of genes enriched for exogenous Mecp2 (FC>2) in NPCs and neurons (M2 group)** | `Neu*_M2_S2S_*`, `NSC*_M2_S2S_*`, `Neu*_S3_*`, `NSC*_S3_*` |

---