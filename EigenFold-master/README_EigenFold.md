
# 🧬 EigenFold (Student Replication)

This is a student replication of **EigenFold: Generative Protein Structure Prediction with Diffusion Models**, originally proposed by Bowen Jing et al. This project reproduces the inference pipeline for protein backbone structure prediction using a pretrained diffusion-based model.

---

## 📚 Paper Reference

> **EigenFold**: [arXiv:2304.02198](https://arxiv.org/abs/2304.02198)  
> Authors: Bowen Jing, Ezra Erives, Peter Pao-Huang, Gabriele Corso, Bonnie Berger, Tommi Jaakkola  
> Year: 2023

---

## 📦 Project Structure

```bash
EigenFold/
├── inference.py              # Main inference script
├── TMscore.exe               # Binary for evaluating structural similarity
├── pretrained_model/         # Contains epoch_7.pt and args.yaml
├── embeddings/               # Generated OmegaFold embeddings (*.npz)
├── structures/               # Ground-truth/reference PDB files
├── splits/                   # Input CSVs with sequence info
└── README.md                 # You're reading this!
```

---

## 🧪 Environment Setup

### 1. Clone the Repository

```bash
git clone https://github.com/<your-username>/eigenfold-student.git
cd eigenfold-student
```

### 2. Create and Activate Virtual Environment

```bash
python -m venv .venv
.venv\Scripts\activate  # On Windows
```

### 3. Install Required Packages

```bash
pip install torch==1.11.0+cu113 -f https://download.pytorch.org/whl/torch_stable.html
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv torch-geometric -f https://data.pyg.org/whl/torch-1.11.0+cu113.html
pip install biopython pandas matplotlib pyyaml
```

---

## 🚀 Inference Pipeline

### Step 1: Generate Embeddings

To generate OmegaFold embeddings for your sequences:

```bash
python make_embeddings.py --out_dir ./embeddings --splits splits/codnas.csv
```

This creates `.npz` files inside the `embeddings/` folder.

---

### Step 2: Run Inference

```bash
python inference.py \
  --model_dir ./pretrained_model \
  --ckpt epoch_7.pt \
  --pdb_dir ./structures \
  --embeddings_dir ./embeddings \
  --embeddings_key name \
  --elbo \
  --num_samples 5 \
  --alpha 1 \
  --beta 3 \
  --elbo_step 0.2 \
  --splits splits/codnas.csv
```

### 🧠 Notes:
- Output `.pdb` files and metrics are saved in `pretrained_model/` with a dynamic experiment folder name.
- Evaluation uses TMscore and optionally lDDT if the binary is installed.

---

## 📊 Evaluation Metrics

- **RMSD** (Root Mean Square Deviation)
- **TM-score** (Topological match)
- **GDT-TS / GDT-HA**
- **lDDT** (if `lddt` binary is present)

---

## 📚 Dataset

We used the **CODNAS** dataset for inference. It includes protein sequences with conformational diversity.

Each CSV in `splits/` must contain:

```csv
name,seqres,path
1xyz,ACDEFGHIKLMNPQRSTVWY,structures/1xyz.pdb
```

---

## 🙋‍♂️ Project Info

- **Project Name:** EigenFold Student Implementation
- **Course:** Data Science Methods and Tools (Spring 2025)
- **Instructor:** Prof. Akash Murthy
- **Institution:** Northeastern University
- **Authors:** Nitish Chowdary K, Rahul Bothra, Aditi Maurya

---

## 🧠 Contributions & Improvements

- Converted `inference.py` for Windows compatibility
- Handled `/dev/null` using `subprocess.DEVNULL`
- Used `tempfile.mkdtemp()` instead of hardcoded temp paths
- Fixed various Torch-Geometric module issues
- Added comments and cleaner logging

---

## 📄 Citation

```
@misc{jing2023eigenfold,
  title={EigenFold: Generative Protein Structure Prediction with Diffusion Models},
  author={Bowen Jing and Ezra Erives and Peter Pao-Huang and Gabriele Corso and Bonnie Berger and Tommi Jaakkola},
  year={2023},
  eprint={2304.02198},
  archivePrefix={arXiv},
  primaryClass={q-bio.BM}
}
```

---


