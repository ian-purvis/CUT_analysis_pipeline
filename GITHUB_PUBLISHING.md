# Publishing This Pipeline to GitHub

Step-by-step guide for first-time GitHub publishers.

---

## Prerequisites

1. **Git** — Install from [git-scm.com](https://git-scm.com/download/win). During setup, choose "Git from the command line and also from 3rd-party software."
2. **GitHub account** — Sign up at [github.com](https://github.com) if needed.

---

## Step 1: Edit the LICENSE (Optional)

Open `LICENSE` and replace the copyright line with your name:

```
Copyright (c) 2025 Your Name
```

---

## Step 2: Create a New Repository on GitHub

1. Go to [github.com](https://github.com) and log in.
2. Click the **+** (top right) → **New repository**.
3. Fill in:
   - **Repository name**: e.g. `cutrun-pipeline` or `CUT-RUN-Tulloch-pipeline`
   - **Description**: e.g. "Snakemake pipeline for CUT&RUN analysis (Tulloch et al. eLife 2025)"
   - **Public**
   - **Do not** check "Add a README" (you already have one)
4. Click **Create repository**.

---

## Step 3: Initialize Git and Push (Command Line)

Open **Git Bash** (or WSL) and run:

```bash
# Navigate to this folder
cd "C:/Users/ophth/OneDrive - The University of Colorado Denver/IanPurvis_OneDrive/OneDrive - The University of Colorado Denver/Brzezinski Lab/Data/Datasets/Github_repo_CUT_analysis_pipeline"

# Initialize git
git init

# Add all files (respects .gitignore)
git add .

# First commit
git commit -m "Initial commit: CUT&RUN pipeline (Tulloch et al. eLife 2025)"

# Add your GitHub repo as remote (replace YOUR_USERNAME and REPO_NAME with yours)
git remote add origin https://github.com/YOUR_USERNAME/REPO_NAME.git

# Push to GitHub
git branch -M main
git push -u origin main
```

**Note:** Git will ask for your GitHub username and password. For password, use a **Personal Access Token** (Settings → Developer settings → Personal access tokens), not your account password.

---

## Step 4: Verify

1. Refresh your repo page on GitHub.
2. Confirm you see: Snakefile, README.md, scripts/, environment.yaml, etc.
3. Confirm you do **not** see: config.yaml, logs/, results/, .snakemake/

---

## What Gets Published (and What Doesn't)

| Published | Not Published (.gitignore) |
|-----------|---------------------------|
| Snakefile, scripts/, README | config.yaml (your local paths) |
| environment.yaml | logs/, results/ |
| run_*.sh, run_from_local.sh | .snakemake/ |
| config.yaml.example | __pycache__/ |

---

## Future Updates

After changing the pipeline:

```bash
git add .
git commit -m "Brief description of changes"
git push
```

---

## Alternative: GitHub Desktop

If you prefer a GUI:

1. Install [GitHub Desktop](https://desktop.github.com).
2. File → Add local repository → select this folder.
3. If it says "not a Git repository," choose "create a repository."
4. Publish to GitHub from the menu.
