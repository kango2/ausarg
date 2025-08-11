
# 🧪 How to install Singularity Containers as Modules on Gadi
---

## 🧱 Overview

This guide walks you through:

- Installing Singularity containers as modules on **Gadi** using **SHPC** (Singularity HPC)
- Pulling containers from Docker Hub or Quay.io
- Locating executables within containers
- Setting up `container.yaml` files
- Installing and updating module aliases for reproducible workflows

---

## 🔧 What is SHPC?

SHPC is a utility that allows you to:

- Install Docker-based containers into Singularity
- Define module-compatible environments on HPC systems (like Gadi)
- Automatically manage metadata, executables, and version tags

---

## 📦 Steps to Install a Singularity Container (e.g. Ensembl-VEP)

### 1. Understand the Docker source

You need to identify:

- The Docker image path (e.g. `ensemblorg/ensembl-vep`)
- The image tags and their SHA256 digests
- The GitHub/tool documentation for command-line tools available inside the container

---

### 2. Create your YAML definition file

Navigate to:

```bash
/g/data/if89/shpcroot/registry/
```

Create a new folder for your tool:

```bash
mkdir ensembl-vep
cd ensembl-vep
```

Then create a file named `container.yaml`:

```yaml
docker: ensemblorg/ensembl-vep
url: https://hub.docker.com/r/ensemblorg/ensembl-vep
maintainer: '@ka6418'
description: VEP (Variant Effect Predictor) predicts the functional effects of genomic variants.

latest:
  release_114.0: sha256:e3dc8fe9303e5246e8483ba81974a2fcef768d2197c9b

tags:
  release_114.0: sha256:e3dc8fe9303e5246e8483ba81974a2fcef768d2197c9b
  release_112.0: sha256:e7612ab7c2923f2b9a78592b939e74874cd29f7494d70
  release_107.0: sha256:f252df581820048dbc90229d38fbc7b94738743bd9eb

```

> **Note:** Make sure the `docker:` path exactly matches the one from DockerHub/Quay.io/your container source.

---

### 3. Install the container

Set up your env

```bash
source /g/data/if89/shpcroot/pyvenv/bin/activate
module load python3/3.10.4
module load singularity
```

To install the default (latest) version:

```bash
shpc install ensembl-vep
```

To install a specific version:

```bash
shpc install ensembl-vep:release_114.0
```

This will:

- Download the Singularity `.sif` file to `/g/data/if89/shpcroot/containers/`
- Build a module entry under your SHPC environment

---

### 4. Find executables in the container

Run the container interactively to locate the paths of binaries:

```bash
singularity exec /path/to/ensembl-vep_*.sif which vep
```

If `which` fails, try:

```bash
singularity exec /path/to/container.sif find / -name 'vep' 2>/dev/null
```

Update the `aliases` section of your YAML file with these full paths.


```yaml
docker: ensemblorg/ensembl-vep
url: https://hub.docker.com/r/ensemblorg/ensembl-vep
maintainer: '@ka6418'
description: VEP (Variant Effect Predictor) predicts the functional effects of genomic variants.

latest:
  release_114.0: sha256:e3dc8fe9303e5246e8483ba81974a2fcef768d2197c9b

tags:
  release_114.0: sha256:e3dc8fe9303e5246e8483ba81974a2fcef768d2197c9b
  release_112.0: sha256:e7612ab7c2923f2b9a78592b939e74874cd29f7494d70
  release_107.0: sha256:f252df581820048dbc90229d38fbc7b94738743bd9eb

aliases:
    vep: /opt/vep/src/ensembl-vep/vep
    variant_recoder: /opt/vep/src/ensembl-vep/variant_recoder
    haplo: /opt/vep/src/ensembl-vep/haplo

```

Then re-run:

```bash
shpc install ensembl-vep
```

---

## ✅ Usage on Gadi

Once installed:

```bash
module load ensemblorg/ensembl-vep/release_114.0
vep --help
```

All commands in `aliases` are now accessible as modules.

If you can't find your modules use ``` module avail ```.
It's also necessary that you add ```module use -a /g/data/if89/shpcroot/modules``` to your bash profile.

---
## How to debug when things don't work

- If your container comes from a different source, chances are one of the containers from that source is already installed. In that case, a yaml file will exist which you can just copy and replace values
- Use manual singularity pull command to debug using the container outside shpc world. There have been cases with programs, such as MitoHifi, where the shpc module results in errors even when correctly installed but the program works perfectly when using it as a .sif image with singularity exec. 
- Containers can be manually pulled like, ``` singularity pull docker://ensemblorg/ensembl-vep:release_114.0 ```. This will save a .sif image in your current directory.
