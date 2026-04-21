# ResPred.AI

**An AI-powered platform for predicting antimicrobial resistance (AMR) and virulence factors directly from protein sequences.**

---

## Overview

Antimicrobial resistance (AMR) is one of the most pressing global health threats of our time. Virulence factors (VFs) compound this threat by increasing disease severity, making infections both harder to treat and more dangerous. Rapid, accurate computational prediction of resistance and virulence is therefore essential for effective surveillance, clinical decision-making, and infection control.

Most existing tools rely on sequence alignment and curated reference databases. While effective for known determinants, these approaches often fail to detect highly divergent, novel, or rapidly evolving resistant and virulent sequences — precisely the variants that pose the greatest surveillance gaps.

**ResPred.AI** addresses this limitation by employing supervised machine learning models trained on curated and enriched protein sequence datasets. By learning sequence-level patterns rather than relying on direct similarity, ResPred.AI can identify resistance and virulence determinants that fall outside the coverage of conventional alignment-based tools.

---

## Features

- **Resistance prediction** — Identifies antimicrobial resistance directly from protein sequences
- **Virulence prediction** — Detects virulence factors and provides functional annotations
- **AMR gene family classification** — Categorizes resistance determinants into established gene families
- **Virulence functional annotation** — Assigns functional roles to predicted virulence factors
- **Beyond alignment** — Captures resistance and virulence patterns in heterogeneous and divergent protein datasets where alignment-based tools fall short
- **Web interface** — Available as a ready-to-use online tool at [mgbio.asab.nust.edu.pk/respred](https://mgbio.asab.nust.edu.pk/respred/)

---

## Live Tool

Access ResPred.AI directly in your browser — no installation required:

🔗 **[https://mgbio.asab.nust.edu.pk/respred/](https://mgbio.asab.nust.edu.pk/respred/)**


## Input Format

ResPred.AI accepts standard FASTA files containing protein sequences:

```
>protein_id_1 [optional description]
MSKTTLVVKDEQTKKKDEEDFFSEKEGVSAEEIKAKMKQLEKDKQALEEQYEKMQKELK...
>protein_id_2
MNFKQILAAVSVALSGFSSIALATDLKTSAEKAIEDLTKQNPSDGSSVYALLQNLKPAK...
```

> **Note:** Input sequences must be amino acid (protein) sequences. Nucleotide sequences are not currently supported.

---

## Methodology

ResPred.AI uses supervised machine learning models trained on curated and enriched protein sequence datasets. The pipeline consists of two main prediction stages:

```
Input FASTA
    │
    ▼
Sequence Preprocessing & Feature Extraction
    │
    ▼
┌───────────────────────────────────────┐
│  Stage 1: Binary Classification       │
│  ├── AMR Predictor                    │
│  └── Virulence Factor Predictor       │
└───────────────────────────────────────┘
    │
    ▼
┌───────────────────────────────────────┐
│  Stage 2: Functional Annotation       │
│  ├── AMR Gene Family Classifier       │
│  └── VF Functional Annotator          │
└───────────────────────────────────────┘
    │
    ▼
Structured Output Reports
```

By capturing sequence-level patterns through machine learning — rather than relying on alignment to curated databases — ResPred.AI extends prediction coverage to divergent and novel determinants that conventional tools miss.

---

## License

This project is licensed under the Apache License 2.0.

---

## Contact

For questions, bug reports, or collaboration inquiries, please open a GitHub Issue or contact the development team.

🌐 [https://mgbio.asab.nust.edu.pk/respred/](https://mgbio.asab.nust.edu.pk/respred/)
