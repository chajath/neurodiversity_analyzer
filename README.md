# 23andMe Neurodiversity SNP Analyzer

Analyzes downloaded 23andMe raw data for genetic variants (SNPs) associated with neurodevelopmental conditions from published genome-wide association studies (GWAS).

## ⚠️ Disclaimer

**This is NOT a diagnostic tool.** These conditions are highly polygenic and multifactorial. No SNP panel can predict or diagnose any of them. Consult a genetic counselor for interpretation.

## Supported Conditions

| Condition | SNPs Checked | Key Genes |
|-----------|-------------|-----------|
| **Autism (ASD)** | 22 | CNTNAP2, CDH9/CDH10, MET, EN2, OXTR, RELN, GABRB3, COMT, BDNF |
| **ADHD** | 22 | DAT1, DRD2, DRD4, DBH, NET, ADRA2A, CLOCK, FOXP2, CDH13 |
| **Dyslexia** | 12 | DCDC2, KIAA0319, DYX1C1, ROBO1, FOXP2, CMIP, ATP2C2 |
| **OCD** | 14 | SLC1A1, HTR2A, GRIN2B, DLGAP1, BTBD3, HTR2C, COMT |
| **Tourette Syndrome** | 8 | FLT3, NTN4, SLITRK1, HDC, TPH2, NRXN1 |
| **Bipolar Disorder** | 10 | CACNA1C, ANK3, FADS2, ODZ4, MTHFR, COMT |
| **Schizophrenia** | 12 | C4A/C4B, MIR137, CACNA1C, ZNF804A, NRG1, ERBB4, DISC1 |

## Requirements

- Python 3.10+
- No external dependencies

## Getting Your 23andMe Raw Data

1. Log in at [23andMe](https://you.23andme.com/)
2. Go to **Settings → 23andMe Data** (or browse to `https://you.23andme.com/tools/data/download/`)
3. Download your raw data file (a `.txt` file)

## Usage

```bash
# Analyze all conditions
python neurodiversity_analyzer.py /path/to/your_23andme_raw_data.txt

# Analyze a specific condition
python neurodiversity_analyzer.py /path/to/data.txt --condition adhd
python neurodiversity_analyzer.py /path/to/data.txt -c ocd
```

Available conditions: `asd`, `adhd`, `dyslexia`, `ocd`, `tourette`, `bipolar`, `scz`

The original ASD-only analyzer is also available:
```bash
python autism_snp_analyzer.py /path/to/your_23andme_raw_data.txt
```

## Output

The report includes:
- Per-condition SNP tables with genotype, risk allele status, and shared-condition tags
- Summary statistics per condition
- Flagged variants with literature references (PMIDs)
- **Cross-condition summary** comparing risk allele proportions across all panels
- **Shared SNP report** showing variants that appear in multiple condition panels

## Limitations

- Individual SNP effect sizes are very small (OR ~1.05–1.3)
- The most impactful variants (rare CNVs, de novo mutations) are **not detectable** from consumer genotyping chips
- Environmental and epigenetic factors are not considered
- Some SNPs may not be present on all 23andMe chip versions
- SNPs shared across conditions reflect common neurodevelopmental biology, not independent risks
