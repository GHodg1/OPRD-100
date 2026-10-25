# OPRD-100 Data Analysis

This directory contains clean (using claude Sonnet 4.5), Jupyter notebooks and associated code for analyzing the OPRD-100 dataset.

## Directory Structure

```
data_analysis/
├── notebooks/           # Jupyter notebooks for analysis
│   ├── dataset_analysis.ipynb           # OPRD-100 dataset analysis
│   └── reaxys_comparison_clean.ipynb    # OPRD-100 vs Reaxys comparison
├── src/                 # Python source files
│   ├── PlottingFunctions.py             # Visualization utilities
│   ├── ReaxysComparison.py              # Reaxys data extraction
│   ├── Stereochemistry.py               # Stereochemistry analysis
│   ├── MetalCatalysts.py               # Metal catalyst identity and role evidence
│   ├── StereoChem.py                    # Stereochemistry utilities
│   ├── reaction_networks.py             # Reaction network generation
│   ├── molecule_presence_checker.py     # PubChem CAS number lookup
│   └── CHEM21_solvent_dict.json         # Solvent classification data
├── output/              # Generated outputs
│   ├── analysis_figures/                # Figures from dataset_analysis.ipynb
│   ├── reaxys_comparison_figures/       # Figures from reaxys_comparison_clean.ipynb
│   └── reaxys_comparison_data/          # CSV/Excel files from Reaxys comparison
├── data/                # Input data files
│   ├── OPRD-100.json                    # Main dataset file
│   └── Reaxys_data/
│       └── reaxys_100_articles_data.xml # Reaxys XML data
└── README.md            # This file
```

## Notebooks

### 1. dataset_analysis.ipynb

Comprehensive analysis of the OPRD-100 dataset covering:
- **Data Locations**: Where reaction data appears in articles (schemes, tables, experimental sections)
- **Yield Types**: Distribution of different yield metrics (isolated, assay, ee, dr, etc.)
- **Yield Distributions**: Statistical analysis of yields by data location
- **Stereochemistry**: Analysis of stereogenic centers gained/lost/preserved
- **Solvent Usage**: Frequency, volume, and CHEM21 environmental classification
- **Reaction Conditions**: Temperature and time distributions
- **Catalyst Usage**: Unique reaction entries per metal, stacked by reported,
  likely, and uncertain catalyst evidence. The analysis reads the repository's
  root `data/OPRD-100.json`; the source data is never modified.

### Metal catalyst classification

`MetalCatalystAnalyzer.analyze(data, output_path)` writes its evidence JSON
directly and returns the same report. RDKit identifies atoms in SMILES;
bounded formula/name tokens and reviewed aliases cover coordination shorthand,
supported catalysts, and common metal names. Acetyl/propyl abbreviations,
supplier codes, lamp wattages, and recognised ferrocene ligands are handled
separately to avoid spurious elements.

The classifier combines explicit catalytic flags, mol%/equivalent loadings
(including loadings recorded under `Concentration`), substrate/catalyst ratios,
and usual reagent roles. Supported-metal wt% is not interpreted as mol%.
**Reported** means an explicit catalytic-amount flag or an explicit positive
substoichiometric loading below 100 mol% (including equivalent amounts and S/C
ratios), after counterion, known noncatalytic-role, and conflict checks.
Thus, `5 mol%` or `0.05 equiv` can establish the reported category without a
`Catalytic Amount` flag. This category describes reported catalytic-scale use,
not mechanistic proof; low loadings can also describe additives.
Counterions and noncatalytic reagents remain in the audit but are not plotted.
Metals in reaction SMILES are recorded separately, never automatically treated
as catalysts. Original identities and amount fields are preserved.

The chart uses `Plot.add_stacked_bar`, counting each dataset entry once per
metal at its strongest classification across all steps: reported > likely >
uncertain. Identical reaction SMILES in different entries remain separate
experimental records. Uncertain bars are candidates, not established catalysts;
the classes are heuristics, not calibrated probabilities. Unresolved names,
missing identities, conflicting loadings, and unknown roles require review.
Formula-token matching is not a complete chemical name parser, and low loading
can also describe an additive, so automated coverage/precision is not guaranteed.

#### Supplier-backed identity registry

`KNOWN_CATALYSTS` in `src/MetalCatalysts.py` is a small, offline, curated registry
of specific compounds and named families, reviewed on 2026-10-25. It is not an
automatically scraped supplier catalogue. Each entry has explicit aliases, its
active metal, an application summary, and a supplier source:

| Identity | Metal | Supplier evidence |
| --- | --- | --- |
| Copper(I) iodide | Cu | [Thermo Fisher: organic synthesis catalyst](https://www.thermofisher.com/order/catalog/product/011606.A1) |
| Scandium(III) triflate | Sc | [Thermo Fisher: Lewis acid catalyst](https://www.thermofisher.com/order/catalog/product/040566.03) |
| Bismuth(III) triflate | Bi | [Thermo Fisher: acylation/substitution catalysis](https://www.thermofisher.com/order/catalog/product/L19687.06) |
| Zinc triflate | Zn | [Thermo Fisher: thioketalization/silylation catalysis](https://www.thermofisher.com/order/catalog/product/L15969.06) |
| Grubbs metathesis catalysts | Ru | [Umicore: product/formula catalogue](https://pmc.umicore.com/en/products/grubbs-catalysts/) |
| Hoveyda-Grubbs metathesis catalysts | Ru | [Umicore: Hoveyda-type catalogue](https://pmc.umicore.com/en/products/grubbs-catalysts/) |
| Wilkinson's catalyst | Rh | [Thermo Fisher: homogeneous hydrogenation](https://www.thermofisher.com/order/catalog/product/010468.06) |
| Lindlar catalyst | Pd | [Thermo Fisher: lead-poisoned Pd on calcium carbonate](https://www.thermofisher.com/order/catalog/product/043172.14) |

The registry additionally covers 21 palladium identities/formulations:

- Precursors: Pd(OAc)2, Pd2(dba)3, PdCl2, and Pd(TFA)2.
- Phosphine complexes: Pd(PPh3)4, PdCl2(PPh3)2, PdCl2(dppf),
  PdCl2(dppf).CH2Cl2, PdCl2(dtbpf), and PdCl2(AmPhos)2.
- Supported catalysts: Pd/C and palladium hydroxide on carbon (Pearlman's).
- Buchwald precatalysts: SPhos/XPhos/RuPhos G2 and G3, tBuXPhos G3, and
  Xantphos G3, including the verified Umicore CX121/122/123/131/132/133/135/136
  trade codes from the [Buchwald catalogue](https://pmc.umicore.com/en/products/buchwald/).
- PEPPSI-IPr, documented in the
  [Sigma-Aldrich catalogue](https://www.sigmaaldrich.com/US/en/product/aldrich/669032).

Each record has its own supplier URL in the registry/evidence JSON. Ferrocene
ligand Fe is retained as ligand metadata, not counted as an Fe catalyst. Solvated
and unsolvated dppf complexes remain separate identities. Incomplete shorthand
such as `Pd(dppf) DCM` or `PdCl2(amphos)` is not silently assigned a specific
supplier compound. Free RuPhos is a ligand, not a ruthenium catalyst.

#### Catalogue/API access checked on 2026-10-25

| Service | Observed access | Use in this analysis |
| --- | --- | --- |
| [Chemspace OpenAPI](https://api.chem-space.com/docs/) | Schema accessible; documented `POST /v5/search/text` returned HTTP 401 without credentials | Catalogue search requires an API key; no authenticated results claimed |
| [Molport OpenAPI](https://docs.molport.com/openapi) | Schema accessible; availability search explicitly requires `X-API-Key` | No authenticated catalogue results claimed |
| [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) | Public CAS/synonym and molecular-formula requests returned HTTP 200 for reviewed records | Identity cross-checks, **not evidence of catalytic use or current supplier stock** |
| Umicore / Thermo Fisher / Sigma-Aldrich | Public supplier product/portfolio pages | Known catalytic applications and trade-name provenance |

Seven Pd records include reviewed PubChem CIDs and identity API URLs alongside
supplier sources. API synonyms are not ingested wholesale: one PEPPSI-IPr CAS
lookup returned a hydrogenated SIPr identity, so that response was rejected;
XPhos G2 CAS lookup returned 404 and has no API identity claim. Their supplier
records remain usable independently. Catalogue identities still need chemical
review; API access alone does not establish a catalytic role.

Aliases match whole normalized names (case/spacing and leading catalytic loading
labels), never arbitrary substrings. The registry can resolve names that do not
spell out their metal, such as Wilkinson's or Grubbs catalysts. It does not
generalize CuI to all copper salts, or infer that a ligand is itself a catalyst.
Supplier recognition supports **likely**, not **reported**, catalysis when the
reaction provides no loading. Explicit noncatalytic flags, stoichiometric loading,
and conflicts take precedence. Known catalytic use does not prove the role in a
particular reaction. Runtime classification requires no network calls.

Evidence JSON schema version 3 records the expanded reported-category criteria
under `method.reported_criteria`. It also includes a registry snapshot under
`method.supplier_registry` and `supplier_evidence` on each assignment (null for
unmatched identities), so the supplier basis remains auditable even when a
reaction-specific flag overrides it. Provenance includes CAS numbers, reviewed
PubChem CIDs/API URLs where available, and ligand-metal metadata. The summary
also reports the number of supplier-matched reagent assignments.

Outputs in `output/analysis_figures/`:

- `metal_catalyst_evidence_OPRD-100.json`: assignments with evidence and indices,
  confidence counts, reaction-component metals, review queue, and unresolved names.
- `catalyst_usage_OPRD-100.png`: confidence-stacked top-10 chart at `figsize=(5, 3)`,
  ranked by total stacked count, without a title or omitted-metal annotation.
  Full counts remain in the evidence JSON. Axes, ticks, and legend use the
  notebook's shared 14-point styling for composite-figure readability.
- `catalyst_reactions_OPRD-100.json`: deduplicated entries with at least one
  reported/likely catalyst; uncertain-only entries are excluded.

Focused regression checks (from the repository root):

```bash
python -m unittest discover -s tests -p 'test_metal_catalysts.py'
```

**Outputs**: 
- Figures saved to `output/analysis_figures/`
- All plots generated from `data/OPRD-100.json`

### 2. reaxys_comparison_clean.ipynb

Comparison between OPRD-100 and Reaxys datasets from the same 100 articles:
- **Data Location Comparison**: Comparison of where data is reported
- **Yield Type Comparison**: Types and distributions across datasets
- **Reaction Networks**: Molecule (node) and reaction (edge) overlap analysis
- **Stereochemistry Conflicts**: Detailed classification of stereochemical differences
- **Molecule Identity**: Classification as intermediates vs reagents (using CAS lookup)
- **Example Networks**: Visualization of reaction networks

**Outputs**:
- Figures saved to `output/reaxys_comparison_figures/`
- Data files (CSV/Excel) saved to `output/reaxys_comparison_data/`

## Installation

### Option 1: Using Conda (Recommended)

The easiest way to install all dependencies, especially rdkit and cairosvg which have system dependencies:

```bash
# Create the conda environment with Python 3.12
conda env create -f environment.yml

# Activate the environment
conda activate oprd100
```
## Running the Notebooks

1. **Activate Environment** (if using conda):
   ```bash
   conda activate oprd100
   ```

2. **Start Jupyter**:
   ```bash
   cd notebooks/
   jupyter notebook
   ```

3. **Run Notebooks**:
   - Open either notebook
   - Run cells sequentially (the notebooks handle path setup automatically)
   - Outputs will be generated in the `output/` directory

## Notes

- The notebooks use relative paths that work from the `notebooks/` directory
- Python modules are automatically imported from `../src/`
- All original analysis logic is preserved without refactoring
- Notebooks include comprehensive markdown documentation for each section

## Data Files

- **OPRD-100.json**: Main dataset containing 3878 reactions from 100 OPRD articles
- **reaxys_100_articles_data.xml**: Reaxys data extracted from the same 100 articles
- **CHEM21_solvent_dict.json**: Solvent classification data for environmental assessment
