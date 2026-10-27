# Submission Instructions

Submit your automatically extracted reaction data to benchmark against the OPRD-100 ground-truth dataset!

## 📋 File Format

Create a JSON file in this directory with the following structure:

```json
{
  "submitter_name": "YourName_or_TeamName",
  "method_description": "Brief description of your extraction method (e.g., LLM-based extraction using GPT-4, Rule-based NLP pipeline, etc.)",
  "repository_url": "https://github.com/your/repo",
  "paper_url": "https://doi.org/your.paper (optional)",
  "contact_email": "your.email@example.com (optional)",
  "reactions": [
    {
        "Reaction": "O=C(O)/C=C/c1ccccc1.Cc1ccc(O)cc1>>Cc1ccc2c(c1)C(c1ccccc1)CC(=O)O2",
        "Reference": "op050024w",
        "Location": {
            "Type": "Scheme",
            "Num": "1"
        },
        "Full text of reaction": null,
        "Amounts": null,
        "Steps": [
            {
                "Step": 1,
                "Yield": {
                    "Isolated": 89.0,
                    "Flow": null,
                    "Conversion wrt SM": null,
                    "Assay Yield": null,
                    "Selectivity": null,
                    "er": null,
                    "er R:S": null,
                    "er S:R": null,
                    "ee": null,
                    "ee (R)": null,
                    "ee (S)": null,
                    "Purity": null,
                    "dr": null,
                    "E:Z": null,
                    "Z:E": null,
                    "Purity type": null,
                    "Mass": null,
                    "other": null,
                    "Quantitative": null,
                    "No reaction": null
                },
                "Reagents": [
                    {
                        "Reagent": "H2SO4",
                        "Amounts": {
                            "Moisture": null,
                            "Concentration": null,
                            "Moles": null,
                            "Mass": null,
                            "Volume": null,
                            "Equivalents": null,
                            "Pressure": null,
                            "Flow Rate": null,
                            "Catalytic Amount": null
                        }
                    }
                ],
                "Solvents": [
                    {
                        "Solvent": null,
                        "Amounts": {
                            "Mass": null,
                            "Volume": null,
                            "ratio": null
                        }
                    }
                ],
                "Time": null,
                "Temperature": "120 - 125 C"
            }
        ]
    }
  ]
}
```

## 📝 Field Descriptions

### Required Fields
- `submitter_name`: Your name or team identifier (will be displayed on leaderboard)
- `method_description`: Brief description of your extraction approach
- `reactions`: Array of reaction data following the OPRD-100 schema

### Optional Fields
- `repository_url`: Optional HTTP(S) link to your extraction/scoring code repository,
  displayed in the leaderboard's **Code repository** column. Omit this field, use
  `null`, or leave it empty (`""`) if you do not wish to share your code. The table
  will show **Not shared**, with no effect on your score or ranking.
- `paper_url`: Link to associated paper or preprint
- `contact_email`: Your contact email

Code does not have to be public to participate. However, submission JSON files,
PRs, and evaluation results in this repository are public: do not include private
repository URLs, access tokens, or confidential contact details. Omitting the code
link keeps the code unshared; it does not make the prediction JSON private.

### Reaction Data Format
Follow the OPRD-100 schema structure. See the [validation_reactions.json](../OPRD-100.json) file for complete examples.

## 🚀 How to Submit

1. **Fork this repository**
2. **Extract reactions** from OPRD-100 papers using your method
3. **Create your submission file** as `data/submissions/your_name.json`
4. **Create a pull request**
   - `git checkout -b your_submission_team_name`
   - `git add data/submissions/your_name.json`
   - `git commit -m "[SUBMISSION] your_submission_team_name"`
   - head over to GitHub to complete the pull request with the following: 
    - Title: `[SUBMISSION] Your Name/Team`
    - Description: Brief overview of your method
5. **Automated evaluation** will:
   - Run scoring against the ground-truth dataset
   - Generate plots in `results/your_name_timestamp/`
   - Update the leaderboard in README.md
   - Comment results on your PR

## ✅ Validation Tips

Before submitting:
- Ensure your JSON is valid (use a JSON validator)
- Match the OPRD-100 schema structure
- Extract from the same papers/locations as the ground-truth dataset
- Include all required fields

## 📊 Scoring Metrics

Your submission will be evaluated using:
- **Reaction SMILES similarity**: Exact string matching
- **Reagent similarity**: Name and amount matching
- **Solvent similarity**: Name and amount matching  
- **Time/Temperature similarity**: Value comparison
- **Yield similarity**: Product and percentage matching

See [example_scoring.ipynb](../notebooks/example_scoring.ipynb) for detailed metric descriptions.

## 🤝 Questions?

Open an issue if you need help with the submission format or have questions about the evaluation process.
