"""Evidence-based catalyst candidates, not a mechanistic assignment of catalysis."""

import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from rdkit import Chem, rdBase


CONFIDENCE_CLASSES = ("reported_catalyst", "likely_catalyst", "uncertain")
CLASS_LABELS = ("Reported catalyst", "Likely catalyst", "Uncertain candidate")
METALS = frozenset(
    "Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr Nb Mo Tc "
    "Ru Rh Pd Ag Cd In Sn Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
    "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf "
    "Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv".split()
)
SALT_CATIONS = frozenset("Li Na K Rb Cs Mg Ca Sr Ba".split())
COMMON_CATALYST_CENTRES = frozenset("Pd Pt Rh Ru Ir Os".split())
METAL_WORDS = {
    "lithium": "Li", "sodium": "Na", "potassium": "K", "cesium": "Cs",
    "caesium": "Cs", "magnesium": "Mg", "calcium": "Ca", "aluminium": "Al",
    "aluminum": "Al", "scandium": "Sc", "titanium": "Ti", "vanadium": "V",
    "chromium": "Cr", "manganese": "Mn", "iron": "Fe", "ferrous": "Fe",
    "ferric": "Fe", "cobalt": "Co", "nickel": "Ni", "copper": "Cu",
    "cuprous": "Cu", "cupric": "Cu", "zinc": "Zn", "zirconium": "Zr",
    "molybdenum": "Mo", "ruthenium": "Ru", "rhodium": "Rh", "palladium": "Pd",
    "silver": "Ag", "tin": "Sn", "cerium": "Ce", "ytterbium": "Yb",
    "tungsten": "W", "osmium": "Os", "iridium": "Ir", "platinum": "Pt",
    "gold": "Au", "mercury": "Hg", "lead": "Pb", "bismuth": "Bi",
}
SYMBOL_PATTERN = re.compile(
    r"(?<![A-Za-z])(" + "|".join(sorted(METALS, key=len, reverse=True)) + r")(?![a-z])"
)
WORD_PATTERN = re.compile(r"\b(" + "|".join(METAL_WORDS) + r")\b", re.I)
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)"
LOADING_PATTERN = re.compile(
    rf"({NUMBER})\s*(mol\s*%|mol\s*percent|equiv(?:alents?)?\.?|eq\.?)"
    r"(?![A-Za-z])", re.I
)
RATIO_PATTERN = re.compile(rf"\bS\s*/\s*C\s*=\s*({NUMBER})", re.I)
STOICHIOMETRIC_REAGENTS = frozenset(
    name.casefold() for name in (
        "Zn", "Mn", "Fe", "iron powder", "CrO3", "CrCl2", "MnO2", "KMnO4",
        "LiAlH4", "Red-Al", "DIBAl-H", "AlMe3", "Bu3SnH", "n-BuLi", "LiHMDS",
        "NaBH4", "NaH", "NaOH", "LiOH", "KOH", "NaOMe", "NaOAc", "K2CO3",
        "Cs2CO3", "AlCl3",
    )
)
NONMETAL_ALIASES = frozenset(
    name.casefold() for name in (
        "DBU", "DMPU", "TBTU", "WSC", "WSCD", "MVK", "PrTMP", "Nh2OH.HCl",
        "RuPhos", "RuPhos ligand", "SPhos", "XPhos", "BrettPhos", "tBuXPhos",
        "Xantphos", "AmPhos", "IPr",
    )
)
SPECTATOR_SALTS = frozenset(("NaCl", "NaBr", "NaI", "KCl", "KBr", "KI"))


@dataclass(frozen=True)
class KnownCatalyst:
    name: str
    metal: str
    aliases: tuple[str, ...]
    source_url: str
    application: str
    cas_number: str | None = None
    pubchem_cid: int | None = None
    ligand_metals: tuple[str, ...] = ()

    def provenance(self) -> dict:
        return {
            "name": self.name, "metal": self.metal, "source_url": self.source_url,
            "application": self.application, "reviewed_on": REGISTRY_REVIEWED_ON,
            "cas_number": self.cas_number, "pubchem_cid": self.pubchem_cid,
            "identity_api_url": (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
                f"{self.pubchem_cid}/synonyms/JSON" if self.pubchem_cid else None
            ),
            "ligand_metals": list(self.ligand_metals),
        }


REGISTRY_REVIEWED_ON = "2026-10-25"
KNOWN_CATALYSTS = (
    KnownCatalyst(
        "Copper(I) iodide", "Cu",
        ("CuI", "copper(I) iodide", "copper (I) iodide", "cuprous iodide"),
        "https://www.thermofisher.com/order/catalog/product/011606.A1",
        "Organic synthesis and coupling catalysis",
    ),
    KnownCatalyst(
        "Scandium(III) triflate", "Sc",
        ("Sc(OTf)3", "scandium triflate", "scandium(III) triflate",
         "scandium(III) trifluoromethanesulfonate"),
        "https://www.thermofisher.com/order/catalog/product/040566.03",
        "Lewis acid catalysis",
    ),
    KnownCatalyst(
        "Bismuth(III) triflate", "Bi",
        ("Bi(OTf)3", "bismuth triflate", "bismuth(III) triflate",
         "bismuth(III) trifluoromethanesulfonate"),
        "https://www.thermofisher.com/order/catalog/product/L19687.06",
        "Lewis acid catalysis, including acylation and substitution",
    ),
    KnownCatalyst(
        "Zinc triflate", "Zn",
        ("Zn(OTf)2", "zinc triflate", "zinc(II) triflate",
         "zinc trifluoromethanesulfonate", "zinc(II) trifluoromethanesulfonate"),
        "https://www.thermofisher.com/order/catalog/product/L15969.06",
        "Lewis acid catalysis, including thioketalization and silylation",
    ),
    KnownCatalyst(
        "Grubbs metathesis catalysts", "Ru",
        ("Grubbs catalyst", "Grubbs' catalyst", "Grubbs I", "Grubbs II",
         "Grubbs catalyst 1st generation", "Grubbs catalyst 2nd generation",
         "first-generation Grubbs catalyst", "second-generation Grubbs catalyst",
         "Grubbs catalyst M102", "Grubbs catalyst M204"),
        "https://pmc.umicore.com/en/products/grubbs-catalysts/",
        "Ruthenium olefin metathesis catalyst family",
    ),
    KnownCatalyst(
        "Hoveyda-Grubbs metathesis catalysts", "Ru",
        ("Hoveyda-Grubbs catalyst", "Hoveyda Grubbs catalyst",
         "Hoveyda-Grubbs I", "Hoveyda-Grubbs II",
         "Hoveyda-Grubbs catalyst M700", "Hoveyda-Grubbs catalyst M720"),
        "https://pmc.umicore.com/en/products/grubbs-catalysts/",
        "Ruthenium olefin metathesis catalyst family",
    ),
    KnownCatalyst(
        "Wilkinson's catalyst", "Rh",
        ("Wilkinson's catalyst", "Wilkinson catalyst", "Wilkinson's",
         "chlorotris(triphenylphosphine)rhodium(I)", "RhCl(PPh3)3",
         "[(C6H5)3P]3RhCl"),
        "https://www.thermofisher.com/order/catalog/product/010468.06",
        "Homogeneous hydrogenation",
    ),
    KnownCatalyst(
        "Lindlar catalyst", "Pd",
        ("Lindlar catalyst", "Lindlar's catalyst", "Lindlar"),
        "https://www.thermofisher.com/order/catalog/product/043172.14",
        "Lead-poisoned palladium on calcium carbonate for selective hydrogenation",
    ),
    KnownCatalyst(
        "Palladium(II) acetate", "Pd",
        ("Pd(OAc)2", "palladium acetate", "palladium diacetate",
         "palladium(II)acetate"),
        "https://pmc.umicore.com/en/products/product-finder/pd-oac-2-3",
        "Palladium precursor for cross-coupling and carbonylation",
        "3375-31-3", 167845,
    ),
    KnownCatalyst(
        "Tris(dibenzylideneacetone)dipalladium(0)", "Pd",
        ("Pd2(dba)3", "Pd2dba3"),
        "https://www.thermofisher.com/order/catalog/product/012760.03",
        "Pd(0) catalyst precursor for cross-coupling",
        "51364-51-3", 9811564,
    ),
    KnownCatalyst(
        "Palladium(II) chloride", "Pd", ("PdCl2", "palladium dichloride"),
        "https://www.thermofisher.com/order/catalog/product/011034.18",
        "Palladium catalyst precursor",
        "7647-10-1",
    ),
    KnownCatalyst(
        "Palladium(II) trifluoroacetate", "Pd", ("Pd(TFA)2",),
        "https://www.thermofisher.com/order/catalog/product/039448.03",
        "Decarboxylation and direct arene cross-coupling",
        "42196-31-6",
    ),
    KnownCatalyst(
        "Tetrakis(triphenylphosphine)palladium(0)", "Pd",
        ("Pd(PPh3)4", "Pd(Ph3P)4", "tetrakis(triphenylphosphine)palladium"),
        "https://www.thermofisher.com/order/catalog/product/010548.09",
        "Palladium phosphine catalyst for cross-coupling",
        "14221-01-3", 11979704,
    ),
    KnownCatalyst(
        "Bis(triphenylphosphine)palladium(II) dichloride", "Pd",
        ("PdCl2(PPh3)2", "Pd(PPh3)2Cl2", "Pd(Ph3P)2Cl2", "PdCl2(Ph3P)2",
         "(Ph3P)2PdCl2", "bis(triphenylphosphine)palladium chloride"),
        "https://www.thermofisher.com/order/catalog/product/299250025",
        "Palladium phosphine cross-coupling precatalyst",
        "13965-03-2", 6102075,
    ),
    KnownCatalyst(
        "PdCl2(dppf)", "Pd", ("Pd(dppf)Cl2", "PdCl2dppf"),
        "https://www.thermofisher.com/order/catalog/product/046503.03",
        "Ferrocene-phosphine palladium cross-coupling precatalyst",
        "72287-26-4", 71310626, ("Fe",),
    ),
    KnownCatalyst(
        "PdCl2(dppf).CH2Cl2", "Pd",
        ("PdCl2.dppf.CH2Cl2", "PdCl2dppf.DCM", "PdCl2(dppf).DCM",
         "Pd(dppf)Cl2.CH2Cl2"),
        "https://www.thermofisher.com/order/catalog/product/348680250",
        "Dichloromethane-solvated dppf palladium cross-coupling precatalyst",
        "95464-05-4", 51341998, ("Fe",),
    ),
    KnownCatalyst(
        "PdCl2(dtbpf)", "Pd", ("Pd(dtbpf)Cl2",),
        "https://www.thermofisher.com/order/catalog/product/044829.03",
        "Ferrocene-phosphine palladium Suzuki cross-coupling precatalyst",
        "95408-45-0", 71310826, ("Fe",),
    ),
    KnownCatalyst(
        "PdCl2(AmPhos)2", "Pd", ("Pd(AmPhos)2Cl2",),
        "https://www.thermofisher.com/order/catalog/product/447890010",
        "Palladium phosphine cross-coupling precatalyst",
        "887919-35-9",
    ),
    KnownCatalyst(
        "Palladium on carbon", "Pd",
        ("Pd/C", "Pd-C", "palladium carbon", "5% Pd/C", "10% Pd/C", "Pd/C 5%"),
        "https://www.thermofisher.com/order/catalog/product/A12012.04",
        "Supported palladium catalyst family for hydrogenation and coupling",
    ),
    KnownCatalyst(
        "Palladium hydroxide on carbon", "Pd",
        ("Pearlman's catalyst", "Pearlman catalyst", "Pd(OH)2/C",
         "Pd(OH)2 on charcoal"),
        "https://www.thermofisher.com/order/catalog/product/042578.04",
        "Supported palladium hydroxide hydrogenolysis catalyst",
        "12135-22-7",
    ),
    KnownCatalyst(
        "SPhos Pd G2", "Pd", ("Pd(SPhos)G2", "SPhos G2", "Umicore CX121", "CX121"),
        "https://pmc.umicore.com/en/products/buchwald/",
        "Second-generation Buchwald cross-coupling precatalyst",
        "1375325-64-6",
    ),
    KnownCatalyst(
        "XPhos Pd G2", "Pd",
        ("Pd(XPhos)G2", "XPhosPd G2", "XPhos G2", "Umicore CX122", "CX122"),
        "https://pmc.umicore.com/en/products/buchwald/",
        "Second-generation Buchwald cross-coupling precatalyst",
        "1310584-14-5",
    ),
    KnownCatalyst(
        "RuPhos Pd G2", "Pd", ("Pd(RuPhos)G2", "RuPhos G2", "Umicore CX123", "CX123"),
        "https://pmc.umicore.com/en/products/buchwald/",
        "Second-generation Buchwald cross-coupling precatalyst",
    ),
    KnownCatalyst(
        "SPhos Pd G3", "Pd", ("Pd(SPhos)G3", "SPhos G3", "Umicore CX131", "CX131"),
        "https://pmc.umicore.com/en/products/buchwald/",
        "Third-generation Buchwald cross-coupling precatalyst",
        "1445085-82-4",
    ),
    KnownCatalyst(
        "XPhos Pd G3", "Pd", ("Pd(XPhos)G3", "XPhos G3", "Umicore CX132", "CX132"),
        "https://pmc.umicore.com/en/products/buchwald/",
        "Third-generation Buchwald cross-coupling precatalyst",
        "1445085-55-1",
    ),
    KnownCatalyst(
        "RuPhos Pd G3", "Pd", ("Pd(RuPhos)G3", "RuPhos G3", "Umicore CX133", "CX133"),
        "https://pmc.umicore.com/en/products/buchwald/",
        "Third-generation Buchwald cross-coupling precatalyst",
        "1445085-77-7",
    ),
    KnownCatalyst(
        "tBuXPhos Pd G3", "Pd",
        ("Pd(tBuXPhos)G3", "tBuXPhos G3", "Umicore CX135", "CX135"),
        "https://pmc.umicore.com/en/products/buchwald/",
        "Third-generation Buchwald amination precatalyst",
        "1447963-75-8",
    ),
    KnownCatalyst(
        "Xantphos Pd G3", "Pd",
        ("Pd(Xantphos)G3", "Xantphos G3", "Umicore CX136", "CX136"),
        "https://pmc.umicore.com/en/products/buchwald/",
        "Third-generation Buchwald cross-coupling precatalyst",
        "1445085-97-1",
    ),
    KnownCatalyst(
        "PEPPSI-IPr catalyst", "Pd", ("PEPPSI-IPr", "Pd-PEPPSI-IPr"),
        "https://www.sigmaaldrich.com/US/en/product/aldrich/669032",
        "Palladium N-heterocyclic carbene cross-coupling precatalyst",
        "905459-27-0",
    ),
)


def _catalyst_key(name: str) -> str:
    name = name.replace("\u2019", "'").strip()
    name = re.sub(rf"^(?:cat\.\s*|{NUMBER}\s*mol\s*%\s*)", "", name, flags=re.I)
    return re.sub(r"\s+", "", name).casefold()


KNOWN_CATALYST_ALIASES = {
    _catalyst_key(alias): catalyst
    for catalyst in KNOWN_CATALYSTS
    for alias in (catalyst.name, *catalyst.aliases)
}


@dataclass(frozen=True)
class MetalIdentity:
    metals: tuple[str, ...]
    ligand_metals: tuple[str, ...]
    method: str
    parsed_structure: bool


@lru_cache(maxsize=None)
def _structure_metals(smiles: str) -> tuple[str, ...] | None:
    # Syntax/atom identity suffices; organometallic valence models need not sanitize.
    params = Chem.SmilesParserParams()
    params.sanitize = False
    params.parseName = False
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles, params)
    if mol is None:
        return None
    return tuple(sorted({atom.GetSymbol() for atom in mol.GetAtoms()} & METALS))


@lru_cache(maxsize=None)
def resolve_metal_identity(name: str) -> MetalIdentity:
    """Resolve atoms first, then bounded formula tokens and reviewed aliases."""
    known = KNOWN_CATALYST_ALIASES.get(_catalyst_key(name))
    if known is not None:
        return MetalIdentity((known.metal,), known.ligand_metals, "supplier_registry_alias", False)
    clean = re.sub(r"\s*\(S\s*/\s*C\s*=.*?\)\s*$", "", name, flags=re.I).strip()
    structure = _structure_metals(clean)
    if structure is not None:
        metals = set(structure)
        method = "smiles_atoms"
    else:
        clean = re.sub(r"\bXPhos(?=Pd)", "XPhos ", clean, flags=re.I)
        # Acetyl/propyl are not actinium/praseodymium. Supplier codes and lamp
        # wattages are not formulae. Keep the original name in the exported audit.
        tokens = re.sub(r"Ac|Pr", " ", clean)
        tokens = re.sub(r"(?<=[a-z])\([IVX]+\)", "", tokens)
        tokens = re.sub(r"\b\d+\s*W\b.*", "", tokens)
        tokens = re.sub(r"\b(?:BASF|Degussa|Johnson Matthey)\b.*?(?=\d+\s*%)",
                        "", tokens)
        tokens = re.split(r"\b(?:CP|E\d+|intermediate|egg|uniform|type)\b",
                          tokens, maxsplit=1)[0]
        metals = set(SYMBOL_PATTERN.findall(tokens))
        # A metal attached to an explicit oxygen atom can occur in formula names.
        metals.update(re.findall(r"O(Zn|Mg|Li|Na|K)(?![a-z])", tokens))
        if "LiAlH4" in name:
            metals.add("Al")
        metals.update(METAL_WORDS[m.group().lower()] for m in WORD_PATTERN.finditer(clean))
        method = "formula_tokens_or_metal_name"
        if name.casefold() in NONMETAL_ALIASES:
            return MetalIdentity((), (), "reviewed_nonmetal_alias", True)
        if name == "Na2Co3":
            return MetalIdentity((), (), "ambiguous_formula_capitalization", False)
    ligand_metals = set()
    if re.search(r"dppf|dtbpf|ferrocen", name, re.I):
        ligand_metals.add("Fe")
        metals.discard("Fe")
        if "palladium" in name.lower():
            metals.add("Pd")
    if name.casefold() == "red-al":
        metals.update(("Na", "Al"))
    return MetalIdentity(tuple(sorted(metals)), tuple(sorted(ligand_metals)),
                         method, structure is not None)


def _loading(name: str, amounts: dict) -> tuple[float | None, list[str], list[str]]:
    values = []
    evidence = []
    issues = []
    for field, raw in [("Reagent", name), *amounts.items()]:
        if raw is None:
            continue
        text = str(raw)
        if re.search(r"\d\s*(?:-|to|/)\s*\d", text, re.I) and re.search(
            r"mol\s*%|equiv|\beq\b", text, re.I
        ):
            issues.append(f"Loading range/fraction requires review in {field}: {text}")
            continue
        for match in LOADING_PATTERN.finditer(text):
            value = float(match[1])
            if not match[2].lower().startswith("mol"):
                value *= 100
            values.append(value)
            evidence.append(f"{field}: {match[0]} ({value:g} mol%)")
        if field == "Equivalents" and re.fullmatch(NUMBER, text.strip()):
            value = float(text) * 100
            values.append(value)
            evidence.append(f"Equivalents: {text} ({value:g} mol%)")
    for match in RATIO_PATTERN.finditer(name):
        ratio = float(match[1])
        if ratio <= 0:
            issues.append("S/C ratio must be positive")
        else:
            values.append(100 / ratio)
            evidence.append(f"S/C = {ratio:g} ({100 / ratio:g} mol%)")
    if any(v <= 0 for v in values):
        issues.append("Catalyst loading must be positive")
    if values and max(values) - min(values) > 1e-6:
        issues.append("Conflicting loading values")
    return (values[0] if values and not issues else None), evidence, issues


def _classify(name: str, metal: str, identity: MetalIdentity, amounts: dict) -> dict:
    loading, evidence, issues = _loading(name, amounts)
    known = KNOWN_CATALYST_ALIASES.get(_catalyst_key(name))
    registry_match = known is not None and known.metal == metal
    supplier_evidence = None
    if known is not None and registry_match:
        supplier_evidence = known.provenance()
        evidence.append(f"Supplier-documented catalyst identity: {known.name}")
    flag = amounts.get("Catalytic Amount")
    reported = flag is True or (isinstance(flag, str) and flag.strip().lower() == "true")
    denied = flag is False or (isinstance(flag, str) and flag.strip().lower() == "false")
    if flag is not None and not reported and not denied:
        issues.append(f"Unrecognised Catalytic Amount: {flag!r}")
    if reported:
        evidence.append("Catalytic Amount explicitly true")
    if denied:
        evidence.append("Catalytic Amount explicitly false")
    other_centres = set(identity.metals) - SALT_CATIONS
    spectator = metal in SALT_CATIONS and (
        bool(other_centres) or name.strip() in SPECTATOR_SALTS
    )
    known_stoichiometric = name.strip().casefold() in STOICHIOMETRIC_REAGENTS
    supported = bool(re.search(
        r"/C(?:\b|aCO3)|on (?:activated )?(?:carbon|charcoal)|Raney|Pd-C", name, re.I
    ))
    if spectator:
        classification = "counterion"
        evidence.append("Reviewed spectator salt or cation accompanying another metal centre")
    elif issues or (reported and loading is not None and loading >= 100):
        classification = "uncertain"
        evidence.extend(issues)
        if reported:
            evidence.append("Catalytic annotation conflicts with stoichiometric loading")
    elif reported:
        classification = "reported_catalyst"
    elif denied or known_stoichiometric or (loading is not None and loading >= 100):
        classification = "noncatalytic_metal_reagent"
        evidence.append("Explicit noncatalytic flag, usual stoichiometric reagent, or >=1 equiv")
    elif metal in SALT_CATIONS:
        classification = "uncertain" if loading is not None else "noncatalytic_metal_reagent"
        evidence.append("Salt/base metal; low loading alone does not establish catalysis")
    elif loading is not None:
        classification = "reported_catalyst"
        evidence.append("Explicit substoichiometric loading reports catalytic-scale use")
    elif registry_match:
        classification = "likely_catalyst"
        evidence.append("Known catalytic use; actual role in this reaction is inferred")
    elif (metal in COMMON_CATALYST_CENTRES
          or (supported and metal in {"Ni", "Cu"})
          or re.search(r"\bcat(?:alyst|alytic)?\.?\b", name, re.I)):
        classification = "likely_catalyst"
        evidence.append("Recognised catalyst family/formulation; loading unavailable")
    else:
        classification = "uncertain"
        evidence.append("Metal identified; catalytic role not established")
    return {
        "metal": metal, "classification": classification,
        "loading_mol_percent": loading, "evidence": evidence,
        "needs_review": bool(issues) or classification == "uncertain",
        "supplier_evidence": supplier_evidence,
    }


class MetalCatalystAnalyzer:
    """Export reagent evidence and unique reaction-metal confidence counts.

    Reported > likely > uncertain is used once per reaction-metal pair, across
    steps. Counterions and noncatalytic reagents remain auditable but unplotted.
    Substrate/product metals never independently generate catalyst assignments.
    """

    def analyze(self, data: list[dict], output_path: str | Path) -> dict:
        assignments = []
        reaction_metals = []
        unresolved = {}
        review_queue = []
        for reaction_index, entry in enumerate(data):
            context = {"reaction_index": reaction_index, "reference": entry.get("Reference")}
            for step_index, step in enumerate(entry.get("Steps") or []):
                for reagent_index, reagent in enumerate(step.get("Reagents") or []):
                    name = reagent.get("Reagent")
                    amounts = reagent.get("Amounts") or {}
                    location = {**context, "step_index": step_index,
                                "step": step.get("Step"), "reagent_index": reagent_index}
                    if not name:
                        review_queue.append({**location, "reason": "Missing reagent identity"})
                        continue
                    identity = resolve_metal_identity(name)
                    if not identity.metals and not identity.parsed_structure:
                        unresolved.setdefault(name, []).append({
                            **location, "amounts": amounts,
                            "identity_method": identity.method,
                        })
                    for metal in identity.metals:
                        result = {**location, "reagent": name, "amounts": amounts,
                                  "identity_method": identity.method,
                                  "ligand_metals": list(identity.ligand_metals),
                                  **_classify(name, metal, identity, amounts)}
                        assignments.append(result)
                        if result["needs_review"]:
                            review_queue.append(result)
            reaction = entry.get("Reaction")
            parts = reaction.split(">") if isinstance(reaction, str) else []
            if len(parts) != 3:
                review_queue.append({**context, "reason": "Missing/invalid reaction SMILES"})
                continue
            for side, smiles in zip(("reactants", "agents", "products"), parts):
                if not smiles:
                    continue
                metals = _structure_metals(smiles)
                if metals is None:
                    review_queue.append({**context, "side": side,
                                         "reason": "Unparsed reaction SMILES side"})
                elif metals:
                    reaction_metals.append({
                        **context, "side": side, "metals": list(metals),
                        "classification": "reaction_component_only",
                        "evidence": "Atom identity only; not used to infer a catalyst",
                    })

        best = {}
        for assignment in assignments:
            category = assignment["classification"]
            if category not in CONFIDENCE_CLASSES:
                continue
            key = (assignment["reaction_index"], assignment["metal"])
            if key not in best or CONFIDENCE_CLASSES.index(category) < CONFIDENCE_CLASSES.index(best[key]):
                best[key] = category
        counts = defaultdict(Counter)
        for (_, metal), category in best.items():
            counts[metal][category] += 1
        metals = sorted(counts, key=lambda m: (-sum(counts[m].values()), m))
        candidate_indices = sorted({index for index, _ in best})
        supported_indices = sorted({index for (index, _), category in best.items()
                                    if category != "uncertain"})
        report = {
            "schema_version": 3,
            "method": {
                "confidence_order": list(CONFIDENCE_CLASSES),
                "counting_unit": "unique dataset entry per metal, strongest evidence across steps",
                "reported_criteria": (
                    "Explicit Catalytic Amount=true or a positive loading below 100 mol% "
                    "(including equivalents and S/C ratios), subject to counterion, "
                    "noncatalytic-role, and conflicting-evidence checks."
                ),
                "supplier_registry": {
                    "reviewed_on": REGISTRY_REVIEWED_ON,
                    "entries": [
                        {**c.provenance(), "aliases": list(c.aliases)}
                        for c in KNOWN_CATALYSTS
                    ],
                },
                "limitations": [
                    "Heuristic role assignments, not calibrated probabilities or mechanistic proof.",
                    "Formula tokens and names are candidates, not full chemical name resolution.",
                    "Missing identities and unparsed names require review; coverage is not complete.",
                    "Substoichiometric loading can describe additives rather than catalysts.",
                    "Metal wt% on a support is not interpreted as mol% loading.",
                    "Supplier recognition establishes known catalytic use, not the role in this reaction.",
                ],
            },
            "summary": {
                "reaction_count": len(data),
                "reactions_with_reported_or_likely_catalysts": len(supported_indices),
                "reactions_with_any_candidate": len(candidate_indices),
                "reaction_metal_pairs": len(best),
                "by_metal": {metal: {c: counts[metal][c] for c in CONFIDENCE_CLASSES}
                             for metal in metals},
                "assignment_classes": dict(Counter(a["classification"] for a in assignments)),
                "supplier_matched_assignments": sum(
                    a["supplier_evidence"] is not None for a in assignments
                ),
                "unresolved_unique_names": len(unresolved),
                "review_items": len(review_queue),
            },
            "assignments": assignments,
            "reaction_component_metals": reaction_metals,
            "review_queue": review_queue,
            "unresolved_reagents": [
                {"reagent": name, "occurrences": locations,
                 "reason": "Not parsed as SMILES and no recognised metal identity"}
                for name, locations in sorted(unresolved.items())
            ],
            "candidate_reaction_indices": candidate_indices,
            "supported_reaction_indices": supported_indices,
        }
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2, allow_nan=False)
            handle.write("\n")
        return report
