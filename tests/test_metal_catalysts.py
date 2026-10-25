"""Focused chemistry regressions and export/counting invariants."""

import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "data_analysis" / "src"))

from MetalCatalysts import (
    CONFIDENCE_CLASSES,
    KNOWN_CATALYSTS,
    MetalCatalystAnalyzer,
    resolve_metal_identity,
)


def reagent(name, **amounts):
    return {"Reagent": name, "Amounts": amounts}


def reaction(*reagents, smiles="CC>>CO"):
    return {"Reference": "test", "Reaction": smiles,
            "Steps": [{"Step": 1, "Reagents": list(reagents)}]}


class MetalCatalystTests(unittest.TestCase):
    def analyze(self, entries):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "nested" / "evidence.json"
            report = MetalCatalystAnalyzer().analyze(entries, path)
            self.assertEqual(json.loads(path.read_text()), report)
            return report

    def test_nickel_catalyst_and_stoichiometric_manganese(self):
        report = self.analyze([reaction(
            reagent("NiCl2.DME", Concentration="5 mol%"),
            reagent("Mn", Equivalents="2 equiv"),
        )])
        ni, mn = report["assignments"]
        self.assertEqual(ni["classification"], "reported_catalyst")
        self.assertEqual(ni["loading_mol_percent"], 5)
        self.assertEqual(mn["classification"], "noncatalytic_metal_reagent")
        self.assertNotIn("Mn", report["summary"]["by_metal"])

    def test_explicit_catalysts_and_string_false(self):
        for name, metal in [("Bi(OTf)3", "Bi"), ("Sc(OTf)3", "Sc"), ("K2CO3", "K")]:
            with self.subTest(name=name):
                report = self.analyze([reaction(reagent(name, **{"Catalytic Amount": "True"}))])
                self.assertEqual(report["summary"]["by_metal"][metal]["reported_catalyst"], 1)
        report = self.analyze([reaction(reagent("Pd/C", **{"Catalytic Amount": "False"}))])
        self.assertEqual(report["assignments"][0]["classification"], "noncatalytic_metal_reagent")

    def test_osmate_counterion_and_spectator_halide(self):
        report = self.analyze([reaction(
            reagent("K2OsO4.H2O"), reagent("NaBr", **{"Catalytic Amount": "True"})
        )])
        classes = {a["metal"]: a["classification"] for a in report["assignments"]}
        self.assertEqual(classes, {"K": "counterion", "Os": "likely_catalyst", "Na": "counterion"})
        self.assertEqual(set(report["summary"]["by_metal"]), {"Os"})

    def test_product_and_reactant_metals_are_context_only(self):
        report = self.analyze([reaction(smiles="O=C([O-])C.[Li+]>>O=C([O-])C.[Li+]")])
        self.assertFalse(report["assignments"])
        self.assertFalse(report["summary"]["by_metal"])
        self.assertEqual(len(report["reaction_component_metals"]), 2)

    def test_reaction_agents_not_silently_used_as_catalysts(self):
        report = self.analyze([reaction(smiles="CC>[Pd]>CO")])
        self.assertFalse(report["summary"]["by_metal"])
        self.assertEqual(report["reaction_component_metals"][0]["side"], "agents")

    def test_ligand_iron_not_active_metal(self):
        for name in ("PdCl2(dppf).CH2Cl2", "PdCl2(dtbpf)",
                     "[1,10-bis(diphenylphosphino)ferrocene]dichloridepalladium-(II)"):
            with self.subTest(name=name):
                identity = resolve_metal_identity(name)
                self.assertEqual(identity.metals, ("Pd",))
                self.assertEqual(identity.ligand_metals, ("Fe",))

    def test_actual_smiles_and_coordination_shorthand(self):
        for name, metals in [
            ("[Na+].[O-][Cu]", ("Cu", "Na")),
            ("[Ir(dtbbpy)(ppy)2]PF6", ("Ir",)),
            ("[Ru(bpy)3]Cl2", ("Ru",)),
            ("Raney nickel", ("Ni",)),
            ("XPhosPd G2", ("Pd",)),
            ("4-FPhOZnBr", ("Zn",)),
            ("i-PrMgCl", ("Mg",)),
            ("LiAlH4", ("Al", "Li")),
        ]:
            with self.subTest(name=name):
                self.assertEqual(resolve_metal_identity(name).metals, metals)

    def test_abbreviations_and_words_are_not_element_substrings(self):
        for name in ("AcOH", "Ac2O", "i-PrOH", "DBU", "DMPU", "TBTU",
                     "WSC", "MVK", "100 W bulb irradiation", "Nh2OH.HCl",
                     "Nitrobenzene", "Morpholine", "phosphorus(V) oxychloride"):
            with self.subTest(name=name):
                self.assertEqual(resolve_metal_identity(name).metals, ())
        self.assertEqual(resolve_metal_identity("Pd(OAc)2").metals, ("Pd",))
        self.assertEqual(resolve_metal_identity("Ti(OiPr)4").metals, ("Ti",))

    def test_supplier_codes_do_not_add_uranium_or_tungsten(self):
        for name in ("BASF 3%Pd(II) CP M/UR Code 00030 type 39",
                     "Degussa 5%Pd(0)/C E196 R/W 5% egg shell type 39"):
            with self.subTest(name=name):
                self.assertEqual(resolve_metal_identity(name).metals, ("Pd",))

    def test_supported_metal_percent_is_not_molar_loading(self):
        report = self.analyze([reaction(reagent("10% Pd/C"))])
        assignment = report["assignments"][0]
        self.assertIsNone(assignment["loading_mol_percent"])
        self.assertEqual(assignment["classification"], "likely_catalyst")

    def test_loading_units_and_substrate_catalyst_ratio(self):
        for name, amounts, expected in [
            ("NiCl2", {"Equivalents": "0.05 equiv"}, 5),
            ("NiCl2", {"Equivalents": 0.05}, 5),
            ("NiCl2", {"Concentration": "5 mol %"}, 5),
            ("0.074 mol% Rh(acac)BDP", {}, 0.074),
            ("[Ru] (S/C = 10000)", {}, 0.01),
        ]:
            with self.subTest(name=name, amounts=amounts):
                assignment = self.analyze([reaction(reagent(name, **amounts))])["assignments"][0]
                self.assertAlmostEqual(assignment["loading_mol_percent"], expected)
                self.assertEqual(assignment["classification"], "reported_catalyst")

    def test_loading_reports_catalytic_scale_with_or_without_supplier_match(self):
        for name in ("CuI", "Pd(OAc)2", "NiCl2.DME", "ZrCl4"):
            for field in ("Concentration", "Moles", "Equivalents"):
                with self.subTest(name=name, field=field):
                    assignment = self.analyze([
                        reaction(reagent(name, **{field: "5 mol%"}))
                    ])["assignments"][0]
                    self.assertEqual(assignment["classification"], "reported_catalyst")
                    self.assertIn("Explicit substoichiometric loading reports catalytic-scale use",
                                  assignment["evidence"])

    def test_loading_does_not_promote_salt_reductant_or_denied_catalyst(self):
        for name, amounts, expected in (
            ("NaCl", {"Concentration": "5 mol%"}, "counterion"),
            ("sodium ascorbate", {"Concentration": "5 mol%"}, "uncertain"),
            ("Mn", {"Equivalents": "0.05 equiv"}, "noncatalytic_metal_reagent"),
            ("CuI", {"Concentration": "5 mol%", "Catalytic Amount": "False"},
             "noncatalytic_metal_reagent"),
            ("CuI", {"Concentration": "100 mol%"}, "noncatalytic_metal_reagent"),
        ):
            with self.subTest(name=name, amounts=amounts):
                assignment = self.analyze([reaction(reagent(name, **amounts))])["assignments"][0]
                self.assertEqual(assignment["classification"], expected)

    def test_conflicting_invalid_and_range_loadings_are_reviewed(self):
        for amounts in (
            {"Concentration": "5 mol%", "Equivalents": "0.2 equiv"},
            {"Concentration": "0 mol%"},
            {"Concentration": "-5 mol%"},
            {"Equivalents": "1/2 equiv"},
            {"Concentration": "5-10 mol%"},
            {"Catalytic Amount": "perhaps"},
            {"Catalytic Amount": "True", "Equivalents": "2 equiv"},
        ):
            with self.subTest(amounts=amounts):
                report = self.analyze([reaction(reagent("Pd(OAc)2", **amounts))])
                self.assertEqual(report["assignments"][0]["classification"], "uncertain")
                self.assertTrue(report["review_queue"])

    def test_strongest_class_once_per_entry_and_metal(self):
        first = reaction(reagent("Pd/C"), reagent("palladium acetate",
                                                **{"Catalytic Amount": "True"}))
        first["Steps"].append({"Step": 2, "Reagents": [reagent("Pd/C")]})
        report = self.analyze([first, reaction(reagent("Pd/C"))])
        self.assertEqual(report["summary"]["by_metal"]["Pd"],
                         {"reported_catalyst": 1, "likely_catalyst": 1, "uncertain": 0})
        self.assertEqual(report["supported_reaction_indices"], [0, 1])
        self.assertEqual(report["summary"]["reaction_metal_pairs"], 2)

    def test_identical_smiles_conditions_remain_distinct_entries(self):
        report = self.analyze([reaction(reagent("Pd/C")), reaction(reagent("Pd/C"))])
        self.assertEqual(report["summary"]["by_metal"]["Pd"]["likely_catalyst"], 2)

    def test_unresolved_missing_and_ambiguous_names_are_visible(self):
        report = self.analyze([
            reaction(reagent("unknown proprietary catalyst"), reagent(None), reagent("Na2Co3"))
        ])
        self.assertEqual(report["summary"]["unresolved_unique_names"], 2)
        self.assertEqual(report["review_queue"][0]["reason"], "Missing reagent identity")
        self.assertFalse(report["summary"]["by_metal"])

    def test_uncertain_metal_not_promoted_to_supported(self):
        report = self.analyze([reaction(reagent("ZrCl4"))])
        self.assertEqual(report["summary"]["by_metal"]["Zr"]["uncertain"], 1)
        self.assertEqual(report["candidate_reaction_indices"], [0])
        self.assertEqual(report["supported_reaction_indices"], [])

    def test_supplier_registry_aliases_resolve_without_element_substrings(self):
        for known in KNOWN_CATALYSTS:
            for alias in (known.name, *known.aliases):
                with self.subTest(alias=alias):
                    identity = resolve_metal_identity(alias)
                    self.assertEqual(identity.metals, (known.metal,))
                    self.assertEqual(identity.method, "supplier_registry_alias")

    def test_supplier_registry_is_likely_not_reported_and_exports_provenance(self):
        report = self.analyze([reaction(reagent(c.name)) for c in KNOWN_CATALYSTS])
        self.assertEqual(len(report["method"]["supplier_registry"]["entries"]),
                         len(KNOWN_CATALYSTS))
        for known, assignment in zip(KNOWN_CATALYSTS, report["assignments"]):
            self.assertEqual(assignment["classification"], "likely_catalyst")
            self.assertEqual(assignment["supplier_evidence"]["source_url"], known.source_url)
            self.assertIsNone(assignment["loading_mol_percent"])
        self.assertEqual(report["schema_version"], 3)

    def test_supplier_evidence_does_not_override_reaction_specific_evidence(self):
        for amounts, expected in [
            ({"Catalytic Amount": "True"}, "reported_catalyst"),
            ({"Catalytic Amount": "False"}, "noncatalytic_metal_reagent"),
            ({"Equivalents": "1 equiv"}, "noncatalytic_metal_reagent"),
            ({"Catalytic Amount": "True", "Equivalents": "2 equiv"}, "uncertain"),
            ({"Concentration": "5-10 mol%"}, "uncertain"),
        ]:
            with self.subTest(amounts=amounts):
                assignment = self.analyze([reaction(reagent("CuI", **amounts))])["assignments"][0]
                self.assertEqual(assignment["classification"], expected)
                self.assertIsNotNone(assignment["supplier_evidence"])

    def test_registry_normalizes_spacing_case_and_loading_but_not_substrings(self):
        for name, metal in [
            (" zinc (II) triflate ", "Zn"), ("Zn(Otf)2", "Zn"),
            ("cat. Sc(OTf)3", "Sc"), ("2 mol% Hoveyda-Grubbs catalyst", "Ru"),
            ("Wilkinson\u2019s catalyst", "Rh"),
        ]:
            with self.subTest(name=name):
                self.assertEqual(resolve_metal_identity(name).metals, (metal,))
        for name in ("Grubbs ligand", "Wilkinson-like ligand", "zinc triflate alternative"):
            with self.subTest(name=name):
                self.assertNotEqual(resolve_metal_identity(name).method, "supplier_registry_alias")

    def test_related_compounds_and_mixtures_not_automatically_promoted(self):
        for name in ("CuBr", "ZnCl2", "CuSO4/Na2SO3"):
            with self.subTest(name=name):
                assignments = self.analyze([reaction(reagent(name))])["assignments"]
                active = [a for a in assignments if a["classification"] != "counterion"]
                self.assertTrue(active)
                self.assertTrue(all(a["classification"] == "uncertain" for a in active))
                self.assertTrue(all(a["supplier_evidence"] is None for a in active))

    def test_palladium_precursors_and_named_precatalysts(self):
        names = (
            "Pd(OAc)2", "Pd2dba3", "Pd(PPh3)4", "(Ph3P)2PdCl2", "PdCl2",
            "Pd(TFA)2", "PdCl2(AmPhos)2", "PdCl2(dppf)", "PdCl2(dppf).CH2Cl2",
            "PdCl2(dtbpf)", "Pd/C", "10% Pd/C", "Pearlman's catalyst",
            "XPhosPd G2", "XPhos G3", "CX122", "RuPhos G2", "RuPhos Pd G3",
            "SPhos G2", "SPhos G3", "tBuXPhos G3", "Xantphos G3", "PEPPSI-IPr",
        )
        report = self.analyze([reaction(reagent(name)) for name in names])
        self.assertEqual(report["summary"]["supplier_matched_assignments"], len(names))
        for assignment in report["assignments"]:
            self.assertEqual(assignment["metal"], "Pd")
            self.assertEqual(assignment["classification"], "likely_catalyst")
            self.assertIsNotNone(assignment["supplier_evidence"])

    def test_palladium_ligand_iron_and_solvates_remain_distinct(self):
        report = self.analyze([reaction(
            reagent("PdCl2(dppf)"), reagent("PdCl2(dppf).CH2Cl2"),
            reagent("PdCl2(dtbpf)"),
        )])
        for assignment in report["assignments"]:
            self.assertEqual(assignment["ligand_metals"], ["Fe"])
            self.assertEqual(assignment["metal"], "Pd")
        evidence = [a["supplier_evidence"] for a in report["assignments"]]
        self.assertEqual(evidence[0]["cas_number"], "72287-26-4")
        self.assertEqual(evidence[1]["cas_number"], "95464-05-4")
        self.assertNotEqual(evidence[0]["pubchem_cid"], evidence[1]["pubchem_cid"])
        self.assertEqual(report["summary"]["reaction_metal_pairs"], 1)

    def test_api_identity_provenance_is_optional_and_not_a_catalysis_claim(self):
        records = {c.name: c for c in KNOWN_CATALYSTS}
        known = records["Palladium(II) acetate"]
        evidence = known.provenance()
        self.assertEqual(evidence["pubchem_cid"], 167845)
        self.assertIn("/cid/167845/synonyms/JSON", evidence["identity_api_url"])
        self.assertEqual(evidence["cas_number"], "3375-31-3")
        self.assertIsNone(records["PEPPSI-IPr catalyst"].provenance()["identity_api_url"])
        self.assertIsNone(records["XPhos Pd G2"].provenance()["pubchem_cid"])
        assignment = self.analyze([reaction(reagent(known.name))])["assignments"][0]
        self.assertEqual(assignment["classification"], "likely_catalyst")

    def test_free_buchwald_ligands_are_not_palladium_or_ruthenium_catalysts(self):
        for name in ("RuPhos", "RuPhos ligand", "ruphos", "XPhos", "SPhos",
                     "tBuXPhos", "Xantphos", "AmPhos", "IPr"):
            with self.subTest(name=name):
                self.assertFalse(resolve_metal_identity(name).metals)
                self.assertFalse(self.analyze([reaction(reagent(name))])["summary"]["by_metal"])
        for name in ("Pd(dppf) DCM", "PdCl2(amphos)", "precatalyst solution"):
            with self.subTest(name=name):
                self.assertNotEqual(resolve_metal_identity(name).method, "supplier_registry_alias")

    def test_registry_aliases_do_not_collide_across_identities(self):
        from MetalCatalysts import _catalyst_key

        seen = {}
        for known in KNOWN_CATALYSTS:
            for alias in (known.name, *known.aliases):
                key = _catalyst_key(alias)
                self.assertIn(seen.get(key), (None, known.name))
                seen[key] = known.name

    def test_empty_input_and_null_fields(self):
        self.assertFalse(self.analyze([])["summary"]["by_metal"])
        report = self.analyze([{"Reaction": None, "Steps": None}])
        self.assertEqual(report["review_queue"][0]["reason"], "Missing/invalid reaction SMILES")

    def test_full_dataset_export_and_counting_invariants(self):
        path = ROOT / "data" / "OPRD-100.json"
        original = path.read_bytes()
        data = json.loads(original)
        report = self.analyze(data)
        self.assertEqual(hashlib.sha256(path.read_bytes()).digest(),
                         hashlib.sha256(original).digest())
        best = {}
        for assignment in report["assignments"]:
            category = assignment["classification"]
            if category in CONFIDENCE_CLASSES:
                key = (assignment["reaction_index"], assignment["metal"])
                rank = CONFIDENCE_CLASSES.index(category)
                best[key] = min(best.get(key, rank), rank)
        for metal, counts in report["summary"]["by_metal"].items():
            self.assertEqual(list(counts), list(CONFIDENCE_CLASSES))
            for rank, category in enumerate(CONFIDENCE_CLASSES):
                self.assertEqual(counts[category], sum(
                    m == metal and r == rank for (_, m), r in best.items()))
        self.assertEqual(sum(sum(c.values()) for c in report["summary"]["by_metal"].values()),
                         len(best))
        self.assertFalse({"Ac", "Pr", "U", "Nh"} & report["summary"]["by_metal"].keys())
        self.assertEqual(report["summary"]["by_metal"]["Os"]["likely_catalyst"], 1)


if __name__ == "__main__":
    unittest.main()
