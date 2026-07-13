"""
This module provides functionality to compare data between manually extracted OPRD-100 entries and generated datasets using the same Schema.
"""
import time
from typing import Optional, Tuple
import json
import os
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity
import re
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.optimize import linear_sum_assignment

# Morgan fingerprint generator (radius 2, 2048 bits), reused across calls. Used by the
# optional graded-Tanimoto reaction-SMILES similarity method.
_MORGAN_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# Valid values for the reaction-SMILES similarity method.
SIMILARITY_METHODS = ("inchikey", "tanimoto")

def json_to_dict(file_path):
    """
    Reads a JSON file and converts it to a dictionary.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")

    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def remove_entry_from_location(location):
        """
        Removes 'entry' from the Location dict, and its corresponding Num.
        """
        if not isinstance(location, dict):
            return location
        types = [t.strip() for t in (location.get("Type") or "").split(",")]
        num_field = location.get("Num") or ""
        nums = [n.strip() for n in num_field.split(";")]
        # Remove 'entry' and its corresponding num
        new_types = []
        new_nums = []
        for t, n in zip(types, nums):
            if t.lower() != "entry":
                new_types.append(t)
                new_nums.append(n)
        location["Type"] = ", ".join(new_types)
        location["Num"] = "; ".join(new_nums)
        return location


class DataComparer():

    def __init__(self, data_path: str, similarity_method: str = "inchikey"):
        """
        Initializes the DataComparer with the ground_truth and validation datasets.

        Args:
            data_path: Path to the validation/submission JSON.
            similarity_method: How reaction SMILES are compared. ``"inchikey"`` (default)
                uses binary InChIKey-set Jaccard of reactants and products; ``"tanimoto"``
                uses graded Morgan-Tanimoto so near-identical molecules (e.g. ethyl vs
                butyl esters, stereo, salt forms) receive partial credit instead of 0.
        """
        if similarity_method not in SIMILARITY_METHODS:
            raise ValueError(
                f"similarity_method must be one of {SIMILARITY_METHODS}, got {similarity_method!r}"
            )
        self.similarity_method = similarity_method
        self.oprd_data = json_to_dict("../data/OPRD-100.json")
        self.val_data = json_to_dict(data_path)

        # check if any of the locations in validation set are not in the OPRD data
        self.validation_pairs = self.generate_location_pairs(self.val_data)
        self.oprd_pairs = self.generate_location_pairs(self.oprd_data)
        # Predicted locations that do not exist in the gold set. These are simply not
        # scored (they cannot be matched to any gold entry) and are reported for
        # transparency. Extra predictions must NOT abort scoring — a submission that
        # contains more reactions than the gold set is expected, not an error.
        self.incorrect_validation_pairs = self.validation_pairs - self.oprd_pairs
        # Only use entries which are present in the validation set
        self.oprd_subset = self.subset_json()
        # Split and sort data by location type -  only considers single location datapoints currently
        self.oprd_experimental_data = self.sort_by_location(self.split_by_location(self.oprd_subset, 'Experimental'))
        self.oprd_scheme_data = self.sort_by_location(self.split_by_location(self.oprd_subset, 'Scheme'))
        self.oprd_table_data = self.sort_by_location(self.split_by_location(self.oprd_subset, 'Table'))
        self.val_experimental_data = self.sort_by_location(self.split_by_location(self.val_data, 'Experimental'))
        self.val_scheme_data = self.sort_by_location(self.split_by_location(self.val_data, 'Scheme'))
        self.val_table_data = self.sort_by_location(self.split_by_location(self.val_data, 'Table'))


    def split_by_location(self, data: list[dict], location_type) -> list[dict]:
        """
        Splits the data by location type.
        """
        if not isinstance(data, list):
            raise ValueError("Input should be a list of dictionaries.")
        if not isinstance(location_type, str):
            raise ValueError("Location type should be a string.")
        split_data = []
        for entry in data:
            if not isinstance(entry, dict):
                raise ValueError("Each entry should be a dictionary.")
            location = entry.get('Location', {})
            if not isinstance(location, dict):
                raise ValueError("Location should be a dictionary.")
            ltype = location.get('Type', '')
            if location_type in ltype: # We will get some crossover entries here

                # Defining a hierarchy for location treatment
                if location_type == 'Scheme':
                    if 'Experimental' in ltype or 'Table' in ltype:
                        continue
                elif location_type == 'Table':
                    if 'Experimental' in ltype:
                        continue
                split_data.append(entry)
        return split_data


    def sort_by_location(self, subset):
        def get_location_key(entry):
            location = entry.get('Location', {})
            location = remove_entry_from_location(location)
            return (
                entry.get('Reference', ''),
                location.get('Type', ''),
                location.get('Num', '')
            )
        return sorted(subset, key=get_location_key)


    def generate_location_pairs(self, data: list[dict]) -> set:
        """
        Generates a set of pairs (Reference, Location) from the val data.
        """
        if not isinstance(data, list):
            raise ValueError("Input should be a list of dictionaries.")
        pairs = set()
        for entry in data:
            if not isinstance(entry, dict):
                raise ValueError("Each entry should be a dictionary.")
            if 'Reference' not in entry or 'Location' not in entry:
                raise ValueError("Each entry should contain 'Reference' and 'Location' keys.")
            location = entry['Location']
            location = remove_entry_from_location(location)
            reference = entry['Reference']
            pairs.add((reference, location.get("Type") or "", location.get("Num") or ""))
        return pairs


    def subset_json(self,) -> list[dict]:
        """
        Subset the JSON data to only include the relevant fields for comparison.
        """
        if not isinstance(self.oprd_data, list) or not isinstance(self.val_data, list):
            raise ValueError("Input should be a list of dictionaries.")
        
        # Find all the Reference - location pairs in the validation data
        pairs = self.generate_location_pairs(self.val_data)
        # Now subset the data to only include the relevant fields
        oprd_subset = []
        for entry in self.oprd_data:
            if not isinstance(entry, dict):
                raise ValueError("Each entry should be a dictionary.")
            if 'Reference' not in entry or 'Location' not in entry:
                raise ValueError("Each entry should contain 'Reference' and 'Location' keys.")
            location = entry['Location']
            location = remove_entry_from_location(location)
            reference = entry['Reference']
            # Check pairs
            if (reference, location.get("Type") or "", location.get("Num") or "") in pairs:
                oprd_subset.append(entry)
            
        return oprd_subset


    def get_entry_counts(self, val_subset: list[dict], oprd_subset: list[dict]) -> dict:
        """
        inputs two subsets of data, one from OPRD and one from val. - best to split by location first or use the attributes of this class.
        Counts the number of reaction entries in each doi-data location pair.
        returns an array with the counts for each pair.
        """
        if not isinstance(oprd_subset, list) or not isinstance(val_subset, list):
            raise ValueError("Both OPRD and val data should be lists of dictionaries.")
        # generate pairs for
        val_pairs = self.generate_location_pairs(val_subset)
        oprd_pairs = self.generate_location_pairs(oprd_subset)
        # Count over the UNION of location pairs so the comparison works even when the
        # submission and gold do not share exactly the same locations (e.g. an automated
        # extractor that misses some locations or adds extras). Pairs absent from one
        # side simply count 0 there. (Human re-extraction has identical pair sets, so
        # this is a no-op for the validation scenario.)
        all_pairs = val_pairs | oprd_pairs
        pair_counts_val = {pair: 0 for pair in all_pairs}
        pair_counts_oprd = {pair: 0 for pair in all_pairs}
        for pair in all_pairs:
            if not isinstance(pair, tuple) or len(pair) != 3:
                raise ValueError("Each pair should be a tuple of (Reference, Type, Num).")
            for entry in val_subset:
                if not isinstance(entry, dict):
                    raise ValueError("Each entry should be a dictionary.")
                if 'Reference' not in entry or 'Location' not in entry:
                    raise ValueError("Each entry should contain 'Reference' and 'Location' keys.")
                location = entry['Location']
                location = remove_entry_from_location(location)
                reference = entry['Reference']
                if (reference, location.get("Type") or "", location.get("Num") or "") == pair:
                    pair_counts_val[pair] += 1
            for entry in oprd_subset:
                if not isinstance(entry, dict):
                    raise ValueError("Each entry should be a dictionary.")
                if 'Reference' not in entry or 'Location' not in entry:
                    raise ValueError("Each entry should contain 'Reference' and 'Location' keys.")
                location = entry['Location']
                location = remove_entry_from_location(location)
                reference = entry['Reference']
                if (reference, location.get("Type") or "", location.get("Num") or "") == pair:
                    pair_counts_oprd[pair] += 1
        return pair_counts_val, pair_counts_oprd


    def compute_linear_sum_assignment(self, cost_matrix: np.ndarray) -> np.ndarray:
        """
        Computes the optimal assignment using the Hungarian algorithm (linear sum assignment).
        Returns the row and column indices of the optimal assignment.
        """
        if not isinstance(cost_matrix, np.ndarray):
            raise ValueError("Cost matrix should be a NumPy array.")
        if cost_matrix.ndim != 2:
            raise ValueError("Cost matrix should be a 2D array.")
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        return row_ind, col_ind


    def find_reaction_matches(self, val_subset: list[dict], oprd_subset: list[dict]) -> list[Tuple[int, int, float]]:
        """
        Match reactions using Hungarian algorithm within identical (Reference, Type, Num) groups.
        Avoids infeasible cost matrices when a group exists on one side only or sizes differ.
        """
        if not isinstance(oprd_subset, list) or not isinstance(val_subset, list):
            raise ValueError("Both OPRD and val data should be lists of dictionaries.")

        def key(entry):
            loc = remove_entry_from_location(entry.get('Location', {}) or {})
            return (entry.get('Reference', ''), loc.get('Type', ''), loc.get('Num', ''))

        # group indices by key
        val_groups: dict[tuple, list[int]] = {}
        for i, e in enumerate(val_subset):
            val_groups.setdefault(key(e), []).append(i)
        oprd_groups: dict[tuple, list[int]] = {}
        for j, e in enumerate(oprd_subset):
            oprd_groups.setdefault(key(e), []).append(j)

        matches: list[Tuple[int, int, float]] = []

        for k in sorted(set(val_groups) & set(oprd_groups)):
            r_idx_list = val_groups[k]
            o_idx_list = oprd_groups[k]

            val_block = [val_subset[i] for i in r_idx_list]
            oprd_block = [oprd_subset[j] for j in o_idx_list]

            # similarity block
            S = self.compare_data(val_block, oprd_block, metric='total_similarity', return_array=True)
            if S.size == 0:
                continue

            # minimize cost = 1 - similarity
            cost = 1.0 - S

            # sanity check: finite in every row/col
            bad_rows = np.where(~np.isfinite(cost).any(axis=1))[0]
            bad_cols = np.where(~np.isfinite(cost).any(axis=0))[0]
            if bad_rows.size or bad_cols.size:
                # skip infeasible submatrix gracefully
                continue

            row_ind, col_ind = self.compute_linear_sum_assignment(cost)

            for r_local, o_local in zip(row_ind, col_ind):
                r_global = r_idx_list[r_local]
                o_global = o_idx_list[o_local]
                matches.append((r_global, o_global, float(S[r_local, o_local])))

        return matches


    def compute_comparison_scores(self, val_subset: list[dict], oprd_subset: list[dict]) -> pd.DataFrame:
        """
        Computes scores for all metrics and returns a dataframe with the results.
        """
        if not isinstance(oprd_subset, list) or not isinstance(val_subset, list):
            raise ValueError("Both OPRD and val data should be lists of dictionaries.")

        matched_entries = self.find_reaction_matches(val_subset, oprd_subset)

        metrics = [
            'reaction_smiles_similarity',
            'reaction_steps_similarity',
            'yield_similarity',
            'reagent_similarity',
            'reagent_name_similarity',
            'reagent_amount_similarity',
            'solvent_similarity',
            'time_similarity',
            'temperature_similarity',
            'total_similarity',
        ]

        # Initialize a dictionary to hold the results
        results = {metric: [] for metric in metrics}

        val_indices = []
        oprd_indices = []
        similarity_scores = []
        for r_idx, o_idx, sim_score in matched_entries:
            val_entry = val_subset[r_idx]
            oprd_entry = oprd_subset[o_idx]
            val_indices.append(r_idx)
            oprd_indices.append(o_idx)
            similarity_scores.append(sim_score)
            comparison = CompareReactionEntries(oprd_entry, val_entry, similarity_method=self.similarity_method)
            for metric in metrics:
                results[metric].append(getattr(comparison, metric, 0))

        # Convert the results to a DataFrame
        df = pd.DataFrame(results)
        # add oprd and val indices and similarity scores
        df['val_index'] = val_indices
        df['oprd_index'] = oprd_indices
    
        return df



    def compare_data(self, val_subset: list[dict], oprd_subset: list[dict], metric: str, return_array: bool = True, return_max: bool = False) -> np.ndarray:
        """
        Compares the OPRD and val data and returns a NumPy array with comparison results.

        valid metrics:
        - 'reaction_smiles_similarity': Similarity of reaction SMILES.
        - 'reaction_steps_similarity': Similarity of reaction steps.
        - 'yield_similarity': Similarity of yields.
        - 'reagent_similarity': Similarity of reagent names and amounts.
        - 'reagent_name_similarity': Similarity of reagent names.
        - 'reagent_amount_similarity': Similarity of reagent amounts (excluded from total sim).
        - 'solvent_similarity': Similarity of solvents (excluded from total sim).
        - 'time_similarity': Similarity of reaction times.
        - 'temperature_similarity': Similarity of reaction temperatures.
        - 'total_similarity': Overall similarity score combining all metrics.
        """
        if not isinstance(oprd_subset, list) or not isinstance(val_subset, list):
            raise ValueError("Both OPRD and val data should be lists of dictionaries.")

        valid_metrics = {
            'reaction_smiles_similarity',
            'reaction_steps_similarity',
            'yield_similarity',
            'reagent_similarity',
            'reagent_name_similarity',
            'reagent_amount_similarity',
            'solvent_similarity',
            'time_similarity',
            'temperature_similarity',
            'total_similarity',
        }
        if metric not in valid_metrics:
            raise ValueError(f"Metric '{metric}' not recognized. Must be one of {valid_metrics}.")

        if return_array:
            comparison_results = []
            for val_entry in val_subset:
                results = []
                for oprd_entry in oprd_subset:
                    comparison = CompareReactionEntries(oprd_entry, val_entry, similarity_method=self.similarity_method)
                    results.append(getattr(comparison, metric, 0))
                comparison_results.append(results)
            return np.array(comparison_results)
        elif return_max:
            max_results = []
            for val_entry in val_subset:
                max_score = 0
                for oprd_entry in oprd_subset:
                    comparison = CompareReactionEntries(oprd_entry, val_entry, similarity_method=self.similarity_method)
                    score = getattr(comparison, metric, 0)
                    if score > max_score:
                        max_score = score
                max_results.append(max_score)
            return np.array(max_results)
        else:
            return None

    def generate_heatmap(
        self,
        comparison_results: np.ndarray,
        val_subset: list[dict],
        oprd_subset: list[dict],
        show_plot: bool = True,
        save_path: str = "heatmap.png",
        title: str = "Comparison Heatmap",
        xlabel: str = "Validation Entries",
        ylabel: str = "OPRD Entries",
        line_col="white",
        figsize: tuple = (6, 6),
        title_fontdict: dict = {'fontsize': 24, 'fontweight': 'bold'},
        label_fontdict: dict = {'fontsize': 20 },
        tick_fontdict: dict = {'fontsize': 18},
        annot=False
    ):
        """
        Generates a heatmap from the comparison results, with grouped reference labels and range bars on axes.
        """
        import seaborn as sns
        import matplotlib.pyplot as plt
        import numpy as np

        def get_reference_groups(subset):
            refs = [entry.get('Reference', '') for entry in subset]
            groups = []
            prev_ref = None
            start = 0
            for i, ref in enumerate(refs + [None]):  # Add sentinel
                if ref != prev_ref:
                    if prev_ref is not None:
                        groups.append((prev_ref, start, i - 1))
                    start = i
                    prev_ref = ref
            return groups

        # Get groups for x and y axes
        x_groups = get_reference_groups(val_subset)
        y_groups = get_reference_groups(oprd_subset)

        # Compute tick positions and labels
        def get_ticks_and_labels(groups):
            ticks = []
            labels = []
            for ref, start, end in groups:
                ticks.append((start + end + 1) / 2)
                labels.append(ref)
            return ticks, labels

        x_ticks, x_labels = get_ticks_and_labels(x_groups)
        y_ticks, y_labels = get_ticks_and_labels(y_groups)

        plt.figure(figsize=figsize)
        ax = sns.heatmap(comparison_results, annot=annot, fmt=".2f", cmap='viridis')
        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize=tick_fontdict.get("fontsize", 14))
        cbar.set_label("Total Similarity Score", fontsize=label_fontdict.get("fontsize", 16), weight=label_fontdict.get("fontweight", "normal"))
        plt.title(title, **title_fontdict)
        plt.xlabel(xlabel, **label_fontdict)
        plt.ylabel(ylabel, **label_fontdict)

        # Set ticks and labels
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_labels, rotation=45, ha='right', **tick_fontdict)
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_labels, rotation=0, **tick_fontdict)

        # Add range bars for x (vertical lines)
        for _, start, end in x_groups:
            ax.vlines([start, end + 1], *ax.get_ylim(), colors=line_col, linestyles='dashed', linewidth=1)

        # Add range bars for y (horizontal lines)
        for _, start, end in y_groups:
            ax.hlines([start, end + 1], *ax.get_xlim(), colors=line_col, linestyles='dashed', linewidth=1)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300)
        if show_plot:
            plt.show()





class CompareReactionEntries():


    def __init__(self, rxn1: dict, rxn2: dict, similarity_method: str = "inchikey"):
        """
        Initializes the DataComparer with paths to the OPRD and val data files.

        Args:
            rxn1: First reaction entry (ground truth).
            rxn2: Second reaction entry (validation/submission).
            similarity_method: ``"inchikey"`` (default) for binary InChIKey-set Jaccard,
                or ``"tanimoto"`` for graded Morgan-Tanimoto reaction-SMILES similarity.
        """
        self.rxn1 = rxn1
        self.rxn2 = rxn2
        self.similarity_method = similarity_method
        self.reaction_smiles_similarity = self.compare_reaction_smiles(rxn1, rxn2)
        self.reaction_steps_similarity = self.compare_reaction_steps(rxn1, rxn2)
        self.yield_similarity = self.compare_yields(rxn1, rxn2)
        self.reagent_name_similarity, self.reagent_amount_similarity = self.compare_reagents_solvents(rxn1, rxn2, key_name="Reagents") 
        self.reagent_similarity = self.reagent_name_similarity * self.reagent_amount_similarity
        self.solvent_name_similarity, self.solvent_amount_similarity = self.compare_reagents_solvents(rxn1, rxn2, key_name="Solvents")
        self.solvent_similarity = self.solvent_name_similarity * self.solvent_amount_similarity
        self.time_similarity = self.compare_times(rxn1, rxn2)
        self.temperature_similarity = self.compare_temperatures(rxn1, rxn2)
        self.total_similarity = (
            self.reaction_smiles_similarity + 
            self.reaction_steps_similarity +
            self.yield_similarity +
            self.reagent_similarity +
            self.solvent_similarity +
            self.time_similarity +
            self.temperature_similarity
        ) / 7
        data_location_match = self.compare_location_data(rxn1, rxn2)
        if not data_location_match:
            self.total_similarity = 0
            self.location_similarity = 0
            self.reaction_smiles_similarity = 0
            self.reaction_steps_similarity = 0
            self.yield_similarity = 0
            self.reagent_similarity = 0
            self.solvent_similarity = 0
            self.time_similarity = 0
            self.temperature_similarity = 0

    def compare_location_data(self, rxn1: dict, rxn2: dict):
        """
        Compares the location data of two reactions.
        Returns a score between 0 and 1, with 1 being a perfect match.
        """
        if not isinstance(rxn1, dict) or not isinstance(rxn2, dict):
            raise ValueError("Both reactions should be dictionaries.")
        # Compare the Location data
        loc1 = rxn1.get('Location', {})
        ref1 = rxn1.get('Reference', '')
        # Remove the 'entry' from the location if it exists
        loc1 = remove_entry_from_location(loc1)
        loc1['Reference'] = ref1
        loc2 = rxn2.get('Location', {})
        ref2 = rxn2.get('Reference', '')
        loc2 = remove_entry_from_location(loc2)
        loc2['Reference'] = ref2
        if not isinstance(loc1, dict) or not isinstance(loc2, dict):
            raise ValueError("Location data should be dictionaries.")
        if loc1 == loc2:
            return True
        return False

    def create_reaction_identifiers(self, data: list[dict]):
        """
        Extracts all molecule information from reaction SMILES and returns sets of InChIKeys for reactants and products.
        """
        r_ids = set()
        p_ids = set()
        for entry in data:
            rxn_smi = entry.get('Reaction', '')
            if rxn_smi:
                rxn = Chem.rdChemReactions.ReactionFromSmarts(rxn_smi, useSmiles=True)
                for mol in rxn.GetReactants():
                    r_ids.add(Chem.MolToInchiKey(mol))
                for mol in rxn.GetProducts():
                    p_ids.add(Chem.MolToInchiKey(mol))
        # sort the identifiers for consistency
        r_ids = sorted(r_ids)
        p_ids = sorted(p_ids)
        return r_ids, p_ids

    def set_similarity_score(self, set1: set, set2: set) -> float:
        """
        Compare similarity of two sets.
        """
        if not set1 and not set2:
            return 1.0  # Both sets empty, consider as perfect match
        intersection = set1 & set2
        union = set1 | set2
        return len(intersection) / len(union)

    def compare_reaction_smiles(self, rxn1: dict, rxn2: dict):
        """
        Compares two reaction entries and returns a number between 0 and 1 with 1 being a perfect match.

        Uses the method selected at construction time: ``"inchikey"`` (binary InChIKey-set
        Jaccard) or ``"tanimoto"`` (graded Morgan-Tanimoto).
        """
        if not isinstance(rxn1, dict) or not isinstance(rxn2, dict):
            raise ValueError("Both reactions should be dictionaries.")
        if getattr(self, "similarity_method", "inchikey") == "tanimoto":
            return self.compare_reaction_smiles_tanimoto(rxn1, rxn2)
        # Compare the InChI keys of the reactants and products
        r1_rids, r1_pids = self.create_reaction_identifiers([rxn1])
        r2_rids, r2_pids = self.create_reaction_identifiers([rxn2])
        reactant_similarity = self.set_similarity_score(set(r1_rids), set(r2_rids))
        product_similarity = self.set_similarity_score(set(r1_pids), set(r2_pids))
        reaction_smiles_similarity = (reactant_similarity + product_similarity) / 2
        return reaction_smiles_similarity

    def _component_fingerprints(self, side_smiles: str) -> list:
        """Morgan fingerprints for each dot-separated molecule on one side of a reaction."""
        fps = []
        for part in (side_smiles or "").split("."):
            part = part.strip()
            if not part:
                continue
            mol = Chem.MolFromSmiles(part)
            if mol is not None:
                fps.append(_MORGAN_GENERATOR.GetFingerprint(mol))
        return fps

    def soft_set_similarity(self, fps_a: list, fps_b: list) -> float:
        """Symmetric graded set similarity: mean of best-match Tanimoto in both directions.

        Reduces to 1.0 for identical molecule sets and to a proportion for partial
        overlap. Empty-vs-empty is 1.0; empty-vs-non-empty is 0.0.
        """
        if not fps_a and not fps_b:
            return 1.0
        if not fps_a or not fps_b:
            return 0.0
        a_best = [max(TanimotoSimilarity(a, b) for b in fps_b) for a in fps_a]
        b_best = [max(TanimotoSimilarity(a, b) for a in fps_a) for b in fps_b]
        return (sum(a_best) + sum(b_best)) / (len(a_best) + len(b_best))

    def compare_reaction_smiles_tanimoto(self, rxn1: dict, rxn2: dict) -> float:
        """Graded reaction-SMILES similarity via Morgan-Tanimoto.

        Mirrors :meth:`compare_reaction_smiles` (averaging reactant and product identity)
        but uses graded Morgan-Tanimoto on each side rather than binary InChIKey-set
        Jaccard, so near-identical molecules receive partial credit instead of 0.
        """
        smi1 = rxn1.get("Reaction", "") or ""
        smi2 = rxn2.get("Reaction", "") or ""
        if ">>" not in smi1 or ">>" not in smi2:
            return 0.0
        r1, p1 = smi1.split(">>", 1)
        r2, p2 = smi2.split(">>", 1)
        reactant_similarity = self.soft_set_similarity(
            self._component_fingerprints(r1), self._component_fingerprints(r2)
        )
        product_similarity = self.soft_set_similarity(
            self._component_fingerprints(p1), self._component_fingerprints(p2)
        )
        return (reactant_similarity + product_similarity) / 2

    def compare_reaction_steps(self, rxn1: dict, rxn2: dict):
        """
        Compare number of steps in two reactions.
        """
        if not isinstance(rxn1, dict) or not isinstance(rxn2, dict):
            raise ValueError("Both reactions should be dictionaries.")
        # get step data
        step1 = rxn1.get('Steps', [])
        step2 = rxn2.get('Steps', [])
        if not isinstance(step1, list) or not isinstance(step2, list):
            raise ValueError("Steps should be a list of dictionaries.")
        if len(step1) != len(step2):
            return 0
        else:
            return 1
    

    def build_entity_amounts(self, steps: list, key_name="Reagents") -> dict[str, dict[str, list]]:
        """
        Get unique reagents/solvents and amounts lists from all steps.
        Splits '+'-concatenated amount strings into separate entries.
        """
        results: dict[str, dict[str, list]] = {}

        for step in steps:
            items = step.get(key_name, [])
            if not isinstance(items, list):
                continue
            for d in items:
                if not isinstance(d, dict):
                    raise ValueError("Each item must be a dict.")
                amounts = d.get("Amounts", {})
                if not isinstance(amounts, dict):
                    amounts = {}
                key = d.get(key_name[:-1])
                if key:
                    bucket = results.setdefault(key, {})
                    # Ensure keys exist
                    for akey in amounts.keys():
                        bucket.setdefault(akey, [])
                    # Append non-empty values (split on '+')
                    for akey, aval in amounts.items():
                        if aval is None:
                            continue
                        # coerce to string; split additive tokens "x + y"
                        parts = str(aval).split("+")
                        for part in parts:
                            sval = part.strip()
                            if not sval:
                                continue
                            bucket[akey].append(sval)
        return results

    def standardize_units(self, value: str) -> str:
        """
        Standardize common unit variations to a canonical form.
        """
        if not isinstance(value, str):
            return value
        s = value.strip().lower()
        s = re.sub(r'\s+', ' ', s)  # Normalize whitespace
        # Common unit standardizations
        replacements = {
            'ml': 'mL',
            'milliliter': 'mL',
            'milliliters': 'mL',
            'l': 'L',
            'liter': 'L',
            'liters': 'L',
            'g': 'g',
            'gram': 'g',
            'grams': 'g',
            'kg': 'kg',
            'kilogram': 'kg',
            'kilograms': 'kg',
            '%': '%',
            'percent': '%',
            'equiv': 'equiv',
            'equivalents': 'equiv',
            'Eq.': 'equiv',
            'eq': 'equiv',
            'eq.': 'equiv',
            'mol': 'mol',
            'mole': 'mol',
            'moles': 'mol',
            'mol%': 'mol%',
            'mole percent': 'mol%',
            'mmol': 'mmol',
            'millimole': 'mmol',
            'millimoles': 'mmol',
            # Add more as needed
        }
        for old, new in replacements.items():
            s = re.sub(r'\b' + re.escape(old) + r'\b', new, s)
        return s

    def parse_mass_to_grams(self, value) -> str | None:
            """
            Parse a mass and return number of grams g as a str.
            Supports: mg, g, ug/µg/μg, kg. If numeric, assume grams.
            Examples:
            "5 mg" -> 0.005
            "1.2 g" -> 1.2
            "250 ug" -> 0.00025
            "2 kg" -> 2000.0
            "1,200 mg" -> 1.2
            """
            if value is None:
                return None
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                return float(value)
            if not isinstance(value, str):
                return None

            s = value.strip()
            if not s:
                return value
            s = s.replace("µ", "u").replace("μ", "u")
            s = re.sub(r'(?<=\d),(?=\d{3}\b)', '', s)

            m = re.search(r'(?P<num>\d+(?:\.\d+)?)\s*(?P<unit>mg|kg|ug|g)\b', s, flags=re.I)
            if not m:
                return value

            num = float(m.group("num"))
            unit = m.group("unit").lower()
            if unit == "g":
                return str(num) + " g"
            if unit == "mg":
                return str(num / 1_000.0) + " g"
            if unit == "ug":
                return str(num / 1_000_000.0) + " g"
            if unit == "kg":
                return str(num * 1_000.0) + " g"
            return value

    def parse_amount_to_moles(self, value) -> str | None:
        """
        Parse an amount and return a canonical string in mol (e.g., '0.0005 mol').
        Supports: umol/µmol/μmol, mmol, mol. If numeric, assume mol.
        Examples:
        "0.5 mmol" -> "0.0005 mol"
        "250 umol" -> "0.00025 mol"
        "1 mol" -> "1.0 mol"
        0.002 -> "0.002 mol"
        """
        if value is None:
            return None
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            return str(float(value)) + " mol"
        if not isinstance(value, str):
            return None

        s = value.strip()
        if not s:
            return value
        s = s.replace("µ", "u").replace("μ", "u")
        s = re.sub(r'(?<=\d),(?=\d{3}\b)', '', s)

        m = re.search(
            r'(?P<num>\d+(?:\.\d+)?)\s*(?P<unit>umol|mmol|mol)\b'
            r'(?!\s*(?:%|percent\b|/|per\b|eq\b|equiv(?:alent)?s?\b|(?:L|kg|g)\s*[-−]?\s*1\b))',
            s,
            flags=re.I
        )
        if not m:
            return value

        num = float(m.group("num"))
        unit = m.group("unit").lower()
        if unit == "mol":
            return f"{num} mol"
        if unit == "mmol":
            return f"{num / 1_000.0} mol"
        if unit == "umol":
            return f"{num / 1_000_000.0} mol"
        return value

    def parse_volume_to_liters(self, value) -> str | None:
            """
            Parse a volume and return liters (L) as float.
            Supports: uL/µL/μL, mL, L. If numeric, assume liters.
            Examples:
            "500 uL" -> 0.0005
            "1.25 mL" -> 0.00125
            "2 L" -> 2.0
            "1,200 uL" -> 0.0012
            0.05 -> 0.05
            """
            if value is None:
                return value
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                return float(value)
            if not isinstance(value, str):
                return value

            s = value.strip()
            if not s:
                return value
            s = s.replace("µ", "u").replace("μ", "u")
            s = re.sub(r'(?<=\d),(?=\d{3}\b)', '', s)

            m = re.search(r'(?P<num>\d+(?:\.\d+)?)\s*(?P<unit>ul|ml|l)\b', s, flags=re.I)
            if not m:
                return value
            
            num = float(m.group("num"))
            unit = m.group("unit").lower()
            if unit == "l":
                return str(num) + " L"
            if unit == "ml":
                return str(num / 1_000.0) + " L"
            if unit == "ul":
                return str(num / 1_000_000.0) + " L"
            return value

    def simplify_entity_amounts(self, amounts:dict)->dict:
        """
        Simplify amounts by parsing and standardizing units.
        Splits '+'-concatenated strings into separate entries.
        Converts to canonical units:
        - Mass -> grams (g)
        - Moles -> moles (mol)
        - Volume -> liters (L)
        """
        simplified = defaultdict(dict)
        for name, amounts_dict in amounts.items():
            for key, value_list in amounts_dict.items():
                simplified[name][key] = value_list
                new_values = []
                for v in value_list:
                    parts = str(v).split("+")
                    for part in parts:
                        part = part.strip()
                        if not part:
                            continue
                        part = self.standardize_units(part)
                        if key.lower() in ['mass']:
                            part = self.parse_mass_to_grams(part)
                        elif key.lower() in ['moles']:
                            part = self.parse_amount_to_moles(part)
                        elif key.lower() in ['volume']:
                            part = self.parse_volume_to_liters(part)
                        if part:
                            new_values.append(part)
                simplified[name][key] = new_values
        return simplified

    def consolidate_entity_amounts(self, simplified: dict, sig_digits: int = 3) -> dict:
        """
        Sum like units within each value list per entity and key.
        Supports canonical tokens produced by simplify_entity_amounts:
        - "<number> g", "<number> mol", "<number> L", "<number> equiv"
        Non-matching tokens are preserved (deduped) and appended.
        """

        def _format_sig(x: float, sig_digits: int = 3) -> str:
            # Format with significant digits, strip trailing zeros and dot
            s = f"{x:.{sig_digits}g}"
            # Normalize 'e' to lowercase and ensure no trailing '.'
            return s

        unit_order = {"g": 0, "mol": 1, "L": 2, "equiv": 3}
        num_unit_re = re.compile(r'([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*(g|mol|l|L|equiv)\b', re.I)

        out = defaultdict(dict)
        for name, amounts_dict in (simplified or {}).items():
            for key, values in (amounts_dict or {}).items():
                if not isinstance(values, list):
                    out[name][key] = values
                    continue

                sums = defaultdict(float)
                others = []
                seen_other = set()

                for v in values:
                    if v is None:
                        continue
                    s = str(v).strip()
                    m = num_unit_re.fullmatch(s)
                    if m:
                        num = float(m.group(1))
                        unit = m.group(2)
                        # normalize unit
                        unit = "L" if unit.lower() == "l" else unit.lower()
                        if unit in ("g", "mol", "equiv") or unit == "L":
                            sums[unit] += num
                        else:
                            if s and s not in seen_other:
                                seen_other.add(s)
                                others.append(s)
                    else:
                        if s and s not in seen_other:
                            seen_other.add(s)
                            others.append(s)

                # Build consolidated values (ordered units first, then others)
                consolidated_values = []
                for unit in sorted(sums.keys(), key=lambda u: unit_order.get(u, 99)):
                    total = sums[unit]
                    if abs(total) > 0:
                        consolidated_values.append(f"{_format_sig(total, sig_digits)} {unit}")

                consolidated_values.extend(others)
                out[name][key] = consolidated_values

        return out

    def standardise_reagents_solvents(self, rxn, key_name)->dict:
        steps = rxn.get('Steps', [])
        reagents = self.build_entity_amounts(steps, key_name=key_name)
        reagents_simplified = self.simplify_entity_amounts(reagents)
        reagents_consolidated = self.consolidate_entity_amounts(reagents_simplified)
        return reagents_consolidated
    

    def compare_reagents_solvents(self, rxn1: dict, rxn2: dict, key_name="Reagents", return_names=False) -> Tuple[float, float]:
        if not isinstance(rxn1, dict) or not isinstance(rxn2, dict):
            raise ValueError("Both reactions should be dictionaries.")
        rgts_1 = self.standardise_reagents_solvents(rxn1, key_name=key_name)
        rgts_2 = self.standardise_reagents_solvents(rxn2, key_name=key_name)

        names_1 = set(i.lower().strip() for i in rgts_1.keys())
        names_2 = set(i.lower().strip() for i in rgts_2.keys())

        # if there are no names in either, return perfect match
        if not names_1 and not names_2:
            if return_names:
                return 1.0, 1.0, set(), set(), set(), set()
            return 1.0, 1.0

        name_score = self.set_similarity_score(names_1, names_2)

        inter = names_1 & names_2
        if not inter and not return_names:
            return name_score, 0.0
        
        diff = names_1 ^ names_2
        if not inter and return_names:
            return name_score, 0.0, names_1, names_2, set(), diff
        

        per_name_scores = []
        for nm in inter:
            a1 = rgts_1.get(nm, {})
            a2 = rgts_2.get(nm, {})
            per_name_scores.append(self.dict_similarity_score(a1, a2))
        amount_score = sum(per_name_scores) / len(per_name_scores) if per_name_scores else 0.0

        if return_names:
            if inter:
                return name_score, amount_score, names_1, names_2, inter, diff
            return name_score, amount_score, names_1, names_2, set(), diff

        return name_score, amount_score





    def _is_missing_value(self, v) -> bool:
        """
        Treat None and empty strings as missing.
        Numeric zeros (0, 0.0, "0", "0.0", etc.) are NOT missing.
        Non-empty strings are considered present.
        Empty containers are considered missing.
        """
        if v is None:
            return True
        if isinstance(v, (list, tuple, set, dict)) and len(v) == 0:
            return True
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            return False  # numeric present, including 0
        if isinstance(v, str):
            s = v.strip().lower()
            if s in {"", "n/a", "na", "none", "null"}:
                return True
            # numeric string (including "0", "0.0") should be present
            try:
                float(s)
                return False
            except Exception:
                return False  # non-empty non-numeric string is present
        return False
    
    def dict_similarity_score(self, d1: dict, d2: dict) -> float:
        """
        General scoring function to compare two dictionaries.
        - Returns 1.0 for exact match.
        - Returns partial score for partial overlap (intersection over union of non-null keys with matching values).
        - Returns 0.0 for no match.
        """
        if not isinstance(d1, dict) or not isinstance(d2, dict):
            raise ValueError("Both inputs should be dictionaries.")
        # Exact match
        if d1 == d2:
            return 1.0
        # Partial match
        non_null_keys1 = {k for k, v in d1.items() if not self._is_missing_value(v)}
        non_null_keys2 = {k for k, v in d2.items() if not self._is_missing_value(v)}
        all_non_null_keys = non_null_keys1.union(non_null_keys2)
        # if all keys are null in both dicts, return 1.0
        if not all_non_null_keys:
            return 1.0
        shared_keys = non_null_keys1.intersection(non_null_keys2)
        matches = [k for k in shared_keys if d1[k] == d2[k]]
        if matches and all_non_null_keys:
            partial_match_score = len(matches) / len(all_non_null_keys)
            return partial_match_score
        return 0.0

    def list_similarity_score(self, list1, list2):
        """
        Compare similarity of two lists.
        """
        if not list1 and not list2:
            return 1.0
        intersection = set(list1) & set(list2)
        union = set(list1) | set(list2)
        return len(intersection) / len(union)

    def ordered_list_similarity_score(self, list1, list2):
        """
        Compare two lists element-wise (order matters).
        Returns the fraction of matching elements (1.0 = perfect match, 0.0 = no match).
        If lists are of different lengths, compares up to the shorter length.
        """
        if not list1 and not list2:
            return 1.0
        min_len = min(len(list1), len(list2))
        if min_len == 0:
            return 0.0
        matches = sum(1 for i in range(min_len) if list1[i] == list2[i])
        return matches / max(len(list1), len(list2))

    def compare_yields(self, rxn1: dict, rxn2: dict, return_dict=False):
        """
        Compare final yields of two reactions.
        """
        if not isinstance(rxn1, dict) or not isinstance(rxn2, dict):
            raise ValueError("Both reactions should be dictionaries.")
        steps1 = rxn1.get('Steps', [])
        final_step1 = steps1[-1] if steps1 else {}
        steps2 = rxn2.get('Steps', [])
        final_step2 = steps2[-1] if steps2 else {}
        yield1 = final_step1.get('Yield', {})
        yield2 = final_step2.get('Yield', {})
        similarity = self.dict_similarity_score(yield1, yield2)
        if return_dict:
            return similarity, yield1, yield2
        else:
            return similarity

    
    def sum_all_times(self, times: list) -> float:
        """
        Sum all time values in a list of time dictionaries.
        """
        total_time = 0.0
        for time in times:
            if time is None:
                continue
            parts = time.split("+")
            for part in parts:
                part = part.strip()
                if not part:
                    continue
                m = re.match(r'(?P<value>\d+(?:\.\d+)?)\s*(?P<unit>\w+)', part)
                if not m:
                    continue
                time_value = m.group('value')
                time_unit = m.group('unit').lower()
                if time_unit in ['m','min', 'mins', 'minute', 'minutes']:
                    total_time += float(time_value)
                elif time_unit in ['h', 'hr', 'hrs', 'hour', 'hours']:
                    total_time += float(time_value) * 60.0
                elif time_unit in ['s', 'sec', 'secs', 'second', 'seconds']:
                    total_time += float(time_value) / 60.0
                elif time_unit in ['day', 'days']:
                    total_time += float(time_value) * 1440.0  # 24*60
                else:
                    # Unknown unit, skip
                    continue
        return total_time

    
    def compare_times(self, rxn1: dict, rxn2: dict):
        """
        Compare times in two reactions.
        """
        if not isinstance(rxn1, dict) or not isinstance(rxn2, dict):
            raise ValueError("Both reactions should be dictionaries.")
        # get all step information
        steps1 = rxn1.get('Steps', [])
        final_step1 = steps1[-1] if steps1 else {}
        steps2 = rxn2.get('Steps', [])
        times1 = []
        for step in steps1:
            if not isinstance(step, dict):
                raise ValueError("Steps should be a list of dictionaries.")
            time = step.get('Time', {})
            times1.append(time)
        total_time1 = self.sum_all_times(times1) # times in minutes
        times2 = []
        for step in steps2:
            if not isinstance(step, dict):
                raise ValueError("Steps should be a list of dictionaries.")
            time = step.get('Time', {})
            times2.append(time)
        total_time2 = self.sum_all_times(times2) # times in minutes
        # If both times are zero, return 1.0
        if total_time1 == 0 and total_time2 == 0:
            return 1.0
        # If one time is zero and the other is not, return 0.0
        if total_time1 == 0 or total_time2 == 0:
            return 0.0
        # compute time score as 1-abs difference over max
        time_score = 1.0 - abs(total_time1 - total_time2) / max(total_time1, total_time2)
        if time_score < 0:
            time_score = 0.0
        return time_score



    def simplify_temperature(self, values: list) -> str | None:
        """
        Take minimum and maximum temperatures and return a simplified string.
        """
        #print(values)
        if not values:
            return None
        # sub reflux, room temperature, ambient, ice bath
        new_values = []
        for v in values:
            if v is None:
                continue
            v = v.lower()
            v = re.sub(r'100', '100', v)
            v = re.sub(r'ambient', '25', v)
            v = re.sub(r'reflux', '3000', v)
            v = re.sub(r'room temperature', '25', v)
            v = re.sub(r'rt\b', '25', v)
            v = re.sub(r'ice bath', '0', v)
            new_values.append(v)
        # extract the numeric part of the temperature values
        temps = re.findall(r'(-?\d+)', ' '.join(new_values))
        if not temps:
            return None
        min_temp = min(temps, key=float)
        max_temp = max(temps, key=float)
        return [min_temp, max_temp]

    def compare_temperatures(self, rxn1: dict, rxn2: dict):
        """
        Compare temperatures in two reactions.
        """
        if not isinstance(rxn1, dict) or not isinstance(rxn2, dict):
            raise ValueError("Both reactions should be dictionaries.")
        # get all step information
        steps1 = rxn1.get('Steps', [])
        final_step1 = steps1[-1] if steps1 else {}
        steps2 = rxn2.get('Steps', [])
        temperatures1 = []
        for step in steps1:
            if not isinstance(step, dict):
                raise ValueError("Steps should be a list of dictionaries.")
            temperature = step.get('Temperature', {})
            temperatures1.append(temperature)
        temperatures2 = []
        for step in steps2:
            if not isinstance(step, dict):
                raise ValueError("Steps should be a list of dictionaries.")
            temperature = step.get('Temperature', {})
            temperatures2.append(temperature)
        if temperatures1 == [] and temperatures2 == []:
            return 1.0
        if temperatures1 == [] or temperatures2 == []:
            return 0.0
        # simplify temperature values
        simplified_temps1 = self.simplify_temperature(temperatures1)
        simplified_temps2 = self.simplify_temperature(temperatures2)
        if simplified_temps1 is None and simplified_temps2 is None:
            return 1.0
        if simplified_temps1 is None or simplified_temps2 is None:
            return 0.0
        if simplified_temps1 == simplified_temps2:
            return 1.0
        else:
            return 0.0


class PubChemResolver:
    """
    Standardize reagent identifiers with PubChem and cache results.
    Supports lookup by written token, name, InChIKey, SMILES, formula, or CID.
    """
    _cache: dict[str, dict] = {}
    _alias_cs: dict[str, str] = {}  # exact match (SMILES, InChIKey, formula, CID str) -> written key
    _alias_ci: dict[str, str] = {}  # case-insensitive (written, standard name) -> written key

    def __init__(self, cache_filename='pubchem_cache.json'):
        self.cache_filename = cache_filename
        self.check_cache_file_exists()

    def _norm_name(self, s: str) -> str:
        return s.strip().lower() if isinstance(s, str) else s

    def _rebuild_alias_indexes(self) -> None:
        self._alias_cs.clear()
        self._alias_ci.clear()
        for written, v in (self._cache or {}).items():
            name = (v.get('iupac_name') or '').strip()
            inchikey = (v.get('inchikey') or '').strip()
            smiles = (v.get('smiles') or '').strip()
            common_name = (v.get('common_name') or '').strip()
            cid = v.get('cid')
            cid_str = str(cid).strip() if cid not in (None, '') else ''

            # Case-insensitive aliases (names)
            for alias in filter(None, {written, name, common_name}):
                self._alias_ci[self._norm_name(alias)] = written

            # Case-sensitive aliases (structure keys)
            for alias in filter(None, {inchikey, smiles, cid_str}):
                self._alias_cs[alias] = written

    def check_cache_file_exists(self):
        """
        Load/create cache JSON and rebuild alias indexes.
        """
        if not os.path.exists(self.cache_filename):
            with open(self.cache_filename, 'w') as f:
                json.dump({}, f)
            self._cache = {}
        else:
            try:
                with open(self.cache_filename, 'r') as f:
                    self._cache = json.load(f)
            except json.JSONDecodeError:
                self._cache = {}
        self._rebuild_alias_indexes()

    def add_to_cache(self, dict_to_add: dict) -> None:
        """
        Merge new entries into cache and persist.
        dict_to_add: {written: {'iupac_name','inchikey','smiles', 'common_name','cid'}}
        """
        self._cache.update(dict_to_add)
        self._rebuild_alias_indexes()
        with open(self.cache_filename, 'w') as f:
            json.dump(self._cache, f, indent=2)

    def lookup(self, query: str) -> dict | None:
        """
        Lookup by any alias:
        - Case-sensitive: SMILES, InChIKey, formula, CID string
        - Case-insensitive: written token, standardized name
        """
        if not isinstance(query, str) or not query.strip():
            return None
        q = query.strip()

        # Try exact-match aliases first (do NOT lowercase)
        written = self._alias_cs.get(q)
        if written:
            return self._cache.get(written)

        # Then try case-insensitive name aliases
        written = self._alias_ci.get(self._norm_name(q))
        if written:
            return self._cache.get(written)

        # Finally, direct key hit (as written)
        return self._cache.get(q)

    def resolve(self, token: str, output_format: str = 'name') -> str | None:
        """
        Resolve token via cache or PubChem. output_format in
        {'iupac_name','inchikey','smiles','common_name','cid'}.
        """
        if not isinstance(token, str) or not token.strip():
            return None

        # 1) Try cache
        rec = self.lookup(token)
        if rec:
            return rec.get(output_format) or None

        # 2) Query PubChem using multiple namespaces
        for ns in ('inchikey', 'smiles', 'formula', 'name', 'cid'):
            try:
                comps = pcp.get_compounds(token, ns)
                if not comps:
                    continue
                c = comps[0]
                rec = {
                    'iupac_name': getattr(c, 'iupac_name', None) or '',
                    'inchikey': getattr(c, 'inchikey', None) or '',
                    'smiles': getattr(c, 'smiles', None) or '',
                    'common_name': getattr(c, 'synonyms', [None])[0] or '',
                    'cid': getattr(c, 'cid', None) or '',
                }
                # Use original token as the written key
                self.add_to_cache({token: rec})
                return rec.get(output_format) or None
            except Exception:
                continue
        # 3) Not found - cache empty record
        self.add_to_cache({token: {'iupac_name':None,
                    'inchikey':None,
                    'smiles': None,
                    'common_name': None,
                    'cid': None}})
        return None

