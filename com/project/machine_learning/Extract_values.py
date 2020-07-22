from pymatgen import MPRester
import pandas as pd
import numpy as np
from pymatgen.core.periodic_table import Element
import json
import math
import os
import sys
from com.project.machine_learning.pymatgen_utils import holder, voronoi_tools
from com.project.examples.gbml import elasticity
import string
from pymatgen.analysis.elasticity import ElasticTensor

from numpy import array

# file to import
IMPORTED_FILE = "elasticElate_ALL_revisionArt_without_Zero.csv"
# file to export
FILTE_TO_EXPORT = 'Extract_Allvalues_descriptors.csv'
# Colonne qui contient la liste des descriptors : IL FAUT RESPECTER L'ORDRE
COLONNE = [
    # mean with power [-4..4] with stdev [0,1] : eGrp, eAtmM, eZ, eRad, eRow, eX, eBlPt, eMlPt :88
    'eGrpHn4A', 'eGrpHn3A', 'eGrpHn2A', 'eGrpHn1A', 'eGrpH0A', 'eGrpH1A', 'eGrpH2A', 'eGrpH3A', 'eGrpH4A', 'eGrpH0S',
    'eGrpH1S',
    'eAtMHn4A', 'eAtMHn3A', 'eAtMHn2A', 'eAtMHn1A', 'eAtMH0A', 'eAtMH1A', 'eAtMH2A', 'eAtMH3A', 'eAtMH4A', 'eAtMH0S',
    'eAtMH1S',
    'eZHn4A', 'eZHn3A', 'eZHn2A', 'eZHn1A', 'eZH0A', 'eZH1A', 'eZH2A', 'eZH3A', 'eZH4A', 'eZH0S',
    'eZH1S',
    'eRadHn4A', 'eRadHn3A', 'eRadHn2A', 'eRadHn1A', 'eRadH0A', 'eRadH1A', 'eRadH2A', 'eRadH3A', 'eRadH4A', 'eRadH0S',
    'eRadH1S',
    'eRowHn4A', 'eRowHn3A', 'eRowHn2A', 'eRowHn1A', 'eRowH0A', 'eRowH1A', 'eRowH2A', 'eRowH3A', 'eRowH4A', 'eRowH0S',
    'eRowH1S',
    'eXHn4A', 'eXHn3A', 'eXHn2A', 'eXHn1A', 'eXH0A', 'eXH1A', 'eXH2A', 'eXH3A', 'eXH4A', 'eXH0S', 'eXH1S',
    'eBlPtHn4A', 'eBlPtHn3A', 'eBlPtHn2A', 'eBlPtHn1A', 'eBlPtH0A', 'eBlPtH1A', 'eBlPtH2A', 'eBlPtH3A', 'eBlPtH4A',
    'eBlPtH0S', 'eBlPtH1S',
    'eMlPtHn4A', 'eMlPtHn3A', 'eMlPtHn2A', 'eMlPtHn1A', 'eMlPtH0A', 'eMlPtH1A', 'eMlPtH2A', 'eMlPtH3A', 'eMlPtH4A',
    'eMlPtH0S', 'eMlPtH1S',

    # autres proprietes:8
    'lvpa', 'cepa', 'cohesive_energy', 'average_electroneg', 'bandGap', 'rho', 'fepa',
    'eah',
    # coordination,bond length,bond angle:33
    'sCoorHn4A', 'sCoorHn3A', 'sCoorHn2A', 'sCoorHn1A', 'sCoorH0A', 'sCoorH1A', 'sCoorH2A', 'sCoorH3A', 'sCoorH4A',
    'sCoorH0S', 'sCoorH1S',
    'sBnLnHn4AH1A', 'sBnLnHn3AH1A', 'sBnLnHn2AH1A', 'sBnLnHn1AH1A', 'sBnLnH0AH1A', 'sBnLnH1AH1A', 'sBnLnH2AH1A',
    'sBnLnH3AH1A', 'sBnLnH4AH1A', 'sBnLnH0SH1A', 'sBnLnH1SH1A',
    'sBnAnHn4AH1A', 'sBnAnHn3AH1A', 'sBnAnHn2AH1A', 'sBnAnHn1AH1A', 'sBnAnH0AH1A', 'sBnAnH1AH1A', 'sBnAnH2AH1A',
    'sBnAnH3AH1A', 'sBnAnH4AH1A', 'sBnAnH0SH1A', 'sBnAnH1SH1A',

    # descriptors of neighbor differences :60
    'sRowADH0AH1A', 'sRowADH1AH1A', 'sRowADH2AH1A', 'sRowADH3AH1A', 'sRowADH4AH1A', 'sRowADH1SH1A', 'sRowSDH1AH1A',
    'sRowSDH2AH1A', 'sRowSDH4AH1A', 'sRowSDH1SH1A', 'sGrpADH0AH1A', 'sGrpADH1AH1A', 'sGrpADH2AH1A', 'sGrpADH3AH1A',
    'sGrpADH4AH1A', 'sGrpADH1SH1A', 'sGrpSDH1AH1A', 'sGrpSDH2AH1A', 'sGrpSDH4AH1A', 'sGrpSDH1SH1A',
    'sAtMADH0AH1A', 'sAtMADH1AH1A', 'sAtMADH2AH1A', 'sAtMADH3AH1A', 'sAtMADH4AH1A', 'sAtMADH1SH1A', 'sAtMSDH1AH1A',
    'sAtMSDH2AH1A', 'sAtMSDH4AH1A', 'sAtMSDH1SH1A',
    'sRadADH0AH1A', 'sRadADH1AH1A', 'sRadADH2AH1A', 'sRadADH3AH1A', 'sRadADH4AH1A', 'sRadADH1SH1A', 'sRadSDH1AH1A',
    'sRadSDH2AH1A', 'sRadSDH4AH1A', 'sRadSDH1SH1A',
    'sXADH0AH1A', 'sXADH1AH1A', 'sXADH2AH1A', 'sXADH3AH1A', 'sXADH4AH1A', 'sXADH1SH1A',
    'sXSDH1AH1A', 'sXSDH2AH1A', 'sXSDH4AH1A', 'sXSDH1SH1A',
    'sZADH0AH1A', 'sZADH1AH1A', 'sZADH2AH1A', 'sZADH3AH1A', 'sZADH4AH1A', 'sZADH1SH1A', 'sZSDH1AH1A',
    'sZSDH2AH1A', 'sZSDH4AH1A', 'sZSDH1SH1A'
]

api = MPRester("eDCEK5m9WVjmajp7e8af")
CAVEAT_AIAB = 'Unable to estimate cohesive energy for material.'
CAVEAT_F_BLOCK = 'Predictions are likely less reliable for materials containing F-block elements.'
CAVEAT_HUBBARD = 'Predictions may be less reliable for materials with non-GGA runs.'
DATAFILE_AIAB = 'data/element_aiab_energy.json'
VERY_SMALL = 1E-12

propsTableauCritere = ['material_id', 'pretty_formula', 'energy', 'energy_per_atom', 'density',
                       'formation_energy_per_atom', 'volume', 'is_hubbard', 'nsites', 'spacegroup', 'band_gap',
                       'e_above_hull', 'structure']
# dictionnaire contenant tous les éléments chimiques du tableau périodique
all_elements = dict()

# material Ids du tableau a remplir dans get_calculated_properties
tableau_material_ids = []


# Recuperer tous les elements du tableau périodique
def get_all_elements():
    for element in Element:
        try:
            all_elements[element.value] = \
                api._make_request("/element/%s/tasks/isolated_atom" % element.value, mp_decode=False)[0]
        except:
            pass
            print("\nProblème de recuperation de l'élément atomique  : " + str(element))


# Import du fichier csv initial
def import_file(fichier):
    return pd.read_csv(fichier, index_col=0)


# Function to handle 'no data' strings
def checkData(field):
    if field is None:
        sys.stderr.write('  Warning: A None field was set to zero in checkData.\n')
        return '0'
    if isinstance(field, str) and field[:7] == 'no data':
        sys.stderr.write('  Warning: A \'no data\' field was set to zero in checkData.\n')
        return '0'
    else:
        return field


def holder_mean_all_powers(values, powers, weights=None, weights_norm=None):
    holder_mean_arrays = []
    for power in powers:
        holder_mean_arrays.append(holder.mean(values, power, weights, weights_norm))
    return holder_mean_arrays


def stdevs_all_powers(values, powers, weights=None, weights_norm=None):
    stdevs_arrays = []
    for power in powers:
        stdevs_arrays.append(holder.stdev(values, power, weights, weights_norm))
    return stdevs_arrays


def build_ascontiguousarray(properties_table):
    ascontiguousarray_output = []

    for entries_by_property in properties_table:
        ascontiguousarray_output.append(entries_by_property)
    return ascontiguousarray_output


def build_check_list(properties_table, num_predictions, nbre_property):
    empty_list = [None] * nbre_property
    for entries_by_property in properties_table:
        if len(entries_by_property) != num_predictions:
            return empty_list


def get_calculated_properties(entries_var):
    matid_list = []
    aiab_problem_list = []
    cohesive_energy_problem_list = []
    properties_table_entries = []
    element_number = 0
    for entry in entries_var:
        properties_table_entry = []
        caveats_str = ''
        aiab_flag = False
        f_block_flag = False
        eElms = []
        eWts = []
        eAIAB = []
        eRow = []
        eGrp = []
        eAtmM = []
        eZ = []
        eRad = []
        eX = []
        eBlPt = []
        eMlPt = []
        average_electroneg = None
        cohesive_energy = None

        # Construct per-element lists for this material
        composition = entry.composition

        for element_key, amount in composition.get_el_amt_dict().items():
            element = Element(element_key)
            eElms.append(element)
            eWts.append(composition.get_atomic_fraction(element))
            aiab_energy = elasticity.get_element_aiab_energy(element_key)  # aiab = atom-in-a-box
            if aiab_energy is None:
                aiab_flag = True
                break
            eAIAB.append(aiab_energy)
            if element.block == 'f':
                f_block_flag = True

            iElms__eBlPt = float(
                ((checkData(element.data['Boiling point'])).replace(' K', '')).replace('(white P) ', ''))
            if iElms__eBlPt == 0:
                if element_key == 'Pa':
                    iElms__eBlPt = 4273
                else:
                    sys.stderr.write(
                        '  Warning: In {} element {} has boiling point of None\n'.format(str(entry.data["material_id"]),
                                                                                         element_key))
            eBlPt.append(float(iElms__eBlPt))

            iElms__eMlPt = float(
                ((checkData(element.data['Melting point'])).replace(' K', '')).replace('(white P) ', ''))
            if iElms__eMlPt == 0:
                sys.stderr.write(
                    '  Warning: In {} element {} has melting point of None\n'.format(str(entry.data["material_id"]),
                                                                                     element_key))
            eMlPt.append(float(iElms__eMlPt))
            eGrp.append(element.group)
            eAtmM.append(element.atomic_mass)
            eZ.append(element.number)
            eRad.append(element.atomic_radius)
            eRow.append(element.row)
            eX.append(element.X)

        # On error, add material to aiab_problem_list and continue with next material
        if aiab_flag:
            aiab_problem_list.append(str(entry.data["material_id"]))
            continue
        # Check caveats
        if bool(entry.parameters['is_hubbard']):
            if len(caveats_str) > 0:
                caveats_str += " "
            caveats_str += CAVEAT_HUBBARD
        if f_block_flag:
            if len(caveats_str) > 0:
                caveats_str += " "
            caveats_str += CAVEAT_F_BLOCK
        # Calculate intermediate weighted averages (WA) for this material
        ewa = np.average(eAIAB, weights=eWts)  # atom-in-a-box energy WA
        try:
            average_electroneg = entry.composition.average_electroneg
        except:
            print("\nErreur calcul average_electroneg de  : " + str(entry.data['material_id']))

        try:
            cohesive_energy = get_cohesive_energy(entry, True)
        except:
            cohesive_energy_problem_list.append(str(entry.data["material_id"]))

        # Construct Voronoi neighbor dictionaries
        structure = entry.data['structure']
        try:
            # TODO change get_voronoi_dicts_3 to get_voronoi_dicts_2 to calculte bond_angles
            # (voronoi_neighbor_sites, voronoi_neighbor_pairs) = voronoi_tools.get_voronoi_dicts(structure)
            (voronoi_neighbor_sites, voronoi_neighbor_pairs) = voronoi_tools.get_voronoi_dicts_3(structure)
            # (voronoi_neighbor_sites, voronoi_neighbor_pairs) = voronoi_tools.get_voronoi_dicts_3(structure)
        except RuntimeError:
            sys.stderr.write('  Error: Voronoi failed for {}\n'.format(str(entry.data["material_id"])))
            continue

        element_number = element_number + 1
        print('  Traitement de l\'élément : {} {}\n'.format(str(entry.data["material_id"]), element_number))
        # Append properties with powers for each material
        # see COLONNE
        for property in [eGrp, eAtmM, eZ, eRad, eRow, eX, eBlPt, eMlPt]:
            properties_table_entry.extend(
                holder_mean_all_powers(property, [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0],
                                       eWts))
            # Append stdevs with powers for each material
            properties_table_entry.extend(stdevs_all_powers(property, [0.0, 1.0], eWts))

        # lvpa
        properties_table_entry.append(math.log10(float(entry.data["volume"]) / float(entry.data["nsites"])))
        # cepa
        properties_table_entry.append(float(entry.energy_per_atom) - ewa)
        # caveats
        # properties_table_entry.append(caveats_str)
        # cohesive energy
        properties_table_entry.append(cohesive_energy)
        properties_table_entry.append(average_electroneg)
        properties_table_entry.append(entry.data["band_gap"])
        properties_table_entry.append(entry.data["density"])
        properties_table_entry.append(entry.data["formation_energy_per_atom"])
        properties_table_entry.append(entry.data["e_above_hull"])

        # Calc and print Voronoi based coordination descriptors:
        #   sCoorHn4A,sCoorHn3A,sCoorHn2A,sCoorHn1A,sCoorH0A,sCoorH1A,sCoorH2A,sCoorH3A,sCoorH4A,sCoorH0S,sCoorH1S
        site_means = voronoi_tools.get_coordinations(voronoi_neighbor_sites)
        properties_table_entry.extend(holder_mean_all_powers(site_means,
                                                             [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]))
        properties_table_entry.extend(stdevs_all_powers(site_means, [0.0, 1.0]))

        # Calc and print Voronoi based bond length descriptors:,
        # sBnLnHn4AH1A, sBnLnHn3AH1A, sBnLnHn2AH1A, sBnLnHn1AH1A, sBnLnH0AH1A, sBnLnH1AH1A, sBnLnH2AH1A, sBnLnH3AH1A, sBnLnH4AH1A, sBnLnH0SH1A, sBnLnH1SH1A
        (site_means, site_stdevs) = voronoi_tools.get_bond_lengths(voronoi_neighbor_sites)
        properties_table_entry.extend(
            holder_mean_all_powers(site_means, [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]))
        properties_table_entry.extend(
            stdevs_all_powers(site_means, [0.0, 1.0]))  # No negative Holder means due to zero valued centered means

        # Calc and print Voronoi based bond angle descriptors:
        #   sBnAnHn4AH1A, sBnAnHn3AH1A, sBnAnHn2AH1A, sBnAnHn1AH1A, sBnAnH0AH1A, sBnAnH1AH1A, sBnAnH2AH1A, sBnAnH3AH1A, sBnAnH4AH1A, sBnAnH0SH1A, sBnAnH1SH1A
        # TODO remove condition to calculte bond_angles
        if len(voronoi_neighbor_pairs) > 0:
            (site_means, site_stdevs) = voronoi_tools.get_bond_angles(voronoi_neighbor_pairs)
            #
            properties_table_entry.extend(
                holder_mean_all_powers(site_means, [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]))
            properties_table_entry.extend(
                stdevs_all_powers(site_means, [0.0, 1.0]))  # No negative Holder means due to zero valued centered means
        else:
            properties_table_entry.extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        # Calculate and print Voronoi based Holder descriptors of neighbor differences:
        # *** With vetting/coding, could add: eBlPt, eMlPt, AIAB ***
        # see COLONNE
        for property in ['row', 'group', 'atomic_mass', 'atomic_radius', 'X', 'Z']:
            (site_means, site_stdevs) = voronoi_tools.get_property_diffs(voronoi_neighbor_sites, property,
                                                                         abs_flag=True)
            properties_table_entry.extend(holder_mean_all_powers(site_means,
                                                                 [0.0, 1.0, 2.0, 3.0,
                                                                  4.0]))  # No negative Holder means due to zero valued diffs
            properties_table_entry.extend(stdevs_all_powers(site_means,
                                                            [
                                                                1.0]))  # No non-positive Holder stdevs due to zero valued diffs

            (site_means, site_stdevs) = voronoi_tools.get_property_diffs(voronoi_neighbor_sites, property,
                                                                         abs_flag=False)
            properties_table_entry.extend(holder_mean_all_powers(site_means,
                                                                 [1.0, 2.0,
                                                                  4.0]))  # No non-positive or cubic Holder means due to negative values
            properties_table_entry.extend(stdevs_all_powers(site_means,
                                                            [
                                                                1.0]))  # No non-positive or cubic Holder means due to negative values

        properties_table_entries.append(properties_table_entry)
        matid_list.append(str(entry.data["material_id"]))

    inverted_properties_table_entries = list(map(list, zip(*properties_table_entries)))
    # Check that at least one valid material was provided
    num_predictions = len(matid_list)
    if num_predictions > 0:
        # Construct descriptor arrays
        build_check_list(inverted_properties_table_entries, num_predictions, len(COLONNE))
        descriptors = np.ascontiguousarray(build_ascontiguousarray(inverted_properties_table_entries), dtype=float)

    print("\nErreur calcul get_cohesive_energy, nombre d'éléments " + str(
        len(cohesive_energy_problem_list)) + "  de  : " + str(cohesive_energy_problem_list))
    print("\nErreur calcul aiab, nombre d'éléments " + str(len(cohesive_energy_problem_list)) + " de  : " + str(
        cohesive_energy_problem_list))

    return pd.DataFrame(descriptors.transpose(), index=matid_list,
                        columns=COLONNE)


def export_additional_properties(data1: pd.DataFrame, data2: pd.DataFrame, fichier):
    # anciennes proprietes du csv avec les nouvelles proprietes calcultées
    mergedDf = data1.merge(data2, left_index=True, right_index=True)
    # les materials qui contiennent au moins une propriete nulle
    df_with_null = mergedDf[mergedDf.isna().any(axis=1)]
    print("\nNombre de propriétés  : " + str(len(list(mergedDf.columns.values))))
    print("\n liste des propriétés : \n")
    print(str(list(mergedDf.columns.values)) + "\n")
    print("\nNombre d'élément avant suppression des lignes contenant null  : " + str(len(mergedDf.index)))
    # mergedDf_filtered = mergedDf[(mergedDf['cohesive_energy'].notnull())]
    # Liste de materials sans ceux qui sont nulls
    mergedDf_filtered = mergedDf.dropna(how='any', axis=0)
    mergedDf_filtered.to_csv(fichier, index=True)
    print("\nNombre d'élément après suppression des lignes contenant null  : " + str(len(mergedDf_filtered.index)))
    print("\nNombre de lignes contenants des null  : " + str(len(df_with_null)) + "\n")
    if len(df_with_null) > 0:
        print(df_with_null)


def get_cohesive_energy(entry, per_atom):
    ebulk = entry.energy / entry.composition.get_integer_formula_and_factor()[1]
    comp_dict = entry.composition.reduced_composition.as_dict()
    isolated_atom_e_sum, n = 0, 0
    for atom in comp_dict.keys():
        e = all_elements.get(atom)
        isolated_atom_e_sum += e['output']["final_energy_per_atom"] * comp_dict[atom]
        n += comp_dict[atom]
    ecoh_per_formula = isolated_atom_e_sum - ebulk
    return ecoh_per_formula / n if per_atom else ecoh_per_formula


# 1- On recupere toutes les lignes du fichiers csv
data_from_cv = import_file(IMPORTED_FILE)

print("\nNombre de tous les éléments dans le fichier csv = {}\n".format(data_from_cv.shape[0]))

# 2- Lister tous les material Ids
Materials_Ids = tuple(list(data_from_cv.index.values))

# 3- Recuperation de tous les entries avec les propriétés correspondantes
entries = api.get_entries({'material_id': {'$in': Materials_Ids}}, property_data=propsTableauCritere)
# materials = api.query(criteria={'material_id': {'$in': Materials_Ids}}, properties=propsTableauCritere)

# 4- Recupération des propriétés des éléments du tableau periodique
get_all_elements()

# 5-Calcul des propriétés à completer au fichier csv
data_additional_prop = get_calculated_properties(entries)
#
# 6- Generation du nouveau fichier csv contenant toutes les proprietes
export_additional_properties(data_from_cv, data_additional_prop, FILTE_TO_EXPORT)
