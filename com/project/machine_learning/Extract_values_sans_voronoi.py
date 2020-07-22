from pymatgen import MPRester
import pandas as pd
import numpy as np
from pymatgen.core.periodic_table import Element
import json
import math
import os
import sys
import string
from pymatgen.analysis.elasticity import ElasticTensor

from numpy import array

COLONNE = ['group-4', 'group-3', 'group-2', 'group-1', 'group-0', 'group1', 'group2', 'group3', 'group4',
           'atomic_mass-4', 'atomic_mass-3', 'atomic_mass-2', 'atomic_mass-1', 'atomic_mass-0', 'atomic_mass1',
           'atomic_mass2', 'atomic_mass3', 'atomic_mass4', 'atomic_number-4', 'atomic_number-3', 'atomic_number-2',
           'atomic_number-1', 'atomic_number-0', 'atomic_number1', 'atomic_number2', 'atomic_number3', 'atomic_number4',
           'atomic_radius-4', 'atomic_radius-3', 'atomic_radius-2', 'atomic_radius-1', 'atomic_radius-0',
           'atomic_radius1', 'atomic_radius2', 'atomic_radius3', 'atomic_radius4', 'row-4', 'row-3', 'row-2', 'row-1',
           'row-0', 'row1', 'row2', 'row3', 'row4', 'x-4', 'x-3', 'x-2', 'x-1', 'x-0', 'x1', 'x2', 'x3', 'x4',
           'eBlPt-4', 'eBlPt-3', 'eBlPt-2', 'eBlPt-1', 'eBlPt-0', 'eBlPt1', 'eBlPt2', 'eBlPt3', 'eBlPt4', 'eMlPt-4',
           'eMlPt-3', 'eMlPt-2', 'eMlPt-1', 'eMlPt-0', 'eMlPt1', 'eMlPt2', 'eMlPt3', 'eMlPt4', 'lvpa', 'cepa',
           'cohesive_energy', 'average_electroneg', 'bandgap', 'density', 'formation_energy-peratom', 'e_above_hull']

api = MPRester("eDCEK5m9WVjmajp7e8af")
CAVEAT_AIAB = 'Unable to estimate cohesive energy for material.'
CAVEAT_F_BLOCK = 'Predictions are likely less reliable for materials containing F-block elements.'
CAVEAT_HUBBARD = 'Predictions may be less reliable for materials with non-GGA runs.'
DATAFILE_AIAB = 'data/element_aiab_energy.json'
VERY_SMALL = 1E-12

propsTableauCritere = ['material_id', 'pretty_formula', 'energy', 'energy_per_atom', 'density',
                       'formation_energy_per_atom', 'volume', 'is_hubbard', 'nsites', 'spacegroup', 'band_gap',
                       'e_above_hull']

all_elements = dict()


# Recuperer tous les elements du tableau périodique
def get_all_elements():
    for element in Element:
        try:
            all_elements[element.value] = \
                api._make_request("/element/%s/tasks/isolated_atom" % element.value, mp_decode=False)[0]
        except:
            pass
            print("\nProblème de recuperation de l'élément atomique  : " + str(element))


# Import du fichier csv initlial
def importer(fichier):
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


# material Ids du tableau a remplir dans get_calculated_properties
tableau_material_ids = []


def holder_mean_all_powers(values, powers, weights, weights_norm=None):
    result_array = []
    for power in powers:
        result_array.append(holder_mean(values, power, weights, weights_norm))

    return result_array


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

    for entry in entries_var:
        properties_table_entry = []
        caveats_str = ''
        aiab_flag = False
        f_block_flag = False
        eElms = []
        weight_list = []
        energy_list = []
        row_list = []
        group_list = []
        atomic_mass_list = []
        atomic_number_list = []
        atomic_radius_list = []
        x_list = []
        eBlPt = []
        eMlPt = []
        average_electroneg = None
        cohesive_energy = None

        # Construct per-element lists for this material
        composition = entry.composition

        for element_key, amount in composition.get_el_amt_dict().items():
            element = Element(element_key)
            eElms.append(element)
            weight_list.append(composition.get_atomic_fraction(element))
            aiab_energy = get_element_aiab_energy(element_key)  # aiab = atom-in-a-box
            if aiab_energy is None:
                aiab_flag = True
                break
            energy_list.append(aiab_energy)
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
            group_list.append(element.group)
            atomic_mass_list.append(element.atomic_mass)
            atomic_number_list.append(element.number)
            atomic_radius_list.append(element.atomic_radius)
            row_list.append(element.row)
            x_list.append(element.X)

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
        ewa = np.average(energy_list, weights=weight_list)  # atom-in-a-box energy WA
        try:
            average_electroneg = entry.composition.average_electroneg
        except:
            print("\nErreur calcul average_electroneg de  : " + str(entry.data['material_id']))

        try:
            cohesive_energy = get_cohesive_energy(entry, True)
        except:
            cohesive_energy_problem_list.append(str(entry.data["material_id"]))

        # print(str(entry.data["material_id"]))

        # Append descriptors for this material to descriptor lists
        for property in [group_list, atomic_mass_list, atomic_number_list, atomic_radius_list, row_list, x_list, eBlPt,
                         eMlPt]:
            result_with_powers = []
            result_with_powers = holder_mean_all_powers(property, [-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0],
                                                        weight_list)
            properties_table_entry.extend(result_with_powers)

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

        properties_table_entries.append(properties_table_entry)
        matid_list.append(str(entry.data["material_id"]))

    inverted_properties_table_entries = list(map(list, zip(*properties_table_entries)))
    # Check that at least one valid material was provided
    num_predictions = len(matid_list)
    if num_predictions > 0:
        # Construct descriptor arrays
        build_check_list(inverted_properties_table_entries, num_predictions, len(COLONNE))
        descriptors = np.ascontiguousarray(build_ascontiguousarray(inverted_properties_table_entries), dtype=float)

    print("\nErreur calcul get_cohesive_energy, nombre d'éléments "+str(len(cohesive_energy_problem_list))+"  de  : " + str(cohesive_energy_problem_list))
    print("\nErreur calcul aiab, nombre d'éléments "+str(len(cohesive_energy_problem_list))+" de  : " + str(cohesive_energy_problem_list))

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


def holder_mean(values, power, weights=None, weights_norm=None):
    """
    Calculate (possibly weighted) Holder (or power) mean
    :param: values: list or array of (typically elemental property) values
    :param: power: the power to which the Holder mean is to be calculated
    :param: weights: option list or array of weights
    :param: weights_norm: optional normalization method to be applied to weights
    :return: holder_mean
    """

    values = np.array(values, dtype=float)
    power = float(power)

    # Make sure weights match length and are normalized
    if weights is None:
        alpha = 1 / len(values)
    else:
        weights = np.array(weights, dtype=float)
        if len(values) != len(weights):
            # warn('Holder.mean returned zero when passed length mis-matched values and weights', UserWarning)
            sys.stderr.write('  Warning: Holder.mean returned zero when passed length mis-matched values and weights\n')
            return 0.0
        if weights_norm is not None:
            if weights_norm == "max" and max(weights) != 1.0:
                weights = weights / max(weights)
            elif weights_norm == "sum" and sum(weights) != 1.0:
                weights = weights / sum(weights)
            else:
                # warn('Holder.mean returned zero when passed unknown weights_norm method', UserWarning)
                sys.stderr.write('  Warning: Holder.mean returned zero when passed unknown weights_norm method\n')
                return 0.0
        alpha = 1 / sum(weights)

    if power == 0.0:  # geometric mean
        if any(abs(value) < VERY_SMALL for value in values):
            return 0.0
        elif any(value < 0 for value in values):
            # warn('Holder.mean returned zero when passed negative value with zero power', UserWarning)
            sys.stderr.write('  Warning: Holder.mean returned zero when passed negative value with zero power\n')
            sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
            return 0.0

        if weights is None:
            return math.exp(alpha * sum(np.log(values)))
        else:
            return math.exp(alpha * sum(weights * np.log(values)))

    elif power == 1.0:  # arithmetic mean
        return np.average(values, weights=weights)

    if any(value < 0 for value in values):
        if power % 2 != 0.0:
            # warn('Holder.mean returned zero when passed negative value with non-even power', UserWarning)
            sys.stderr.write('  Warning: Holder.mean returned zero when passed negative value with non-even power\n')
            sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
            return 0.0

    if weights is None:
        return pow(alpha * sum(np.power(values, power)), 1 / power)
    else:
        return pow(alpha * sum(weights * np.power(values, power)), 1 / power)


def stdev(values, power, weights=None, weights_norm=None):
    values = array(values, dtype=float)
    power = float(power)

    # Check for single value in values
    if len(values) is 1:
        return 0.0

    # Make sure weights match length and are normalized
    if weights is None:
        beta = 1 / (len(values) - 1)
    else:
        weights = array(weights, dtype=float)
        if len(values) != len(weights):
            # warn('Holder.stdev returned zero when passed length mis-matched values and weights', UserWarning)
            sys.stderr.write(
                '  Warning: Holder.stdev returned zero when passed length mis-matched values and weights\n')
            return 0.0
        if weights_norm is not None:
            if weights_norm == "max" and max(weights) != 1.0:
                weights = weights / max(weights)
            elif weights_norm == "sum" and sum(weights) != 1.0:
                weights = weights / sum(weights)
            else:
                # warn('Holder.stdev returned zero when passed unknown weights_norm method', UserWarning)
                sys.stderr.write('  Warning: Holder.stdev returned zero when passed unknown weights_norm method\n')
                return 0.0
        alpha = sum(weights)  # Note: Alpha is defined differently here than in mean function!
        beta = alpha / (alpha ** 2 - sum(np.power(weights, 2)))
    # TODO
    # holder_mean = mean(values, power, weights=weights)

    if power == 0.0:  # geometric stdev (unbiased estimate)
        if any(value <= 0 for value in values):
            # warn('Holder.stdev returned zero when passed non-positive value with zero power', UserWarning)
            sys.stderr.write('  Warning: Holder.stdev returned zero when passed non-positive value with zero power\n')
            sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
            return 0.0
        if abs(holder_mean) < VERY_SMALL:
            # warn('Holder.stdev returned zero when passed values with near-zero mean with zero power', UserWarning)
            sys.stderr.write(
                '  Warning: Holder.stdev returned zero when passed values with near-zero mean with zero power\n')
            sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
            return 0.0
        if weights is None:
            return math.exp(math.sqrt(beta * sum(np.power(np.log(values / holder_mean), 2))))
        else:
            return math.exp(math.sqrt(beta * sum(weights * np.power(np.log(values / holder_mean), 2))))

    holder_mean_centered_values = values - holder_mean

    if power == 1.0:  # arithmetic stdev (unbiased estimate)
        if weights is None:
            return math.sqrt(beta * sum(np.power(holder_mean_centered_values, 2)))
        else:
            return math.sqrt(beta * sum(weights * np.power(holder_mean_centered_values, 2)))

    # else:
    #   #warn('Holder.stdev returned zero when passed power other than 0 or 1', UserWarning)
    #   sys.stderr.write('  Warning: Holder.stdev returned zero when passed power other than 0 or 1\n')
    #   sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
    #   return 0.0

    # if sum(np.abs(holder_mean_centered_values)) < VERY_SMALL:
    #   return 0.0

    if power < 0.0:
        if any(abs(centered_value) < VERY_SMALL for centered_value in holder_mean_centered_values):
            # warn('Holder.stdev returned zero when passed near-zero value with negative power', UserWarning)
            sys.stderr.write('  Warning: Holder.stdev returned zero when passed near-zero value with negative power\n')
            sys.stderr.write('    values = {:s}  power = {:f}\n'.format(holder_mean_centered_values, power))
            return 0.0

    if any(centered_value < 0 for centered_value in holder_mean_centered_values):
        if power % 1 != 0.0:
            # warn('Holder.stdev returned zero when passed non-positive value with non-integer power', UserWarning)
            sys.stderr.write(
                '  Warning: Holder.stdev returned zero when passed negative value with non-integer power\n')
            sys.stderr.write('    values = {:s}  power = {:f}\n'.format(holder_mean_centered_values, power))
            return 0.0

    if weights is None:
        # sys.stderr.write('    values = {:s}  power = {:f}\n'.format(holder_mean_centered_values, power))
        return pow(beta * sum(np.power(holder_mean_centered_values, 2 * power)), 1 / 2 / power)
    else:
        # sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(holder_mean_centered_values, weights, power))
        return pow(beta * sum(weights * np.power(holder_mean_centered_values, 2 * power)), 1 / 2 / power)


def get_element_aiab_energy(element):
    """
    Lookup atom-in-a-box (aiab) energy for specified element.
    Used to estimate cohesive energy of a compound from the compound's VASP energy.
    The elemental atom-in-a-box energies were provided by Wei Chen and Maarten de Jong.
    Returns atom-in-a-box energy for specified element.
    :param element:
    :return: aiab_energy
    """

    element_aiab_energy_dict = None

    try:
        with open(os.path.join(os.path.dirname(__file__), DATAFILE_AIAB), 'r') as json_file:
            element_aiab_energy_dict = json.load(json_file)

    finally:
        if element_aiab_energy_dict is None:
            return None

        object = element_aiab_energy_dict.get(element)
        if object is not None:
            return object[0]

        return None


# 1- On recupere toutes les lignes du fichiers csv
data_from_cv = importer("elasticElate_ALL_revisionArt_without_Zero.csv")

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
export_additional_properties(data_from_cv, data_additional_prop, 'Extract_Allvalues_descriptors_sans_voronoi.csv')
