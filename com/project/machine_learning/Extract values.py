from pymatgen import MPRester
import pandas as pd
import numpy as np
from pymatgen.core.periodic_table import Element
import json
import math
import os
import sys
from numpy import array
import string
from pymatgen.analysis.elasticity import ElasticTensor


api = MPRester("eDCEK5m9WVjmajp7e8af")
CAVEAT_AIAB = 'Unable to estimate cohesive energy for material.'
CAVEAT_F_BLOCK = 'Predictions are likely less reliable for materials containing F-block elements.'
CAVEAT_HUBBARD = 'Predictions may be less reliable for materials with non-GGA runs.'
DATAFILE_AIAB = 'data/element_aiab_energy.json'
VERY_SMALL = 1E-12

propsTableauCritere = ['material_id', 'pretty_formula', 'energy', 'energy_per_atom', 'density',
                       'formation_energy_per_atom', 'volume', 'is_hubbard', 'nsites', 'spacegroup', 'band_gap', 'e_above_hull']

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
def checkData( field ):
  if field is None:
    sys.stderr.write('  Warning: A None field was set to zero in checkData.\n')
    return '0'
  if isinstance(field, str) and field[:7] == 'no data':
    sys.stderr.write('  Warning: A \'no data\' field was set to zero in checkData.\n')
    return '0'
  else:
    return field

# material Ids du tableau a remlir dans get_calculated_properties
tableau_material_ids = []


def get_calculated_properties(entries_var):
    lvpa_list = []
    cepa_list = []
    group1_list = []
    group2_list = []
    group3_list = []
    group4_list = []
    group0_list = []
    groupNEG1_list = []
    groupNEG2_list = []
    groupNEG3_list = []
    groupNEG4_list = []
    atomicmass1_list = []
    atomicmass2_list = []
    atomicmass3_list = []
    atomicmass4_list = []
    atomicmass0_list = []
    atomicmassNEG1_list = []
    atomicmassNEG2_list = []
    atomicmassNEG3_list = []
    atomicmassNEG4_list = []
    atomicnumber1_list = []
    atomicnumber2_list = []
    atomicnumber3_list = []
    atomicnumber4_list = []
    atomicnumber0_list = []
    atomicnumberNEG1_list = []
    atomicnumberNEG2_list = []
    atomicnumberNEG3_list = []
    atomicnumberNEG4_list = []
    atomicradius1_list = []
    atomicradius2_list = []
    atomicradius3_list = []
    atomicradius4_list = []
    atomicradius0_list = []
    atomicradiusNEG1_list = []
    atomicradiusNEG2_list = []
    atomicradiusNEG3_list = []
    atomicradiusNEG4_list = []
    eBlPt1_list = []
    eBlPt2_list = []
    eBlPt3_list = []
    eBlPt4_list = []
    eBlPt0_list = []
    eBlPtNEG1_list = []
    eBlPtNEG2_list = []
    eBlPtNEG3_list = []
    eBlPtNEG4_list = []
    eMlPt1_list = []
    eMlPt2_list = []
    eMlPt3_list = []
    eMlPt4_list = []
    eMlPt0_list = []
    eMlPtNEG1_list = []
    eMlPtNEG2_list = []
    eMlPtNEG3_list = []
    eMlPtNEG4_list = []
    rowH1A_list = []
    rowH2A_list = []
    rowH3A_list = []
    rowH4A_list = []
    rowH0A_list = []
    rowHn1A_list = []
    rowHn2A_list = []
    rowHn3A_list = []
    rowHn4A_list = []
    xH4A_list = []
    xH3A_list = []
    xH2A_list = []
    xH1A_list = []
    xH0A_list = []
    xHn4A_list = []
    xHn3A_list = []
    xHn2A_list = []
    xHn1A_list = []
    matid_list = []
    caveats_list = []
    aiab_problem_list = []
    cohesive_energy_list = []
    cohesive_energy_problem_list = []
    average_electroneg_list = []
    bandgap_list = []
    density_list = []
    formation_energie_peratom_list = []
    energie_above_hull_list = []
    Atomicmass_0D = []
    Atomicmass_1D = []
    Atomicnumber_0D = []
    Atomicnumber_1D = []
    Atomicradius_0D = []
    Atomicradius_1D = []
    BoilingTem_0D = []
    BoilingTem_1D = []
    Electro_0D =[]
    Electro_1D =[]
    Groupnumber_0D =[]
    Groupnumber_1D = []
    MelingTem_0D = []
    MelingTem_1D = []
    rownumber_0D = []
    rownumber_1D = []
    for entry in entries_var:
        caveats_str = ''
        aiab_flag = False
        f_block_flag = False
        eElms = []
        eBlPt = []
        eMlPt = []
        weight_list = []
        energy_list = []
        row_list = []
        group_list = []
        atomicmass_list = []
        atomicnumber_list = []
        atomicradius_list = []
        x_list = []
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

            iElms__eBlPt = float(((checkData(element.data['Boiling point'])).replace(' K', '')).replace('(white P) ', ''))
            if iElms__eBlPt == 0:
                if element_key == 'Pa':
                    iElms__eBlPt = 4273
                else:
                    sys.stderr.write('  Warning: In {} element has boiling point of None\n'.format(element))
            eBlPt.append(float(iElms__eBlPt))

            iElms__eMlPt = float(((checkData(element.data['Melting point'])).replace(' K', '')).replace('(white P) ', ''))
            if iElms__eMlPt == 0:
                sys.stderr.write('  Warning: In {} element has melting point of None\n'.format(element))
            eMlPt.append(float(iElms__eMlPt))

            row_list.append(element.row)
            group_list.append(element.group)
            atomicmass_list.append(element.atomic_mass)
            atomicnumber_list.append(element.number)
            atomicradius_list.append(element.atomic_radius)
            x_list.append(element.X)

        # On error, add material to aiab_problem_list and continue with next material
        if aiab_flag:
            aiab_problem_list.append(str(entry.data["material_id"]))
            continue
        # Check caveats
        if bool(entry.parameters['is_hubbard']):
            if len(caveats_str) > 0: caveats_str += " "
            caveats_str += CAVEAT_HUBBARD
        if f_block_flag:
            if len(caveats_str) > 0: caveats_str += " "
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
        lvpa_list.append(math.log10(float(entry.data["volume"]) / float(entry.data["nsites"])))
        cepa_list.append(float(entry.energy_per_atom) - ewa)
        group1_list.append(holder_mean(group_list, 1.0, weights=weight_list))
        group2_list.append(holder_mean(group_list, 2.0, weights=weight_list))
        group3_list.append(holder_mean(group_list, 3.0, weights=weight_list))
        group4_list.append(holder_mean(group_list, 4.0, weights=weight_list))
        group0_list.append(holder_mean(group_list, 0.0, weights=weight_list))
        groupNEG1_list.append(holder_mean(group_list, -1.0, weights=weight_list))
        groupNEG2_list.append(holder_mean(group_list, -2.0, weights=weight_list))
        groupNEG3_list.append(holder_mean(group_list, -3.0, weights=weight_list))
        groupNEG4_list.append(holder_mean(group_list, -4.0, weights=weight_list))
        atomicmass1_list.append(holder_mean(atomicmass_list, 1.0, weights=weight_list))
        atomicmass2_list.append(holder_mean(atomicmass_list, 2.0, weights=weight_list))
        atomicmass3_list.append(holder_mean(atomicmass_list, 3.0, weights=weight_list))
        atomicmass4_list.append(holder_mean(atomicmass_list, 4.0, weights=weight_list))
        atomicmass0_list.append(holder_mean(atomicmass_list, 0.0, weights=weight_list))
        atomicmassNEG1_list.append(holder_mean(atomicmass_list, -1.0, weights=weight_list))
        atomicmassNEG2_list.append(holder_mean(atomicmass_list, -2.0, weights=weight_list))
        atomicmassNEG3_list.append(holder_mean(atomicmass_list, -3.0, weights=weight_list))
        atomicmassNEG4_list.append(holder_mean(atomicmass_list, -4.0, weights=weight_list))
        atomicnumber1_list.append(holder_mean(atomicnumber_list, 1.0, weights=weight_list))
        atomicnumber2_list.append(holder_mean(atomicnumber_list, 2.0, weights=weight_list))
        atomicnumber3_list.append(holder_mean(atomicnumber_list, 3.0, weights=weight_list))
        atomicnumber4_list.append(holder_mean(atomicnumber_list, 4.0, weights=weight_list))
        atomicnumber0_list.append(holder_mean(atomicnumber_list, 0.0, weights=weight_list))
        atomicnumberNEG1_list.append(holder_mean(atomicnumber_list, -1.0, weights=weight_list))
        atomicnumberNEG2_list.append(holder_mean(atomicnumber_list, -2.0, weights=weight_list))
        atomicnumberNEG3_list.append(holder_mean(atomicnumber_list, -3.0, weights=weight_list))
        atomicnumberNEG4_list.append(holder_mean(atomicnumber_list, -4.0, weights=weight_list))
        atomicradius1_list.append(holder_mean(atomicradius_list, 1.0, weights=weight_list))
        atomicradius2_list.append(holder_mean(atomicradius_list, 2.0, weights=weight_list))
        atomicradius3_list.append(holder_mean(atomicradius_list, 3.0, weights=weight_list))
        atomicradius4_list.append(holder_mean(atomicradius_list, 4.0, weights=weight_list))
        atomicradius0_list.append(holder_mean(atomicradius_list, 0.0, weights=weight_list))
        atomicradiusNEG1_list.append(holder_mean(atomicradius_list, -1.0, weights=weight_list))
        atomicradiusNEG2_list.append(holder_mean(atomicradius_list, -2.0, weights=weight_list))
        atomicradiusNEG3_list.append(holder_mean(atomicradius_list, -3.0, weights=weight_list))
        atomicradiusNEG4_list.append(holder_mean(atomicradius_list, -4.0, weights=weight_list))
        eBlPt1_list.append(holder_mean(eBlPt, 1.0, weights=weight_list))
        eBlPt2_list.append(holder_mean(eBlPt, 2.0, weights=weight_list))
        eBlPt3_list.append(holder_mean(eBlPt, 3.0, weights=weight_list))
        eBlPt4_list.append(holder_mean(eBlPt, 4.0, weights=weight_list))
        eBlPt0_list.append(holder_mean(eBlPt, 0.0, weights=weight_list))
        eBlPtNEG1_list.append(holder_mean(eBlPt, -1.0, weights=weight_list))
        eBlPtNEG2_list.append(holder_mean(eBlPt, -2.0, weights=weight_list))
        eBlPtNEG3_list.append(holder_mean(eBlPt, -3.0, weights=weight_list))
        eBlPtNEG4_list.append(holder_mean(eBlPt, -4.0, weights=weight_list))
        eMlPt1_list.append(holder_mean(eMlPt, 1.0, weights=weight_list))
        eMlPt2_list.append(holder_mean(eMlPt, 2.0, weights=weight_list))
        eMlPt3_list.append(holder_mean(eMlPt, 3.0, weights=weight_list))
        eMlPt4_list.append(holder_mean(eMlPt, 4.0, weights=weight_list))
        eMlPt0_list.append(holder_mean(eMlPt, 0.0, weights=weight_list))
        eMlPtNEG1_list.append(holder_mean(eMlPt, -1.0, weights=weight_list))
        eMlPtNEG2_list.append(holder_mean(eMlPt, -2.0, weights=weight_list))
        eMlPtNEG3_list.append(holder_mean(eMlPt, -3.0, weights=weight_list))
        eMlPtNEG4_list.append(holder_mean(eMlPt, -4.0, weights=weight_list))
        rowH1A_list.append(holder_mean(row_list, 1.0, weights=weight_list))
        rowH2A_list.append(holder_mean(row_list, 2.0, weights=weight_list))
        rowH3A_list.append(holder_mean(row_list, 3.0, weights=weight_list))
        rowH4A_list.append(holder_mean(row_list, 4.0, weights=weight_list))
        rowH0A_list.append(holder_mean(row_list, 0.0, weights=weight_list))
        rowHn1A_list.append(holder_mean(row_list, -1.0, weights=weight_list))
        rowHn2A_list.append(holder_mean(row_list, -2.0, weights=weight_list))
        rowHn3A_list.append(holder_mean(row_list, -3.0, weights=weight_list))
        rowHn4A_list.append(holder_mean(row_list, -4.0, weights=weight_list))
        xH4A_list.append(holder_mean(x_list, 4.0, weights=weight_list))
        xH3A_list.append(holder_mean(x_list, 3.0, weights=weight_list))
        xH2A_list.append(holder_mean(x_list, 2.0, weights=weight_list))
        xH1A_list.append(holder_mean(x_list, 1.0, weights=weight_list))
        xH0A_list.append(holder_mean(x_list, 0.0, weights=weight_list))
        xHn4A_list.append(holder_mean(x_list, -4.0, weights=weight_list))
        xHn3A_list.append(holder_mean(x_list, -3.0, weights=weight_list))
        xHn2A_list.append(holder_mean(x_list, -2.0, weights=weight_list))
        xHn1A_list.append(holder_mean(x_list, -1.0, weights=weight_list))
        Atomicmass_0D.append(stdev(atomicmass_list, 0.0, weights=weight_list))
        Atomicmass_1D.append(stdev(atomicmass_list, 1.0, weights=weight_list))
        Atomicnumber_0D.append(stdev(atomicnumber_list, 0.0, weights=weight_list))
        Atomicnumber_1D.append(stdev(atomicnumber_list, 1.0, weights=weight_list))
        Atomicradius_0D.append(stdev(atomicradius_list, 0.0, weights=weight_list))
        Atomicradius_1D.append(stdev(atomicradius_list, 1.0, weights=weight_list))
        BoilingTem_0D.append(stdev(eBlPt, 0.0, weights=weight_list))
        BoilingTem_1D.append(stdev(eBlPt, 1.0, weights=weight_list))
        Electro_0D.append(stdev(x_list, 0.0, weights=weight_list))
        Electro_1D.append(stdev(x_list, 1.0, weights=weight_list))
        Groupnumber_0D.append(stdev(group_list, 0.0, weights=weight_list))
        Groupnumber_1D.append(stdev(group_list, 1.0, weights=weight_list))
        MelingTem_0D.append(stdev(eMlPt, 0.0, weights=weight_list))
        MelingTem_1D.append(stdev(eMlPt, 1.0, weights=weight_list))
        rownumber_0D.append(stdev(row_list, 0.0, weights=weight_list))
        rownumber_1D.append(stdev(row_list, 1.0, weights=weight_list))
        caveats_list.append(caveats_str)
        cohesive_energy_list.append(cohesive_energy)
        average_electroneg_list.append(average_electroneg)
        bandgap_list.append(entry.data["band_gap"])
        density_list.append(entry.data["density"])
        matid_list.append(str(entry.data["material_id"]))
        formation_energie_peratom_list.append(entry.data["formation_energy_per_atom"])
        energie_above_hull_list.append(entry.data["e_above_hull"])

    # Check that at least one valid material was provided
    num_predictions = len(matid_list)
    if num_predictions > 0:
        # Construct descriptor arrays
        if (len(lvpa_list) != num_predictions or len(cepa_list) != num_predictions or
                len(rowH1A_list) != num_predictions or len(rowHn3A_list) != num_predictions or
                len(xH4A_list) != num_predictions or len(xHn4A_list) != num_predictions or
                len(cohesive_energy_list) != num_predictions or len(average_electroneg_list) != num_predictions):
            return (None, None, None, None, None, None)
        descriptors = np.ascontiguousarray(
            [lvpa_list, cepa_list, group1_list, group2_list, group3_list, group4_list, group0_list, groupNEG1_list, groupNEG2_list,
             groupNEG3_list, groupNEG4_list, atomicmass1_list, atomicmass2_list, atomicmass3_list, atomicmass4_list, atomicmass0_list,
             atomicmassNEG1_list, atomicmassNEG2_list, atomicmassNEG3_list, atomicmassNEG4_list, atomicnumber1_list, atomicnumber2_list, atomicnumber3_list, atomicnumber4_list, atomicnumber0_list,
             atomicnumberNEG1_list, atomicnumberNEG2_list, atomicnumberNEG3_list, atomicnumberNEG4_list, atomicradius1_list, atomicradius2_list, atomicradius3_list,
             atomicradius4_list, atomicradius0_list, atomicradiusNEG1_list, atomicradiusNEG2_list, atomicradiusNEG3_list, atomicradiusNEG4_list, eBlPt1_list, eBlPt2_list,
             eBlPt3_list, eBlPt4_list, eBlPt0_list, eBlPtNEG1_list, eBlPtNEG2_list, eBlPtNEG3_list, eBlPtNEG4_list, eMlPt1_list, eMlPt2_list, eMlPt3_list, eMlPt4_list,
             eMlPt0_list, eMlPtNEG1_list, eMlPtNEG2_list, eMlPtNEG3_list, eMlPtNEG4_list,
             rowH1A_list, rowH2A_list, rowH3A_list, rowH4A_list, rowH0A_list, rowHn1A_list, rowHn2A_list, rowHn3A_list,
             rowHn4A_list, xH4A_list, xH3A_list, xH2A_list, xH1A_list, xH0A_list, xHn4A_list, xH3A_list, xH2A_list, xH1A_list, Atomicmass_0D, Atomicmass_1D,
             Atomicnumber_0D, Atomicnumber_1D, Atomicradius_0D, Atomicradius_1D, BoilingTem_0D, BoilingTem_1D, Electro_0D, Electro_1D, Groupnumber_0D, Groupnumber_1D,
             rownumber_0D, rownumber_1D, MelingTem_0D, MelingTem_1D,
             cohesive_energy_list, average_electroneg_list, bandgap_list, density_list, formation_energie_peratom_list, energie_above_hull_list] ,
            dtype=float)

    print("\nErreur calcul get_cohesive_energy de  : " + str(cohesive_energy_problem_list))
    print("\nErreur calcul aiab de  : " + str(cohesive_energy_problem_list))

    return pd.DataFrame(descriptors.transpose(), index=matid_list,
                        columns=['lvpa', 'cepa', 'group1','group2', 'group3', 'group4', 'group0', 'group-1', 'group-2', 'group-3', 'group-4',
                                 'atomic_mass1', 'atomic_mass2', 'atomic_mass3', 'atomic_mass4', 'atomic_mass0', 'atomic_mass-1', 'atomic_mass-2',
                                 'atomic_mass-3', 'atomic_mass-4',  'atomic_number1', 'atomic_number2', 'atomic_number3', 'atomic_number4', 'atomic_number0', 'atomic_number-1',
                                 'atomic_number-2','atomic_number-3', 'atomic_number-4','atomicRadius1', 'atomicRadius2', 'atomicRadius3', 'atomicRadius4', 'atomicRadius0', 'atomicRadius-1',
                                 'atomicRadius-2', 'atomicRadius-3', 'atomicRadius-4', 'eBlPt1', 'eBlPt2', 'eBlPt3', 'eBlPt4', 'eBlPt0', 'eBlPt-1', 'eBlPt-2',
                                 'eBlPt-3', 'eBlPt-4', 'eMlPt1', 'eMlPt2', 'eMlPt3', 'eMlPt4', 'eMlPt0', 'eMlPt-1', 'eMlPt-2',
                                 'eMlPt-3', 'eMlPt-4', 'rowH1A', 'rowH2A', 'rowH3A', 'rowH4A', 'rowH0A',
                                 'rowHn1A', 'rowHn2A', 'rowHn3A', 'rowHn4A', 'xH4A', 'xH3A', 'xH2A', 'xH1A', 'xH0A','xHn4A',
                                 'xHn3A', 'xHn2A', 'xHn1A', 'Atomicmass_0D', 'Atomicmass_1D', 'Atomicnumber_0D', 'Atomicnumber_1D', 'Atomicradius_0D',
                                 'Atomicradius_1D','BoilingTem_0D', 'BoilingTem_1D', 'Electro_0D', 'Electro_1D', 'Groupnumber_0D', 'Groupnumber_1D',
                                 'rownumber_0D', 'rownumber_1D', 'MelingTem_0D', 'MelingTem_1D', 'cohesive_energy',
                                 'average_electroneg', 'bandgap', 'density', 'formation_energy-peratom', 'e_above_hull'])

def export_additional_properties(data1, data2, fichier):
    mergedDf = data1.merge(data2, left_index=True, right_index=True)
    mergedDf_filtered = mergedDf[(mergedDf['cohesive_energy'].notnull())]
    mergedDf_filtered.to_csv(fichier, index=True)


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
            sys.stderr.write('  Warning: Holder.stdev returned zero when passed length mis-matched values and weights\n')
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

    holdermean = holder_mean(values, power, weights=weights)

    if power == 0.0:  # geometric stdev (unbiased estimate)
        if any(value <= 0 for value in values):
            # warn('Holder.stdev returned zero when passed non-positive value with zero power', UserWarning)
            sys.stderr.write('  Warning: Holder.stdev returned zero when passed non-positive value with zero power\n')
            sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
            return 0.0
        if abs(holdermean) < VERY_SMALL:
            # warn('Holder.stdev returned zero when passed values with near-zero mean with zero power', UserWarning)
            sys.stderr.write('  Warning: Holder.stdev returned zero when passed values with near-zero mean with zero power\n')
            sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
            return 0.0
        if weights is None:
            return math.exp(math.sqrt(beta * sum(np.power(np.log(values / holdermean), 2))))
        else:
            return math.exp(math.sqrt(beta * sum(weights * np.power(np.log(values / holdermean), 2))))

    holder_mean_centered_values = values - holdermean

    if power == 1.0:  # arithmetic stdev (unbiased estimate)
        if weights is None:
            return math.sqrt(beta * sum(np.power(holder_mean_centered_values, 2)))
        else:
            return math.sqrt(beta * sum(weights * np.power(holder_mean_centered_values, 2)))

    else:
      #warn('Holder.stdev returned zero when passed power other than 0 or 1', UserWarning)
      sys.stderr.write('  Warning: Holder.stdev returned zero when passed power other than 0 or 1\n')
      sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
      return 0.0

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
            sys.stderr.write('  Warning: Holder.stdev returned zero when passed negative value with non-integer power\n')
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
export_additional_properties(data_from_cv, data_additional_prop, 'Extract_Allvalues_descriptors.csv')
