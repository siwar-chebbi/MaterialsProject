from __future__ import division, unicode_literals

import math
import numpy as np
from numpy import array
from pyhull.voronoi import VoronoiTess
from pymatgen.analysis.local_env import NearNeighbors
from pymatgen.analysis.structure_analyzer import solid_angle
from pymatgen.analysis.structure_analyzer import VoronoiNN

from com.project.machine_learning.pymatgen_utils.structure_analyser_bis import VoronoiCoordFinder
from pymatgen.core.periodic_table import Element
from pymatgen.util.coord import get_angle

# Local variables
VORONOI_CUTOFF_RADIUS = 15.0  # default was 10.0
VORONOI_CUTOFF_WEIGHT_FRACTION = 0.03  # 0.04  # 1/25


# Warning: This script requires an edit to get_voronoi_polyhedra() in VoronoiCoordFinder in pymatgen.analysis.structure_analyzer
# Line 70 is:   center.coords, self.default_cutoff)
# Line 70 was:  center.coords, VoronoiCoordFinder.default_cutoff)

def get_voronoi_dicts(structure, cutoff_radius=VORONOI_CUTOFF_RADIUS,
                      cutoff_weight_fraction=VORONOI_CUTOFF_WEIGHT_FRACTION):
    voronoi_finder = VoronoiCoordFinder(structure)

    voronoi_finder.cutoff = cutoff_radius

    voronoi_finder_2 = VoronoiNN()
    voronoi_finder_2.cutoff = cutoff_radius

    all_neighbor_sites = {}
    all_neighbor_pairs = {}

    for center_site_index in range(0, len(structure.sites)):
        neighbor_sites_dict = {}
        neighbor_pairs_list = []
        # print "center_site_index: ", center_site_index

        # Construct initial voronoi polyhedra using cutoff_radius
        neighborhood_sites = [structure[center_site_index]]  # begin list with center site

        polyhedra_sites = voronoi_finder.get_voronoi_polyhedra(center_site_index)

        polyhedra_sites_2 = voronoi_finder_2.get_voronoi_polyhedra(structure, center_site_index)

        cutoff_weight = cutoff_weight_fraction * sum(polyhedra_sites.values())
        for neighbor_site, neighbor_weight in polyhedra_sites.items():
            if neighbor_weight > cutoff_weight:
                neighborhood_sites.append(neighbor_site)

        # Re-tesselate only with sites that also meet cutoff_weight_fraction criteria
        voronoi_input_coords = [site.coords for site in neighborhood_sites]
        try:
            voronoi_tess = VoronoiTess(voronoi_input_coords)
        except IndexError:
            raise RuntimeError('VoronoiTess: details unknown')

        local_solid_angles = []
        center_coords = structure[center_site_index].coords
        n_neighborhood_sites = len(neighborhood_sites)
        n_voronoi_vertices = len(voronoi_tess.vertices)
        neighbor_lookup = np.zeros(shape=(n_voronoi_vertices, n_voronoi_vertices, 2), dtype=np.int8)

        # Loop over neighbor sites (which surrond center_site) to:
        # - calculate solid angles
        # - construct neighbor_lookup array (which maps each voronoi edge to two adjacent neighbors)
        for neighbor_index in range(1, n_neighborhood_sites):
            facet_vertex_list = voronoi_tess.ridges[(0, neighbor_index)]  # 0 = center site
            if 0 in facet_vertex_list:
                raise RuntimeError("Pathological structure: voronoi facet includes vertex at infinity")

            # Calculate solid angle of facet between center site and this neighbor site
            facet_coords = [voronoi_tess.vertices[i] for i in facet_vertex_list]
            local_solid_angles.append(solid_angle(center_coords, facet_coords))

            # Add this voronoi ridge to neighbor_lookup array, to identify adjacent voronoi facets
            n_vertices = len(facet_vertex_list)
            for vertex1_index in range(0, n_vertices):
                vertex2_index = vertex1_index + 1
                if vertex2_index == n_vertices:
                    vertex2_index = 0  # wrap vertex2_index

                # Properly order vertex indicies to only use upper triangle of neighbor_lookup array
                (low_vert_index, high_vert_index) = (facet_vertex_list[vertex1_index], facet_vertex_list[vertex2_index])
                if low_vert_index > high_vert_index:
                    (low_vert_index, high_vert_index) = (high_vert_index, low_vert_index)

                # Store adjacent voronoi neighbors on different levels of neighbor_lookup, at location defined
                # by common voronoi edge running between low_vert_index (row) and high_vert_index (column).
                if neighbor_lookup[low_vert_index, high_vert_index, 0] == 0:
                    neighbor_lookup[low_vert_index, high_vert_index, 0] = neighbor_index  # first neighbor
                else:
                    neighbor_lookup[low_vert_index, high_vert_index, 1] = neighbor_index  # second neighbor

        # Loop over neighbor sites again to construct neighbor_sites_dict with solid angle weights
        max_local_solid_angle = max(local_solid_angles)
        for neighbor_index in range(1, n_neighborhood_sites):
            # Note: neighborhood_sites (which starts with center) is one longer than local_solid_angles
            neighbor_sites_dict[neighborhood_sites[neighbor_index]] = \
                local_solid_angles[neighbor_index - 1] / max_local_solid_angle

        # Loop through upper triangle of neighbor_lookup and build a list of tuples of neighbor pair sites
        for low_vert_index in range(1, n_voronoi_vertices):
            for high_vert_index in range(low_vert_index + 1, n_voronoi_vertices):
                neighbor1_index = neighbor_lookup[low_vert_index, high_vert_index, 0]
                neighbor2_index = neighbor_lookup[low_vert_index, high_vert_index, 1]
                if neighbor1_index > 0 and neighbor2_index > 0:  # Both are non-zero
                    neighbor_pairs_list.append(
                        (neighborhood_sites[neighbor1_index], neighborhood_sites[neighbor2_index]))
                elif neighbor1_index != 0 or neighbor2_index != 0:  # Only one is non-zero
                    raise RuntimeError("Unexpected neighbor_lookup matrix: non-paired edge")

        # Append per site values
        all_neighbor_sites[structure.sites[center_site_index]] = neighbor_sites_dict
        all_neighbor_pairs[structure.sites[center_site_index]] = neighbor_pairs_list

    return (all_neighbor_sites, all_neighbor_pairs)


def get_voronoi_dicts_2(structure, cutoff_radius=VORONOI_CUTOFF_RADIUS,
                        cutoff_weight_fraction=VORONOI_CUTOFF_WEIGHT_FRACTION):
    voronoi_finder_2 = VoronoiNN()
    voronoi_finder_2.cutoff = cutoff_radius

    all_neighbor_sites = {}
    all_neighbor_pairs = {}

    for center_site_index in range(0, len(structure.sites)):
        neighbor_sites_dict = {}
        neighbor_pairs_list = []

        neighborhood_sites = [structure[center_site_index]]

        sites_property = voronoi_finder_2.get_nn_info(structure, center_site_index)

        cutoff_weight = cutoff_weight_fraction * sum(
            voronoi_finder_2.get_weights_of_nn_sites(structure, center_site_index))

        for neighbor_index in range(0, len(sites_property)):
            neighbor_weight = voronoi_finder_2.get_nn_info(structure, center_site_index)[neighbor_index]['weight']
            neighbor_site = voronoi_finder_2.get_nn_info(structure, center_site_index)[neighbor_index]['site']
            if neighbor_weight > cutoff_weight:
                neighborhood_sites.append(neighbor_site)
                neighbor_sites_dict[neighbor_site] = neighbor_weight

        # Re-tesselate only with sites that also meet cutoff_weight_fraction criteria
        voronoi_input_coords = [site.coords for site in neighborhood_sites]
        try:
            voronoi_tess = VoronoiTess(voronoi_input_coords)
        except IndexError:
            raise RuntimeError('VoronoiTess: details unknown')

        n_neighborhood_sites = len(neighborhood_sites)
        n_voronoi_vertices = len(voronoi_tess.vertices)
        neighbor_lookup = np.zeros(shape=(n_voronoi_vertices, n_voronoi_vertices, 2), dtype=np.int8)

        # Loop over neighbor sites (which surrond center_site) to:
        # - calculate solid angles
        # - construct neighbor_lookup array (which maps each voronoi edge to two adjacent neighbors)
        for neighbor_index in range(1, n_neighborhood_sites):
            facet_vertex_list = voronoi_tess.ridges[(0, neighbor_index)]  # 0 = center site
            if 0 in facet_vertex_list:
                raise RuntimeError("Pathological structure: voronoi facet includes vertex at infinity")

            # Add this voronoi ridge to neighbor_lookup array, to identify adjacent voronoi facets
            n_vertices = len(facet_vertex_list)
            for vertex1_index in range(0, n_vertices):
                vertex2_index = vertex1_index + 1
                if vertex2_index == n_vertices:
                    vertex2_index = 0  # wrap vertex2_index

                # Properly order vertex indicies to only use upper triangle of neighbor_lookup array
                (low_vert_index, high_vert_index) = (facet_vertex_list[vertex1_index], facet_vertex_list[vertex2_index])
                if low_vert_index > high_vert_index:
                    (low_vert_index, high_vert_index) = (high_vert_index, low_vert_index)

                # Store adjacent voronoi neighbors on different levels of neighbor_lookup, at location defined
                # by common voronoi edge running between low_vert_index (row) and high_vert_index (column).
                if neighbor_lookup[low_vert_index, high_vert_index, 0] == 0:
                    neighbor_lookup[low_vert_index, high_vert_index, 0] = neighbor_index  # first neighbor
                else:
                    neighbor_lookup[low_vert_index, high_vert_index, 1] = neighbor_index  # second neighbor

        # Loop through upper triangle of neighbor_lookup and build a list of tuples of neighbor pair sites
        for low_vert_index in range(1, n_voronoi_vertices):
            for high_vert_index in range(low_vert_index + 1, n_voronoi_vertices):
                neighbor1_index = neighbor_lookup[low_vert_index, high_vert_index, 0]
                neighbor2_index = neighbor_lookup[low_vert_index, high_vert_index, 1]
                if neighbor1_index > 0 and neighbor2_index > 0:  # Both are non-zero
                    neighbor_pairs_list.append(
                        (neighborhood_sites[neighbor1_index], neighborhood_sites[neighbor2_index]))
                elif neighbor1_index != 0 or neighbor2_index != 0:  # Only one is non-zero
                    raise RuntimeError("Unexpected neighbor_lookup matrix: non-paired edge")

        # Append per site values
        all_neighbor_sites[structure.sites[center_site_index]] = neighbor_sites_dict
        all_neighbor_pairs[structure.sites[center_site_index]] = neighbor_pairs_list

    return (all_neighbor_sites, all_neighbor_pairs)



def get_voronoi_dicts_3(structure, cutoff_radius=VORONOI_CUTOFF_RADIUS,
                        cutoff_weight_fraction=VORONOI_CUTOFF_WEIGHT_FRACTION):
    voronoi_finder_2 = VoronoiNN()
    voronoi_finder_2.cutoff = cutoff_radius

    all_neighbor_sites = {}
    all_neighbor_pairs = {}

    for center_site_index in range(0, len(structure.sites)):
        neighbor_sites_dict = {}
        neighbor_pairs_list = []

        neighborhood_sites = [structure[center_site_index]]

        sites_property = voronoi_finder_2.get_nn_info(structure, center_site_index)

        cutoff_weight = cutoff_weight_fraction * sum(
            voronoi_finder_2.get_weights_of_nn_sites(structure, center_site_index))

        for neighbor_index in range(0, len(sites_property)):
            neighbor_weight = voronoi_finder_2.get_nn_info(structure, center_site_index)[neighbor_index]['weight']
            neighbor_site = voronoi_finder_2.get_nn_info(structure, center_site_index)[neighbor_index]['site']
            if neighbor_weight > cutoff_weight:
                neighborhood_sites.append(neighbor_site)
                neighbor_sites_dict[neighbor_site] = neighbor_weight

        # Append per site values
        all_neighbor_sites[structure.sites[center_site_index]] = neighbor_sites_dict
        #all_neighbor_pairs[structure.sites[center_site_index]] = neighbor_pairs_list

    return (all_neighbor_sites, all_neighbor_pairs)



def get_coordinations(voronoi_neighbor_sites):
    coordinations = []
    for neighbor_sites in voronoi_neighbor_sites.values():
        coordinations.append(sum(neighbor_sites.values()))
    return coordinations


def get_bond_lengths(voronoi_neighbor_sites):
    bond_lengths_avg = []
    bond_lengths_std = []
    for center_site, neighbor_sites in voronoi_neighbor_sites.items():
        local_bond_lengths = []
        for neighbor_site in neighbor_sites:
            local_bond_lengths.append(math.sqrt(np.sum((neighbor_site.coords - center_site.coords) ** 2.0)))
        bond_lengths_avg.append(np.mean(local_bond_lengths))
        bond_lengths_std.append(np.std(local_bond_lengths))  # sample stdev (for population, set ddof=1)
    return (bond_lengths_avg, bond_lengths_std)


def get_property_diffs(voronoi_neighbor_sites, property, abs_flag=False):
    property_diff_avg = []
    property_diff_std = []
    for center_site, neighbor_sites in voronoi_neighbor_sites.items():
        local_property_diffs = []
        for neighbor_site in neighbor_sites:
            property_diff = float(getattr(center_site.specie, property) - getattr(neighbor_site.specie, property))
            if abs_flag:
                property_diff = abs(property_diff)
            local_property_diffs.append(property_diff)
        property_diff_avg.append(np.mean(local_property_diffs))
        property_diff_std.append(np.std(local_property_diffs))  # sample stdev (for population, set ddof=1)
    return (property_diff_avg, property_diff_std)


def get_bond_angles(voronoi_neighbor_pairs):
    bond_angles_avg = []
    bond_angles_std = []
    for center_site, neighbor_sites in voronoi_neighbor_pairs.items():
        local_bond_angles = []
        for neighbor1_site, neighbor2_site in neighbor_sites:
            neighbor1_vector = neighbor1_site.coords - center_site.coords
            neighbor2_vector = neighbor2_site.coords - center_site.coords
            local_bond_angles.append(get_angle(neighbor1_vector, neighbor2_vector, units="degrees"))
        bond_angles_avg.append(np.mean(local_bond_angles))
        bond_angles_std.append(np.std(local_bond_angles))  # sample stdev (for population, set ddof=1)
    return (bond_angles_avg, bond_angles_std)
