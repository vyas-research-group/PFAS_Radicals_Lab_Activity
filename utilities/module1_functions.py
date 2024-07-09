import cclib.io as cc
from cclib.parser import utils as ccp
import numpy as np


# ----------------------------
# students write this function
# ----------------------------
def calculate_bond_length(atom1_xyz, atom2_xyz) -> float:
    bond_length = np.linalg.norm(atom2_xyz - atom1_xyz)
    return bond_length

# ----------------------------
# students write this function
# ----------------------------
def calculate_bond_angle(atom1_xyz, atom2_xyz, atom3_xyz) -> float:
    bond1_unit_vector = (atom1_xyz - atom2_xyz) / np.linalg.norm(atom1_xyz - atom2_xyz)
    bond2_unit_vector = (atom3_xyz - atom2_xyz) / np.linalg.norm(atom3_xyz - atom2_xyz)

    bond_angle_radians = np.arccos(np.dot(bond1_unit_vector, bond2_unit_vector))
    bond_angle_degrees = np.rad2deg(bond_angle_radians)

    return bond_angle_degrees

def get_average_bond_lengths_and_angles(optimized_geometry) -> "tuple[float, float]":
    carbon_xyz_position = optimized_geometry[0]
    heteroatom1_xyz_position = optimized_geometry[1]
    heteroatom2_xyz_position = optimized_geometry[2]
    heteroatom3_xyz_position = optimized_geometry[3]

    # bond lengths
    bond1_length = calculate_bond_length(heteroatom1_xyz_position, carbon_xyz_position)
    bond2_length = calculate_bond_length(heteroatom2_xyz_position, carbon_xyz_position)
    bond3_length = calculate_bond_length(heteroatom3_xyz_position, carbon_xyz_position)
    bond_lengths = [bond1_length, bond2_length, bond3_length]
    average_bond_length = np.average(bond_lengths)
    bond_length_range = np.max(bond_lengths) - np.min(bond_lengths)
    assert bond_length_range < 0.1 , "The range in bond lengths is greater than 0.1 angstrom. Check the structure with a molecular visualizer." 

    # bond angles
    bond_angle1 = calculate_bond_angle(heteroatom1_xyz_position, carbon_xyz_position, heteroatom2_xyz_position)
    bond_angle2 = calculate_bond_angle(heteroatom2_xyz_position, carbon_xyz_position, heteroatom3_xyz_position)
    bond_angle3 = calculate_bond_angle(heteroatom3_xyz_position, carbon_xyz_position, heteroatom1_xyz_position)
    bond_angles = [bond_angle1, bond_angle2, bond_angle3]
    average_bond_angle = np.average(bond_angles)
    bond_angle_range = np.max(bond_angles) - np.min(bond_angles)
    assert bond_angle_range < 0.05 , "The range in bond angles is greater than 0.05 degrees. Check the structure with a molecular visualizer." 

    return average_bond_length, average_bond_angle

# ----------------------------
# students write this function
# ----------------------------
def calculate_force_constant(IR_frequency_wavenumbers : float, reduced_mass_amu : float) -> float:
    # define conversion factors
    kg_per_amu = 1.660538921E-27  # [kg/amu]
    speed_of_light = 2.99792458E10 # [cm/s]
    # convert values
    reduced_mass_kg = reduced_mass_amu * kg_per_amu # [kg]
    IR_frequency_Hz = IR_frequency_wavenumbers * speed_of_light # [Hz]
    # compute force constant
    force_constant = (2 * np.pi * IR_frequency_Hz)**2 * reduced_mass_kg # [kg/s^2]
    return force_constant

# ----------------------------
# students write this function
# ----------------------------
def get_partial_charges(outfile) -> list:
    carbon_partial_charge = outfile.atomcharges['mulliken'][0] # [-]
    heteroatom_partial_charge = np.average(outfile.atomcharges['mulliken'][1:4]) # [-]
    return [carbon_partial_charge, heteroatom_partial_charge]

def get_frontier_orbital_energies(outfile, heteroatom) -> list:
    # get the orbital number for the frontier orbitals
    n_homo = min(outfile.homos)
    n_somo = max(outfile.homos)
    n_lumo = n_somo + 1

    # convert the energies to an array
    moenergies = np.array(outfile.moenergies)

    # assign the mo energies to variables
    # note that the energies are normalized to a single electron
    homo_energy = sum(moenergies[:, n_homo])/2 # [hartree]
    somo_a_energy = moenergies[0, n_somo] # [hartree]
    somo_b_energy = moenergies[1, n_somo] # [hartree]
    lumo_energy = sum(moenergies[:, n_lumo])/2 # [hartree]

    lumo_plus_1_energy = sum(moenergies[:, n_lumo+1])/2 # [hartree]
    lumo_plus_2_energy = sum(moenergies[:, n_lumo+1])/2 # [hartree]

    # get the orbital energies for orbitals below the HOMO
    core_orbital_energies = []
    if heteroatom == 'F':
        for index in range(1,15):
            energy = sum(moenergies[:, n_homo-index])/2 # [hartree]
            core_orbital_energies.append(energy)

    elif heteroatom == 'H':
        for index in range(1,15):
            if index <= 3:
                # the CH3 radical has fewer MOs than CF3
                # but we want the data to be the same shape for both molecules
                energy = sum(moenergies[:, n_homo-index])/2 # [hartree]
            else:
                energy = 0
            core_orbital_energies.append(energy)

    datalist = [homo_energy, somo_a_energy, somo_b_energy, lumo_energy, \
                lumo_plus_1_energy, lumo_plus_2_energy, *core_orbital_energies]

    return datalist


def parse_outfile(outfile_path : str) -> "tuple[list, list]":

    outfile = cc.ccread(outfile_path)

    # get the heteroatom used in the calculation
    heteroatoms = {"9" : "F", "1" : "H"}
    heteroatom = heteroatoms[str(outfile.atomnos[1])]

    # get the SCF energy in hartree
    SCF_energy = ccp.convertor(outfile.scfenergies[-1], "eV", "hartree") # [hartree]

    # compute the bond length and angle
    optimized_geometry = outfile.atomcoords[-1]
    bond_length, bond_angle = get_average_bond_lengths_and_angles(optimized_geometry) # [angstroms] , [degrees]

    # get the IR frequency of the C-X bond stretch
    IR_frequency = outfile.vibfreqs[-2] # [wavenumbers]

    # compute the force constant of the C-X bond stretch
    carbon_mass = 12.011 # [amu]
    heteroatom_mass = {'H' : 1.00784, 'F' : 18.998403} # [amu]
    reduced_mass = heteroatom_mass[heteroatom] * carbon_mass / (heteroatom_mass[heteroatom] + carbon_mass) # [amu]
    force_constant = calculate_force_constant(IR_frequency, reduced_mass) # [kg/s^2]
    
    # get the partial charges
    partial_charges = get_partial_charges(outfile)

    # get the frontier orbital energies
    frontier_orbital_energies = get_frontier_orbital_energies(outfile, heteroatom)

    # compile the data to return
    data_names = ['heteroatom', 'bond_angle[deg]', 'scf_energy[hartree]', 'bond_length[angstrom]', 'ir_frequency[wavenumbers]', 'force_constant[mg/s^2]', \
                  'carbon_partial_charge[-]', 'heteroatom_partial_charge[-]', 'homo_energy[hartree]', 'somo-a_energy[hartree]', 'somo-b_energy[hartree]', 'lumo_energy[hartree]', \
                  'lumo+1_energy[hartree]', 'lumo+2_energy[hartree]', 'homo-1_energy[hartree]', 'homo-2_energy[hartree]', 'homo-3_energy[hartree]', \
                  'homo-4_energy[hartree]', 'homo-5_energy[hartree]', 'homo-6_energy[hartree]', 'homo-7_energy[hartree]', 'homo-8_energy[hartree]', \
                  'homo-9_energy[hartree]', 'homo-10_energy[hartree]', 'homo-11_energy[hartree]', 'homo-12_energy[hartree]', 'homo-13_energy[hartree]', \
                  'homo-14_energy[hartree]', 'homo-15_energy[hartree]']
    
    data = [heteroatom, bond_angle, SCF_energy, bond_length, IR_frequency, force_constant, *partial_charges, *frontier_orbital_energies]
    
    return data_names, data
