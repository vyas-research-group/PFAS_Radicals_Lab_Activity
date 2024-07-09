import numpy as np
from math import isclose
import os
from . import module1_functions as m1
import pandas as pd

# ----------------------------
# students write this function
# ----------------------------
def calculate_z_displacement(bond_angle_deg : float) -> float:
    bond_angle_rad = np.deg2rad(bond_angle_deg)
    z_displacement =  np.sqrt((0.5 * np.sqrt(3) / np.sin(bond_angle_rad / 2))**2 - 1)
    return z_displacement

# ----------------------------
# students write this function
# ----------------------------
def scale_bond_lengths(bond_length : float, carbon_xyz_position, heteroatom_xyz_position_list : list) -> list:
    new_bond_vectors = []
    for heteroatom_position in heteroatom_xyz_position_list:
            bond_vector = heteroatom_position - carbon_xyz_position
            bond_unit_vector =  bond_vector / np.linalg.norm(bond_vector)
            new_bond_vectors.append(bond_unit_vector * bond_length)
    return new_bond_vectors

def generate_geometry(bond_angle_deg : float, bond_length : float) -> list:

    # generate geometry using simple scaled diagram
    z_displacement = calculate_z_displacement(bond_angle_deg)
    assert z_displacement >= 0, "The z displacement is either nan or less than zero. Check that you are not taking the square root of a negative number."

    heteroatom_xyz_positions = [np.array([1, 0, 0]), np.array([-np.cos(np.deg2rad(60)), np.sin(np.deg2rad(60)), 0]), np.array([-np.cos(np.deg2rad(60)), -np.sin(np.deg2rad(60)), 0])]
    carbon_xyz_position = np.array([0, 0, z_displacement])

    # adjust bond lengths to the equilibrium lengths
    heteroatom_xyz_positions = scale_bond_lengths(bond_length, carbon_xyz_position, heteroatom_xyz_positions)
    assert False not in check_bond_angles(bond_angle_deg, heteroatom_xyz_positions), "There actual angle is not equal to the input angle.\n\t\t \
    Check the z displacement calculation and ensure you are not changing the angles when scaling bond lengths."

    return heteroatom_xyz_positions

def make_input_file_text(heteroatom : str, bond_angle_deg : float, bond_length : float) -> str:

    # generate coordinates
    heteroatom_xyz_positions = generate_geometry(bond_angle_deg, bond_length)

    # format coordinates
    formatted_lines = []
    for heteroatom_coordinates in heteroatom_xyz_positions:
        atom_line = heteroatom + " " + " ".join(([ "{:0.10f}".format(coordinate) for coordinate in heteroatom_coordinates ]))
        formatted_lines.append(atom_line)
    heteroatom_section = "\n".join(formatted_lines)

    angle = "{:0.1f}".format(bond_angle_deg)

    input_file_text = f"""
    # C{heteroatom}3 radical {angle} degrees fixed
    ! UKS TightSCF wB97x-D3 def2-TZVPD xyzfile opt freq
    %geom Constraints
            {{ A * 0 * C }}
        end
    end
    * xyz 0 2
    C 0.00 0.00 0.00
    {heteroatom_section}
    *
    """
    input_file_text = format_text_indentation(input_file_text)
    print(input_file_text)
    return input_file_text

def write_input_file(heteroatom : str, bond_angle_deg : float, bond_length : float):

    filename = f"C{heteroatom}3_radical_{bond_angle_deg}_degrees.inp"
    input_file_text = make_input_file_text(heteroatom, bond_angle_deg, bond_length)

    with open(filename, "w") as f:
        f.write(input_file_text)

# ----------------------------
# students write this function
# ----------------------------
def write_coordinate_scan_input_files(heteroatom, bond_length, low : float, high : float, step : float):
    scan_angles = np.arange(low, high, step)
    for angle in scan_angles:
        write_input_file(heteroatom, angle, bond_length)

def parse_outfiles_from_folder(folderpath : str) -> list:
    folder_data = []
    files = os.listdir(folderpath)
    outfiles = list(filter(lambda ext: ".out" in ext, files))
    for file in outfiles:
        print(f"NOW PARSING {file}")
        filepath = os.path.join(folderpath, file)
        data_names, data = m1.parse_outfile(filepath)
        file_data = [file] + data
        folder_data.append(file_data)
    print("PARSING COMPLETE")
    folder_data = [['filename', *data_names]] + folder_data
    return folder_data

def write_data_to_csv(folder_data : list):
    dataframe = pd.DataFrame(folder_data[1:], columns=folder_data[0])
    dataframe.to_csv("./summary_data.csv", index=False)
 
# -----------------------------
# other custom helper functions
# -----------------------------
def check_bond_angles(bond_angle_deg : float, bond_vectors : float) -> list:
    n = len(bond_vectors)
    bond_vectors = bond_vectors + bond_vectors
    results = []
    for bond_index in range(n):
        dot_product = np.dot(bond_vectors[bond_index], bond_vectors[bond_index + 1])
        norms = np.linalg.norm(bond_vectors[bond_index]) * np.linalg.norm(bond_vectors[bond_index + 1])
        actual_bond_angle = np.rad2deg(np.arccos(dot_product / norms))
        results.append(isclose(bond_angle_deg, actual_bond_angle, abs_tol=4))
    return results

def format_text_indentation(text : str) -> str:
    newlines= []
    for line in text.split("\n"):
        if line[:4] == ' '*4 :
            newlines.append(line[4:])
        else:
            newlines.append(line)
    return "\n".join(newlines)
