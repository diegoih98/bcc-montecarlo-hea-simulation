import numpy as np
import pandas as pd
import random
import os
from math import sqrt

def create_composition_folder(elements, molar_percentage):
    folder_name = "".join([f"{el}{perc:.2f}" for el, perc in zip(elements, molar_percentage)])
    os.makedirs(folder_name, exist_ok=True)
    return folder_name

# Prompt user for composition input
print("Enter elements and their molar percentages in the format: Element1 Percentage1 Element2 Percentage2 ...")
print("Example: Al 40 Fe 60")
user_input = input("Composition: ").split()

num_points = int(input("Enter the number of points in the lattice (recommended: 12): "))


elements = user_input[::2]  # Extract element symbols
molar_percentage = list(map(float, user_input[1::2]))  # Extract percentages

# Validate input
if len(elements) != len(molar_percentage) or not np.isclose(sum(molar_percentage), 100):
    raise ValueError("Invalid input. Ensure percentages sum to 100.")

# Convert to decimal fraction
molar_percentage = [p / 100 for p in molar_percentage]

# Create output folder
output_folder = create_composition_folder(elements, molar_percentage)

def save_results(output_folder, LRO, element_LRO):
    np.savetxt(f'{output_folder}/LRO.txt', LRO, fmt='%f', delimiter=',')
    for element, data in element_LRO.items():
        np.savetxt(f'{output_folder}/{element}.txt', data, fmt='%f', delimiter=',')


def create_bcc_lattice(a, num_points):
    
    a1 = np.array([a,0.,0.])
    a2 = np.array([0.,a,0.])
    a3 = np.array([0.,0.,a])

    Lattice_alpha = [ [0.,0.,0.] ]
    for n3 in range(2*num_points):
        for n2 in range(num_points):
            for n1 in range(num_points):
                    R = n1*a1 + n2*a2 + n3*a3              
                    Lattice_alpha = np.append(Lattice_alpha, [R],0)   
                    
    Lattice_alpha = np.delete(Lattice_alpha, 0,0)

    Lattice_beta = [ [0.,0.,0.] ]
    for n3 in range(2*num_points):
        for n2 in range(num_points):
            for n1 in range(num_points):
                    R = (0.5 + n1)*a1 + (0.5 + n2)*a2 + (0.5 + n3)*a3
                    Lattice_beta = np.append(Lattice_beta, [R],0)                
    Lattice_beta = np.delete(Lattice_beta, 0,0) 

# Create an empty array to store the result
    # Create an empty array to store the result
    row_alpha, col_alpha = np.shape(Lattice_alpha)
    row_beta, col_beta = np.shape(Lattice_beta)
    assert col_alpha == col_beta, 'number of cols should be same'
    Lattice = np.ravel([Lattice_alpha,Lattice_beta],order="F").reshape(col_alpha,row_alpha+row_beta).T
    
    return np.array(Lattice)  #,beta_index

def rearrange_by_pairs(arr):
    if len(arr) % 2 == 0:
        reshaped_arr = arr.reshape(-1, 2)

        # Shuffle the pairs along axis 0
        np.random.shuffle(reshaped_arr)

        # Flatten the reshaped array back to 1D and append the last element
        result = reshaped_arr.flatten()

    else:
        # If the array is odd-sized, shuffle the array except for the last element
        #print("Array length is odd, leaving the last element in its original position.")
        last_element = arr[-1]
        arr = arr[:-1]
        reshaped_arr = arr.reshape(-1, 2)
        np.random.shuffle(reshaped_arr)
        result = np.concatenate([reshaped_arr.flatten(), [last_element]])

    return result

def create_hea_bcc(elements, bcc_points,molar_percentage ):
    num_points = len(bcc_points)
    num_elements = len(elements)
    molar_percentage=np.array(molar_percentage)   
        #populate the lattice
        #hea_elements = np.random.choice(elements, num_points, molar_percentage)

        # Shuffle the elements to emulate disorder
        #np.random.shuffle(hea_elements)

        #start_list = ['A', 'B', 'C', 'D']

        # Count of repetitions for each element
    repetitions = molar_percentage*num_points
    
    # Check if the sum matches
    repetitions = np.round(repetitions).astype(int)
    while np.sum(repetitions) != len(bcc_points):
        # Choose a random index
        random_index = np.random.randint(0, len(repetitions))

        # Add or subtract 1 to the random element
        repetitions[random_index] += 1 if np.sum(repetitions) < len(bcc_points) else -1
        # Creating the repeated array
    hea_elements = np.repeat(elements, repetitions)
    
    hea_elements = rearrange_by_pairs(hea_elements)
    
    
    return hea_elements


def build_structured_neighbor_list(num_points):
    """
    Build a neighbor list for a structured BCC lattice with alternating α and β sites,
    using full 3D periodic boundary conditions (wrap in x, y, and z).


    Returns:
    - neighbor_list: list where neighbor_list[i] contains indices of the 8 nearest neighbors of atom i
    """
    neighbor_list = []
    total_layers = 2 * num_points # z has 2*num_points layers in your construction


    # Map (i, j, k, sublattice) to global index; sublattice: 0=alpha, 1=beta
    def index(i, j, k, sublattice):
        return 2 * (i + num_points * (j + num_points * k)) + sublattice


    for k in range(total_layers):
        for j in range(num_points):
            for i in range(num_points):
                for sublattice in [0, 1]: # 0: alpha, 1: beta
                    neighbors = []


                    # Each site has 8 NN, all on the opposite sublattice in BCC.
                    # Use wrapping (modulo) to enforce PBC in every direction.
                    for dx in [0, -1]:
                        for dy in [0, -1]:
                            for dz in [0, -1]:
                                ni = (i + dx) % num_points
                                nj = (j + dy) % num_points
                                nk = (k + dz) % total_layers
                                neighbor_idx = index(ni, nj, nk, 1 - sublattice)
                                neighbors.append(neighbor_idx)


                    # With PBC, every site has exactly 8 neighbors.
                    neighbor_list.append(neighbors)


    return neighbor_list


def find_nearest_neighbors(lattice, selected_atom_index, num_points):
    """
    Drop-in replacement for finding nearest neighbors in a structured BCC lattice.

    Inputs:
    - lattice: output from create_bcc_lattice()
    - selected_atom_index: integer index in lattice
    - num_points: number of unit cells along x and y (used to generate the lattice)

    Returns:
    - nearest_neighbors: positions of the 8 nearest neighbors (shape: [8, 3])
    - nearest_neighbor_indices: indices of those neighbors in the lattice
    """

    total_layers = 2 * num_points  # Matches create_bcc_lattice logic

    def reverse_index(idx):
        cell_index = idx // 2
        sublattice = idx % 2
        i = cell_index % num_points
        j = (cell_index // num_points) % num_points
        k = (cell_index // (num_points * num_points)) % total_layers
        return i, j, k, sublattice

    def index(i, j, k, sublattice):
        return 2 * (i + num_points * (j + num_points * k)) + sublattice

    i, j, k, sublattice = reverse_index(selected_atom_index)
    neighbor_indices = []

    for dx in [0, -1]:
        for dy in [0, -1]:
            for dz in [0, -1]:
                ni = (i + dx) % num_points
                nj = (j + dy) % num_points
                nk = (k + dz) % total_layers
                neighbor_idx = index(ni, nj, nk, 1 - sublattice)
                neighbor_indices.append(neighbor_idx)


    neighbor_positions = lattice[neighbor_indices]

    return neighbor_positions, neighbor_indices




def ordering(lattice, hea_elements, elements, c):
    LRO = np.zeros(len(elements))
    n_alpha = np.zeros(len(elements))
    n_beta = np.zeros(len(elements))
    
    for i in range(len(elements)):
        n_alpha[i] = np.count_nonzero(hea_elements[::2] == elements[i])
        n_beta[i] = np.count_nonzero(hea_elements[1::2] == elements[i])
        LRO[i] = (n_alpha[i] - n_beta[i]) / (n_alpha[i] + n_beta[i])
    
    c = np.array(c)
    Total_LRO = sqrt(np.sum(10 * c * LRO * LRO))
    
    return (Total_LRO,) + tuple(LRO)

def energy(atom_l, atom_m, nn_l, nn_m, v, hea_elements):
    delta_H = 0.0

    # Before swap: atom_l with its neighbors
    for neighbor in nn_l:
        if neighbor == atom_m:
            continue
        delta_H -= v[hea_elements[atom_l]][hea_elements[neighbor]]

    # Before swap: atom_m with its neighbors
    for neighbor in nn_m:
        if neighbor == atom_l:
            continue
        delta_H -= v[hea_elements[atom_m]][hea_elements[neighbor]]

    # After swap: atom_l becomes atom_m
    for neighbor in nn_l:
        if neighbor == atom_m:
            continue
        delta_H += v[hea_elements[atom_m]][hea_elements[neighbor]]

    # After swap: atom_m becomes atom_l
    for neighbor in nn_m:
        if neighbor == atom_l:
            continue
        delta_H += v[hea_elements[atom_l]][hea_elements[neighbor]]

    return delta_H

        
def MC_precomputed(bcc_points, hea_elements, elements, molar_percentage, vij_df, sample, T, neighbor_list):
    for _ in range(sample):
        atom_l = random.randint(0, len(hea_elements) - 1)
        nn_l_indices = neighbor_list[atom_l]
        if len(nn_l_indices) == 0:
            continue
        atom_m = np.random.choice(nn_l_indices)
        nn_m_indices = neighbor_list[atom_m]

        H = energy(atom_l, atom_m, nn_l_indices, nn_m_indices, vij_df, hea_elements)

        if H >= 0 and np.exp(-H / (8.617333262145E-5 * T)) < random.uniform(0.0, 1.0):
            continue

        hea_elements[atom_l], hea_elements[atom_m] = hea_elements[atom_m], hea_elements[atom_l]

    return (hea_elements,) + ordering(bcc_points, hea_elements, elements, molar_percentage)




        



# Load DataFrame from Excel
v = pd.read_excel("corrected_vij_matrix.xlsx", index_col=0)

# Convert to nested dictionary
vij_df= {
    row: {col: v.loc[row, col] for col in v.columns}
    for row in v.index
}
    
# SET-UP

k = 8.617333262145E-5  # Boltzmann constant in eV/K
T = list(range(1800, 0, -1))  # Temperature in Kelvin
sample = 10000  # Number of swaps per temperature
a = 1.0  # Lattice constant

# Create lattice and initial configuration
bcc_points = create_bcc_lattice(a, num_points)
hea_elements = create_hea_bcc(elements, bcc_points, molar_percentage)
neighbor_list = build_structured_neighbor_list(num_points)

# Prepare storage for LRO tracking
LRO = np.zeros(len(T))
element_LRO = {el: np.zeros(len(T)) for el in elements}

# Run MC simulation over temperatures
for j in range(len(T)):
    print(f"T = {T[j]} K")
    results = MC_precomputed(bcc_points, hea_elements, elements, molar_percentage,
                             vij_df, sample, T[j], neighbor_list)
    hea_elements, LRO[j] = results[0], results[1]
    for i, el in enumerate(elements):
        element_LRO[el][j] = results[i + 2]


    
  
    


save_results(output_folder, LRO, element_LRO)



