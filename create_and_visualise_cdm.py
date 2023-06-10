import pandas as pd
import copy
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import scipy.spatial.distance

def read_data(cell_data_name):
    cell_data = pd.read_csv(cell_data_name, sep='\t')
    cell_data.rename(columns={cell_data.columns[0]: 'cellID'}, inplace=True)

    return cell_data

def adjacency_matrix(data_dict, cells):
    n = len(cells)
    adjacency_matrix = np.zeros(shape=(n, n))

    for i in range(0, n - 1):
        for j in range(i + 1, n):
            adjacency_matrix[i][j] = scipy.spatial.distance.euclidean(data_dict[cells[i]], data_dict[cells[j]])
            adjacency_matrix[j][i] = adjacency_matrix[i][j]

    return adjacency_matrix

def reduce_gene_expr(cell_data, num_of_pca_components):
    gene_exp = cell_data.copy()
    gene_exp.drop('cellID', axis=1, inplace=True)
    gene_exp.drop('x', axis=1, inplace=True)
    gene_exp.drop('y', axis=1, inplace=True)

    scaling = StandardScaler()
    scaling.fit(gene_exp)
    scaled_data = scaling.transform(gene_exp)

    pca = PCA(n_components=num_of_pca_components)
    reduced_data = pca.fit_transform(scaled_data)
    pca.explained_variance_ratio_.cumsum()

    return reduced_data

def normalize(adjacency_matrix):
    max_value = np.max(adjacency_matrix)

    adjacency_matrix /= max_value

    return adjacency_matrix

def normalize_adjacency_matrices(cell_data, num_of_pca_componenets):
    cell_dict = {}
    cells = []

    for index, row in cell_data.iterrows():
        cell_dict[row['cellID']] = (row['x'], row['y'])
        cells.append(row['cellID'])

    adjacency_matrix_1 = adjacency_matrix(cell_dict, cells)

    reduced_data = reduce_gene_expr(cell_data, num_of_pca_componenets)
    reduced_data_dict = {}
    for i in range(0, len(reduced_data)):
        reduced_data_dict[cells[i]] = reduced_data[i]

    adjacency_matrix_2 = adjacency_matrix(reduced_data_dict, cells)

    norm_adj_matrix_1 = normalize(adjacency_matrix_1)
    norm_adj_matrix_2 = normalize(adjacency_matrix_2)

    return norm_adj_matrix_1, norm_adj_matrix_2

def cdm_1(norm_adj_matrix_1, norm_adj_matrix_2):
    CDM = norm_adj_matrix_1 - norm_adj_matrix_2
    upper_indices = np.triu_indices(CDM.shape[0], k=1)
    upper_values = CDM[upper_indices]

    figure, ax = plt.subplots(nrows=1, ncols=1)
    figure.dpi = 100
    figure.set_figheight(10)
    figure.set_figwidth(16)
    _ = ax.hist(upper_values, bins=1000)

def create_and_visualise_cdm(data_set, num_of_pca_componenets):
    start_time = datetime.now()
    print(f"====== Creating and visualising CDM for '{data_set}'")

    cell_data = read_data(data_set)
    norm_adj_matrix_1, norm_adj_matrix_2 = normalize_adjacency_matrices(cell_data, num_of_pca_componenets)
    cdm_1(norm_adj_matrix_1, norm_adj_matrix_2)

    print(f"====== job completed in: {datetime.now() - start_time}")
    