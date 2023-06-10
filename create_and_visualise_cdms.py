"""

"""
import argparse

from create_cell_data import create_cell_data
from create_and_visualise_cdm import create_and_visualise_cdm

def create_and_visualise_cdms(data_sets, num_of_pca_componenets):
    for i in len(data_sets):
        data_set = data_sets[i]
        num_marker = "" if i == 0 else "_" + str(i)
        cell_data_name = "cell_data" + num_marker

        create_cell_data(data_set, cell_data_name)
        create_and_visualise_cdm(data_set, num_of_pca_componenets)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Comparing CDM values for different sets of data')
    parser.add_argument('-ds', '--data_sets', default='', type=str, nargs='+',
                        help='data sets on which CDM will be computed and visualised')
    parser.add_argument('-m', '--metric', default='euclidian', type=str, 
                        help='metric used to create CDM matrix')
    parser.add_argument('-pca', '--pca_componenets', default=50, choices=range(len(4800) + 1), type=int, 
                        help='number of PCA gene componenets')
    
    args = parser.parse_args()

    create_and_visualise_cdms(args.data_sets, args.pca_componenets)

