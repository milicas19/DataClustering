import pandas as pd
from datetime import datetime
import os

def create_cell_data(data_set, cell_data_name):
    if os.path.exists(cell_data_name):
        print(f"====== Cell data for data set '{data_set}' already exists")
    else:
        start_time = datetime.now()
        print(f"====== Creating cell data for data set '{data_set}'")

        gene_exp_data = pd.read_csv(data_set, sep = '\t')

        columns = set(gene_exp_data['geneID'].values)
        columns.add('x')
        columns.add('y')

        rows = set(gene_exp_data['cell'].values)

        cell_data = pd.DataFrame(0, index=list(rows), columns=list(columns))

        for index, row in gene_exp_data.iterrows():
            if cell_data.loc[row['cell'], ['x']].values[0] == 0:
                cell_data.loc[row['cell'], ['x']] = row['x']
                cell_data.loc[row['cell'], ['y']] = row['y']
            cell_data.loc[row['cell'], [row['geneID']]] = row['MIDCounts']

        cell_data.to_csv(cell_data_name, sep='\t')

        print(f"====== job completed in: {datetime.now() - start_time}")
