import pandas as pd

gene_exp_data = pd.read_csv('gene_expression.tsv',sep = '\t')
print(gene_exp_data)

columns = set(gene_exp_data['geneID'].values)
columns.add('x')
columns.add('y')
print(len(columns))

rows = set(gene_exp_data['cell'].values)
print(len(rows))

cell_data = pd.DataFrame(0, index=list(rows), columns = list(columns))
print(cell_data)

for index, row in gene_exp_data.iterrows():
  print(index)
  if cell_data.loc[row['cell'], ['x']].values[0] == 0:
    cell_data.loc[row['cell'], ['x']] = row['x']
    cell_data.loc[row['cell'], ['y']] = row['y']
  cell_data.loc[row['cell'], [row['geneID']]] = row['MIDCounts']

print(cell_data)

cell_data.to_csv('cell_data.tsv', sep='\t')
