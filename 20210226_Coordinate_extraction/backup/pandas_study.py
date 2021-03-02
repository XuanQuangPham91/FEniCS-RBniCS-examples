import pandas as pd

df = pd.DataFrame({
    'name': ['Raphael', 'Donatello'],
    'mask': ['red', 'purple'],
    'weapon': ['sai', 'bo staff']
})
df.to_csv(index=False)

compression_opts = dict(method='zip', archive_name='out.csv')
df.to_csv('out.zip', index=False, compression=compression_opts)
