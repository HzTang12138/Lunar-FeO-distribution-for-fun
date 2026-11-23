## Data Process
from pathlib import Path
import numpy as np
import pandas as pd

base = Path(__file__).resolve().parent
raw = base / 'lpgrs_high1_elem_abundance_2deg.tab'
if not raw.exists():
    raw = base.parent / 'lpgrs_high1_elem_abundance_2deg.tab'
output = base / 'processed_feo_grid.csv'
cols = ['PIXEL_INDEX', 'MIN_LAT', 'MAX_LAT', 'MIN_LON', 'MAX_LON',
        'AM', 'NEUTRON_DEN', 'W_MGO', 'W_AL2O3', 'W_SIO2', 'W_CAO',
        'W_TIO2', 'W_FEO', 'W_K', 'W_TH', 'W_U']
all_cols = cols + [f'ERROR_MATRIX_TERM_{i:02d}' for i in range(16,61)]
df = pd.read_csv(raw, sep=r'\s+', header=None, names=all_cols)
assert df.shape == (11306, 61), f'Unexpected shape: {df.shape}'
assert df['PIXEL_INDEX'].is_unique
assert df[all_cols].notna().all().all()
assert (df['MAX_LAT'] > df['MIN_LAT']).all() and (df['MAX_LON'] > df['MIN_LON']).all()
site_lat, site_lon = -45.4446, 177.5991
site = df[(df.MIN_LAT <= site_lat) & (site_lat < df.MAX_LAT) &
          (df.MIN_LON <= site_lon) & (site_lon < df.MAX_LON)]
assert len(site) == 1, f'Expected one landing grid, got {len(site)}'
prepared = df[cols[:5] + ['W_FEO']].copy()
prepared['FeO_wt_percent'] = prepared['W_FEO'] * 100.0
prepared['CENTER_LAT'] = (prepared.MIN_LAT + prepared.MAX_LAT) / 2
prepared['CENTER_LON'] = (prepared.MIN_LON + prepared.MAX_LON) / 2
prepared['CE4_LANDING_PIXEL'] = prepared.PIXEL_INDEX.eq(int(site.iloc[0].PIXEL_INDEX))
prepared['ZERO_FEO_REVIEW'] = prepared.W_FEO.eq(0)
prepared.to_csv(output, index=False, float_format='%.9g')
print('Raw table:', df.shape)
print('FeO range (wt.%):', prepared.FeO_wt_percent.min(), prepared.FeO_wt_percent.max())
print('Zero FeO values (flagged, not removed):', int(prepared.ZERO_FEO_REVIEW.sum()))
print('Landing grid:', prepared.loc[prepared.CE4_LANDING_PIXEL].to_string(index=False))
print('Saved:', output)
