## global analysis
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator

folder = Path(__file__).resolve().parent
input_path = folder.parent / 'Lunar_FeO_Stage2' / 'processed_feo_grid.csv'
if not input_path.exists():
    input_path = folder / 'processed_feo_grid.csv'
d = pd.read_csv(input_path)
assert len(d) == 11306
assert d.PIXEL_INDEX.is_unique and d.FeO_wt_percent.notna().all()
assert (d.MAX_LAT > d.MIN_LAT).all() and (d.MAX_LON > d.MIN_LON).all()

# Exactly integrate geographic bin area on sphere; radius cancels for normalized means.
weights = np.deg2rad(d.MAX_LON-d.MIN_LON) * (np.sin(np.deg2rad(d.MAX_LAT)) - np.sin(np.deg2rad(d.MIN_LAT)))
assert np.isclose(weights.sum(), 4*np.pi, rtol=0, atol=1e-10), 'Bins must cover full sphere exactly once'
fe = d.FeO_wt_percent.to_numpy()
summary = [
    ('n_grid_cells', len(d), 'cells'),
    ('sphere_area_coverage_fraction', float(weights.sum()/(4*np.pi)), 'fraction'),
    ('unweighted_grid_mean', float(fe.mean()), 'wt.%'),
    ('area_weighted_global_mean', float(np.average(fe, weights=weights)), 'wt.%'),
    ('minimum_all_records', float(fe.min()), 'wt.%'),
    ('maximum_all_records', float(fe.max()), 'wt.%'),
    ('median_grid_value', float(np.median(fe)), 'wt.%'),
    ('grid_p05', float(np.percentile(fe,5)), 'wt.%'),
    ('grid_p95', float(np.percentile(fe,95)), 'wt.%'),
    ('zero_valued_records_flagged', int((fe==0).sum()), 'cells'),
    ('grid_mean_excluding_zero_sensitivity_only', float(fe[fe!=0].mean()), 'wt.%'),
    ('area_weighted_mean_excluding_zero_sensitivity_only', float(np.average(fe[fe!=0],weights=weights[fe!=0])), 'wt.%'),
]
ce4 = d.loc[d.CE4_LANDING_PIXEL]
assert len(ce4) == 1
summary.append(('CE4_landing_cell_FEO', float(ce4.iloc[0].FeO_wt_percent), 'wt.%'))
summary.append(('CE4_landing_cell_pixel_index', int(ce4.iloc[0].PIXEL_INDEX), 'index'))
pd.DataFrame(summary, columns=['metric','value','unit']).to_csv(folder/'global_statistics.csv',index=False)

# Global rectangular map, with each cell drawn from original boundaries.
vertices = []
for a,b,c,e in d[['MIN_LON','MAX_LON','MIN_LAT','MAX_LAT']].itertuples(index=False,name=None):
    vertices.append([(a,c),(b,c),(b,e),(a,e)])
fig,ax=plt.subplots(figsize=(13.2,6.8))
coll=PolyCollection(vertices, array=fe, cmap='viridis', norm=Normalize(vmin=0,vmax=30),edgecolors='none',rasterized=True)
ax.add_collection(coll)
ax.set_xlim(-180,180); ax.set_ylim(-90,90);ax.set_aspect('equal', adjustable='box')
ax.xaxis.set_major_locator(MultipleLocator(60));ax.yaxis.set_major_locator(MultipleLocator(30))
ax.set_xlabel('Longitude (degrees east; -180 to +180)');ax.set_ylabel('Latitude (degrees)')
ax.set_title('Global lunar FeO abundance | Lunar Prospector GRS 2-degree product',pad=12)
ax.plot(177.5991,-45.4446,marker='o',mfc='none',mec='white',mew=1.7,ms=10,zorder=4)
ax.annotate('Chang\'e-4 landing region',xy=(177.5991,-45.4446),xytext=(38,-78),color='white',fontsize=10,
            arrowprops={'arrowstyle':'->','color':'white','lw':1.2},ha='left',va='bottom',zorder=5)
cb=fig.colorbar(coll,ax=ax,orientation='horizontal',pad=0.12,fraction=0.055,shrink=0.86)
cb.set_label('FeO oxide-equivalent abundance (wt.%)')
fig.text(0.5,0.035,'Each rectangle represents an original PDS bin; no interpolation. The landing marker is a location, not rover-measured FeO.',ha='center',fontsize=9,color='#444444')
fig.savefig(folder/'figure1_global_feo_map.png',dpi=220,bbox_inches='tight');fig.savefig(folder/'figure1_global_feo_map.pdf',bbox_inches='tight');plt.close(fig)

# Area-weighted FeO histogram: fractional surface area in each abundance interval.
fig,ax=plt.subplots(figsize=(9,5.1))
bins=np.arange(0,32,2)
area_frac=weights/weights.sum()*100
ax.hist(fe,bins=bins,weights=area_frac,color='#477a89',edgecolor='white')
mean=float(np.average(fe,weights=weights))
ax.axvline(mean,color='#b84435',ls='--',lw=2,label=f'Area-weighted mean = {mean:.2f} wt.%')
ax.set_xlabel('FeO oxide-equivalent abundance (wt.%)')
ax.set_ylabel('Fraction of lunar surface area (%)')
ax.set_title('Area-weighted global distribution of lunar FeO')
ax.legend(frameon=False)
ax.grid(axis='y',alpha=0.2)
fig.tight_layout()
fig.savefig(folder/'figure2_area_weighted_histogram.png',dpi=220)
fig.savefig(folder/'figure2_area_weighted_histogram.pdf')
plt.close(fig)
print(pd.DataFrame(summary,columns=['metric','value','unit']).to_string(index=False))
print('Output:',folder)
