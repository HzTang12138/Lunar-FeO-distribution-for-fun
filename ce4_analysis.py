"""Stage 4 — Chang'e-4 regional FeO analysis from Stage 2 processed NASA PDS table.
Run: python stage4_ce4_analysis.py
Dependencies: numpy, pandas, matplotlib.
The neighborhood uses center-based great-circle distance on a sphere, not
partial-cell clipping. Values refer to orbital grid abundance, not rover VNIS.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator

out = Path(__file__).resolve().parent
source = out.parent/'Lunar_FeO_Stage3'/'processed_feo_grid.csv'
if not source.exists(): source=out/'processed_feo_grid.csv'
d = pd.read_csv(source)
assert len(d)==11306
LON0,LAT0=177.5991,-45.4446
lat0,lon0=np.deg2rad([LAT0,LON0])
lat=np.deg2rad(d.CENTER_LAT.to_numpy()); lon=np.deg2rad(d.CENTER_LON.to_numpy())
cosd=np.sin(lat0)*np.sin(lat)+np.cos(lat0)*np.cos(lat)*np.cos(lon-lon0)
d['angular_distance_deg']=np.rad2deg(np.arccos(np.clip(cosd,-1,1)))
w=np.deg2rad(d.MAX_LON-d.MIN_LON)*(np.sin(np.deg2rad(d.MAX_LAT))-np.sin(np.deg2rad(d.MIN_LAT)))
d['solid_angle_sr']=w
fe=d.FeO_wt_percent.to_numpy()
global_mean=np.average(fe,weights=w)
land=d.loc[(d.MIN_LAT<=LAT0)&(LAT0<d.MAX_LAT)&(d.MIN_LON<=LON0)&(LON0<d.MAX_LON)]
assert len(land)==1 and int(land.iloc[0].PIXEL_INDEX)==1716
single=float(land.iloc[0].FeO_wt_percent)
rows=[{'sample':'Global (all cells)','selection':'Whole Moon','n_cells':len(d),'area_weighted_mean_wt_percent':global_mean,
       'minimum_wt_percent':fe.min(),'maximum_wt_percent':fe.max(),'note':'Area-weighted full sphere'}]
for radius in [5,10,15]:
    cut=d[d.angular_distance_deg<=radius]
    a=np.average(cut.FeO_wt_percent, weights=cut.solid_angle_sr)
    rows.append({'sample':f'CE4 neighborhood {radius} deg','selection':f'Bin centers within {radius} degrees angular distance',
                 'n_cells':len(cut),'area_weighted_mean_wt_percent':a,
                 'minimum_wt_percent':cut.FeO_wt_percent.min(), 'maximum_wt_percent':cut.FeO_wt_percent.max(),
                 'note':'Selection is by grid center; boundary bins not fractionally clipped'})
rows.append({'sample':'CE4 landing bin','selection':'Grid containing exact landing coordinate','n_cells':1,
             'area_weighted_mean_wt_percent':single,'minimum_wt_percent':single,'maximum_wt_percent':single,
             'note':'Orbital 2-degree bin, not a Yutu-2 in-situ measurement'})
stats=pd.DataFrame(rows)
stats.to_csv(out/'ce4_comparison_statistics.csv',index=False)
region=d[d.angular_distance_deg<=10].copy()
region.to_csv(out/'ce4_neighborhood_10deg_cells.csv',index=False)

# Scientific local map: shift longitude edges continuously about 180-degree seam.
# Include all native bin rectangles intersecting the viewport, without resampling.
view=d[(d.MAX_LAT>=-63)&(d.MIN_LAT<=-28)].copy()
polys=[];vals=[]
for rec in view.itertuples(index=False):
    mid=(rec.MIN_LON+rec.MAX_LON)/2
    shift=360*round((LON0-mid)/360)
    lo,hi=rec.MIN_LON+shift,rec.MAX_LON+shift
    if hi<152 or lo>204: continue
    polys.append([(lo,rec.MIN_LAT),(hi,rec.MIN_LAT),(hi,rec.MAX_LAT),(lo,rec.MAX_LAT)])
    vals.append(rec.FeO_wt_percent)
fig,ax=plt.subplots(figsize=(10,6.3))
coll=PolyCollection(polys,array=np.asarray(vals),cmap='viridis',norm=Normalize(vmin=0,vmax=25),edgecolors='none',rasterized=True)
ax.add_collection(coll)
ax.plot(LON0,LAT0,marker='*',color='white',markeredgecolor='black',markeredgewidth=.8,markersize=16,zorder=8,label="Chang'e-4 landing site")
ax.set_xlim(152,204);ax.set_ylim(-63,-28);ax.set_aspect('equal',adjustable='box')
ax.set_xticks([155,165,175,185,195,200]);ax.set_xticklabels(['155°E','165°E','175°E','175°W','165°W','160°W'])
ax.set_yticks([-60,-55,-50,-45,-40,-35,-30]);ax.set_yticklabels(['60°S','55°S','50°S','45°S','40°S','35°S','30°S'])
ax.set_xlabel('Longitude (wrapped continuously at 180°)');ax.set_ylabel('Latitude')
ax.set_title("FeO abundance around Chang'e-4 | Lunar Prospector GRS",pad=14)
# 10-degree angular neighborhood boundary on sphere for visual guide.
angles=np.linspace(0,2*np.pi,600); rr=np.deg2rad(10)
latc=np.arcsin(np.sin(lat0)*np.cos(rr)+np.cos(lat0)*np.sin(rr)*np.cos(angles))
lonc=lon0+np.arctan2(np.sin(angles)*np.sin(rr)*np.cos(lat0),np.cos(rr)-np.sin(lat0)*np.sin(latc))
long=np.rad2deg(lonc); long=LON0+((long-LON0+180)%360-180)
ax.plot(long,np.rad2deg(latc),ls='--',color='white',lw=1.6,zorder=7,label='10° angular distance')
ax.legend(loc='upper right',framealpha=.94)
cbar=fig.colorbar(coll,ax=ax,pad=.015,shrink=.85); cbar.set_label('FeO abundance (wt.%)')
fig.text(.5,.025,'Colors: original orbital-grid bins; dashed circle: 10° angular neighborhood. No interpolation or rover measurements.',ha='center',fontsize=8.7)
fig.tight_layout(rect=(0,.04,1,1))
fig.savefig(out/'figure3_ce4_regional_feo_map.png',dpi=240,bbox_inches='tight')
fig.savefig(out/'figure3_ce4_regional_feo_map.pdf',bbox_inches='tight')
plt.close(fig)

# Plain category comparison. Caveat: local samples overlap and are not statistically independent.
labels=['Global\n(area-weighted)','CE-4 landing\nbin','CE-4 ±5°\nneighborhood','CE-4 ±10°\nneighborhood','CE-4 ±15°\nneighborhood']
values=[global_mean,single,*[float(stats.loc[stats['sample']==f'CE4 neighborhood {r} deg','area_weighted_mean_wt_percent'].iloc[0]) for r in (5,10,15)]]
fig,ax=plt.subplots(figsize=(9.6,5.1))
bars=ax.bar(labels,values,color=['#788b9a','#b66a4f','#5e9d96','#447e86','#7da5ac'])
for bar,v in zip(bars,values):
    ax.text(bar.get_x()+bar.get_width()/2,v+.12,f'{v:.2f}',ha='center',fontsize=10)
ax.set_ylim(0,max(values)+1.5)
ax.set_ylabel('FeO abundance (wt.%)')
ax.set_title("Global versus Chang'e-4 orbital-grid FeO abundance")
ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
fig.text(.5,.01,'Neighborhood radius is angular distance from landing site to grid center. All measurements are from Lunar Prospector GRS.',ha='center',fontsize=8.4)
fig.tight_layout(rect=(0,.04,1,1))
fig.savefig(out/'figure4_ce4_global_comparison.png',dpi=240,bbox_inches='tight')
fig.savefig(out/'figure4_ce4_global_comparison.pdf',bbox_inches='tight')
plt.close(fig)

print(stats.to_string(index=False))
