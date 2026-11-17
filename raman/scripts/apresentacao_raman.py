"""Create paired thesis figures from saved spectra and exploratory fits.
Both panels use the same spectrum. The thesis selection is read from JSON.
This script reads saved fits and does not refit spectra.
"""
from pathlib import Path
import json,csv,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
PACKAGE_ROOT=Path(__file__).resolve().parents[1]
SOURCE=None
OUTPUT=None
TITLES={
'nfes_ca':'Carbonized CA (NFES)',
'ffes_ca':'Carbonized CA (FFES)',
'ffes_ca_cb10':'Carbonized CA + 10% carbon black (FFES)',
'ffes_ca_cb20':'Carbonized CA + 20% carbon black (FFES)',
'ffes_ca_graphite':'Carbonized CA + 0.2% graphite (FFES)',
'ffes_novolak':'Carbonized Novolak (FFES)',
'ffes_novolak_cb15':'Carbonized Novolak + 15% carbon black (FFES)',
'cb_reference':'Carbon black reference'}
COLORS=['#8c6bb1','#d95f02','#1b9e77','#2774ae','#cc79a7']
BANDS=['D4','D1','D3','G','D2']
def plot(r):
 a=np.load(SOURCE/(r['id']+'.npz'));norm=max(a['fit_y']);x=a['fit_x'];y=a['fit_y']/norm
 fig=plt.figure(figsize=(10.4,4.9),layout='constrained')
 gs=fig.add_gridspec(2,2,height_ratios=[3.4,1])
 left=fig.add_subplot(gs[:,0]);right=fig.add_subplot(gs[0,1]);res=fig.add_subplot(gs[1,1])
 fig.suptitle(TITLES[r['group']],fontsize=14,fontweight='bold')
 left.plot(x,y,color='black',lw=.9)
 for key,name,col,level in [('observed_D_cm1','D','#d95f02',1.16),('observed_G_cm1','G','#2774ae',1.08)]:
  pos=r[key];left.axvline(pos,color=col,lw=1,ls='--')
  left.annotate(f'{name}: {pos:.1f}'+r' cm$^{-1}$',(pos,level),xytext=(-6,0),textcoords='offset points',ha='right',va='top',fontsize=11,color=col)
 left.set_ylim(-.1,1.24);left.set_title('(a) Baseline-corrected spectrum',fontsize=12)
 left.set_xlabel(r'Raman shift (cm$^{-1}$)');left.set_ylabel('Normalized intensity (a.u.)')
 right.plot(x,y,color='.5',lw=.7,label='Measured')
 right.plot(x,a['fitted']/norm,color='black',lw=1.25,label='Total fit')
 for band,curve,c in zip(BANDS,a['components'],COLORS):
  right.plot(x,curve/norm,color=c,lw=1.1,label=band);right.fill_between(x,0,curve/norm,color=c,alpha=.13)
 right.set_title('(b) Five-band fit and band areas',fontsize=12);right.set_ylabel('Normalized intensity (a.u.)')
 right.set_ylim(-.12,1.45);right.legend(loc='upper left',ncol=3,fontsize=9.5,frameon=False,columnspacing=.8,handlelength=1.5)
 right.tick_params(labelbottom=False)
 res.plot(x,a['residual']/norm,color='.3',lw=.7);res.axhline(0,color='.6',lw=.7);res.set_ylabel('Residual');res.set_xlabel(r'Raman shift (cm$^{-1}$)')
 for ax in [left,right,res]:
  ax.set_xlim(1000,1800);ax.spines[['top','right']].set_visible(False)
  ax.set_xticks([1000,1200,1400,1600,1800])
 for ext in ['png','svg']:fig.savefig(OUTPUT/(r['id']+'_pair.'+ext),dpi=300)
 plt.close(fig)
def main():
    import argparse
    from datetime import datetime
    global SOURCE, OUTPUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=PACKAGE_ROOT / "saved_fits")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--selection", type=Path, default=PACKAGE_ROOT / "presentation_selection.json")
    args = parser.parse_args()
    SOURCE = args.results
    OUTPUT = args.output or PACKAGE_ROOT / "results" / ("figures_" + datetime.now().strftime("%Y%m%d_%H%M%S"))
    metadata = json.loads((SOURCE / "resultados.json").read_text(encoding="utf-8"))
    results = metadata["results"]
    frozen_selection = json.loads(args.selection.read_text(encoding="utf-8"))
    present_groups = list(dict.fromkeys(result["group"] for result in results))
    selection = {}
    summary = []
    OUTPUT.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({"font.size":12,"axes.labelsize":12,"xtick.labelsize":10,"ytick.labelsize":10})
    for group in present_groups:
        candidates = [result for result in results if result["group"] == group]
        selected_id = frozen_selection[group]["selected"]
        if selected_id not in {result["id"] for result in candidates}:
            raise ValueError("Missing selected spectrum: " + selected_id)
        ordered = sorted(candidates, key=lambda result: result["id"] != selected_id)
        selection[group] = {
            "selected": selected_id,
            "additional": [result["id"] for result in ordered[1:]],
            "title": TITLES[group],
        }
        for result in ordered:
            plot(result)
            summary.append({
                "material": TITLES[group], "id": result["id"],
                "location": "main" if result["id"] == selected_id else "appendix",
                "height_ratio": result["height_ratio"], "area_ratio": result["area_ratio"],
                "d3_area_fraction": result["d3_area_fraction"],
                "warnings": "; ".join(result["warnings"]),
            })
        print(group, "selected", selected_id, flush=True)
    (OUTPUT / "selection.json").write_text(json.dumps(selection, indent=2), encoding="utf-8")
    with (OUTPUT / "presentation_summary.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    (OUTPUT / "source_results_sha256.txt").write_text(
        hashlib.sha256((SOURCE / "resultados.json").read_bytes()).hexdigest(), encoding="ascii")
    print("FIGURES_DIRECTORY=" + str(OUTPUT.resolve()), flush=True)

if __name__ == "__main__":
    main()
