"""Adapted COPY of Baseline Correction + vertical lines.py (see originais/).
Use --results <folder produced by raman_integracao.py>.
Lines mark OBSERVED envelope maxima, not fixed literature positions.
The displayed curve is unsmoothed, baseline-corrected experimental data.
"""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from raman_integracao import configure_plots

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',required=True,type=Path);args=parser.parse_args()
    results=json.loads((args.results/'resultados.json').read_text(encoding='utf-8'))['results']
    configure_plots()
    for group in dict.fromkeys(r['group'] for r in results):
        spectra=[r for r in results if r['group']==group];dest=args.results/(group+'_picos.png')
        if dest.exists():raise FileExistsError('Refusing to overwrite '+str(dest))
        fig,axes=plt.subplots(1,len(spectra),figsize=(5.3*len(spectra),3.9),squeeze=False,layout='constrained')
        for j,r in enumerate(spectra):
            a=np.load(args.results/(r['id']+'.npz'));x=a['x'];y=a['corrected'];sel=(x>=1000)&(x<=1800)
            ax=axes[0,j];ax.plot(x[sel],y[sel]/max(y[sel]),color='black',lw=1)
            for key,name,color,level in [('observed_D_cm1','D envelope','#d95f02',1.14),('observed_G_cm1','G envelope','#2774ae',1.03)]:
                pos=r[key];ax.axvline(pos,color=color,ls='--',lw=1)
                ax.text(pos-8,level,f'{name}\n{pos:.1f}'+r' cm$^{-1}$',ha='right',va='top',color=color,fontsize=10)
            ax.set_title(f"({chr(97+j)}) {r['filename']}");ax.set_xlabel(r'Raman shift (cm$^{-1}$)')
            ax.set_ylabel('Normalized corrected intensity (a.u.)');ax.set_xlim(1000,1800);ax.set_ylim(-.10,1.24)
        fig.savefig(dest);fig.savefig(args.results/(group+'_picos.svg'));plt.close(fig)
        print(dest,flush=True)

if __name__=='__main__':main()
