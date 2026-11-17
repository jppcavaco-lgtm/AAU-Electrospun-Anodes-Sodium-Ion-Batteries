"""Portable copy of the five-band analysis used for the thesis figures.
Run from the package root with: python scripts/raman_integracao.py
Fit parameters are exploratory outputs; see Raman_Scripts_Documentation.txt.
Five-band fitting: Sadezky et al. (2005), doi:10.1016/j.carbon.2005.02.018.
Numerical bounds/windows are documented analysis choices, not paper prescriptions.
Never modifies original spectra. Creates a NEW output directory on every run.
"""
import argparse, csv, hashlib, json, platform
from datetime import datetime
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.optimize import least_squares
from scipy.signal import savgol_filter

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
BANDS=['D4','D1','D3','G','D2']
COLORS=['#8c6bb1','#d95f02','#1b9e77','#2774ae','#cc79a7']
SETTINGS=dict(baseline_windows=[[700,900],[1750,1800]],fit_window=[1000,1800],
              smooth_before_fit=False,multistart=16,random_seed=20261029,
              peak_detection_savgol=dict(window=11,order=3),
              centre_bounds=[[1100,1270],[1300,1400],[1450,1550],[1550,1605],[1605,1655]],
              fwhm_bounds=[[40,450],[40,350],[60,400],[15,180],[10,180]])

def sample_list(manifest_path):
    """Resolve the packaged two-column spectra listed in the manifest."""
    samples = json.loads(manifest_path.read_text(encoding="utf-8"))
    seen = set()
    for sample in samples:
        if sample["id"] in seen:
            raise ValueError("Duplicate spectrum ID: " + sample["id"])
        seen.add(sample["id"])
        source = Path(sample["source"])
        if not source.is_absolute():
            source = manifest_path.parent / source
        if not source.is_file():
            raise FileNotFoundError(source)
        sample["source"] = str(source.resolve())
    return samples


def lorentzian(x,height,centre,fwhm):
    return height/(1+4*((x-centre)/fwhm)**2)

def gaussian(x,height,centre,fwhm):
    return height*np.exp(-4*np.log(2)*((x-centre)/fwhm)**2)

def components(x,p):
    return np.array([(gaussian if name=='D3' else lorentzian)(x,*p[3*i:3*i+3])
                     for i,name in enumerate(BANDS)])

def model(x,p):
    return components(x,p).sum(axis=0)+p[-1]

def areas(p):
    # Analytic areas of complete fitted band functions (no background).
    return np.array([p[3*i]*p[3*i+2]*(np.sqrt(np.pi/(4*np.log(2))) if name=='D3' else np.pi/2)
                     for i,name in enumerate(BANDS)])

def ratios(p):
    a=areas(p)
    return np.array([p[3]/p[9],a[1]/a[3],a[2]/a.sum()])

def observed_maximum(x,y,lo,hi):
    # Filtering only locates a noise-robust observed envelope maximum; never fits filtered data.
    sm=savgol_filter(y,11,3)
    idx=np.flatnonzero((x>=lo)&(x<=hi)); k=idx[np.argmax(sm[idx])]
    return float(x[k])

def analyze_and_plot(sample):
    path=Path(sample['source']); a=np.loadtxt(path)
    if a.ndim!=2 or a.shape[1]!=2 or not np.isfinite(a).all():raise ValueError('Invalid spectrum '+str(path))
    a=a[np.argsort(a[:,0])];x,raw=a.T
    if np.any(np.diff(x)<=0):raise ValueError('Repeated Raman shifts')
    mask=np.zeros(len(x),dtype=bool)
    for lo,hi in SETTINGS['baseline_windows']:mask|=(x>=lo)&(x<=hi)
    coeff=np.polyfit(x[mask],raw[mask],1);baseline=np.polyval(coeff,x);corrected=raw-baseline
    sel=(x>=1000)&(x<=1800);xf,yf=x[sel],corrected[sel]
    scale=float(max(yf));yn=yf/scale
    lower=[];upper=[]
    for c,w in zip(SETTINGS['centre_bounds'],SETTINGS['fwhm_bounds']):
        lower += [0,c[0],w[0]];upper += [3,c[1],w[1]]
    lower=np.array(lower+[-.5]);upper=np.array(upper+[.5])
    rng=np.random.default_rng(SETTINGS['random_seed']);fits=[]
    for attempt in range(SETTINGS['multistart']):
        p0=[]
        for h,c,w in zip([.12,.75,.25,.65,.20],[1180,1350,1500,1580,1620],[170,155,170,70,65]):
            p0 += [h,c,w] if attempt==0 else [rng.uniform(.04,.85),c+rng.normal(0,12),w*rng.uniform(.5,1.7)]
        p0=np.clip(p0+[0],lower+1e-6,upper-1e-6)
        f=least_squares(lambda p:model(xf,p)-yn,p0,bounds=(lower,upper),x_scale='jac',
                        max_nfev=2500,ftol=1e-9,xtol=1e-9,gtol=1e-9)
        if f.success and np.isfinite(f.fun).all():fits.append(f)
    if not fits:raise RuntimeError('No fit converged for '+sample['id'])
    fits.sort(key=lambda f:float(f.fun@f.fun));f=fits[0];rss=float(f.fun@f.fun)
    p=f.x.copy();amp=list(range(0,15,3))+[15];p[amp]*=scale
    covn=np.linalg.pinv(f.jac.T@f.jac)*rss/(len(xf)-len(p))
    transform=np.ones(len(p));transform[amp]=scale;cov=covn*np.outer(transform,transform)
    pe=np.sqrt(np.maximum(0,np.diag(cov)))
    rj=np.zeros((3,len(p)))
    for j in range(len(p)):
        step=1e-5*max(1,abs(p[j]));d=np.zeros(len(p));d[j]=step
        rj[:,j]=(ratios(p+d)-ratios(p-d))/(2*step)
    re=np.sqrt(np.maximum(0,np.diag(rj@cov@rj.T)))
    near=[q for q in fits if q.fun@q.fun<=rss*1.01+1e-14]
    nr=np.array([ratios(q.x) for q in near]);warnings=[];bounded=[]
    for i,name in enumerate(BANDS):
        for off,field in [(1,'centre'),(2,'FWHM')]:
            j=3*i+off;frac=(f.x[j]-lower[j])/(upper[j]-lower[j])
            if frac<.01 or frac>.99:bounded.append(name+'_'+field)
        if f.x[3*i]<.005:bounded.append(name+'_weak')
    if bounded:warnings.append('bound_or_weak_component:'+','.join(bounded))
    singular=np.linalg.svd(f.jac,compute_uv=False);condition=float(singular[0]/max(singular[-1],1e-15))
    if condition>1e6:warnings.append('ill_conditioned_fit')
    if np.ptp(nr[:,1])>.1*ratios(p)[1]:warnings.append('area_ratio_multistart_sensitive')
    residual=yf-model(xf,p);rat=ratios(p);aa=areas(p)
    bands=[dict(band=name,shape='Gaussian' if name=='D3' else 'Lorentzian',
       centre_cm1=float(p[3*i+1]),centre_se=float(pe[3*i+1]),height=float(p[3*i]),
       height_se=float(pe[3*i]),fwhm_cm1=float(p[3*i+2]),fwhm_se=float(pe[3*i+2]),area=float(aa[i]))
       for i,name in enumerate(BANDS)]
    result=dict(**sample,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
       baseline_slope=float(coeff[0]),baseline_intercept=float(coeff[1]),
       observed_D_cm1=observed_maximum(x,corrected,1250,1420),observed_G_cm1=observed_maximum(x,corrected,1520,1660),
       height_ratio=float(rat[0]),area_ratio=float(rat[1]),d3_area_fraction=float(rat[2]),
       height_ratio_se=float(re[0]),area_ratio_se=float(re[1]),d3_area_fraction_se=float(re[2]),
       rmse=float(np.sqrt(np.mean(residual**2))),r_squared=float(1-np.sum(residual**2)/np.sum((yf-yf.mean())**2)),
       jacobian_condition=condition,successful_starts=len(fits),near_optimal_starts=len(near),
       near_optimal_area_ratio_range=[float(min(nr[:,1])),float(max(nr[:,1]))],warnings=warnings,bands=bands)
    arrays=dict(x=x,raw=raw,baseline=baseline,corrected=corrected,fit_x=xf,fit_y=yf,
                fitted=model(xf,p),components=components(xf,p),residual=residual,parameters=p,covariance=cov)
    return result,arrays

def configure_plots():
    plt.rcParams.update({'font.size':11,'axes.labelsize':12,'axes.spines.top':False,
                        'axes.spines.right':False,'savefig.dpi':240,'legend.fontsize':9})

def plot_group(group,rs,arrays,out):
    n=len(rs);fig,axes=plt.subplots(2,n,figsize=(5.3*n,6.1),squeeze=False,layout='constrained',
                                   gridspec_kw={'height_ratios':[3.2,1]})
    for j,r in enumerate(rs):
        a=arrays[r['id']];norm=max(a['fit_y']);ax,res=axes[:,j]
        ax.plot(a['fit_x'],a['fit_y']/norm,color='.40',lw=.85,label='Measured, baseline corrected')
        ax.plot(a['fit_x'],a['fitted']/norm,color='black',lw=1.4,label='Five-band fit')
        for i,(name,color) in enumerate(zip(BANDS,COLORS)):
            curve=a['components'][i]/norm;ax.plot(a['fit_x'],curve,color=color,lw=1.2,label=name)
            ax.fill_between(a['fit_x'],0,curve,color=color,alpha=.12)
        ax.set_title(f"({chr(97+j)}) {r['filename']}");ax.set_xlim(1000,1800);ax.set_ylim(bottom=-.12)
        ax.set_ylabel('Normalized intensity (a.u.)');ax.legend(loc='upper left',ncol=2,frameon=False)
        res.plot(a['fit_x'],a['residual']/norm,color='.25',lw=.75);res.axhline(0,color='.6',lw=.7)
        res.set_xlim(1000,1800);res.set_ylabel('Residual');res.set_xlabel(r'Raman shift (cm$^{-1}$)')
    fig.savefig(out/(group+'_integracao.png'));fig.savefig(out/(group+'_integracao.svg'));plt.close(fig)

def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8-sig') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path)
    parser.add_argument('--manifest',type=Path,default=PACKAGE_ROOT/'samples.json')
    parser.add_argument('--group');args=parser.parse_args()
    out=args.output or PACKAGE_ROOT/'results'/('resultados_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
    out.mkdir(parents=True,exist_ok=False);configure_plots();results=[];arrays={}
    selected=[s for s in sample_list(args.manifest) if not args.group or s['group']==args.group]
    if not selected:raise ValueError('No matching sample group')
    for s in selected:
        print('Fitting '+s['id']+' '+s['filename'],flush=True)
        r,a=analyze_and_plot(s);results.append(r);arrays[r['id']]=a
        np.savez_compressed(out/(r['id']+'.npz'),**a)
        print(f"  H_D1/H_G={r['height_ratio']:.4f}; A_D1/A_G={r['area_ratio']:.4f}; R2={r['r_squared']:.5f}; {r['warnings']}",flush=True)
    metadata=dict(settings=SETTINGS,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                  matplotlib=matplotlib.__version__,results=results)
    (out/'resultados.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    flat=[]
    for r in results:
        row={k:v for k,v in r.items() if k not in ('bands','warnings','near_optimal_area_ratio_range')}
        row['warnings']='; '.join(r['warnings']);row['near_optimal_area_ratio_min'],row['near_optimal_area_ratio_max']=r['near_optimal_area_ratio_range']
        for band in r['bands']:
            for k in ['centre_cm1','height','fwhm_cm1','area']:row[band['band']+'_'+k]=band[k]
        flat.append(row)
    write_csv(out/'resumo_raman.csv',flat)
    write_csv(out/'parametros_bandas.csv',[dict(spectrum_id=r['id'],filename=r['filename'],**b) for r in results for b in r['bands']])
    for g in dict.fromkeys(r['group'] for r in results):plot_group(g,[r for r in results if r['group']==g],arrays,out)
    print('RESULTS_DIRECTORY='+str(out.resolve()),flush=True)

if __name__=='__main__':main()
