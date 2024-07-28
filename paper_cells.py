# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 13:19:57 2024

This file is used to document and provide means to reproduce
run and also debug cells if needed.

@author: simo nikula
"""
# %% rectangles
import matplotlib.pyplot as plt
w=100
for h in (100,80,60,30,10,2):
    gPlotDone=False
    for ec_in_h in (10,):
        ms=h/1000./10/(ec_in_h*ec_in_h)
        section=None
        runfile('primitive.py',#noqa
          args=f"""-A -W={w/1000} -H={h/1000}
          --mesh_size={ms}
          --primitive=rectangle""")
        if not gPlotDone:
            pdf=plt.figure()
            fn=section.gfn(section.default_filename(".pdf","geometry"))
            section.geometry.plot_geometry(
                labels=()
                ,title=f"solid rectangle {w}x{h} mm"
                ,cp=False
                ,legend=False
                ,filename=fn)
            print(f'Wrote {fn}')
            gPlotDone=True
        (fig,ax)=section.contour_warping_values(levels=51,title='')
        fn=section.gfn(section.default_filename(".pdf","contour"))
        plt.savefig(fn)
        print(f'Wrote {fn}')
        plt.show();
        (fig,ax)=section.plot_warping_values(title='')
        fn=section.gfn(section.default_filename(".pdf","3d"))
        plt.savefig(fn,bbox_inches='tight')
        print(f'Wrote {fn}')
        plt.show();
        section.write_json()
        section.write_warping_csv()
        section.write_triangles_csv()
        section.write_warping_gltf()
        section.write_warping_ply()
# %% swxy
import matplotlib.pyplot as plt
import numpy as np
xv=[]
wv=[]
iv=[]
w=100
awc=pow(w/1000,3)/144
aic=w/1000/3
for h in np.geomspace(w/40,w,num=25):
    ec_in_h=10
    ms=1e-6*h*w/(ec_in_h*ec_in_h)
    runfile('primitive.py',#noqa
      args=f"""-A -W={w/1000} -H={h/1000}
      --mesh_size={ms}
      --primitive=rectangle""")
    xv.append(w/h)
    dc=pow((h/1000),3)
    av=awc*dc
    nv=section.get_gamma()
    wv.append(nv/av)
    av=aic*dc
    nv=section.get_j()
    iv.append(nv/av)    
fig, ax = plt.subplots()
ax.set(xlim=(1,max(xv)),ylim=(0,1))
#ax.set_title('Ratio of actual and thin wall theoretical values')
ax.set_xlabel(r'$l/t$')
ax.set_ylabel(r'Ratio')
ax.plot(xv,wv,'k',label=r'$I_{\omega}$')
ax.plot(xv,iv,'r',label=r'$I_{t}$')
legend = ax.legend(loc='lower center', shadow=True, fontsize='x-large')
legend.get_frame().set_facecolor('C0')
for i in (-1, -3):
    ax.text(xv[i], wv[i], f" ({xv[i]:2.2},{wv[i]:2.2})", ha="left")
fn='gen/swxy.pdf'
plt.savefig(fn)
print(f'Wrote {fn}')
plt.show()        
# %% U and SHS
import matplotlib.pyplot as plt
import time
for p in ("rhs",): # "rhs","u"
    match p:#noqa
        case "rhs":
            script="primitive"
            primitive=f"--primitive {p}"
            h=150
            w=h
            t=8
            ms=1e-4
        case "u":
            script="cold-formed-u"
            primitive=""
            h=100
            w=50
            t=4
            ms=1e-4
    for r in ("r",): # "s","r"
        if r=="s":
            n_r_s=(0,)
        else:
            n_r_s=range(8,)
        for n_r in n_r_s:
            ts=time.time()        
            runfile(f'{script}.py',#noqa
              args=f"""-A -W={w/1000} -H={h/1000} --thickness={t/1000}
              --mesh_size={ms}
              {primitive} --n_r={n_r}""")
            elapsed=time.time()-ts
            uc=section.default_filename("","solve")
            print(f'{uc} took {elapsed:.3f} seconds')
            pdf=plt.figure()
            fn=section.gfn(section.default_filename(".pdf","geometry"))
            section.geometry.plot_geometry(
                labels=()
                ,title=f"{p}-{r} {w}x{h}x{t} mm"
                ,cp=False
                ,legend=False
                ,filename=fn)
            print(f'Wrote {fn}')
            (fig,ax)=section.contour_warping_values(levels=51, title='')
            fn=section.gfn(section.default_filename(".pdf","contour"))
            plt.savefig(fn)
            print(f'Wrote {fn}')
            plt.show();
            (fig,ax)=section.plot_warping_values(title='')
            fn=section.gfn(section.default_filename(".pdf","3d"))
            plt.savefig(fn,bbox_inches='tight')
            print(f'Wrote {fn}')
            plt.show();
            section.write_json()
            section.write_warping_csv()
            section.write_triangles_csv()
            section.write_warping_gltf()
            section.write_warping_ply()
# %% shs-n_r
"""
Experiment on n_r option of rhs (points for describing rounding)
and mesh size.
"""
import matplotlib.pyplot as plt
import numpy as np
import time
import gc
p="rhs"
script="primitive"
primitive=f"--primitive {p}"
h=150
w=h
t=8
msa=["0.0001","0.00001"] # 1e-4,1e-5
mss=len(msa)
xv=list(range(2,9,1))+list(range(10,61,5))#range(2,9,1)+(range(10,61,5)
n_rs=len(xv)
wv=np.zeros((mss,n_rs))
st=np.zeros((mss,n_rs))
for msi in range(mss):
    nri=0
    for n_r in xv:
        gc.collect()
        ts=time.time()        
        runfile(f'{script}.py',#noqa
          args=f"""-A -W={w/1000} -H={h/1000} --thickness={t/1000}
          --mesh_size={msa[msi]}
          {primitive} --n_r={n_r}""")
        elapsed=time.time()-ts
        uc=section.default_filename("","solve")
        print(f'{uc} took {elapsed:.3f} seconds')
        section.plot_mesh(
            title=f'Mesh with mesh_size={msa[msi]} and n_r={n_r}',
                          materials=False)
        nv=section.get_gamma()
        wv[msi,nri]=nv  
        st[msi,nri]=elapsed
        # for table {tab:shs-values-rounded}
        print(f'''Section-Properties({2*t},{msa[msi]},{n_r},\
{section.num_nodes},{len(section.elements)})\
 & {nv*1e12:.1f} \\(10^{{-12}}\\)\\\\''')
        nri=nri+1
for pic in range(2):
    match pic:
        case 0:
            axv=xv[:8]
        case 1:
            axv=xv[5:]
    fig, ax = plt.subplots()
    ax.set_xticks(axv,axv)
    ax.set_xlabel(r'n_r')
    ax.set_ylabel(r'$I_{\omega} [m^6]$')
    for msi in range(mss):
        ms=msa[msi]
        match pic:
            case 0:
                awv=wv[msi,:8]
            case 1:
                awv=wv[msi,5:]
        ax.plot(axv,awv,label=f'ms={ms}')
    legend = ax.legend(loc='lower center', shadow=True, fontsize='x-large')
    fn=f'gen/shs-n_r-iw-{pic}.pdf'
    plt.savefig(fn)
    print(f'Wrote {fn}')
    plt.show()        
fig, ax = plt.subplots()
#ax.set(xlim=(1,max(xv)),ylim=(0,max(st)))
ax.set_xlabel(r'n_r')
ax.set_ylabel(r'solve time [s]')
for msi in range(mss):
    ms=msa[msi]
    ax.plot(xv,st[msi],label=f'ms={ms}')
legend = ax.legend(loc='upper left', shadow=True, fontsize='x-large')
fn='gen/shs-n_r-time.pdf'
plt.savefig(fn)
print(f'Wrote {fn}')
plt.show()        
# %% gbtul
"""
Experiment on n_r option with GBTUL.
"""
import matplotlib.pyplot as plt
import numpy as np
import time
import gc
for p in ("u",): # "rhs","u"
    match p:#noqa
        case "rhs":
            script="primitive"
            primitive=f"--primitive {p}"
            h=150
            w=h
            t=8
            xv=list(range(1,5,1))
        case "u":
            script="cold-formed-u"
            primitive=""
            h=100
            w=50
            t=4
            xv=list(range(1,8,1))
    n_rs=len(xv)
    wv=np.zeros(n_rs)
    st=np.zeros(n_rs)
    nri=0
    for n_r in xv:
        gc.collect()
        ts=time.time()      
        runfile(f'{script}.py',#noqa
          args=f"""--gbtul -W={w/1000} -H={h/1000} --thickness={t/1000}
          {primitive} --n_r={n_r}""") 
        elapsed=time.time()-ts
        if not section.gbtul:
            raise Exception("Gbtul failed")
        uc=section.default_filename("","gbtul")
        print(f'{uc} took {elapsed:.3f} seconds')
        nv=section.gbtul.gamma
        wv[nri]=nv  
        st[nri]=elapsed
        (fig,ax)=section.plot_gbt()
        ax.set_title(label=f"Iw from GBTUL={nv:.3}")
        plt.show()
        # for table {tab:shs-values-rounded}
        print(f'''GBTUL({n_r})\
 & {nv*1e12:.1f} \\(10^{{-12}}\\)\\\\''')
        nri=nri+1
    for pic in range(1):
        match pic:
            case 0:
                axv=xv[:8]
            case 1:
                axv=xv[5:]
        fig, ax = plt.subplots()
        ax.set_xticks(axv,axv)
        ax.set_xlabel(r'n_r')
        ax.set_ylabel(r'$I_{\omega} [m^6]$')
        match pic:
            case 0:
                awv=wv[:8]
            case 1:
                awv=wv[5:]
        ax.plot(axv,awv)
        fn=f'gen/{p}-n_r-iw-gbtul-{pic}.pdf'
        plt.savefig(fn)
        print(f'Wrote {fn}')
        plt.show()        
    fig, ax = plt.subplots()
    ax.set_xlabel(f'n_r for corner of {p}')
    ax.set_ylabel(r'solve time [s]')
    ax.plot(xv,st)
    fn=f'gen/{p}-n_r-gbtul-time.pdf'
    plt.savefig(fn)
    print(f'Wrote {fn}')
    plt.show()        
