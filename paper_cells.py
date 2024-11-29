# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 13:19:57 2024

This file is used to document and provide means to reproduce
run and also debug cells if needed.

@author: simo nikula
"""
# %% common config
plot_stress_vector=True
plt_pause=0.5
# %% close plot windows
import matplotlib.pyplot as plt
plt.close('all')
# %% rectangles
import matplotlib.pyplot as plt
w=100
for h in (100,80,60,30,10): #100,80,60,30,10
    gPlotDone=False
    for ec_in_h in (10,): #5,10,20
        ms=h/1000./10/(ec_in_h*ec_in_h)
        section=None
        runfile('primitive.py',#noqa
          args=f"""-A -W={w/1000} -H={h/1000}
          --mesh_size={ms}
          --primitive=rectangle""")
        section.log_write=True
        if not gPlotDone:
            ax=section.geometry.plot_geometry(
                labels=()
                ,title=f"solid rectangle {w}x{h} mm"
                ,num="geometry",clear=True
                ,cp=False
                ,legend=False)
            ax.set_xlabel("x [m]")
            ax.set_ylabel("y [m]")
            fn=section.gfn(section.default_filename(".pdf","geometry"))
            plt.savefig(fn)
            if section.log_write:
                print(f'Wrote {fn}')
            gPlotDone=True
        (fig,ax)=section.contour_warping_values(levels=51,title='',
                label='Warping [m$^2$]',
                num='contour',clear=True)
        fn=section.gfn(section.default_filename(".pdf","contour"))
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        plt.savefig(fn,bbox_inches='tight')
        if section.log_write:
            print(f'Wrote {fn}')
        plt.pause(plt_pause)
        (fig,ax)=section.plot_warping_values(title='',
                num='contour3d',clear=True)
        ax.set_xlabel("\nx [m]",linespacing=2.6)
        ax.set_ylabel("y [m]")
        ax.set_zlabel("\nwarping [m$^2$]",linespacing=1.6)
        fn=section.gfn(section.default_filename(".pdf","3d"))
        plt.savefig(fn,bbox_inches='tight',pad_inches=0.3)
        if section.log_write:
            print(f'Wrote {fn}')
        plt.pause(plt_pause);
        if plot_stress_vector:
            stress = section.calculate_stress(mzz=0.001)
            sv_ax=stress.plot_stress_vector(stress="mzz_zxy", 
                                      cmap="viridis",
                                      fmt="{x:.3f}",
                                      normalize=False,
                                      num='stress_vector',
                                      clear=True
                                      )
            sv_ax.set_xlabel("x [m]")
            sv_ax.set_ylabel("y [m]")
            fn=section.gfn(section.default_filename(".pdf","stress-vector"))
            plt.savefig(fn,bbox_inches='tight')
            if section.log_write:
                print(f'Wrote {fn}')
        plt.pause(plt_pause)
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
fig, ax = plt.subplots(num='xy',clear=True)
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
plt.pause(plt_pause)        
# %% box_girder
import matplotlib.pyplot as plt
import numpy as np
import time
import math
import sympy
import concurrent.futures
import sectionproperties.pre.library.steel_sections as steel_sections
import simo.dev
import types
b=150
h0=150
t0=8
d=t0
A=b*h0-(b-2*d)*(h0-2*t0)
# Keep A and d constant
print(f'A={1e-6*A:.3g} d={1e-3*d:.1g}')
w=b
plot_it=False
plot_geometry=False
print_each_point=False
x_ticks=50
w_s=sympy.symbols('w')
def save_plot(fig,ax,pdf_name):
    fig.show()
    fn=f'gen/{pdf_name}.pdf'
    try:
        fig.savefig(fn)
        print(f'Wrote {fn}')
    except PermissionError as e:
        print(f'Write of {fn} failed due to {e}')
fig_iw, ax_iw = plt.subplots(num='Iw',clear=True)
ax_iw.set_yscale("log", base=10)
if plot_geometry:
    fig_g, ax_g = plt.subplots(num='Geometry',clear=True)
ax_iw.set_xlabel(r'$h/w$')
ax_iw.set_ylabel(r'$I_w [m^6]$')
if plot_it:
    fig_it, ax_it = plt.subplots(num='It',clear=True)
    ax_it.set_xlabel(r'$h/w$')
    ax_it.set_ylabel(r'$I_t$')
def run_solve(t,h,w,index):
    geometry = steel_sections.box_girder_section(h/1000,
        w/1000,w/1000,
        t/1000,t/1000,d/1000)
    ms=1e-6*h*w/100
    geometry.create_mesh(mesh_sizes=[ms])
    section=simo.dev.DevSection(geometry)
    section.log_write=False
    args=types.SimpleNamespace()
    args.primitive=simo.dev.BOX
    args.width=round(w/1000,3)
    args.height=round(h/1000,3)
    args.thickness=round(t/1000,3)
    args.web_thickness=round(d/1000,3)
    args.gen='gen'
    args.z_scale=0.5
    section.set_args(args)
    section.calculate_geometric_properties()
    section.calculate_warping_properties()
    section.write_json()
    section.write_warping_gltf()
    return (t,h,w,index,section)
if not print_each_point:
    print("d/t=", end="")
# waiting for better concurrency
with concurrent.futures.ThreadPoolExecutor() as executor:
    for t in np.linspace(0.75*t0,1.5*t0,num=8):
        fs=[]
        xv=np.zeros(x_ticks)
        wv=np.zeros(x_ticks)
        iv=np.zeros(x_ticks)
        index=0
        d_over_t=d/t
        if d_over_t>0.95:
            label=rf'$\frac{{d}}{{t}}={d_over_t:.3g}$'
        else:
            label=rf'$\frac{{d}}{{t}}={d_over_t:.2g}$'
        h1=math.sqrt(d_over_t)*h0
        dh=h1/x_ticks
        if print_each_point:
            print(f"\nd/t={d_over_t:.4g}, t={t:.4g}, h1={h1:.4g}, dh={dh:.4g}")
        else:
            print(f"{d_over_t:.4g}", end=" ")
        hi=0    
        def getW(h):
            sol=sympy.solve(A-w_s*h+(w_s-2*d)*(h-2*t),w_s)
            return float(sol[0].evalf())
        for h in np.concatenate(
                (np.linspace(0.5*h1,h1,num=x_ticks//2),
                 np.linspace(h1+dh,1.2*h1,num=x_ticks//2))):
            hi=hi+1
            h=int(round(h,0))
            w=round(getW(h),0)
            if print_each_point:
                print(f"h={h}, w={w}, h/w={h/w:.4g}, hi={hi}")
            if w<2*d:
                print(f'loop ended as w={w:.3g} < 2d={2*d:.3g}')
                break
            fs.append(executor.submit(run_solve,t,h,w,index))
            index=index+1
        for future in concurrent.futures.as_completed(fs):
            (t,h,w,index,section)=future.result()
            xv[index]=h/w
            wv[index]=section.get_gamma()
            iv[index]=section.get_j()
            if print_each_point:
                print(f"Iw(t={t:.3g},h={h:.3g},w={w:.3g})={wv[index]:.3g}")
            if plot_geometry:
                section.geometry.plot_geometry(num='Geometry',
                labels=(),
                cp=False,
                legend=False,
                title=f'{section.default_filename("","geometry")}',
                clear=True)
                plt.pause(plt_pause)
        ax_iw.plot(xv,wv,label=label)
        ax_iw.legend(loc='upper right', shadow=True, fontsize='x-large')    
        fig_iw.show()
        plt.pause(plt_pause)
        if plot_it:
            ax_it.plot(xv,iv,label=label)
            ax_it.legend(loc='upper right', shadow=True, fontsize='x-large')    
            plt.pause(plt_pause)    
print()
save_plot(fig_iw,ax_iw,'girder_iw')
if plot_it:
    save_plot(fig_it,ax_it,'girder_it')
# %% U, SHS and RHS
import matplotlib.pyplot as plt
import time#noqa
for p in ("u","shs","rhs"): # "shs","u","rhs"
    match p:#noqa
        case "rhs":
            script="primitive"
            primitive="--primitive rhs"
            h=150
            w=120
            t=8
            ms=1e-4
        case "shs":
            script="primitive"
            primitive="--primitive rhs"
            h=150
            w=h
            t=8
            ms=1e-5 # 1e-4 -> 264/792, 1e-5 -> 1600/1716
        case "u":
            script="cold-formed-u"
            primitive=""
            h=100
            w=50
            t=4
            ms=1e-5# 1e-5 -> 383/495
        case _:
            raise RuntimeError(f"""profile {p} is not supported.
Check spelling or add support""")
    for r in ("s","r"): # "s","r"
        if r=="s":
            n_r_s=(0,)
        else:
            n_r_s=range(8,9)
        for n_r in n_r_s:
            ts=time.time()        
            runfile(f'{script}.py',#noqa
              args=f"""-A -W={w/1000} -H={h/1000} --thickness={t/1000}
              --mesh_size={ms}
              {primitive} --n_r={n_r}""")
            elapsed=time.time()-ts
            section.log_write=True
            uc=section.default_filename("","solve")
            print(f'{uc} took {elapsed:.3f} seconds')
            fn=section.gfn(section.default_filename(".pdf","geometry"))        
            ax=section.geometry.plot_geometry(num='geometry',clear=True,
                labels=()
                ,title=f"{p}-{r} {w}x{h}x{t} mm"
                ,cp=False
                ,legend=False)
            ax.set_xlabel("x [m]")
            ax.set_ylabel("y [m]")
            plt.savefig(fn,bbox_inches='tight')
            if section.log_write:
                print(f'Wrote {fn}')
            (fig,ax)=section.contour_warping_values(
                label='Warping [m$^2$]',
                num='contour',clear=True,levels=51, title='')
            fn=section.gfn(section.default_filename(".pdf","contour"))
            ax.set_xlabel("x [m]")
            ax.set_ylabel("y [m]")
            plt.savefig(fn,bbox_inches='tight')
            if section.log_write:
                print(f'Wrote {fn}')
            plt.pause(plt_pause);
            (fig,ax)=section.plot_warping_values(
                num='contour3d',clear=True,title='')
            fn=section.gfn(section.default_filename(".pdf","3d"))
            ax.set_xlabel("\nx [m]",linespacing=1)
            ax.set_ylabel("\ny [m]",linespacing=2.6)
            ax.set_zlabel("\nwarping [m$^2$]",linespacing=1.6)
            plt.savefig(fn,bbox_inches='tight',pad_inches=0.3)
            if section.log_write:
                print(f'Wrote {fn}')
            if plot_stress_vector:
                stress = section.calculate_stress(mzz=0.001)
                sv_ax=stress.plot_stress_vector(stress="mzz_zxy", 
                                          cmap="viridis",
                                          fmt="{x:.3f}",
                                          normalize=False,
                                          num='stress_vector',
                                          clear=True
                                          )
                sv_ax.set_xlabel("x [m]")
                sv_ax.set_ylabel("y [m]")
                fn=section.gfn(
                    section.default_filename(".pdf","stress-vector"))
                plt.savefig(fn,bbox_inches='tight')
                if section.log_write:
                    print(f'Wrote {fn}')                
            plt.pause(plt_pause)
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
for num in ['geometry','contour','contour3d','stress_vector']:
    plt.close(num)
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
                          materials=False
                          ,num='mesh',clear=True)
        plt.pause(plt_pause)
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
    fig, ax = plt.subplots(num='n_r/Iw',clear=True)
    ax.set_xticks(axv,axv)
    ax.set_xlabel(r'n_r')
    ax.set_ylabel(r'$I_{w} [m^6]$')
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
    if section.log_write:
        print(f'Wrote {fn}')
    plt.pause(plt_pause)
fig, ax = plt.subplots(num='solveTime',clear=True)
#ax.set(xlim=(1,max(xv)),ylim=(0,max(st)))
ax.set_xlabel(r'n_r')
ax.set_ylabel(r'solve time [s]')
for msi in range(mss):
    ms=msa[msi]
    ax.plot(xv,st[msi],label=f'ms={ms}')
legend = ax.legend(loc='upper left', shadow=True, fontsize='x-large')
fn='gen/shs-n_r-time.pdf'
plt.savefig(fn)
if section.log_write:
    print(f'Wrote {fn}')
plt.pause(plt_pause)  
# %% gbtul
"""
Experiment on n_r option with GBTUL.
"""
import matplotlib.pyplot as plt
import numpy as np
import time
import os.path
for dir in (".","..","../.."):
    gbt_path=dir+"/gbtul/GBT/GBT.exe"
    if os.path.exists(gbt_path):
        break
if not os.path.exists(gbt_path):
    raise Exception("GBT not found")
for p in ("shs","u","rhs"): # "shs","u","rhs"
    match p:#noqa
        case "shs":
            script="primitive"
            primitive="--primitive rhs"
            h=150
            w=h
            t=8
            xv=list(range(1,5,1))
        case "rhs":
            script="primitive"
            primitive="--primitive rhs"
            h=150
            w=120
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
        ts=time.time()      
        runfile(f'{script}.py',#noqa
          args=f"""--gbtul -W={w/1000} -H={h/1000} --thickness={t/1000}
          --gbt {gbt_path}
          {primitive} --n_r={n_r}""") 
        elapsed=time.time()-ts
        if not section.gbtul:
            raise Exception("Gbtul failed")
        uc=section.default_filename("","gbtul")
        print(f'{uc} took {elapsed:.3f} seconds')
        nv=section.gbtul.gamma
        wv[nri]=nv  
        st[nri]=elapsed
        (fig,ax)=section.plot_gbt(num="geometry",clear=True)
        ax.set_title(label=f"Iw from GBTUL={nv:.3}")
        plt.pause(plt_pause)
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
        fig, ax = plt.subplots(num="n_r/Iw",clear=True)
        ax.set_xticks(axv,axv)
        ax.set_xlabel(r'n_r')
        ax.set_ylabel(r'$I_{w} [m^6]$')
        match pic:
            case 0:
                awv=wv[:8]
            case 1:
                awv=wv[5:]
        ax.plot(axv,awv)
        fn=f'gen/{p}-n_r-iw-gbtul-{pic}.pdf'
        plt.savefig(fn)
        if section.log_write:
            print(f'Wrote {fn}')
        plt.pause(plt_pause)        
    fig, ax = plt.subplots(num="solveTime",clear=True)
    ax.set_xlabel(f'n_r for corner of {p}')
    ax.set_ylabel(r'solve time [s]')
    ax.plot(xv,st)
    fn=f'gen/{p}-n_r-gbtul-time.pdf'
    plt.savefig(fn)
    if section.log_write:
        print(f'Wrote {fn}')
    plt.pause(plt_pause)  
