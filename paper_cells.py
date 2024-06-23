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
for p in ("rhs",):
    for r in ("s","r"):
        if r=="s":
            n_r=0
        else:
            n_r=8
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
        runfile(f'{script}.py',#noqa
          args=f"""-A -W={w/1000} -H={h/1000} --thickness={t/1000}
          --mesh_size={ms}
          {primitive} --n_r={n_r}""")
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

