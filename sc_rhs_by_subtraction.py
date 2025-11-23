# -*- coding: utf-8 -*-
"""
sharp corner RHS as box and by subtraction
@author: simo nikula
"""
# %% common config
import matplotlib.pyplot as plt
import numpy as np
import simo.dev
import types
import time
import sectionproperties.pre.library.steel_sections as steel_sections
import sectionproperties.pre.library.primitive_sections as primitive_sections
do_plots=False
plot_stress_vector=False
plt_pause=0.5
write_files=False
write_table_line=True
log_parts=True
#E=210e9#
E=1e9#Something soft, e.g. some plastic
nu=0.3
G=E/(2*(1+nu))
def kc(It,Iw):
    return np.sqrt((G*It)/(E*Iw))
def theta(T,It,Iw,L,x):
    if Iw==0.:
        return T/(G*It)*x
    k=kc(It,Iw)
    c0=T/(k*G*It)
    if np.tanh(k*L)==1:
        y=c0*(np.exp(-k*x)-1+k*x)
    else:
        y=c0*((np.tanh(k*L)*(np.cosh(k*x)-1))-np.sinh(k*x)+k*x)
    return y
# %% close plot windows
plt.close('all')
# %% configure 
W=160
H=80
T=10
d=H/1000
b_t=W/1000
b_b=W/1000
t_ft=T/1000
t_fb=T/1000
t_w=T/1000
gPlotDone=False
br_A=(W-T)/1000*(H-T)/1000
br_j=4*np.pow(br_A,2)/((2*(W-T)+2*(H-T))/T)
print(f"J from Bredt's formula {1e6*br_j:.5g} e-6")
# %% various ways
if not do_plots:
    plt.close('all')
for ec_in_h in range(10,40,5): #4,10,14,15,30,35,40
    geometry = steel_sections.box_girder_section(
        d=d,b_t=b_t,b_b=b_b,t_ft=t_ft,t_fb=t_fb,t_w=t_w)
    ob_geometry=primitive_sections.rectangular_section(d=d,b=b_t)
    ib_geometry=primitive_sections.rectangular_section(d=d-t_ft-t_fb,
                                             b=b_t-2*t_w)
    ms=np.pow(d/ec_in_h,2)
    ga=[geometry,]#ob_geometry,ib_geometry
    sa=[]
    for i in range(len(ga)):
        start=time.time()
        go=ga[i];
        go.create_mesh(mesh_sizes=[ms])
        section=simo.dev.DevSection(go)
        section.log_write=True
        args=types.SimpleNamespace()
        match i:
            case 0:
                args.primitive=simo.dev.BOX
                args.width=b_t
                args.height=d
            case 1:
                args.primitive=simo.dev.RECTANGLE
                args.width=b_t
                args.height=d
            case 2:
                args.primitive=simo.dev.RECTANGLE
                args.width=b_t-2*t_w
                args.height=d-t_ft-t_fb
        args.thickness=t_ft
        args.web_thickness=t_w
        args.gen='gen'
        args.z_scale=0.5
        section.set_args(args)
        section.calculate_geometric_properties()
        section.calculate_warping_properties()
        end=time.time()
        sa.append(section)
        if do_plots: 
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
        if write_files:
            section.write_json()
            section.write_warping_csv()
            section.write_triangles_csv()
            section.write_warping_gltf()    
    if write_table_line:
        if len(sa)>1:
            #print("|{0}|{1}|{2:.3g}|{3:.3g}|{4:.3g}|{5:.3g}|".#
            print("{0}&{1}&{2:.3g}&{3:.3g}&{4:.3g}&{5:.3g} \\\\".
                  format(ec_in_h,
              sa[0].mesh_nodes.shape[0],
              1e9*sa[0].get_gamma(),
              1e6*sa[0].get_j(),
              1e9*(sa[1].get_gamma()-sa[2].get_gamma()),
              1e6*(sa[1].get_j()-sa[2].get_j())))
        else:
            #print("|{0}|{1}|{2:.4g}|{3:.4g}|{4:.4g}|".
            print("{0} & {1} & {2:.4g} & {3:.4g} & {4:.2g} \\\\".
                  format(ec_in_h,
              sa[0].mesh_nodes.shape[0],
              1e6*sa[0].get_j(),
              1e9*sa[0].get_gamma(),
              (end-start)
              ))
    if log_parts and len(sa)>1:
        print("|{0}|{1}|{2:.3g}|{3:.3g}|{4:.3g}|{5:.3g}|".format(
          sa[1].mesh_nodes.shape[0],
          sa[2].mesh_nodes.shape[0],
          1e9*sa[1].get_gamma(),
          1e6*sa[1].get_j(),
          1e9*sa[2].get_gamma(),
          1e6*sa[2].get_j()))
# %% selected results
L=2
T=1e3
It=sa[0].get_j()
Iw=sa[0].get_gamma()
k=kc(It,Iw)
print(f"k={k:.3g}")
theta_v=theta(T,br_j,0,L,L)
print(f"It={br_j:.4}, theta={theta_v*180/np.pi:.4g}°")
theta_v=theta(T,It,0,L,L)
print(f"It={It:.4}, Iw=0, theta={theta_v*180/np.pi:.4g}°")
theta_v=theta(T,It,Iw,L,L)
print(f"It={It:.4}, Iw={Iw:.4g}, theta={theta_v*180/np.pi:.4g}°")
It=10.17e-6
Iw=1.388e-9
theta_v=theta(T,It,Iw,L,L)
print(f"It={It:.4}, Iw={Iw:.4g}, theta={theta_v*180/np.pi:.4g}°")
# %% bilinear
T=1e3
def beam_theory_to_latex_table(title,It,Iw):
    theta_v=theta(T,It,Iw,L,L)
    print(f"Warping torsion {title} & {It*1e6:.2f} & "
          f"{Iw*1e9:.3f} & {theta_v*180/np.pi:.4g} \\\\"
          )
beam_theory_to_latex_table("$120 \\times 60$",
                           10.724e-6,2.5394e-9)
beam_theory_to_latex_table("$280 \\times 140$",
                           10.603e-6,2.4938e-9)