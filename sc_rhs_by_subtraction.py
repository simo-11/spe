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
import sectionproperties.pre.library.steel_sections as steel_sections
import sectionproperties.pre.library.primitive_sections as primitive_sections
do_plots=False
plot_stress_vector=False
plt_pause=0.5
write_files=False
write_table_line=True
log_parts=True
# %% close plot windows
plt.close('all')
# %% configure 
import matplotlib.pyplot as plt
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
for ec_in_h in (4,10,14,15,30,35,40): #
    geometry = steel_sections.box_girder_section(
        d=d,b_t=b_t,b_b=b_b,t_ft=t_ft,t_fb=t_fb,t_w=t_w)
    ob_geometry=primitive_sections.rectangular_section(d=d,b=b_t)
    ib_geometry=primitive_sections.rectangular_section(d=d-t_ft-t_fb,
                                             b=b_t-2*t_w)
    ms=np.pow(d/ec_in_h,2)
    ga=[geometry,ob_geometry,ib_geometry]
    sa=[]
    for go in ga:
        go.create_mesh(mesh_sizes=[ms])
        section=simo.dev.DevSection(go)
        section.log_write=False
        args=types.SimpleNamespace()
        args.primitive=simo.dev.BOX
        args.width=b_t
        args.height=d
        args.thickness=t_ft
        args.web_thickness=t_w
        args.gen='gen'
        args.z_scale=0.5
        section.set_args(args)
        section.calculate_geometric_properties()
        section.calculate_warping_properties()
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
        print("|{0}|{1}|{2:.3g}|{3:.3g}|{4:.3g}|{5:.3g}|".format(ec_in_h,
          sa[0].mesh_nodes.shape[0],
          1e9*sa[0].get_gamma(),
          1e6*sa[0].get_j(),
          1e9*(sa[1].get_gamma()-sa[2].get_gamma()),
          1e6*(sa[1].get_j()-sa[2].get_j())))
    if log_parts:
        print("|{0}|{1}|{2:.3g}|{3:.3g}|{4:.3g}|{5:.3g}|".format(
          sa[1].mesh_nodes.shape[0],
          sa[2].mesh_nodes.shape[0],
          1e9*sa[1].get_gamma(),
          1e6*sa[1].get_j(),
          1e9*sa[2].get_gamma(),
          1e6*sa[2].get_j()))