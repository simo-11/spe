# -*- coding: utf-8 -*-
"""

@author: simo nikula
"""
# %% common config
do_plots=False
plot_stress_vector=False
plt_pause=0.5
write_files=False
write_table_line=True
# %% close plot windows
import matplotlib.pyplot as plt
plt.close('all')
# %% princeton
import matplotlib.pyplot as plt
w=3.2024
h=12.377
gPlotDone=False
if not do_plots:
    plt.close('all')
for ec_in_h in (10,14,15,17,20,23,25,30,35,40,50,100): #
    ms=h/1000./10/(ec_in_h*ec_in_h)
    section=None
    runfile('primitive.py',#noqa
      args=f"""-A -W={w/1000} -H={h/1000}
      --mesh_size={ms}
      --primitive=rectangle""")
    section.log_write=True
    if not gPlotDone and do_plots:
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
        print("|{0}|{1}|{2:.3g}|{3:.3g}|\n".format(ec_in_h,
          section.mesh_nodes.shape[0],
          1e18*section.get_gamma(),
          1e12*section.get_j()))