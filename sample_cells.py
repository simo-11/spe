# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 19:42:59 2023

This file is used to speedup testing of various fragments of code
known as cells e.g. in Spyder which provides means to run and debug cells.
Matlab has quite similar feature for matlab-scripts but uses
%% as cell separator


@author: snikula
"""
# %% rectangle
runfile('primitive.py',#noqa
  args="""-A -W 0.085 -H 0.118 --mesh_size=0.00003 --primitive=rectangle""")
# %% rectangle 1000-1000
runfile('primitive.py',#noqa
  args="""-A -W=1 -H=1 --mesh_size=0.003 --primitive=rectangle""")
# %% rectangle 100-100-0.0003
runfile('primitive.py',#noqa
  args="""-A -W=0.1 -H=0.1 --mesh_size=3e-4 --primitive=rectangle""")
# %% rectangle 100-100-0.0001
runfile('primitive.py',#noqa
  args="""-A -W=0.1 -H=0.1 --mesh_size=1e-4 --primitive=rectangle""")
# %% rectangle 100-100-0.00003
runfile('primitive.py',#noqa
  args="""-A -W=0.1 -H=0.1 --mesh_size=3e-5 --primitive=rectangle""")
# %% rectangle 100-10
runfile('primitive.py',#noqa
  args="""-A -W=0.1 -H=0.01 --mesh_size=1e-5 --primitive=rectangle""")
# %% rectangle 100-4
runfile('primitive.py',#noqa
  args="""-A -W=0.004 -H=0.1 --mesh_size=1e-5 --primitive=rectangle""")
# %% circle
runfile('primitive.py',#noqa
  args="""-A --diameter 1 --mesh_size=0.05 --primitive=circular""")
# %% rhs
runfile('primitive.py',#noqa
  args="""-A -W 1 -H 1 --thickness=0.03 --mesh_size=0.05 --primitive=rhs
  --radius=0 --n_r=0""")
# %% Timo's rhs
# On the distortion and warping of cantilever beams with hollow section
runfile('primitive.py',#noqa
  args="""-A -W 0.1 -H 0.1 --thickness=0.004 --mesh_size=1e-5
  --primitive=rhs --n_r=0""")
# %% Timo's rhs with rounder corners
runfile('primitive.py',#noqa
  args="""-A -W 0.1 -H 0.1 --thickness=0.004 --mesh_size=1e-6
  --primitive=rhs --n_r=8""")
# %% rsh 200x100x4
runfile('primitive.py',#noqa
  args="""-A -W 0.1 -H 0.2 --thickness=0.004 --mesh_size=1e-6
  --primitive=rhs --n_r=8""")
# %% chs
runfile('primitive.py',#noqa
  args="""-A --diameter 1 --thickness 0.03 --mesh_size=0.2 --primitive=chs
   --count=105""")
# %% cold-formed-u
runfile('cold-formed-u.py',#noqa
  args="""-A -W 0.05 -H 0.1 --thickness=0.004 --n_r=0 --mesh_size=0.05""")
# %% plot_geometry
section.geometry.plot_geometry()#noqa
# %% plot_mesh
import matplotlib.pyplot as plt
s=section#noqa
w=s.args.width
h=s.args.height
alpha=0.5
mask=None
x=s._mesh_nodes[:,0]
y=s._mesh_nodes[:,1]
t=s._mesh_elements[:, 0:3]
fig, ax = plt.subplots()
s.set_box_aspect(ax)
ax.triplot(x,y,t,lw=0.5,color="black",alpha=alpha,mask=mask)
axins = ax.inset_axes([0.1, 0.1, 0.7, 0.7])
s.set_box_aspect(axins)
axins.triplot(x,y,t,lw=0.5,color="black",alpha=alpha,mask=mask)
# subregion of the original image
if min(x)<0 and min(y)<0:
    x1, x2, y1, y2 = -0.1*w, 0.05*w, 0.36*h, 0.51*h
else:
    x1, x2, y1, y2 = 0.45*w, 0.56*w, 0.9*h, 1.01*h
els=s.find_elements_in_region(x1,x2,y1,y2)
print(els)
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticklabels([])
axins.set_yticklabels([])
ax.indicate_inset_zoom(axins, edgecolor="black")
plt.show()
# %% plot_warping_values
section.plot_warping_values()#noqa
# %% contour_warping_values
section.contour_warping_values(levels=51)#noqa
# %% write_warping_csv
section.write_warping_csv()#noqa
# %% write_warping_gltf
section.write_warping_gltf()
# %% torsion stress plots from upstream
"""
Note to myself, figure details on how warping function values (at nodes)
are used with shape function derivates (at gauss integration points)
to get correct shear stresses.
https://sectionproperties.readthedocs.io/en/latest/rst/theory.html
refers to usage of smoothing matrix.
fea.py:element_stress
              sig_zxy_mzz_gp[i, :] = (
                  self.material.elastic_modulus
                  * Mzz
                  / j
                  * (B.dot(omega) - np.array([Ny, -Nx]))

"""
stress_post = section.calculate_stress(Mzz=1e6)#noqa
#ax_v=stress_post.plot_vector_mzz_zxy()
ax_c_xy=stress_post.plot_stress_mzz_zxy(normalize=False)
#ax_c_x=stress_post.plot_stress_mzz_zx()
#ax_c_y=stress_post.plot_stress_mzz_zy()
# %% analytic for circular, diameter=1: 5.09 MPa
# section-properties: 5.16 MPa (+ 1 %)
import math
d=1
r=d/2
Mzz=1e6
It=math.pi*math.pow(d,4)/32
tau=Mzz/It*r
print("d={0}, It={1:.3g}, Mzz={2:.3g}, tau={3:.3g}".format(d,It,Mzz,tau))
# %% analytic for chs, outer diameter=1, t=0.03: 23.2 MPa, It=0.0215
# section-properties(n=32), 965 elements: 24.3 MPa (+ 5 %), It=0.0212 (- 1 %)
# section-properties(n=64): 128-1023 elemenets 23.8 MPa (+ 3 %), It=0.0215
# section-properties(n=128): 256-512 elements 23.2 MPa It=0.0215
import math
d=1
t=0.03
d2=d-2*t
r=d/2
Mzz=1e6
It=math.pi*(math.pow(d,4)-math.pow(d2,4))/32
tau=Mzz/It*r
print("d={0}, It={1:.3g}, Mzz={2:.3g}, tau={3:.3g}".format(d,It,Mzz,tau))
# %% inset axes for contour_warping_values
import matplotlib.pyplot as plt
levels=51
fn=None
title=None
s=section#noqa
w=s.args.width#noqa
h=s.args.height#noqa
if s.is_shs():#noqa
    fn='warping-of-shs.pdf'
    title='Warping distribution of SHS using section properties'
(fig,ax)=s.contour_warping_values(levels=levels)#noqa
if title:
    ax.set_title(title)
axins = ax.inset_axes([0.1, 0.1, 0.7, 0.7])
s.set_box_aspect(axins)
x=s._mesh_nodes[:,0]
y=s._mesh_nodes[:,1]
triangles=s.get_triangles()#noqa
z=s.section_props.omega#noqa
trictr = axins.tricontourf(x, y, triangles, z,levels=levels)
# subregion of the original image
if s.args.primitive in ['circular','chs']:
    x1, x2, y1, y2 = -0.1*w, 0.05*w, 0.36*h, 0.51*h
else:
     x1, x2, y1, y2 = 0.92*w, 1.0*w, 0.92*h, 1.0*h
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticklabels([])
axins.set_yticklabels([])
ax.indicate_inset_zoom(axins, edgecolor="black")
if fn!=None:
    plt.savefig(fn)
    print('Wrote {0}'.format(fn))
plt.show()
