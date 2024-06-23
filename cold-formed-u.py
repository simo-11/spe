r"""

Cold-formed-u-section
--------------

Mesh is refined until relative change of torsion and warping constants
is not more than rtol
"""
import math
import argparse
#from sectionproperties.analysis.section import Section
import simo.dev
import numpy as np
from shapely.geometry import Polygon
import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.pre as pre
from sectionproperties.pre.library.utils import draw_radius
import matplotlib.pyplot as plt
def u_section(
    d: float,
    b: float,
    t: float,
    r_out: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a U section (typical of cold-formed steel)
    with the bottom left corner at the
    origin *(0, 0)*, with depth/height *d*, width *b*, thickness *t*
    and outer radius *r_out*,
    using *n_r* points to construct the radius.
    If the outer radius is less than the thickness,
    the inner radius is set to zero.
    Code is based on steel_sections.cee_section

    :param float d: Depth/Height
    :param float b: Width
    :param float t: Thickness
    :param float r_out: Outer radius
    :param int n_r: Number of points discretising the outer radius
    :param Optional[sectionproperties.pre.pre.Material]:
        Material to associate with this geometry

    """
    points = []
    # calculate internal radius
    r_in = max(r_out - t, 0)
    # construct the outer bottom left radius
    points += draw_radius([r_out, r_out], r_out, np.pi, n_r)
    # bottom right corner
    points.append([b,0])
    points.append([b,t])
    # construct the inner bottom left radius
    points += draw_radius([t + r_in, t + r_in], r_in, 1.5 * np.pi, n_r, False)
    # construct the inner top left radius
    points += draw_radius([t + r_in, d - t - r_in], r_in, np.pi, n_r, False)
    # top right corner
    points.append([b,d-t])
    points.append([b,d])
    # construct the outer top left radius
    points += draw_radius([r_out, d - r_out], r_out, 0.5 * np.pi, n_r)
    polygon = Polygon(points)
    return geometry.Geometry(polygon, material)

parser = argparse.ArgumentParser(description=
    ('Calculate section properties for cold-formed U-section.'),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
simo.dev.add_common_arguments(parser)
args = parser.parse_args()
args.primitive=simo.dev.COLD_FORMED_U
simo.dev.check_arguments(parser,args)
if args.n_r>0:
    if args.radius<args.thickness:
        args.radius=2*args.thickness
    if args.n_r==1:
       args.n_r=2
else:
    args.radius=0
args.title=("""{5}: width={0:.5g}, height={1:.5g},
USection: thickness={2:.5g}, outer radius={3:.5g}, n_r={4}""".
      format(args.width, args.height,args.thickness,
             args.radius,args.n_r,args.section_type))
bending=args.bending
frame_analysis=args.frame_analysis
geometry = u_section(args.height, args.width,
                              args.thickness, args.radius, args.n_r)
if args.plot_section:
    fig, axes = plt.subplots()
    fig.set_size_inches(4,4)
    axes.set_aspect("equal", anchor="C")
    axes.set_title('U-{0:g}x{1:g}x{2:g}'.
                   format(args.height,args.width,args.thickness))
    # plot outline
    for (f) in geometry.facets:
        axes.plot(
            [geometry.points[f[0]][0], geometry.points[f[1]][0]],
            [geometry.points[f[0]][1], geometry.points[f[1]][1]],
            'k-',
            )
    axes.set_xticks([0,args.width])
    axes.set_yticks([0,args.height])
    t=args.thickness
    r_in=args.radius-t
    n_r=args.n_r
    if n_r>0:
        ai=n_r//2-1
        ap=draw_radius([t + r_in, t + r_in], r_in, 1.5 * np.pi, n_r, False)
        axes.annotate('r={0:.5g}'.format(args.thickness),
                      xycoords='data',
                      xy=(ap[ai][0],ap[ai][1]),
                      xytext=(0.25*args.width,0.25*args.height),
                      arrowprops=dict(arrowstyle='->')
                      )
    fn='gen/USection-{0:g}x{1:g}x{2:g}.pdf'.format(*tuple([f * 1000 for f in
                       (args.height,args.width,args.thickness)]));
    plt.tight_layout()
    plt.savefig(fn);
    print("Saved {0}".format(fn))
    plt.show()
if args.plot_geometry:
    geometry.plot_geometry()
a=geometry.calculate_area()
it0=a
iw0=a
ms=math.pow(6*args.thickness,2)
vertices0=0 # sometimes requesting smaller mesh size generates same mesh
it_num=0
while simo.dev.run(args):
    ms=0.82*ms
    if args.mesh_size:
        ms=args.mesh_size
    geometry.create_mesh(mesh_sizes=[ms])
    vertices=geometry.mesh.get('vertices').size
    if vertices0==vertices:
        continue
    vertices0=vertices
    it_num=it_num+1
    section = simo.dev.DevSection(geometry)
    section.set_args(args,it_num)
    if args.plot_mesh:
        section.plot_mesh()
    section.calculate_geometric_properties()
    if bending:
        print(("A = {0:.3g}, Ixx = {2:.3g}, Iyy = {1:.3g}, Ixy = {3:.3g}")
              .format(section.get_area(),*section.get_ic()))
        print(("Centroid: ({0:.3g},{1:.3g})".format(*section.get_c()))
              .format(section.get_area(),*section.get_ic()))
        bending=False
    if frame_analysis:
        (area, ixx, iyy, ixy, it, phi)=section.calculate_frame_properties()
        print(("f: A = {0:.3g}, Ixx = {2:.3g}, Iyy = {1:.3g}, "+
               "Ixy = {3:.3g}, J = {4:.3g}")
              .format(area,ixx,iyy,ixy,it))
        iwDiff=0
        iw=0
    else:
        section.calculate_warping_properties()
        it = section.get_j()
        if math.isnan(it):
            continue
        iw = section.get_gamma()
        iwDiff=abs((iw-iw0)/iw0)
        print(("It = {0:.3g}, Iw = {1:.3g}, k(steel) = {2:.2f}")
              .format(it,iw,math.sqrt(it/(2.6*iw))))
        print("Shear center: ({0:.3g},{1:.3g})".format(*section.get_sc()))
        if args.plot_warping_values:
            section.plot_warping_values()
        if args.write_warping_csv:
            section.write_warping_csv()
        if args.write_triangles_csv:
            section.write_triangles_csv()
    itDiff=abs((it-it0)/it0)
    if section.done(ms,itDiff,iwDiff):
        break
    else:
        it0=it
        iw0=iw
    section.plot_centroids()
