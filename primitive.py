r"""

Primirive
--------------

Calculate section properties of
 * rectangle(default) based on width(W) and height(H)
 * circular based on diameter(D) and count(N)
Mesh is refined until relative change of torsion and warping constants
is not more than rtol unless mesh_size is given
"""
import math
import argparse
import sectionproperties.pre.library.primitive_sections as sections
import sectionproperties.pre.library.steel_sections as steel_sections
#from sectionproperties.analysis.section import Section
import simo.dev
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--primitive", help="type of primitive",
                    default=simo.dev.RECTANGLE,
                    choices=simo.dev.PRIMITIVE_CHOICES)
simo.dev.add_common_arguments(parser)
args = parser.parse_args()
simo.dev.check_arguments(parser,args)
bending=args.bending
frame_analysis=args.frame_analysis
if args.primitive==simo.dev.RECTANGLE:
    args.title=("{2}: width = {0:.5g} and height = {1:.5g}".
      format(args.width, args.height,args.primitive))
    geometry = sections.rectangular_section(args.height, args.width)
elif args.primitive==simo.dev.RHS:
    if args.n_r>0:
        if args.radius<args.thickness:
            args.radius=2*args.thickness
        if args.n_r==1:
           args.n_r=2
    else:
        args.radius=0
    args.title=("""{2}: width={0:.5g}, height={1:.5g},
 thickness={3:.5g}, outer radius={4:.5g}, n_r={5}""".
       format(args.width, args.height,args.primitive,args.thickness,
              args.radius,args.n_r))
    geometry = steel_sections.rectangular_hollow_section(args.width,
        args.height,args.thickness,args.radius,args.n_r)
elif args.primitive==simo.dev.CIRCULAR:
    args.title=("{2}: diameter = {0:.5g} and count = {1}".
      format(args.diameter, args.count,args.primitive))
    geometry = sections.circular_section(args.diameter, args.count)
    args.width=args.diameter
    args.height=args.diameter
elif args.primitive==simo.dev.CHS:
    args.title=("""{2}: outer diameter={0:.5g},
thickness={1:.5g} and count={3}""".
       format(args.diameter, args.thickness, args.primitive, args.count))
    geometry = steel_sections.circular_hollow_section(args.diameter,
        args.thickness,args.count)
    args.width=args.diameter
    args.height=args.diameter
if args.plot_geometry:
    geometry.plot_geometry()
a=geometry.calculate_area()
it0=a
iw0=a
ms=math.pow(min(args.width,args.height)/2,2)
vertices0=0 # sometimes requesting smaller mesh size generates same mesh
section=None
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
    section = simo.dev.DevSection(geometry, time_info=args.time_info)
    section.set_args(args,it_num)
    if args.plot_mesh:
        section.plot_mesh()
    section.calculate_geometric_properties()
    if bending:
        print(("A = {0:.3g}, Ixx = {2:.3g}, Iyy = {1:.3g}, Ixy = {3:.3g}")
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
