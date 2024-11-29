# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 16:36:57 2022

@author: simo.nikula@gmail.com
"""
#import argparse
from sectionproperties.analysis.section import Section
import pygltflib#noqa
import numpy as np
import csv
import math
import json
import subprocess
import os.path
import matplotlib.pyplot as plt
import sectionproperties.pre.library.utils as sp_utils
from sectionproperties.post.post import SectionProperties
from plyfile import PlyData, PlyElement
import types

consts = types.SimpleNamespace()

RECTANGLE='rectangle'
CIRCULAR='circular'
RHS='rhs'
CHS='chs'
BOX='box'
BOX_GIRDER='box_girder'
PRIMITIVE_CHOICES=[RECTANGLE,CIRCULAR,RHS,CHS,BOX,BOX_GIRDER]
COLD_FORMED_U='cold-formed-u'
consts.RHS=RHS
consts.COLD_FORMED_U=COLD_FORMED_U
class DevSection(Section):
    def __init__(self,*args, **kwargs):
        super().__init__(*args, **kwargs)
        self.triangles=None
        self.args=None
        self.gbtul=None
        self.log_write=True
    """
    """
    def get_triangles(self):
        if not isinstance(self.triangles,np.ndarray):
            ma=self.mesh_elements
            ne=len(ma)
            nt=4*ne
            ti=0
            triangles=np.empty([nt, 3],dtype='uint16')
            for i in range(0,ne):
                me=ma[i]
                triangles[ti,0]=me[0]
                triangles[ti,1]=me[3]
                triangles[ti,2]=me[5]
                ti+=1
                triangles[ti,0]=me[3]
                triangles[ti,1]=me[1]
                triangles[ti,2]=me[4]
                ti+=1
                triangles[ti,0]=me[3]
                triangles[ti,1]=me[4]
                triangles[ti,2]=me[5]
                ti+=1
                triangles[ti,0]=me[5]
                triangles[ti,1]=me[4]
                triangles[ti,2]=me[2]
                ti+=1
            self.triangles=triangles
        return self.triangles

    def find_elements_in_region(self,x1,x2,y1,y2):
        els=[]
        for el in self.elements:
            add=True
            for i in range(0,3):
                x=el.coords[0][i]
                if x<x1 or x>x2:
                    add=False
                    break
                y=el.coords[1][i]
                if y<y1 or y>y2:
                    add=False
                    break
            if add:
                els.append(el)
        return els

    def set_args(self,args,it_num=0):
        self.args=args
        if it_num==1:
            print(self.args.title)

    def get_box_aspect(self):
        if self.args is None:
            return (4,4,2)
        z=min(self.args.width,self.args.height)*self.args.z_scale
        return (self.args.width,self.args.height,z)

    def plot_warping_values(self,title=None,cmap=None,
                            **fig_kw:any):
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"},**fig_kw)
        ax.set_box_aspect(self.get_box_aspect())
        x=self.mesh_nodes[:,0]
        y=self.mesh_nodes[:,1]
        z=self.section_props.omega
        triangles=self.get_triangles()
        if cmap==None:
            cmap=plt.cm.seismic
        ax.plot_trisurf(x, y,triangles, z, cmap=cmap)
        if title==None:
            title=('{2}, {0} nodes, {1} elements'.format
            (self.num_nodes,len(self.elements),self.args.title))
        ax.set_title(title)
        xticks=np.linspace(min(x),max(x),3)
        ax.set_xticks(xticks)
        yticks=np.linspace(min(y),max(y),3)
        ax.set_yticks(yticks)
        zticks=np.linspace(min(self.section_props.omega),
                           max(self.section_props.omega),3)
        ax.set_zticks(zticks,labels=
                      ['','0',f'{zticks[2]:2.2}'])
        return (fig,ax)

    def contour_warping_values(self, title=None,
                               label="Warping",
                               levels=None,cmap=None,
                               **fig_kw:any):
        fig, ax = plt.subplots(**fig_kw)
        self.set_box_aspect(ax)
        x=self.mesh_nodes[:,0]
        y=self.mesh_nodes[:,1]
        triangles=self.get_triangles()
        z=self.section_props.omega
        if cmap==None:
            cmap=plt.cm.seismic
        trictr = ax.tricontourf(x, y, triangles, z,levels=levels,cmap=cmap)
        fig.colorbar(trictr, label=label, format="%.4g")
        if title==None:
            title=('{2}\n{0} nodes, {1} elements'.format
            (self.num_nodes,len(self.elements),self.args.title))
        ax.set_title(title)
        return (fig,ax)

    def set_box_aspect(self,ax):
        ax.set_box_aspect(self.args.height/self.args.width);

    def get_k(self,nu: float=0.3):
        return math.sqrt(self.get_j()/((2*(1+nu))*self.get_gamma()))

    def gfn(self,fn):
        """
        Provides fn prefixed with value of gen-parameter
        """
        return(self.args.gen+'/'+fn)

    def default_filename(self,suffix,use_case):
        words=[use_case]
        words.append(self.args.primitive)
        if self.uses_diameter():
            words.append("{0:g}-{1}".format(
                1000*self.args.diameter,self.args.count))
        else:
            words.append("{0:g}-{1:g}".format
                         (1000*self.args.height,1000*self.args.width))
        if self.uses_t():
            words.append("{0:g}".format(1000*self.args.thickness))
        if self.uses_web_thickness():
            words.append("{0:g}".format(1000*self.args.web_thickness))
        if self.uses_r():
            words.append("{0:g}".format(1000*self.args.radius))
        if self.uses_n_r():
            words.append("{0}".format(self.args.n_r))
        if not self.section_props.omega is None:
            words.append('{0}'.format(len(self.section_props.omega)))
        return '-'.join(words)+suffix

    def uses_diameter(self):
        return self.args.primitive in (CIRCULAR,CHS)

    def uses_web_thickness(self):
        return self.args.primitive in (BOX)

    def uses_t(self):
        return self.args.primitive in (RHS,CHS,COLD_FORMED_U,BOX)

    def uses_r(self):
        return self.args.primitive in (RHS,COLD_FORMED_U)

    def uses_n_r(self):
        return self.args.primitive in (RHS,COLD_FORMED_U)

    def wrote(self,name):
        if self.log_write:
            print(f"Wrote {name}")

    def write_json(self,fn=None):
        if fn==None:
            fn=self.default_filename(suffix='.json',use_case='results')
        with open(self.gfn(fn), 'w', newline='') as file:
            data={
                'area':self.get_area(),
                'c':self.get_c(),
                'gamma':self.get_gamma(),
                'ic':self.get_ic(),
                'ig':self.get_ig(),
                'ip':self.get_ip(),
                'j':self.get_j(),
                'q':self.get_j(),
                'sc':self.get_sc(),
                'z':self.get_z(),
            }
            print(json.dumps(data, indent=2), file=file)
            self.wrote(f"{fn}")

    def write_warping_csv(self,fn=None):
        if fn==None:
            fn=self.default_filename(suffix='.csv',use_case='warping')
        x=self.mesh_nodes[:,0]
        y=self.mesh_nodes[:,1]
        z=self.section_props.omega
        rows=np.empty([len(x),3],dtype=float)
        rows[:,0]=x
        rows[:,1]=y
        rows[:,2]=z
        with open(self.gfn(fn), 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['x','y','w'])
            writer.writerows(rows)
        self.wrote("{0}".format(fn))

    def write_triangles_csv(self,fn=None):
        if fn==None:
            fn=self.default_filename(suffix='.csv',use_case='triangles')
        rows=self.get_triangles();
        with open(self.gfn(fn), 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['f','s','t'])
            writer.writerows(rows)
        self.wrote("{0}".format(fn))
        
    def write_warping_gltf(self,fn=None):
        if fn==None:
            fn=self.default_filename(suffix='.glb',use_case='warping')
        ps=len(self.mesh_nodes)
        points=np.empty((ps,3),dtype="float32")
        points[:,0]=self.mesh_nodes[:,0]
        points[:,1]=self.mesh_nodes[:,1]
        points[:,2]=self.section_props.omega
        triangles=self.get_triangles()
        triangles_binary_blob = triangles.flatten().tobytes()
        points_binary_blob = points.tobytes()
        n_times=31
        d=n_times-1
        times=np.empty(n_times,dtype="float32")
        scales=np.ones((n_times,3),dtype="float32")
        scaler=0.4*self.get_box_aspect()[2]/max(self.section_props.omega)
        for i in range(0,n_times):
            times[i]=i
            scales[i,2]=math.sin(times[i]/d*2*math.pi)*scaler
        flat_scales=np.ones((n_times,3),dtype="float32")
        flat_scales[:,2]=0
        times_blob=times.tobytes()
        scales_blob=scales.tobytes()
        gltf = pygltflib.GLTF2(
        scene=0,
        scenes=[pygltflib.Scene(nodes=[0])],
        nodes=[pygltflib.Node(children=[1,2])
               ,pygltflib.Node(mesh=0,name="cross-section with warping")
               ,pygltflib.Node(mesh=1,name="undeformed cross-section",
                               scale=[1,1,0])
               ],
        meshes=[
            pygltflib.Mesh(
                primitives=[
                    pygltflib.Primitive(
                        attributes=pygltflib.Attributes(POSITION=1),
                        indices=0,
                        material=0
                    )])
            ,pygltflib.Mesh(
                primitives=[
                    pygltflib.Primitive(
                        attributes=pygltflib.Attributes(POSITION=1),
                        indices=0,
                        material=1
                    )])
        ],
        materials=[
            pygltflib.Material(pbrMetallicRoughness=
                               pygltflib.PbrMetallicRoughness(
                                   baseColorFactor=[0.2,0.2,0.2,0.85]),
                               doubleSided=True,
                               alphaCutoff=None,
                               name='cross-section with warping',
                               alphaMode='BLEND')
            ,pygltflib.Material(pbrMetallicRoughness=
                               pygltflib.PbrMetallicRoughness(
                                   baseColorFactor=[0.5,0.1,0.1,0.]),
                               doubleSided=True,
                               alphaCutoff=None,
                               name='undeformed cross-section',
                               alphaMode='OPAQUE')
        ],
        accessors=[
            pygltflib.Accessor(
                bufferView=0,
                componentType=pygltflib.UNSIGNED_SHORT,
                count=triangles.size,
                type=pygltflib.SCALAR,
                max=[int(triangles.max())],
                min=[int(triangles.min())],
            ),
            pygltflib.Accessor(
                bufferView=1,
                componentType=pygltflib.FLOAT,
                count=len(points),
                type=pygltflib.VEC3,
                max=points.max(axis=0).tolist(),
                min=points.min(axis=0).tolist(),
            ),
            pygltflib.Accessor(
                bufferView=2,
                componentType=pygltflib.FLOAT,
                count=n_times,
                type=pygltflib.SCALAR,
                max=[times.max().item()],
                min=[0],
            ),
            pygltflib.Accessor(
                bufferView=3,
                componentType=pygltflib.FLOAT,
                count=n_times,
                type=pygltflib.VEC3,
                max=scales.max(axis=0).tolist(),
                min=scales.min(axis=0).tolist(),
            ),
        ],
        bufferViews=[
            pygltflib.BufferView(
                buffer=0,
                byteLength=len(triangles_binary_blob),
                target=pygltflib.ELEMENT_ARRAY_BUFFER,
                name='triangles',
            ),
            pygltflib.BufferView(
                buffer=0,
                byteLength=len(points_binary_blob),
                byteOffset=len(triangles_binary_blob),
                target=pygltflib.ARRAY_BUFFER,
                name='points',
            ),
            pygltflib.BufferView(
                buffer=0,
                byteLength=len(times_blob),
                byteOffset=len(triangles_binary_blob)
                    +len(points_binary_blob)
                ,name='times',
            ),
            pygltflib.BufferView(
                buffer=0,
                byteLength=len(scales_blob),
                byteOffset=len(triangles_binary_blob)
                    +len(points_binary_blob)
                    +len(times_blob)
                ,name='scales',
            ),
        ],
        buffers=[
            pygltflib.Buffer(byteLength=len(triangles_binary_blob)
                             +len(points_binary_blob)
                             +len(times_blob)
                             +len(scales_blob)
                             )
        ],
        animations=[
            pygltflib.Animation(name="warping",
                                channels=[pygltflib.AnimationChannel(
                                    sampler=0,
                                    target=pygltflib.AnimationChannelTarget(
                                        node=0,path='scale')
                                    )],
                                samplers=[pygltflib.AnimationSampler(
                                    input=2,output=3),
                                    ]),
        ]
        )
        gltf.set_binary_blob(triangles_binary_blob + points_binary_blob+
                             times_blob+scales_blob)
        #
        #
        gfn=self.gfn(fn)
        gltf.save(gfn)
        self.wrote("{0}".format(fn))

    def write_warping_ply(self,fn=None):
        if fn==None:
            fn=self.default_filename(suffix='.ply',use_case='warping')
        scaler=0.4*self.get_box_aspect()[2]/max(self.section_props.omega)
        ps=len(self.mesh_nodes)
        data=[(self.mesh_nodes[i,0],
               self.mesh_nodes[i,1],
               scaler*self.section_props.omega[i]) for i in range(ps)]
        points=np.array(data,dtype=[('x', 'f4'), ('y', 'f4'),
                            ('z', 'f4')])
        face=self.get_triangles()
        face_data=[(([face[i,0],
               face[i,1],
               face[i,2]],)) for i in range(len(face))]
        # use int as e.g. gigamesh reader does not support shorts
        faces=np.array(face_data,dtype=[('vertex_indices', 'int', (3,))])
        ply_nodes = PlyElement.describe(points, 'vertex')
        ply_elements = PlyElement.describe(faces,'face',
                                 comments=['simo.dev.write_warping_ply'])
        
        gfn=self.gfn(fn)
        PlyData([ply_nodes,ply_elements], text=False).write(gfn);
        self.wrote("{0}".format(fn))

    def init_gbt_u(self):
        b=self.args.width
        d=self.args.height
        t=self.args.thickness
        ht=t/2
        r=max(self.args.radius-ht,0)
        n_r=self.args.n_r
        points=[]
        points.append((b,ht))
        # construct the inner bottom left radius
        points += sp_utils.draw_radius([ht + r, ht + r], 
                              r, 1.5 * np.pi, n_r, False)
        # construct the inner top left radius
        points += sp_utils.draw_radius([ht + r, d - ht - r], 
                              r, np.pi, n_r, False)
        # top right corner
        points.append((b,d-ht))
        self.gbt_points=points
        self.gbt_elements=[]
        mi=len(points)
        for i in range(mi-1):
            j=i+1
            e={'i':i,'j':j,'m':1,'inc':0,'t':t}
            self.gbt_elements.append(e)

    def init_gbt_rhs(self):
        b=self.args.width
        d=self.args.height
        t=self.args.thickness
        ht=t/2
        r=max(self.args.radius-ht,0)
        n_r=self.args.n_r
        # adapted based on steel_sections.rectangular_hollow_section
        # construct the center line points
        points=[]
        points += sp_utils.draw_radius((ht + r, ht + r), r, np.pi, n_r)
        points += sp_utils.draw_radius(
            (b - ht - r, ht + r), r, 1.5 * np.pi, n_r
        )
        points += sp_utils.draw_radius((b - ht - r, d - ht - r), r, 0, n_r)
        points += sp_utils.draw_radius(
            (ht + r, d - ht - r), r, 0.5 * np.pi, n_r
        )
        self.gbt_points=points
        self.gbt_elements=[]
        mi=len(points)
        for i in range(mi):
            if i==mi-1:
                j=0
            else:
                j=i+1
            e={'i':i,'j':j,'m':1,'inc':0,'t':t}
            self.gbt_elements.append(e)

    def plot_gbt(self,**fig_kw:any):
        fig, ax = plt.subplots(**fig_kw)
        for f in self.gbt_elements:
            i_node=self.gbt_points[f['i']]
            j_node=self.gbt_points[f['j']]               
            ax.plot(
                [i_node[0], j_node[0]],
                [i_node[1], j_node[1]],
                "ko-",
                markersize=2,
                linewidth=1.5,
            )
        ax.set_aspect('equal',adjustable='box')
        return (fig,ax)   
            
    def write_gbtul(self):
        match self.args.primitive:#noqa
            case consts.RHS:
                self.init_gbt_rhs()
            case consts.COLD_FORMED_U:
                self.init_gbt_u()
            case _:
                raise Exception(f'''{self.args.primitive} \
 is not supported''')
        with open(self.gfn('inp1.txt'), 'w') as f:
            f.write(f"{len(self.gbt_points)}\n")
            for p in self.gbt_points:
                f.write(f"{p[0]:.4} {p[1]:.4}\n")                
            f.write(f"{len(self.gbt_elements)}\n")
            for e in self.gbt_elements:
                f.write((f"{e['i']+1} {e['j']+1} "
                          f"{e['m']} {e['inc']} {e['t']}\n"))
            f.write("1\n") # material, currently fixed steel
            E=21e10
            nu=0.3
            rho=7800
            f.write(f"{E:.4} {E:.4} {nu} {nu} {E/(2*(1+nu)):.4} {rho}\n")
            #Length-Distributed Elastic Supports and Additional Masses
            #Distributed along longitudinal strips
            f.write("0\n0\n")
                
    def read_gbtul(self):
        na=np.fromfile(self.gfn('out1.txt'),sep=' ')
        self.gbtul=SectionProperties(area=na[0],
                                    gamma=na[3],
                                    j=na[4]
                                    )
        
    def run_gbtul(self):
        if not os.path.isfile(self.args.gbt):
            raise Exception(f'GBT is not available at {self.args.gbt}')
        sp_args=[self.args.gbt,"1"]
        od=self.gfn("/Output_Files")
        if not os.path.isdir(od):
            os.makedirs(od)
        try:
            self.write_gbtul()
            cps=subprocess.run(sp_args,
                       cwd=self.args.gen,
                       capture_output=True,
                       check=True)
            with open(self.gfn(f'gbtout-{self.args.n_r}.txt'), 'w') as f:
                f.write(cps.stdout.decode())
            with open(self.gfn(f'gbterr-{self.args.n_r}.txt'), 'w') as f:
                f.write(cps.stderr.decode())
            self.read_gbtul()
        except subprocess.CalledProcessError as e:
            print(f"""Running gbtul solver {self.args.gbt} failed
stdout={e.stdout.decode()}
stderr={e.stderr.decode()}             
""")
            raise
        except BaseException:
            print(f"""Running gbtul solver {self.args.gbt} failed""")
            raise

    def done(self,ms,itDiff,iwDiff):
        if self.args.mesh_size:
            print(("meshSize = {0:.3g}, {1} nodes, {2} elements")
              .format(ms,self.num_nodes,len(self.elements)))
        else:
            print(("meshSize = {0:.3g}, {3} nodes, {4} elements, "+
                 "itDiff = {1:.3g}, iwDiff = {2:.3g}")
              .format(ms,itDiff,iwDiff,
                      self.num_nodes,len(self.elements)))
        return self.args.mesh_size or (
            itDiff<self.args.rtol and iwDiff<self.args.rtol)

    def is_shs(self):
        return (self.args.primitive=='rhs'
                and self.args.width==self.args.height)

def add_common_arguments(parser):
    parser.add_argument("--rtol", help="relative tolerance",
                        default=1e-2,type=float)
    parser.add_argument("--mesh_size",
                        help="fixed meshSize, rtol is not used",
                        type=float)
    parser.add_argument("-M","--plot_mesh", help="plot each mesh",
                        action="store_true")
    parser.add_argument("-G","--plot_geometry", help="plot geometry",
                        action="store_true")
    parser.add_argument("-P","--plot_section", help="Plot section",
                    action="store_true")
    parser.add_argument("-B","--bending",
                        help="show bending related constants",
                        action="store_true")
    parser.add_argument("-F","--frame_analysis",
                        help="show frame analysis results",
                        action="store_true")
    parser.add_argument("-A","--run_analysis",
                        help="run analysis",
                        action="store_true")
    parser.add_argument("--time_info",
                        help="show detailed info for computation",
                        action="store_true")
    parser.add_argument("--plot_warping_values",
                        help="plot warping values for each iteration",
                        action="store_true")
    parser.add_argument("--z_scale",
        help="scaling of z in plot_warping_values and write_warping_gltf",
        default=0.5,type=float)
    parser.add_argument("--write_warping_csv",
                        help="write warping values for each iteration",
                        action="store_true")
    parser.add_argument("--write_triangles_csv",
                        help="write triangles for each iteration",
                        action="store_true")
    parser.add_argument("-W","--width", help="width",
                        default=1,type=float)
    parser.add_argument("-H","--height", help="height",
                        default=1,type=float)
    parser.add_argument("-T","--thickness", help="thickness",
                        default=0.004,type=float)
    parser.add_argument("--web_thickness", help="web thickness",
                        default=0.004,type=float)
    parser.add_argument("-R","--radius", help="""outer radius,
    if < thickness, 2*thickness is used""",
                        default=0,type=float)
    parser.add_argument("--n_r", help="number of points in radius, 0 or >1",
                        default=4,type=int)
    parser.add_argument("-D","--diameter", help="diameter",
                        default=1,type=float)
    parser.add_argument("-N","--count", help="count of points for diameter",
                        default=32,type=int)
    parser.add_argument("--title")
    parser.add_argument("--section_type")
    parser.add_argument("--gen", help="""directory for generated files""",
                        default="gen")
    parser.add_argument("--gbtul", help="""run gbtul""",
                        action="store_true")
    parser.add_argument("--gbt",help="location of GBT solver", 
                        default="../../gbtul/GBT/GBT.exe")

def check_arguments(parser,args):
    if (not args.plot_section and not args.plot_geometry
        and not args.run_analysis and not args.mesh_size
        and not args.plot_warping_values
        and not args.gbtul
        ):
        parser.print_help()
    if(not os.path.exists(args.gen)):
        os.makedirs(args.gen)
        print(f"created directory {args.gen} for generated files")

def run(args):
    return (args.run_analysis or args.mesh_size 
           or args.plot_warping_values or args.gbtul)
