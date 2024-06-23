# Topics

## Provide warping function information
 * Plot
 * Numerical values
 * Function fitting approximation

### Analysis
Warping values are available in section.section_props.omega\[el.node_ids\] where el in section.mesh_elements.
coordinates in section.mesh_nodes.
See 
 * assemble_sc_warping_integrals in section.py
 * shear_warping_integrals in fea.py

### Changes
Samples in [primitive.py](primitive.py) and [cold-formed-u.py](cold-formed-u.py).

Run using [sample-cells.py](sample-cells.py).
In Spyder set working directory to section-properties\examples\11-simo and press CTRL-Enter on wanted cell.


Common parts in [dev.py](simo/dev.py)

#### matplotlib.pyplot plotting
Usable but a bit limited 3D-plotting can be created e.g. using plot_trisurf in matplotlib.pyplot.

```
runfile('C:/github/section-properties/examples/11-simo/primitive.py',args='--primitive=rectangle -H 1 -W 1 --plot_warping_values -A --rtol=0.001')
Rectangle: width = 1 and height = 1, rtol=0.001
It = 0.146, Iw = 0.000174
meshSize = 0.205, 25 nodes, 8 elements, itDiff = 0.854, iwDiff = 1
It = 0.143, Iw = 0.000156
meshSize = 0.113, 41 nodes, 16 elements, itDiff = 0.0204, iwDiff = 0.102
It = 0.141, Iw = 0.000133
meshSize = 0.0623, 73 nodes, 28 elements, itDiff = 0.012, iwDiff = 0.146
It = 0.141, Iw = 0.000132
meshSize = 0.0419, 85 nodes, 34 elements, itDiff = 0.000743, iwDiff = 0.00728
It = 0.141, Iw = 0.000136
meshSize = 0.0344, 101 nodes, 42 elements, itDiff = 0.0005, iwDiff = 0.0316
It = 0.141, Iw = 0.000138
meshSize = 0.0282, 145 nodes, 64 elements, itDiff = 0.00121, iwDiff = 0.00952
It = 0.141, Iw = 0.000135
meshSize = 0.0231, 163 nodes, 70 elements, itDiff = 0.000712, iwDiff = 0.0179
It = 0.141, Iw = 0.000136
meshSize = 0.0189, 194 nodes, 85 elements, itDiff = 1.15e-05, iwDiff = 0.00513
It = 0.141, Iw = 0.000134
meshSize = 0.0155, 236 nodes, 103 elements, itDiff = 0.000415, iwDiff = 0.0108
It = 0.141, Iw = 0.000135
meshSize = 0.0127, 283 nodes, 126 elements, itDiff = 2.25e-05, iwDiff = 0.00604
It = 0.141, Iw = 0.000135
meshSize = 0.0104, 338 nodes, 153 elements, itDiff = 0.000111, iwDiff = 0.00332
It = 0.141, Iw = 0.000135
meshSize = 0.00857, 401 nodes, 184 elements, itDiff = 0.000111, iwDiff = 0.000247
```
![image](https://user-images.githubusercontent.com/1210784/218270574-bd96f8fc-b424-4c5c-8935-d090531b7151.png)
..
![image](https://user-images.githubusercontent.com/1210784/218270602-a13d46d0-77de-4c03-9ace-ee667b0cd00b.png)
..
![image](https://user-images.githubusercontent.com/1210784/218270622-84012158-5339-4198-954f-d338a41e3022.png)

#### write coordinates, warping and triangles to csv 
 * --write_warping_csv
   * coordinates for nodes x, y 
   * warping value w 
 * --write_triangles_csv
   * triangles for visualization
   * each parabolic (six noded) triangle is subdivided to four three noded triangles 
   * indexes are 0-based
```
runfile('C:/github/section-properties/examples/11-simo/cold-formed-u.py',args='-A -B --write_warping_csv --write_triangles_csv')
Cold-formed-U: width = 0.05, height = 0.1,
thickness= 0.004, outer radius=0.008, n_r=4
rtol=0.001
A = 0.000744, Ixx = 1.81e-07, Iyy = 1.12e-06, Ixy = 2.12e-22
Centroid: (0.0144,0.05)
It = 3.91e-09, Iw = 2.75e-10
Shear center: (-0.0166,0.05)
Wrote USection-100x50x4-8-4-249.csv
Wrote USection-tri-100x50x4-8-4-249.csv
meshSize = 0.025, 249 nodes, 86 elements, itDiff = 1, iwDiff = 1
It = 3.91e-09, Iw = 2.75e-10
Shear center: (-0.0166,0.05)
Wrote USection-100x50x4-8-4-267.csv
Wrote USection-tri-100x50x4-8-4-267.csv
```
##### Usage in matlab for interpolation of scattered data
```
t=readtable('C:/github/section-properties/examples/11-simo/USection-100x50x4-8-4-267.csv');
F = scatteredInterpolant(t.x,t.y,t.w);
F([0.05 0.05 0.05],[0 0.002 0.004])
ans =
   -0.0032   -0.0033   -0.0034
```

##### Usage in matlab e.g. for plotting
```
t=readtable('C:/github/section-properties/examples/11-simo/USection-100x50x4-8-4-267.csv');
t2=readtable('C:/github/section-properties/examples/11-simo/USection-tri-100x50x4-8-4-267.csv');
T=[t2.f, t2.s, t2.t]+1; % +1 to get proper indexes to coordinates in t
TO=triangulation(T,t.x,t.y,t.w);
trisurf(TO)
axis equal;
```
Or using show_section
```
show_section('USection','100x50x4-8-4-267');
```
Use daspect command in matlab to scale data for viewing e.g. instead of default [1 1 1] due to axis equal
```
daspect([1 1 0.2])
```
![image](https://user-images.githubusercontent.com/1210784/181523464-fcfa25a5-6be3-4f1d-bf58-70d1d6397418.png)
..
![image](https://user-images.githubusercontent.com/1210784/181523584-a583d1db-8aca-4058-a63b-55e7a41bbb1d.png)

## Examples

  * [rectangle](rectangle.md)
