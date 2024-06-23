# Rectangle

## Square 100-100

### In spyder

Create data files in csv format for warping functions
```
runcell('rectangle 100-100', 'C:/Users/simon/github/section-properties/examples/11-simo/sample_cells.py')
rectangle: width = 0.1 and height = 0.1
It = 1.43e-05, Iw = 1.3e-10, k(steel) = 205.28
meshSize = 0.001, 41 nodes, 16 elements
runcell('write_warping_csv', 'C:/Users/simon/github/section-properties/examples/11-simo/sample_cells.py')
Wrote warping-rectangle-100-100-41.csv

rectangle: width = 0.1 and height = 0.1
It = 1.41e-05, Iw = 1.35e-10, k(steel) = 200.51
meshSize = 0.0001, 357 nodes, 162 elements
runcell('write_warping_csv', 'C:/Users/simon/github/section-properties/examples/11-simo/sample_cells.py')
Wrote warping-rectangle-100-100-357.csv
```
### In matlab 

Create surface fit to data provided by csv files using Curve fitting toolbox
 * cubicinterp - Cubic spline interpolation
 * poly44 - polynomial surfaces where i is the degree in x and j is the degree in y. The maximum for both i and j is five. For rectangle 44 gives quite good results.
 * anonymous function - see fittype for details
   * sinhs (sinh series) test process by using known analytical solution https://www.researchgate.net/publication/361446204_Efficient_modeling_and_order_reduction_of_new_3D_beam_elements_with_warping_via_absolute_nodal_coordinate_formulation
     (or https://en.wikiversity.org/wiki/Warping_functions#Example_3:_Rectangular_Cylinder) for few lowest values of n  

using 41 node function as input to cubicinterp produces reasonable approximation and recalculated value of Iw is about 8 % too high
```
> fprintf("%s\n",pwd)
C:\Users\simon\github\section-properties\examples\11-simo
>> dir gen/*.csv
warping-rectangle-100-100-41.csv
>> t41=readtable('gen/warping-rectangle-100-100-41.csv');
>> c41=fit([t41.x t41.y],t41.w,'cubicinterp');
>> plot(c41,[t41.x t41.y],t41.w);
```
![image](https://github.com/simo-11/section-properties/assets/1210784/7a02dc7d-5467-40ac-988a-1167b797ca06)
```
>> fun41=@(x,y)c41(x,y).^2;
>> fprintf("%.3g\n",integral2(fun41,0.,0.1,0,0.1))
1.46e-10
```
Using 357 nodes reduces error in Iw to less than 1 %.
```
>> t2=readtable('gen/warping-rectangle-100-100-357.csv');
>> f2=fit([t2.x t2.y],t2.w,'cubicinterp');
>> plot(f2,[t2.x t2.y],t2.w);
>> w2=@(x,y)f2(x,y).^2;
>> fprintf("%.3g\n",integral2(w2,0.,0.1,0,0.1))
1.36e-10
```
Create mex_function using [w_100_100_41.m](w_100_100_41.m) is not supported
```
x=0.1;y=0.025;codegen w_100_100_41.m -args {x,y}
Undefined function or variable 'f1'. For code generation, all variables must be fully
defined before use.
>> x=0.1;y=0.025;codegen w_100_100_41.m -args {f1,x,y}
Function input at args{1} does not have a valid type.

Caused by:
    Type conversion failed at <obj>(1).expr.
        Class function_handle is not supported by coder.Type.
```
### Upper level script
Reads all files selected by rectangle dimensions and performs cubicinterp and poly44 fits.
```
>> c100=testRectangle(height=100,width=100)
>> c100{1}
 struct with fields:

                  t: [357×3 table]
               file: "C:\Users\simon\github\section-properties\examples\11-simo\gen/warping-rectangle-100-100-357.csv"
        cubicinterp: 1.3570e-10
    cubicinterp_fit: [1×1 sfit]
             poly44: 1.3591e-10
         poly44_fit: [1×1 sfit]
>> figure(101)
>> plot(c100{1}.poly44_fit,[c100{1}.t.x c100{1}.t.y],c100{1}.t.w);
```
![image](https://github.com/simo-11/section-properties/assets/1210784/0d7dee0b-3db8-4541-8883-41deed27b57b)

Checking parameters for models providing them.
Comments after --
```
>> w=testRectangle(debugLevel=2,models=["poly44"],plot=0);
file=warping-rectangle-100-100-357.csv
x values: 0 - 0.1
y values: 0 - 0.1
w values: -0.000366 - 0.000366
model=poly44
     Linear model Poly44:
     f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y 
                    + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 
                    + p13*x*y^3 + p04*y^4
     Coefficients (with 95% confidence bounds):
       p00 =   1.644e-07  (-6.632e-06, 6.961e-06) -- not important <0.001*w  
       p10 =    -0.03778  (-0.03834, -0.03721) -- important, transforms origin
       p01 =     0.03776  (0.0372, 0.03833) -- important, transforms origin
       p20 =       1.134  (1.115, 1.153) -- important 
       p11 =  -0.0001405  (-0.0147, 0.01442) -- not important <0.003*w
       p02 =      -1.133  (-1.152, -1.114) -- important
       p30 =      -7.564  (-7.828, -7.3) -- important
       p21 =      -22.67  (-22.86, -22.47) -- important
       p12 =       22.67  (22.47, 22.86) -- important
       p03 =       7.548  (7.282, 7.815) -- important
       p40 =     0.04689  (-1.232, 1.325) -- not important
       p31 =       151.1  (150, 152.2) -- important
       p22 =     0.02697  (-1.039, 1.093) -- not important
       p13 =      -151.1  (-152.2, -150) -- important
       p04 =     0.03326  (-1.261, 1.327) -- not important
```

### sinh series with various numbers of n
![image](https://github.com/simo-11/section-properties/assets/1210784/ce6fd20c-cd56-4385-8753-74d2c95f5f45)

![image](https://github.com/simo-11/section-properties/assets/1210784/fb15f3f0-07b0-4815-a59e-6a73236e3940)

![image](https://github.com/simo-11/section-properties/assets/1210784/e42e64af-48b2-412f-b412-cf58422712e1)

## 100-10
Iw should be about 6.6e-12

### Spyder
```
runcell('rectangle 100-10', 'C:/Users/simon/github/section-properties/examples/11-simo/sample_cells.py')
Reloaded modules: simo, simo.dev
rectangle: width = 0.1 and height = 0.01
It = 3.13e-08, Iw = 6.64e-12, k(steel) = 42.54
meshSize = 1e-05, 375 nodes, 160 elements
runcell('write_warping_csv', 'C:/Users/simon/github/section-properties/examples/11-simo/sample_cells.py')
Wrote warping-rectangle-10-100-375.csv
```
### Matlab
```
>> r=testRectangle(height=10,width=100,plot=1);
```
![image](https://github.com/simo-11/section-properties/assets/1210784/8ac719b9-badd-45c2-b8c1-9966f7061771)

```
>> w=testRectangle(debugLevel=2,models=["poly44"],plot=0,height=10,width=100);
file=warping-rectangle-10-100-375.csv
x values: 0 - 0.1
y values: 0 - 0.01
w values: -0.000226 - 0.000226
model=poly44
     Linear model Poly44:
     f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y 
                    + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 
                    + p13*x*y^3 + p04*y^4
     Coefficients (with 95% confidence bounds):
       p00 =  -0.0002339  (-0.0002361, -0.0002317) -- important
       p10 =    0.003616  (0.003414, 0.003819) -- important
       p01 =     0.04727  (0.04568, 0.04885) -- important
       p20 =     0.03242  (0.02553, 0.0393) -- important 
       p11 =     -0.7427  (-0.7869, -0.6985) -- important
       p02 =     -0.1243  (-0.6633, 0.4147)
       p30 =     -0.2292  (-0.3249, -0.1336)
       p21 =      -6.124  (-6.757, -5.491) -- important
       p12 =       3.039  (-2.694, 8.772)
       p03 =       7.493  (-68.29, 83.27)
       p40 =      0.1197  (-0.3439, 0.5834)
       p31 =        40.7  (37.15, 44.25) -- important
       p22 =       1.631  (-29.47, 32.73)
       p13 =      -215.7  (-526.4, 95.02)
       p04 =       100.2  (-3572, 3772)
           sse: 1.7980e-09
       rsquare: 0.9995
           dfe: 360
    adjrsquare: 0.9994
          rmse: 2.2348e-06

Iw=6.62e-12
```

## 100-4
Iw should be about 4.4e-13

### Spyder

```
runcell('rectangle 100-4', 'C:/Users/simon/github/section-properties/examples/11-simo/sample_cells.py')
rectangle: width = 0.004 and height = 0.1
It = 2.09e-09, Iw = 4.42e-13, k(steel) = 42.66
meshSize = 1e-05, 171 nodes, 56 elements
runcell('write_warping_csv', 'C:/Users/simon/github/section-properties/examples/11-simo/sample_cells.py')
Wrote warping-rectangle-100-4-171.csv

runcell('rectangle 100-4', 'C:/Users/simon/github/section-properties/examples/11-simo/sample_cells.py')
rectangle: width = 0.004 and height = 0.1
It = 2.08e-09, Iw = 4.41e-13, k(steel) = 42.58
meshSize = 1e-06, 1397 nodes, 620 elements
runcell('write_warping_csv', 'C:/Users/simon/github/section-properties/examples/11-simo/sample_cells.py')
Wrote warping-rectangle-100-4-1397.csv
```
### matlab
```
>> w=testRectangle(debugLevel=2,models=["poly44"],plot=0,height=100,width=4);
file=warping-rectangle-100-4-1397.csv
x values: 0 - 0.004
y values: 0 - 0.1
w values: -9.62e-05 - 9.62e-05
model=poly44
     Linear model Poly44:
     f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y 
                    + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 
                    + p13*x*y^3 + p04*y^4
     Coefficients (with 95% confidence bounds):
       p00 =   9.869e-05  (9.853e-05, 9.885e-05) -- important
       p10 =    -0.04944  (-0.04973, -0.04914) -- important, transforms origin
       p01 =   -0.001887  (-0.001902, -0.001872) -- important, transforms origin
       p20 =     0.05794  (-0.197, 0.3129) -- not important
       p11 =      0.9445  (0.9359, 0.9531) -- important
       p02 =   -0.002523  (-0.00304, -0.002006) -- important
       p30 =       2.479  (-87.79, 92.75) -- not important
       p21 =       -2.29  (-5.106, 0.5266) -- not important
       p12 =       1.353  (1.229, 1.476) -- important
       p03 =     0.01494  (0.007673, 0.02221) -- not important
       p40 =       -2769  (-1.379e+04, 8253) -- not important
       p31 =       394.7  (7.439, 781.9) -- not important
       p22 =      0.1674  (-15.76, 16.1) -- not important
       p13 =      -8.992  (-9.699, -8.286) -- important
       p04 =     0.01662  (-0.01898, 0.05222) -- not important
           sse: 1.5201e-10
       rsquare: 0.9999
           dfe: 1382
    adjrsquare: 0.9999
          rmse: 3.3165e-07

Iw=4.4e-13
file=warping-rectangle-100-4-171.csv
x values: 0 - 0.004
y values: 0 - 0.1
w values: -9.6e-05 - 9.59e-05
model=poly44
     Linear model Poly44:
     f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y 
                    + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 
                    + p13*x*y^3 + p04*y^4
     Coefficients (with 95% confidence bounds):
       p00 =   9.991e-05  (9.886e-05, 0.000101)
       p10 =       4e+11  (2.392e+11, 5.608e+11)
       p01 =   -0.001938  (-0.00204, -0.001837)
       p20 =  -2.356e+14  (-3.25e+14, -1.462e+14)
       p11 =  -1.163e+12  (-4.14e+12, 1.815e+12)
       p02 =   -0.001633  (-0.005214, 0.001948)
       p30 =   1.684e+15  (-1.382e+16, 1.719e+16)
       p21 =   8.719e+14  (-1.361e+15, 3.105e+15)
       p12 =       0.804  (0.1193, 1.489)
       p03 =     0.01259  (-0.03783, 0.063)
       p40 =   8.053e+18  (3.213e+18, 1.289e+19)
       p31 =  -1.453e+17  (-5.175e+17, 2.269e+17)
       p22 =     -0.8057  (-81.37, 79.76)
       p13 =      -5.904  (-9.805, -2.003)
       p04 =   -0.004945  (-0.2506, 0.2407)
           sse: 9.0319e-11
       rsquare: 0.9998
           dfe: 156
    adjrsquare: 0.9997
          rmse: 7.6090e-07

Iw=6.96e+12
Model poly44 for file warping-rectangle-100-4-171.csv rejected, Iw=6.96e+12>1.47e-11
```
![image](https://github.com/simo-11/section-properties/assets/1210784/7dd79516-126c-4dee-be22-dfe763fa7476)

![image](https://github.com/simo-11/section-properties/assets/1210784/49bf0827-c37c-4ad0-ad57-badac2ce67d6)
