# Introduction
Profile used in Princeton beam experiment is common reference case for ANCF implementations e.g. [https://www.sciencedirect.com/science/article/pii/S0997753822002546].
In this context various interpolations are compared.

# Results from section-properties
Source is [princeton.py](./princeton.py) and results are summarized below
rectangle of 12.377 x 3.2024 is used (too many digits, but this is how they are usually given)

|ec_in_h|Number of nodes|Iw 10E-18 m^6|It 10E-12 m^4|
|--|---------------|-----|------|
|10|24|334|118|
|14|27|336|118|
|15|38|337|117|
|17|46|335|116|
|20|57|335|116|
|23|65|327|114|
|25|90|328|114|
|30|114|328|114|
|35|148|328|114|
|40|204|328|113|
|50|294|328|113|
|100|1060|328|113|

# Analytical results
[princeton.m](./princeton.m) cell  analytical solid rectangle which uses [iw_rect.m](./iw_rect.m)
```
Analytical Iw for 3.2024x12.377, for thin 4.32E-16
Iw=3.28E-16 (n=6) took 1.4912 ms
Iw=3.27E-16 (n=5) took 0.4313 ms
Iw=3.28E-16 (n=4) took 0.3262 ms
Iw=3.28E-16 (n=3) took 0.2764 ms
Iw=3.28E-16 (n=2) took 0.2899 ms
Iw=3.29E-16 (n=1) took 0.2419 ms
Iw=3.33E-16 (n=0) took 0.1278 ms
```

# interpolate based on section properties results
```
file=warping-rectangle-12.377-3.2024-24.csv
model=cubicinterp-integral2, Iw=3.39e-16, cub took 1.56 ms
model=linearinterp-integral2, Iw=3.33e-16, cub took 0.718 ms
file=warping-rectangle-12.377-3.2024-57.csv
model=cubicinterp-integral2, Iw=3.33e-16, cub took 1.88 ms
model=linearinterp-integral2, Iw=3.3e-16, cub took 0.74 ms
file=warping-rectangle-12.377-3.2024-114.csv
model=cubicinterp-integral2, Iw=3.28e-16, cub took 2.24 ms
model=linearinterp-integral2, Iw=3.28e-16, cub took 0.944 ms
file=warping-rectangle-12.377-3.2024-294.csv
model=cubicinterp-integral2, Iw=3.28e-16, cub took 4.46 ms
model=linearinterp-integral2, Iw=3.27e-16, cub took 1.32 ms
file=warping-rectangle-12.377-3.2024-1060.csv
model=cubicinterp-integral2, Iw=3.28e-16, cub took 9.77 ms
model=linearinterp-integral2, Iw=3.28e-16, cub took 1.12 ms
```

# write warping results based on interpolation f
uses [save_interpolation_results.m](./save_interpolation_results.m)
```
aved gen/cubicinterp-12.377-3.2024-24-10x10.xlsx, max error percent=9.6
Saved gen/cubicinterp-12.377-3.2024-24-100x100.xlsx, max error percent=9.3
Saved gen/linearinterp-12.377-3.2024-24-10x10.xlsx, max error percent=8.7
Saved gen/linearinterp-12.377-3.2024-24-100x100.xlsx, max error percent=8.4
Saved gen/cubicinterp-12.377-3.2024-57-10x10.xlsx, max error percent=5.3
Saved gen/cubicinterp-12.377-3.2024-57-100x100.xlsx, max error percent=5.6
Saved gen/linearinterp-12.377-3.2024-57-10x10.xlsx, max error percent=5.4
Saved gen/linearinterp-12.377-3.2024-57-100x100.xlsx, max error percent=6.5
Saved gen/cubicinterp-12.377-3.2024-114-10x10.xlsx, max error percent=2.3
Saved gen/cubicinterp-12.377-3.2024-114-100x100.xlsx, max error percent=2.4
Saved gen/linearinterp-12.377-3.2024-114-10x10.xlsx, max error percent=2.8
Saved gen/linearinterp-12.377-3.2024-114-100x100.xlsx, max error percent=3.2
Saved gen/cubicinterp-12.377-3.2024-294-10x10.xlsx, max error percent=1.2
Saved gen/cubicinterp-12.377-3.2024-294-100x100.xlsx, max error percent=1.5
Saved gen/linearinterp-12.377-3.2024-294-10x10.xlsx, max error percent=1.2
Saved gen/linearinterp-12.377-3.2024-294-100x100.xlsx, max error percent=2.2
Saved gen/cubicinterp-12.377-3.2024-1060-10x10.xlsx, max error percent=0.73
Saved gen/cubicinterp-12.377-3.2024-1060-100x100.xlsx, max error percent=0.72
Saved gen/linearinterp-12.377-3.2024-1060-10x10.xlsx, max error percent=0.73
Saved gen/linearinterp-12.377-3.2024-1060-100x100.xlsx, max error percent=0.73
```


