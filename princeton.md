# Introduction
Profile used in Princeton beam experiment is common reference case for ANCF implementations e.g. [https://www.sciencedirect.com/science/article/pii/S0997753822002546].
In this context various interpolations are compared.

# Results from section-properties
Source is [princeton.py](./princeton.py) and results are summarized below
rectangle of 12.377 x 3.2024 is used (too many digits, but this is how they are usually given)

|ec_in_h|Number of nodes|Iw 10E-18 m^6|It 10E-12 m^4|
|--|---------------|-----|------|
|10|24|334|118|
|20|57|335|116|
|30|114|328|114|
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


# Results from interpolation
