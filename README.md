# multi-ridges-detection
This is the codes for multiple curves extraction algorithm
- The usage of MultiCurveExt can be found in 'MultiCurveExt_demo.m' in the 2-ridges and 3-ridges folders, respectively.
- The original C code are 'CurveMultiExt_init_2curves.c' in the 2-ridges folder and 'CurveMultiExt_init.c' in the 3-ridges folder.

# Iterative warping
- The main algorithm, Warping, is in 'iterWarping.m'
- Contains another curve extraction version: 'CurveExt_withFund.c'. This is the fundamental-informed version  of MultiCurveExt.
