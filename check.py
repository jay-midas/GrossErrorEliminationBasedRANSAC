#!/usr/bin/python3

import numpy as np
import sys
d = '''
-0.458007
-1.0589
326.578
-0.531866
-0.713735
263.334
-0.00146522
-0.00317854
1
'''
d = '''
0.0447704 -0.131353 308.551
-0.124555 0.850661 42.9577
-0.000533284 -0.000154999 1.13586
'''
if __name__ == "__main__":
    x = d.split()
    H = np.array([
        [float(x[0]), float(x[1]), float(x[2])],
        [float(x[3]), float(x[4]), float(x[5])],
        [float(x[6]), float(x[7]), float(x[8])]
    ])
    print(H)
    l = np.array([[float(sys.argv[1])], [float(sys.argv[2])], [float(sys.argv[3])]])
    r = np.dot(H, l)
    print(r)