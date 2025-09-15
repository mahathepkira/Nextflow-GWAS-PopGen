#!/usr/bin/env python
import sys

filename = sys.argv[1]
n        = int(sys.argv[2])
fd = open(filename)
line = fd.readline().strip().split(':')
order_attr = line[1].strip().split(',')
num  = int(line[2].strip()) + 1
x = order_attr[:n]
#for i in range(1, n + 1) :
#    x += str(i) + ','
x.append(str(num))
print(','.join(x))

