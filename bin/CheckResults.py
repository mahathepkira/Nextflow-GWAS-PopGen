#!/usr/bin/env python
import sys, os

if __name__=='__main__' :
    filename = sys.argv[1]
    fd = open(filename)
    lines = fd.readlines()
    begin_read = 0
    correct = 0.0
    for line in lines :
        if line.strip() == '=== Error on training data ===' :
            begin_read = 1
        elif begin_read :
            if line.strip() != '' :
                if line.strip().split()[0] == 'Correctly' :
                    correct = '%.2f %s' % (float(line.strip().split()[-2]),'%')
    fd.close()
    print(filename, correct)
