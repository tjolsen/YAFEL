#! /usr/bin/python

import sys

if(len(sys.argv) < 4):
    print "Must specify number of truss sections, length of section, and output filename"
    print "Usage: ./trussGen <N> <len> <fname>"
    sys.exit()

N = int(sys.argv[1])
L = float(sys.argv[2])
fname = sys.argv[3]

if(N < 1):
    print "Must have at least 1 section"
    sys.exit()

if(L <= 0):
    print "Must have positive segment length"
    sys.exit()

print "You have specified %d sections of length %f" % (N, L)
print "The gmsh file will be written to %s" % fname

fd = open(fname, 'w')

ptcount = 1
line = 'Point(1) = {0,0,0,'+str(10*L)+'};\n'
fd.write(line)

for n in range(N):
    line = 'Point('+str(ptcount+1)+') = {'+str((n+1)*L)+','+str(0)+',0,'+str(10*L)+'};\n'
    fd.write(line)
    line = 'Point('+str(ptcount+2)+') = {'+str((n+1)*L)+','+str(L)+',0,'+str(10*L)+'};\n'
    fd.write(line)
    ptcount = ptcount + 2

line = 'Point('+str(ptcount+1)+') = {'+str((N+1)*L)+','+str(0)+',0,'+str(10*L)+'};\n'
fd.write(line)


line = 'Line(1) = {1,2};\n'
fd.write(line)
line = 'Line(2) = {1,3};\n'
fd.write(line)
line = 'Line(3) = {2,3};\n'
fd.write(line)

linecount = 3
for n in range(N-1):
    line = 'Line('+str(linecount+1)+') = {'+str(2*(n+1))+','+str(2*(n+2))+'};\n'
    fd.write(line)
    line = 'Line('+str(linecount+2)+') = {'+str(2*(n+1))+','+str(2*(n+2)+1)+'};\n'
    fd.write(line)
    line = 'Line('+str(linecount+3)+') = {'+str(2*(n+1)+1)+','+str(2*(n+2)+1)+'};\n'
    fd.write(line)
    line = 'Line('+str(linecount+4)+') = {'+str(2*(n+2))+','+str(2*(n+2)+1)+'};\n'
    fd.write(line)
    linecount = linecount + 4


line = 'Line('+str(linecount+1)+') = {'+str(ptcount-1)+','+str(ptcount+1)+'};\n'
fd.write(line)
line = 'Line('+str(linecount+2)+') = {'+str(ptcount)+','+str(ptcount+1)+'};\n'
fd.write(line)

fd.write('Physical Point(1) = {1};\n')
fd.write('Physical Point(2) = {'+str(ptcount+1)+'};\n')

line='{'
for i in range(1,linecount+2):
    line = line + str(i) + ','

line = line + str(linecount+2) + '};'

fd.write('Physical Line(3) = ' + line)

fd.close()
