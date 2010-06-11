#!/usr/bin/env python2.6

#
# NOTE: initially, tried to loop over (a[i], b[j], c[k]) tripple.
#
#a=[1,2,3,4,5]
#b=['a','b','c','d']
#c=['!','@','#']
#for i in range(len(a)*len(b)*len(c)):
#	print a[(i/(len(b)*len(c)))%len(a)], b[(i/len(c))%len(b)], c[i%len(c)]

# TODO
# 1. execption handling when arg_name_l does not match arg_val_l
#

import glob, os, sys, getopt
from subprocess import *
execfile("def.py")

arg_name_l=[
	"a",
	"b",
	"c",
	"d"
]

arg_val_l={
	"a":[1,2,3],
	"b":['x','y','z'],
	"c":['@','#'],
	"d":range_float(0.0, 0.4, 0.2)
}

#
# part_len: partial length required to compute the below
# 	a[(i/(len(b)*len(c)))%len(a)], b[(i/len(c))%len(b)], c[i%len(c)]
#
part_len=[]
for i in range(len(arg_name_l)):
	part_len.insert(i,1)
	for j in range(1+i, len(arg_name_l)):
		part_len[i] = part_len[i] * len(arg_val_l[arg_name_l[j]])
#print part_len

#
# Utility Functions
#
def f_idx(i, j):
	return (i/part_len[j])%len(arg_val_l[arg_name_l[j]])
#
# the total length of loop
#
loop_len = 1
for n in arg_name_l: loop_len = loop_len * len(arg_val_l[n])
#print loop_len

#
# Set up arguments for the big loop
#
arg=[]
for i in range(loop_len):
	arg.append({})
	for j in range(len(arg_name_l)):
		arg[i][arg_name_l[j]] = arg_val_l[arg_name_l[j]][f_idx(i,j)]
#
# Main loop: do something interesting inside the loop
#
for i in range(loop_len):
	print "{0}: ".format(i),
	for n in arg_name_l: print "{0} = ".format(n), arg[i][n], ", ",
	print
