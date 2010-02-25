#!/usr/bin/env python2.6

import glob, os, sys, getopt
from subprocess import *

def usage():
	print "Usage: edp.py -c -i -s source_pattern"
	print "  -i: to turn on interactive mode. Default off."
	print "  -c: to print the number of matches. Default on"

def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "cis:")
	except getopt.GetoptError, err:
		print str(err)
		usage()
		sys.exit(2)
	print_count = 0
	source_pattern = "*.exr"
	interactive_mode = 0
	for o, a in opts:
		if o == "-s":
			source_pattern = a
		elif o == "-c":
			print_count = 1
		elif o == "-i":
			interactive_mode = 1
	list=glob.glob(source_pattern)
	list_len= len(list)
	print "The number files are {0}.".format(list_len)
	if print_count == 1:
		sys.exit(2)
	if source_pattern == "*.exr":
		k = raw_input("Would you like to display all the files in this folder? Press 'y' to continue. ")
		if k != 'y': sys.exit(2)
		else: print "Exiting."
	count = 0;
	for source in list:
		count = count + 1
		Popen("exrdisplay {0} >& /dev/null &".format(source), 
			stdout=PIPE, shell=True).stdout.read()
		print "[{0}/{1}]: {2}".format(count, list_len, source)
		if interactive_mode == 1:
			raw_input("Press any key to continue:")

if __name__ == "__main__":
    main()
