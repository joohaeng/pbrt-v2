#!/usr/bin/env python2.6

import glob, os, sys, getopt

def usage():
	print "Usage: rename.py -s source_pattern -t target_pattern -x extension"

def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "s:t:x:")
	except getopt.GetoptError, err:
		print str(err)
		usage()
		sys.exit(2)
	for o, a in opts:
		if o == "-s":
			source_pattern = a
		elif o == "-t":
			target_pattern = a
		elif o == "-x":
			extension = a
	list=glob.glob("*" + source_pattern + "*." + extension)
	for source in list:
		target = source.replace(source_pattern, target_pattern)
		os.system("mv " + source + " " + target)

if __name__ == "__main__":
    main()
