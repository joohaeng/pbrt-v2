#!/usr/bin/env python2.6

#
# Joo-Haeng Lee ( joohaeng at gmail dot com )
# 
# 2009-12-10
#
# Copyright 2ME
#
import glob, os, sys, getopt

def usage():
	print "Usage: exrtotiff.py -s source_pattern"
	print "- (ex) exrtotiff.py -s "

def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "s:")
	except getopt.GetoptError, err:
		print str(err)
		usage()
		sys.exit(2)
	source_pattern = ''
	for o, a in opts:
		if o == "-s":
			source_pattern = a
	if source_pattern == '':
		usage()
		sys.exit(2)
	list=glob.glob("*" + source_pattern + "*.exr")
	for source in list:
		target = source.replace(".exr", ".tif")
		print ">>> Converting {0} to {1}....".format(source, target),
		os.system("exrtotiff " + source + " " + target)
		print "Done!"

if __name__ == "__main__":
    main()
