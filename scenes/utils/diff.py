#! /usr/bin/env python2.6

#
# Author: Joo-Haeng Lee
# Date: 2009-11-08
#
from subprocess import *

ls_list = 1, 4, 16, 64
ds_list = 1, 4, 16, 64

# cmd
cmd = "exrdiff"
file_o = "diff-cfg_3-ls_{0}-ds_{1}" 						
file_1 = "layered-ior_1.5-d_3-mat_arbi-bt_metal-br_0.01-cr_0-cn_white-an_blue-ls_64-mp_2-ds_64-smp_0-bo_false-cfg_3-nb_10-ra_0-dbg_9-L_ufz.exr"
file_2 = "layered-ior_1.5-d_3-mat_arbi-bt_metal-br_0.01-cr_0-cn_white-an_blue-ls_{0}-mp_2-ds_{1}-smp_0-bo_false-cfg_3-nb_10-ra_0-dbg_9-L_ufz.exr"

for ls in ls_list:
	for ds in ds_list:
		#
		# diff
		#
		file_out = file_o.format(ls,ds)
		arg = "-o {0}.exr {1} {2}".format(file_out, file_1, file_2.format(ls,ds))
		result = Popen("{0} {1}".format(cmd, arg), stdout=PIPE, shell=True).stdout.read()
		print result
		result = Popen("echo \"{0}\" > {1}.txt".format(result, file_out), stdout=PIPE, shell=True).stdout.read()

if raw_input("Press 'y' to display all the diff images: ") <> 'y':
	print "Exit"
else:
	for ls in ls_list:
		for ds in ds_list:
			file_out = file_o.format(ls,ds)
			Popen("exrdisplay -1 {0}.exr".format(file_out), shell=True).stdout
