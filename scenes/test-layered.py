#! /usr/bin/env python2.6

#
# Author: Joo-Haeng Lee
# Date: 2009-11-08
#
################################################################################
#
# Usage:
# - materials: 
#	(1) named_materials, 
#	(2) ["arbi"] arbitrary
#
################################################################################
execfile("utils/def.py")
execfile("materials/materials.py")
################################################################################

arg_scene = "test-ball"

arg_name_l = [
	"ior",
	"thickness",
	"material",
	"base_type", # matte, measured, metal, mirror, plastic, brdf_merl (bm)
	"base_roughness",
	"coat_roughness",
	"coat_name",
	"absorption_name",
	"light_nsamples",
	"metropolis_mp",
	"metropolis_ds",
	"mfnormal",
	"sampling_method",
	"baseonly",
	"configuration",
	"nbundles",
	"rotate_angle",
	"tir",
	"brdf_merl",
]

arg_name_short_l = {
	"ior":"ior",
	"thickness":"d",
	"material":"mat",
	"base_type":"bt",
	"base_roughness":"br",
	"coat_roughness":"cr",
	"coat_name":"cn",
	"absorption_name":"an",
	"light_nsamples":"ls",
	"metropolis_mp":"mp",
	"metropolis_ds":"ds",
	"mfnormal":"mf",
	"sampling_method":"smp",
	"baseonly":"bo",
	"configuration":"cfg",
	"nbundles":"nb",
	"rotate_angle":"ra",
	"tir":"tir",
	"brdf_merl":"bm",
}

def usage():
	print 'Usage: test-layered.py -a "arg.py"'

import getopt, sys, os
try: 
	opts, args = getopt.getopt(sys.argv[1:], "a:")
except getopt.GetoptError, err:
	print str(err)
	usage()
	sys.exit(2)
arg_file = 'arg_1.py'
for o, a in opts:
	if o == "-a":
		arg_file = a
print '>>> {0} will be used as a argument file.'.format(arg_file)
execfile(arg_file)

#
task = "layered"
cmd = "pbrt-" + task
scene_file = arg_scene + ".pbrt"

#
from subprocess import *
from datetime import datetime, timedelta

d_mat = d_named_materials

d={}

#
# LightSample.pbrt
#
chunk = "LightSample"
content = """\
"""
d[chunk] = content

#
# Metropolis.pbrt
#
chunk = "Metropolis"
content = """\
Renderer "metropolis" 
	"integer maxdepth" 7\
	"bool bidirectional" ["true"]\
	"bool dodirectseparately" ["true"]\
"""
d[chunk] = content

#
# Material.pbrt
#
chunk = "Material"
content = """\
Material "layered"\
"""
d[chunk] = content

#
# Film.pbrt
#
chunk = "Film"
content = """\
Film "image"\
"""
content += '\n\t' + '"integer xresolution" {0}'.format(arg_xresolution)
content += '\n\t' + '"integer yresolution" {0}'.format(arg_yresolution)
d[chunk] = content


#
# part_len: partial length required to compute the below
# 	a[(i/(len(b)*len(c)))%len(a)], b[(i/len(c))%len(b)], c[i%len(c)]
#
part_len=[]
for i in range(len(arg_name_l)):
	part_len.insert(i,1)
	for j in range(1+i, len(arg_name_l)):
		if arg_name_l[j] in arg_val_l:
			part_len[i] = part_len[i] * len(arg_val_l[arg_name_l[j]])
#print part_len

#
# Utility Functions
#
def f_exponent(r):
	if r >= 0.00001: 
		return 1/r
	else:
		return 10000

def formatTD(td):
	hours 	=	td.seconds // 3600
	minutes =	(td.seconds % 3600) // 60
	seconds =	td.seconds % 60
	return '%s days %s:%s:%s' % (td.days, hours, minutes, seconds)

def f_idx(i, j):
	return (i/part_len[j])%len(arg_val_l[arg_name_l[j]])

def f_pattern_name():
	str = ''
	str += '"{0}/{1}'.format(arg_scene, task)
	for n in arg_name_l: 
		if n in arg_val_l:
			if len(arg_val_l[n]) > 1:
				str += '-{0}_{1}'.format(arg_name_short_l[n], "*")
			else:
				str += '-{0}_{1}'.format(arg_name_short_l[n], arg_val_l[n][0])
	if arg_cmt != "":
		str += '-{0}'.format(arg_cmt)
	str += '.exr"'
	return str

def f_filename(argi):
	str = ''
	str += '"{0}/{1}'.format(arg_scene, task)
	for n in arg_name_l: 
		if n in arg_val_l:
			str += '-{0}_{1}'.format(arg_name_short_l[n], argi[n])
	if arg_cmt != "":
		str += '-{0}'.format(arg_cmt)
	str += '.exr"'
	return str

def f_filename_2(argi):
	return 	'"{0}/{0}-{1}-bo_{2}-tir_{3}-d_{4}-ior_{5}-mf_{6}-{7}-smp_{8}-mp_{9}-ds_{10}-{11}.exr"'.format(
		arg_scene, 
		task, 
		arg_baseonly, 
		arg_tir, 
		argi["thickness"], 
		argi["ior"], 
		arg_mfnormal, 
		argi["material"], 
		arg_smp,
		arg_mp,
		arg_ds,
		arg_cmt,
	)

#
# the total length of loop
#
loop_len = 1
for n in arg_name_l: 
	if n in arg_val_l:
		loop_len = loop_len * len(arg_val_l[n])
#print loop_len

#
# Set up arguments for the big loop
#
arg=[]
for i in range(loop_len):
	arg.append({})
	for j in range(len(arg_name_l)):
		if arg_name_l[j] in arg_val_l:
			arg[i][arg_name_l[j]] = arg_val_l[arg_name_l[j]][f_idx(i,j)]

#
# Main loop: do something interesting inside the loop
#
print '>>> Total {0} tasks.'.format(loop_len)

print ">>> Common Parameters:"
for n in arg_name_l: 
	if n in arg_val_l:
		if len(arg_val_l[n]) == 1:
			print "{0} = ".format(n), arg[i][n], ", ",
print
k = raw_input(">>> Would you like to continue? Press 'y' to continue: ")
if k != 'y': 
	print ">>> Exiting..."
	sys.exit(2)

for i in range(loop_len):
	print "[{0}/{1}]: ".format(i+1, loop_len),
	for n in arg_name_l: 
		if n in arg_val_l:
			if len(arg_val_l[n]) > 1:
				print "{0} = ".format(n), arg[i][n], ", ",
	print

k = raw_input(">>> Would you like to continue? Press 'y' to continue: ")
if k != 'y': 
	print ">>> Exiting..."
	sys.exit(2)
else: 
	print ">>> Starting..."

#
# Big Loop
#
for i in range(loop_len):
	argi = arg[i]
	#
	# Named Material
	#
	content_var = ""
	if argi["material"] == "arbi":
		content_var += '\n' + make_named_material(
			"mf_coat", 
			"plastic", 
			argi["coat_roughness"], 
			d_named_material_param[argi["coat_name"]+"_kd"], 
			d_named_material_param[argi["coat_name"]+"_ks"], 
			)
		content_var += '\n\t' + '"float ior" {0}'.format(argi["ior"])
		if argi["base_type"] == "measured": 
			content_var += '\n' + make_named_material_measured(
				"mf_base",
				argi["base_type"],
				argi["brdf_merl"]
				)
		else: 
			content_var += '\n' + make_named_material(
				"mf_base",
				argi["base_type"],
				argi["base_roughness"],
				arg_base_color
				)
	#
	# Material
	#
	chunk = "Material"
	content_var += '\n' + d[chunk]
	content_var += '\n' + d_mat[argi["material"]]
	content_var += '\n\t' + '"float exponent" {0}'.format(f_exponent(argi["sampling_method"]))
	content_var += '\n\t' + '"integer samplingmethod" {0}'.format(argi["sampling_method"])
	content_var += '\n\t' + '"integer configuration" {0}'.format(argi["configuration"])
	content_var += '\n\t' + '"integer nbundles" {0}'.format(argi["nbundles"])
	content_var += '\n\t' + '"float ior" {0}'.format(argi["ior"])
	content_var += '\n\t' + '"float thickness" {0}'.format(argi["thickness"])
	if "mfnormal" in arg_val_l:
		content_var += '\n\t' + '"bool mfnormal" "{0}"'.format(argi["mfnormal"])
	if "tir" in arg_val_l:
		content_var += '\n\t' + '"bool tir" "{0}"'.format(argi["tir"])
	content_var += '\n\t' + '"bool baseonly" "{0}"'.format(argi["baseonly"])
	content_var += '\n\t' + '"texture absorption" "a_coat_{0}"'.format(argi["absorption_name"])
	content_var = content_var.replace("\"","\\\"")
	Popen("echo \"{0}\" > {1}/{2}.pbrt".format(content_var, arg_scene, chunk),\
		stdout=PIPE, shell=True).stdout.read()
	#
	# LightSample
	#
	chunk = "LightSample"
	content_var = d[chunk]
	if "light_nsamples" in arg_val_l:
		content_var += '\t' + '"integer nsamples" {0}'.format(argi["light_nsamples"])
	content_var = content_var.replace("\"","\\\"")
	Popen("echo \"{0}\" > {1}/{2}.pbrt".format(content_var, "lights", chunk),\
		stdout=PIPE, shell=True).stdout.read()
	#
	# Metropolis
	#
	chunk = "Metropolis"
	content_var = d[chunk]
	content_var += '\n\t' + '"integer samplesperpixel" {0}'.format(argi["metropolis_mp"])
	content_var += '\n\t' + '"integer directsamples" {0}'.format(argi["metropolis_ds"])
	content_var = content_var.replace("\"","\\\"")
	Popen("echo \"{0}\" > {1}/{2}.pbrt".format(content_var, arg_scene, chunk),\
		stdout=PIPE, shell=True).stdout.read()
	#
	# Sampler
	#
	chunk = "Sampler"
	content_var = 'Sampler "bestcandidate" "integer pixelsamples" {0}'.format(argi["metropolis_ds"])
	content_var = content_var.replace("\"","\\\"")
	Popen("echo \"{0}\" > {1}/{2}.pbrt".format(content_var, arg_scene, chunk),\
	stdout=PIPE, shell=True).stdout.read()
	
	#
	# Film
	#
	chunk = "Film"
	content_var = ""
	content_var += d[chunk]
	filename = f_filename(argi)
	content_var += '\n\t' + '"string filename" ' + filename
	content_var = content_var.replace("\"","\\\"")
	#print content_var
	Popen("echo \"{0}\" > {1}/{2}.pbrt".format(content_var, arg_scene, chunk),\
		 stdout=PIPE, shell=True).stdout.read()
	#
	# Xform (Transform)
	#
	chunk = "Xform"
	content_var = ""
	if "rotate_angle" in arg_val_l:
		content_var += 'Rotate {0} 0 0 1'.format(argi["rotate_angle"])
	content_var = content_var.replace("\"","\\\"")
	#print content_var
	print Popen("echo \"{0}\" > {1}/{2}.pbrt".format(content_var, arg_scene, chunk),\
		 stdout=PIPE, shell=True).stdout.read()
	#
	# Execute
	#
	print '>>>[{0}/{1}] {2}'.format(i+1,loop_len,filename)
	if arg_skip_render == 0:
		if arg_overwrite == 1 or (not(os.path.exists(filename.replace("\"","")))):
			t1 = datetime.now()
			print ">>> Started:", 	t1.strftime("%A, %d. %B %Y %H:%M:%S")
			Popen("{0} {1}".format(cmd, scene_file), stdout=PIPE, shell=True).stdout.read()
			t2 = datetime.now()
			print ">>> Finished:", 	t2.strftime("%A, %d. %B %Y %H:%M:%S")
			print ">>> Duration:", 	formatTD(t2-t1)
		else:
			print ">>> Skipped not to overwrite."
	else:
		print ">>> Skipped: arg_skip_render == 1"

print ">>> Generated images: " + f_pattern_name()
k = raw_input(">>> Would you like to display generated images? Press 'y' to continue: ")
if k == 'y': 
	count = 0
	for i in range(loop_len):
		count = count + 1
		argi = arg[i]
		fname = f_filename(argi)
		Popen('exrdisplay {0} >& /dev/null &'.format(fname), stdout=PIPE, shell=True).stdout.read()
		print "[{0}/{1}]: {2}".format(count, loop_len, fname)
		raw_input("Press any key to continue:")
