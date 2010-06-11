def make_texture_float_const(name, value):
	str = ''
	str += 'Texture "{0}" "float" "constant" "float value" {1}'.format(name, value)
	return str

def pbrt_vector(vec):
	str = '['
	for v in vec: str += '{0} '.format(v)
	str += ']'
	return str

bm_root = "/Users/joo-haenglee/Dev/BRDFs/"
bm_ext = ".binary"
def make_named_material_measured(name, type, bm_name):
	str = ''
	str += '\n' + 'MakeNamedMaterial "{0}" "string type" "{1}"'.format(name, type)
	str += '\n\t' + '"string filename" "{0}"'.format(bm_root + bm_name + bm_ext)
	return str

def make_named_material(name, type, roughness, kd=[0.1, 0.45, 0.35], ks=[0.7, 0.7, 0.7]):
	str = ''
	str += make_texture_float_const("tex_"+name,roughness)
	str += '\n' + 'MakeNamedMaterial "{0}" "string type" "{1}"'.format(name, type)
	if type == "plastic":
		str += '\n\t' + '"texture roughness" "{0}"'.format("tex_"+name)
		str += '\n\t' + '"color Kd" ' + pbrt_vector(kd)
		str += '\n\t' + '"color Ks" ' + pbrt_vector(ks)
	elif type == "matte":
		str += '\n\t' + '"float sigma" 0'
		str += '\n\t' + '"color Kd" ' + pbrt_vector(kd)
	elif type == "metal":
		str += '\n\t' + '"texture roughness" "{0}"'.format("tex_"+name)

	return str

named_materials = [
	"mwpf", "mwpg", "mwpr",
	"mfpf", "mfpg", "mfpr",
	"mgpf", "mgpg", "mgpr",
	"mrpf", "mrpg", "mrpr",
	"pg2fpf", "pg2fpg", "pg2fpr",
	"pg2gpf", "pg2gpg", "pg2gpr",
	"pg2rpf", "pg2rpg", "pg2rpr",
]

#
# This is usually used in plastic coat
#
d_named_material_param = {\
	# white
	"white_kd":[0.0, 0.0, 0.0],
	"white_ks":[1.0, 1.0, 1.0],
	# green2
	"green2_kd":[0.1, 0.45, 0.35],
	"green2_ks":[0.7, 0.7, 0.7],
	# blue2
	"blue2_kd":[0.1, 0.35, 0.45],
	"blue2_ks":[0.7, 0.7, 0.7],
}
d_named_materials = {\
	"arbi":"""\
	"string namedmaterial1" "mf_coat" 
	"string namedmaterial2" "mf_base"\
""",\
	"mwpf":"""\
	"string namedmaterial1" "plastic_coat_flat" 
	"string namedmaterial2" "matte_white"\
""",\
	"mwpg":"""\
	"string namedmaterial1" "plastic_coat_glossy" 
	"string namedmaterial2" "matte_white"\
""",\
	"mwpr":"""\
	"string namedmaterial1" "plastic_coat_rough" 
	"string namedmaterial2" "matte_white"\
""",\
	"mfpf":"""\
	"string namedmaterial1" "plastic_coat_flat" 
	"string namedmaterial2" "metal_flat"\
""",\
	"mfpg":"""\
	"string namedmaterial1" "plastic_coat_glossy" 
	"string namedmaterial2" "metal_flat"\
""",\
	"mfpr":"""\
	"string namedmaterial1" "plastic_coat_rough" 
	"string namedmaterial2" "metal_flat"\
""",\
	"mgpf":"""\
	"string namedmaterial1" "plastic_coat_flat" 
	"string namedmaterial2" "metal_glossy"\
""",\
	"mgpg":"""\
	"string namedmaterial1" "plastic_coat_glossy" 
	"string namedmaterial2" "metal_glossy"\
""",\
	"mgpr":"""\
	"string namedmaterial1" "plastic_coat_rough" 
	"string namedmaterial2" "metal_glossy"\
""",\
	"mrpf":"""\
	"string namedmaterial1" "plastic_coat_flat" 
	"string namedmaterial2" "metal_rough"\
""",\
	"mrpg":"""\
	"string namedmaterial1" "plastic_coat_glossy" 
	"string namedmaterial2" "metal_rough"\
""",\
	"mrpr":"""\
	"string namedmaterial1" "plastic_coat_rough" 
	"string namedmaterial2" "metal_rough"\
""",\
	"pg2fpf":"""\
	"string namedmaterial1" "plastic_coat_flat" 
	"string namedmaterial2" "plastic_green2_flat"\
""",\
	"pg2fpg":"""\
	"string namedmaterial1" "plastic_coat_glossy" 
	"string namedmaterial2" "plastic_green2_flat"\
""",\
	"pg2fpr":"""\
	"string namedmaterial1" "plastic_coat_rough" 
	"string namedmaterial2" "plastic_green2_flat"\
""",\
	"pg2gpf":"""\
	"string namedmaterial1" "plastic_coat_flat" 
	"string namedmaterial2" "plastic_green2_glossy"\
""",\
	"pg2gpg":"""\
	"string namedmaterial1" "plastic_coat_glossy" 
	"string namedmaterial2" "plastic_green2_glossy"\
""",\
	"pg2gpr":"""\
	"string namedmaterial1" "plastic_coat_rough" 
	"string namedmaterial2" "plastic_green2_glossy"\
""",\
	"pg2rpf":"""\
	"string namedmaterial1" "plastic_coat_flat" 
	"string namedmaterial2" "plastic_green2_rough"\
""",\
	"pg2rpg":"""\
	"string namedmaterial1" "plastic_coat_glossy" 
	"string namedmaterial2" "plastic_green2_rough"\
""",\
	"pg2rpr":"""\
	"string namedmaterial1" "plastic_coat_rough" 
	"string namedmaterial2" "plastic_green2_rough"\
"""
}

brdf_merl_list = [
"alum-bronze",
"alumina-oxide",
"aluminium",
"aventurnine",
"beige-fabric",
"black-fabric",
"black-obsidian",
"black-oxidized-steel",
"black-phenolic",
"black-soft-plastic",
"blue-acrylic",
"blue-fabric",
"blue-metallic-paint",
"blue-metallic-paint2",
"blue-rubber",
"brass",
"cherry-235",
"chrome-steel",
"chrome",
"colonial-maple-223",
"color-changing-paint1",
"color-changing-paint2",
"color-changing-paint3",
"dark-blue-paint",
"dark-red-paint",
"dark-specular-fabric",
"delrin",
"fruitwood-241",
"gold-metallic-paint",
"gold-metallic-paint2",
"gold-metallic-paint3",
"gold-paint",
"gray-plastic",
"grease-covered-steel",
"green-acrylic",
"green-fabric",
"green-latex",
"green-metallic-paint",
"green-metallic-paint2",
"green-plastic",
"hematite",
"ipswich-pine-221",
"light-brown-fabric",
"light-red-paint",
"maroon-plastic",
"natural-209",
"neoprene-rubber",
"nickel",
"nylon",
"orange-paint",
"pearl-paint",
"pickled-oak-260",
"pink-fabric",
"pink-fabric2",
"pink-felt",
"pink-jasper",
"pink-plastic",
"polyethylene",
"polyurethane-foam",
"pure-rubber",
"purple-paint",
"pvc",
"red-fabric",
"red-fabric2",
"red-metallic-paint",
"red-phenolic",
"red-plastic",
"red-specular-plastic",
"silicon-nitrade",
"silver-metallic-paint",
"silver-metallic-paint2",
"silver-paint",
"special-walnut-224",
"specular-black-phenolic",
"specular-blue-phenolic",
"specular-green-phenolic",
"specular-maroon-phenolic",
"specular-orange-phenolic",
"specular-red-phenolic",
"specular-violet-phenolic",
"specular-white-phenolic",
"specular-yellow-phenolic",
"ss440",
"steel",
"teflon",
"tungsten-carbide",
"two-layer-gold",
"two-layer-silver",
"violet-acrylic",
"violet-rubber",
"white-acrylic",
"white-diffuse-bball",
"white-fabric",
"white-fabric2",
"white-marble",
"white-paint",
"yellow-matte-plastic",
"yellow-paint",
"yellow-phenolic",
"yellow-plastic",
]
