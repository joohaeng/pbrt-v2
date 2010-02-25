arg_scene = "test-happy"

arg_val_l = {
	"ior":		[1.5], #range_float(1.5, 3.6, 0.5),
	"thickness":[0], #range_float(1.0, 5.1, 1.0),
	"material":["arbi"], #["arbi"],
	"base_type":["measured"], #["metal", "matte", "plastic", "measured"],
	"base_roughness":[1], #[0, 0.0001, 0.001, 0.01, 0.1], #[1.0],
	"coat_roughness":[0.01, 0.0001], #[.001, 0.01, 0.1],
	"coat_name":["white"], # white, blue2, green2
	"absorption_name":["blue","red_2", "white_1"], #["white","white_1", "red","red_2","green","green_2","blue","blue_2"], # white, red, red_2, green, green_2, blue, blue_2
	"metropolis_mp":[2], #[8,16,32,64,128,256,512,1024,2048],
	"metropolis_ds":[512],
	"mfnormal":			["true"],
	"sampling_method":	[0],
	"baseonly":			["false"],
	"configuration":	[3],
	"nbundles":			[1], #range_float(10, 101, 20), #[10, 100],
	"rotate_angle":		[0], #range_float(0,72,2),
	"tir":				["false"],
	"brdf_merl":		["gold-metallic-paint"], brdf_merl_list,
}

arg_cmt = "dbg_11-L_ufz" # dbg4: LayeredBxDF::Pdf(), 
arg_skip_render = 0 # 1 to skip render
arg_overwrite 	= 1 # 1 to overwrite
arg_xresolution = 256
arg_yresolution = 256
