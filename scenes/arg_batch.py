arg_val_l = {
	"ior":		[1.5, 3.0], #range_float(1.5, 3.6, 0.5),
	"thickness":[1.0, 3.0], #range_float(1.0, 5.1, 1.0),
	"material":["arbi"],
	"base_type":["metal"], #["metal", "matte", "plastic"],
	"base_roughness":[0, 0.001, 0.01, 0.1], #[0, 0.0001, 0.001, 0.01, 0.1], #[1.0],
	"coat_roughness":[0, 0.001, 0.01, 0.1], #[1.0],
	"coat_name":["white"], # white, blue2, green2
	"absorption_name":["blue"], #["white","red","red_2","green","green_2","blue","blue_2"], # white, red, red_2, green, green_2, blue, blue_2
	"metropolis_mp":[16], #[8,16,32,64,128,256,512,1024,2048],
	"metropolis_ds":[256],
	"mfnormal":			["true"],
	"sampling_method":	[0],
	"baseonly":			["false"],
	"configuration":	[1],
	"nbundles":			[1], #range_float(10, 101, 20), #[10, 100],
}

arg_cmt = "dbg_4-L_ufz" # dbg4: LayeredBxDF::Pdf(), 
arg_skip_render = 0 # 1 to skip render
arg_overwrite 	= 1 # 1 to overwrite
arg_xresolution = 1024
arg_yresolution = 1024

arg_tir = "off"
