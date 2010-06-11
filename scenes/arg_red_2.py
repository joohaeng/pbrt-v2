arg_val_l = {
	"ior":		[1.5, 3.5], #range_float(1.5, 3.6, 0.5),
	"thickness":[0.1, 0.3, 1.0], #range_float(1.0, 5.1, 1.0),
	"material":["arbi"],
	"base_type":["metal"], #["metal", "matte", "plastic"],
	"base_roughness":[0.0001], #[1, 0.1, 0.01, 0.001, 0.0001], #[0, 0.0001, 0.001, 0.01, 0.1], #[1.0],
	"coat_roughness":[0.01, 0.0001],
	"coat_name":["white"], # white, blue2, green2
	"absorption_name":["red_2"], #["white","white_1", "red","red_2","green","green_2","blue","blue_2"], # white, red, red_2, green, green_2, blue, blue_2
	"metropolis_mp":[2], #[8,16,32,64,128,256,512,1024,2048],
	"metropolis_ds":[2048],
	"mfnormal":			["true"],
	"sampling_method":	[0],
	"baseonly":			["false"],
	"configuration":	[3],
	"nbundles":			[1], #range_float(10, 101, 20), #[10, 100],
	"rotate_angle":		[0], #range_float(0,72,2),
	"tir":				["false"],
}

arg_cmt = "dbg_12-L_ufz" # 12: G(wor,wir,wh)
arg_skip_render = 0 # 1 to skip render
arg_overwrite 	= 0 # 1 to overwrite
arg_xresolution = 256
arg_yresolution = 256

arg_tir = "off"
