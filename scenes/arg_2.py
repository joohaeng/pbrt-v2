arg_val_l = {
	"ior":		[1.5], #range_float(1.5, 3.6, 0.5),
	"thickness":[3], #range_float(1.0, 5.1, 1.0),
	"material":["arbi"],
	"base_type":["metal"], #["metal", "matte", "plastic"],
	"base_roughness":[0], #[1, 0.1, 0.01, 0.001, 0.0001], #[0, 0.0001, 0.001, 0.01, 0.1], #[1.0],
	"coat_roughness":[0],
	"coat_name":["white"], # white, blue2, green2
	"absorption_name":["red_2"], #["white","white_1", "red","red_2","green","green_2","blue","blue_2"], # white, red, red_2, green, green_2, blue, blue_2
	"light_nsamples":[16], #[8,16,32,64,128,256,512,1024,2048],
	"metropolis_mp":[2], #[8,16,32,64,128,256,512,1024,2048],
	"metropolis_ds":[16],
	"mfnormal":			["true"],
	"sampling_method":	[0,1,2,3],
	"baseonly":			["true"],
	"configuration":	[4],
	"nbundles":			[32], #range_float(10, 101, 20), #[10, 100],
	"tir":				["false"],
}

arg_cmt = "dbg_17-L_ufz" # 14: seeding again in RNG
arg_skip_render = 0 # 1 to skip render
arg_overwrite 	= 1 # 1 to overwrite
arg_xresolution = 256
arg_yresolution = 256
arg_base_color = [1,1,1] #[0.5, 0.01, 0.1] #red_2 [0.5, 0.01, 0.1]
