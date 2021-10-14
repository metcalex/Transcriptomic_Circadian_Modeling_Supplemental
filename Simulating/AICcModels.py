import numpy as np

def cosine_1_term_model(t,offset, peak_time, peak_height):
	return peak_height*(1+np.cos((t-peak_time)/24*2*np.pi))+offset

def cosine_2_term_model(t,offset, peak_time_1, peak_height_1, peak_time_2, peak_height_2):
	return peak_height_1*(1+np.cos((t-peak_time_1)/24*2*np.pi)) + peak_height_2*(1+np.cos((t-peak_time_2)/24*4*np.pi)) + offset

def square_model(t,offset, peak_time, peak_height):
	t_new=(t-peak_time)%24
	return offset+peak_height*np.piecewise(t_new,
		[t_new<=12, t_new>12],
		[1,0])

def square_width_model(t, offset, peak_time, peak_height, width):
	t_new=(t-peak_time)%24
	return offset+peak_height*np.piecewise(t_new,
		[t_new<=width, t_new>width],
		[1,0])

def hat(t,offset, peak_time, peak_height, half_width):
	t_new=(t-(peak_time-half_width))%24
	return offset+np.piecewise(t_new,
		[t_new<=half_width, (t_new<=2*half_width) & (t_new>half_width), t_new>2*half_width],
		[lambda t_new: t_new/half_width*peak_height,lambda t_new: (2*half_width-t_new)/(half_width)*peak_height,0])


def fixed_decay(t,offset, cliff_time, peak_height, char_time):
	t_new=(t-cliff_time)%24
	return offset+np.piecewise(t_new,
		[t_new<=12, t_new>12],
		[lambda t_new: peak_height*(np.exp(-char_time*t_new)-np.exp(-12*char_time)), lambda t_new: peak_height*(1-np.exp(-char_time*(t_new-12)))])

def variable_decay(t, offset, cliff_time, peak_height, char_time, decay_time):
	t_new=(t-cliff_time)%24
	return offset+np.piecewise(t_new,
		[t_new<=decay_time, t_new> decay_time],
		[lambda t_new: peak_height*(np.exp(-char_time*t_new)), lambda t_new: peak_height*((1-np.exp(-char_time*decay_time/(24-decay_time)*(t_new-decay_time)))+np.exp(-char_time*decay_time))])

def constant_model(t,offset):
	return offset+t*0

def Kronecker_delta(t,offset,t_peak,height):
	return np.piecewise(t,
		[abs(t-t_peak)>=0.1, abs(t-t_peak)<0.1],[offset,height])

def hat_mult_model(t, peak_time1, peak_height1, peak_time2, peak_height2):
	t_new1=(t-(peak_time1-12))%24
	t_new2=(t-(peak_time2-12))%24
	return np.piecewise(t_new1,
		[t_new1<=12, (t_new1<=24) & (t_new1>12), t_new1>24],
		[lambda t_new1: t_new1/12*peak_height1,lambda t_new1: (24-t_new1)/12*peak_height1,0])*np.piecewise(t_new2,
			[t_new2<=12, (t_new2<=24) & (t_new2>12), t_new2>24],
			[lambda t_new2: t_new2/12*peak_height2,lambda t_new2: (24-t_new2)/12*peak_height2,0])

def cos_mult_model(t, peak_time1, peak_height1, peak_time2, peak_height2):
	t_new1=(t-(peak_time1-12))%24
	t_new2=(t-(peak_time2-12))%24
	return  peak_height1/2*(np.cos((t-peak_time1)/24*2*np.pi)+1)* peak_height2/2*(np.cos((t-peak_time2)/24*2*np.pi)+1)
