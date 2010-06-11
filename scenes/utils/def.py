def range_float(start, end, step):
	v = start
	c=0
	l=[]
	while v < end:
		l.append(v)
		v = v+step
	return l
