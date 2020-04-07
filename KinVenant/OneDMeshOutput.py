print "1D Mesh"
x_0 = 0.;
x_N = 200;
NumEl = 20;
NumNodes = NumEl+1;
hx = (x_N - x_0)/NumEl;
hz = 0.01;
x = [];
y = [];
z = [];
x.append(x_0);
y.append(0);
z.append(0);

for i in range(1,NumEl+1):
	newx = x[i-1]+hx
	newy = 0;
	newz = 0;
	#newz = z[i-1]+hz;
	x.append(newx);
	y.append(newy);
	z.append(newz);

print NumEl

for i in range(0,NumEl+1):
	print "{iteration} \t{xcor} \t {ycor} \t {zcor} \t\t {bval} \t\t 0.0 \t 0.0 \t 0.0".format(iteration= i, xcor=x[i], ycor=y[i], zcor=z[i], bval=1)
