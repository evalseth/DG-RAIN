NumKinRoutes = 1
print NumKinRoutes
x_0 = 2.5;
x_N = 2.5;
y_0 = 600.0;
#y_0 = 600*0.3048;
y_N = 0;
S0 = 0.0016;
z_0 = -S0*(y_0);

nval = 0.016825;
#nval = 0.025;
NumEl = 800;
NumNodes = NumEl+1;
hx = (x_N - x_0)/NumEl;
hy = (y_N - y_0)/NumEl;
x = [];
y = [];
z = [];
n = [];
x.append(x_0);
y.append(y_0);
z.append(z_0);
n.append(nval);

for i in range(1,NumEl+1):
	newx = x[i-1]+hx
	newy = y[i-1]+hy;
	newz = z[i-1] - 0.0016*hy;
	x.append(newx);
	y.append(newy);
	z.append(newz);
	n.append(nval);

print NumEl

for i in range(0,NumEl+1):
	print "{iteration} \t{xcor} \t {ycor} \t {zcor} \t\t {nval}".format(iteration= i, xcor=x[i], ycor=y[i], zcor=z[i], nval=n[i])

print " "
print "0 0 0"

#x_0 = 42.;
#x_N = 42;
#y_0 = 100;
#y_N = 0;
#z_0 = -0.25;
#nval = 0.3;
#NumEl = 25;
#NumNodes = NumEl+1;
#hx = (x_N - x_0)/NumEl;
#hy = (y_N - y_0)/NumEl;
#hz = 0.01;
#x = [];
#y = [];
#z = [];
#n = [];
#x.append(x_0);
#y.append(y_0);
#z.append(z_0);
#n.append(nval);
#
#for i in range(1,NumEl+1):
#	newx = x[i-1]+hx
#	newy = y[i-1]+hy;
#	newz = z[i-1]+hz;
#	x.append(newx);
#	y.append(newy);
#	z.append(newz);
#	n.append(nval);
#
#print NumEl
#
#for i in range(0,NumEl+1):
#	print "{iteration} \t{xcor} \t {ycor} \t {zcor} \t\t {nval}".format(iteration= i, xcor=x[i], ycor=y[i], zcor=z[i], nval=n[i])
#

