from math import cos
from math import pi

print "1D Mesh"
x_0 = 0.;
x_N = 1000;
NumEl = 200;
NumNodes = NumEl+1;
h = (x_N - x_0)/NumEl;
x = [];
y = [];
z = [];
b = [];
x.append(x_0);
y.append(0);
b.append(5);

for i in range(1,NumEl+1):
	newx = x[i-1]+h
	newy = 0;
	x.append(newx);
	y.append(newy);
	newb = 0;
	if newx <= 100:
		newb = 5
	elif newx < 400:
		newb = 5-(5-3.587)*cos((newx-250)/300*pi)**2
	else:
		newb = 5
	
	b.append(newb);

print NumNodes

for i in range(0,NumEl+1):
	#newz  = -max(0,0.2-0.05*pow(x[i]-10,2));
	#a = 1;
	#h0 = 0.5;
	#L = 4;
	#newz = -h0*((x[i]-L/2)*(x[i]-L/2) - 1);
	#z.append(newz);
	z.append(0);
	print "{iteration} \t{xcor} \t {ycor} \t {zcor} \t\t {bval} \t\t 0.0 \t 0.0 \t 0.0".format(iteration= i, xcor=x[i], ycor=y[i], zcor=z[i], bval=b[i])
