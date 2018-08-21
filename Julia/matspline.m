% spline.m
x = cumsum(rand(10,1));
x = (x-x(1))*pi/(x(10)-x(1));
y = sin(x);
xx = linspace(0, pi, 97);
yy = spline(x,y,xx);
yyy = sin(xx);
plot(x,y,'o',xx,yy,'.',xx,yyy,'-')