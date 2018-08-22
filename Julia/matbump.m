% matbump.m
x = [0.0, 0.1, 0.2, 0.3, 0.35, 0.55, 0.65, 0.75];
y = [0.0, 0.01, 0.02, 0.03, 0.5, 0.51, 0.52, 0.53];
xx = linspace(0.0,0.75,400);
yy = spline(x,y,xx);
plot(x,y,'o',xx,yy,'-')
title("spline")
figure()
yy = pchip(x,y,xx);
plot(x,y,'o',xx,yy,'-')
title("pchip")
