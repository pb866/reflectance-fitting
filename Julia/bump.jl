using Interp
using PyPlot

function tbump()
    x = [0.0, 0.1, 0.2, 0.3, 0.35, 0.55, 0.65, 0.75]
    y = [0.0, 0.01, 0.02, 0.03, 0.5, 0.51, 0.52, 0.53]
    cs = CubicSpline(x,y)
    figure()
    v = range(0.0,stop=0.75,length=100);
    yy = [interp(cs,vv) for vv in v]
    plot(x,y,"o",v,yy)
end

tbump()
