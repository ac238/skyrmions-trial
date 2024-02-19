
N=50; % num of pts in each direction

% create bounds of graph
xlow=0;
ylow=0;
xhigh=10;
yhigh=10;
[xx,yy]=meshgrid(linspace(xlow,xhigh,N),linspace(ylow,yhigh,N));

% apply function of position
uu=sin(xx);
vv=5*cos(yy);

% plot
quiver(xx,yy,uu,vv)
axis([0 10 0 10])