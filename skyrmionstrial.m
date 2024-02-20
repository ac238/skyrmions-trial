

% create bounds of graph
N=50;
xlow=-5;
ylow=-5;
xhigh=5;
yhigh=5;
[xx,yy]=meshgrid(linspace(xlow,xhigh,N),linspace(ylow,yhigh,N));

% initialize skyrmion according to QHMF 35
n=4;
z_0=0;
lambda=2;
for x_index = 1:N
    for y_index = 1:N
        omega=((xx(x_index,y_index)+yy(x_index,y_index)*1i)/lambda)^n;
        m1_init(x_index,y_index)=4*real(omega)/((abs(omega))^2+4);
        m2_init(x_index,y_index)=4*imag(omega)/((abs(omega))^2+4);
    end
end

% plot
quiver(xx,yy,m1_init,m2_init)
axis([xlow xhigh ylow yhigh])
