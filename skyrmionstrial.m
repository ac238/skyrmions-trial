

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
for x_ind = 1:N
    for y_ind = 1:N
        omega=((xx(x_ind,y_ind)+yy(x_ind,y_ind)*1i)/lambda)^n;
        m1_init(x_ind,y_ind)=4*real(omega)/((abs(omega))^2+4);
        m2_init(x_ind,y_ind)=4*imag(omega)/((abs(omega))^2+4);
    end
end

% plot
quiver(xx,yy,m1_init,m2_init)
axis([xlow xhigh ylow yhigh])

% dynamics
t=0;
dt=0.1;
while t<10
    for x_ind = 1:N
        for y_ind = 1:N
            % differential equations
            m1(x_ind,y_ind) = m1_init(x_ind,y_ind) - m2_init(x_ind,y_ind)*dt;
            m2(x_ind,y_ind) = m2_init(x_ind,y_ind) + m1_init(x_ind,y_ind)*dt;
        end
    end
    
    % plot
    quiver(xx,yy,m1,m2)
    drawnow
    
    % reset for new loop
    m1_init=m1;
    m2_init=m2;
    t=t+dt;
end


