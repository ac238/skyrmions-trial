

% create bounds of graph
N=50;
xlow=-50;
ylow=-50;
xhigh=50;
yhigh=50;
[xx,yy]=meshgrid(linspace(xlow,xhigh,N),linspace(ylow,yhigh,N));

% initialize skyrmion according to QHMF 35
n=-1;
z_0=0;
lambda=10;
omega = ((xx + yy*1i - z_0)/lambda).^n;
m1_init=4*real(omega)./((abs(omega)).^2.+4);
m2_init=4*imag(omega)./((abs(omega)).^2.+4);
m3_init = sqrt(1-(m1_init).^2-(m2_init).^2);

% plot
quiver(xx,yy,m1_init,m2_init)
axis([xlow xhigh ylow yhigh])

% dynamics
landau = 0; % broken
zeeman = 1;

t=0;
t_final=1;
dt=0.001;
dx=(xhigh-xlow)/(N-1);
dy=(yhigh-ylow)/(N-1);
kappa = 0.01;
gmuB = 10;
while t<t_final

    m1 = m1_init;
    m2 = m2_init;
    m3 = sqrt(1 - (m1_init).^2 - (m2_init).^2);

    % Zeeman term (RK4)
    m1_k1 = zeeman*gmuB*m2_init;
    m2_k1 = -zeeman*gmuB*m1_init;
    m1_k2 = zeeman*gmuB*(m2_init+dt/2*m2_k1);
    m2_k2 = -zeeman*gmuB*(m1_init+dt/2*m1_k1);
    m1_k3 = zeeman*gmuB*(m2_init+dt/2*m2_k2);
    m2_k3 = -zeeman*gmuB*(m1_init+dt/2*m1_k2);
    m1_k4 = zeeman*gmuB*(m2_init+dt*m2_k3);
    m2_k4 = -zeeman*gmuB*(m1_init+dt*m1_k3);
    m1 = m1 + dt/6*(m1_k1 + 2*m1_k2 + 2*m1_k3 + m1_k4);
    m2 = m2 + dt/6*(m2_k1 + 2*m2_k2 + 2*m2_k3 + m2_k4);


    %for x_i = 1:N
        %for y_i = 1:N
            % Landau-Lifshitz term (broken)
            %m1(x_i,y_i) = m1(x_i,y_i) + dt*( landau*kappa*( m2_init(x_i,y_i)*(m3_init(mod(x_i-2,N)+1,y_i)-2*m3_init(x_i,y_i)+m3_init(mod(x_i,N)+1,y_i))/dx^2 + (m3_init(x_i,mod(y_i-2,N)+1)-2*m3_init(x_i,y_i)+m3_init(x_i,mod(y_i,N)+1))/dx^2) );
            %m1(x_i,y_i) = m1(x_i,y_i) - dt*( landau*kappa*( m3_init(x_i,y_i)*(m2_init(mod(x_i-2,N)+1,y_i)-2*m2_init(x_i,y_i)+m2_init(mod(x_i,N)+1,y_i))/dx^2 + (m2_init(x_i,mod(y_i-2,N)+1)-2*m2_init(x_i,y_i)+m2_init(x_i,mod(y_i,N)+1))/dx^2) );
            %m2(x_i,y_i) = m2(x_i,y_i) + dt*( landau*kappa*( m3_init(x_i,y_i)*(m1_init(mod(x_i-2,N)+1,y_i)-2*m1_init(x_i,y_i)+m1_init(mod(x_i,N)+1,y_i))/dx^2 + (m1_init(x_i,mod(y_i-2,N)+1)-2*m1_init(x_i,y_i)+m1_init(x_i,mod(y_i,N)+1))/dx^2) );
            %m2(x_i,y_i) = m2(x_i,y_i) - dt*( landau*kappa*( m1_init(x_i,y_i)*(m3_init(mod(x_i-2,N)+1,y_i)-2*m3_init(x_i,y_i)+m3_init(mod(x_i,N)+1,y_i))/dx^2 + (m3_init(x_i,mod(y_i-2,N)+1)-2*m3_init(x_i,y_i)+m3_init(x_i,mod(y_i,N)+1))/dx^2) );
            % try open boundary condition, not periodic
        %end
    %end

    % Pontryagin density
    m3 = sqrt(1 - (m1).^2 - (m2).^2);
    triple = m1(1:N-1,1:N-1).*m2(2:N,1:N-1).*m3(1:N-1,2:N) + m2(1:N-1,1:N-1).*m3(2:N,1:N-1).*m1(1:N-1,2:N) + m3(1:N-1,1:N-1).*m1(2:N,1:N-1).*m2(1:N-1,2:N) - m1(1:N-1,1:N-1).*m3(2:N,1:N-1).*m2(1:N-1,2:N) - m3(1:N-1,1:N-1).*m2(2:N,1:N-1).*m1(1:N-1,2:N) - m2(1:N-1,1:N-1).*m1(2:N,1:N-1).*m3(1:N-1,2:N) ;
    


    % topological charge of previous frame
    Q=sum(sum(triple))/(4*pi*dx*dy)
    
    % plot
    quiver(xx,yy,m1,m2)
    drawnow
    
    % reset for new loop
    m1_init=m1;
    m2_init=m2;
    m3_init = sqrt(1-(m1_init).^2-(m2_init).^2);
    t=t+dt;
end


