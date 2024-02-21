

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
for x_i = 1:N
    for y_i = 1:N
        omega=((xx(x_i,y_i)+yy(x_i,y_i)*1i)/lambda)^n;
        m1_init(x_i,y_i)=4*real(omega)/((abs(omega))^2+4);
        m2_init(x_i,y_i)=4*imag(omega)/((abs(omega))^2+4);
        m3_init(x_i,y_i) = sqrt(1-(m1_init(x_i,y_i))^2-(m2_init(x_i,y_i))^2);
    end
end

% plot
quiver(xx,yy,m1_init,m2_init)
axis([xlow xhigh ylow yhigh])

% dynamics
landau = 0; % broken
zeeman = 1;

t=0;
t_final=1;
dt=0.001;
dx=(xhigh-xlow)/N;
dy=(yhigh-ylow)/N;
kappa = 0.01;
gmuB = 10;
while t<t_final
    for x_i = 1:N
        for y_i = 1:N
            m1(x_i,y_i) = m1_init(x_i,y_i);
            m2(x_i,y_i) = m2_init(x_i,y_i);

            % Landau-Lifshitz term (broken)
            m1(x_i,y_i) = m1(x_i,y_i) + dt*( landau*kappa*( m2_init(x_i,y_i)*(m3_init(mod(x_i-2,N)+1,y_i)-2*m3_init(x_i,y_i)+m3_init(mod(x_i,N)+1,y_i))/dx^2 + (m3_init(x_i,mod(y_i-2,N)+1)-2*m3_init(x_i,y_i)+m3_init(x_i,mod(y_i,N)+1))/dx^2) );
            m1(x_i,y_i) = m1(x_i,y_i) - dt*( landau*kappa*( m3_init(x_i,y_i)*(m2_init(mod(x_i-2,N)+1,y_i)-2*m2_init(x_i,y_i)+m2_init(mod(x_i,N)+1,y_i))/dx^2 + (m2_init(x_i,mod(y_i-2,N)+1)-2*m2_init(x_i,y_i)+m2_init(x_i,mod(y_i,N)+1))/dx^2) );
            m2(x_i,y_i) = m2(x_i,y_i) + dt*( landau*kappa*( m3_init(x_i,y_i)*(m1_init(mod(x_i-2,N)+1,y_i)-2*m1_init(x_i,y_i)+m1_init(mod(x_i,N)+1,y_i))/dx^2 + (m1_init(x_i,mod(y_i-2,N)+1)-2*m1_init(x_i,y_i)+m1_init(x_i,mod(y_i,N)+1))/dx^2) );
            m2(x_i,y_i) = m2(x_i,y_i) - dt*( landau*kappa*( m1_init(x_i,y_i)*(m3_init(mod(x_i-2,N)+1,y_i)-2*m3_init(x_i,y_i)+m3_init(mod(x_i,N)+1,y_i))/dx^2 + (m3_init(x_i,mod(y_i-2,N)+1)-2*m3_init(x_i,y_i)+m3_init(x_i,mod(y_i,N)+1))/dx^2) );
            
            % Zeeman term
            m1_Heun = m1(x_i,y_i) + dt*zeeman*gmuB*m2_init(x_i,y_i);
            m2_Heun = m2(x_i,y_i) - dt*zeeman*gmuB*m1_init(x_i,y_i);
            m1(x_i,y_i) = m1(x_i,y_i) + dt*zeeman*gmuB/2*(m2_init(x_i,y_i)+m2_Heun);
            m2(x_i,y_i) = m2(x_i,y_i) - dt*zeeman*gmuB/2*(m1_init(x_i,y_i)+m1_Heun);

            % Pontryagin density (triple product, not solid angle)
            mvec = [m1_init(x_i,y_i) m2_init(x_i,y_i) m3_init(x_i,y_i)];
            mvec_x = [m1_init(mod(x_i,N)+1,y_i) m2_init(mod(x_i,N)+1,y_i) m3_init(mod(x_i,N)+1,y_i)];
            mvec_y = [m1_init(x_i,mod(y_i,N)+1) m2_init(x_i,mod(y_i,N)+1) m3_init(x_i,mod(y_i,N)+1)];
            triple(x_i,y_i) = dot(mvec,cross(mvec_x,mvec_y));
            denom = norm(mvec)*norm(mvec_x)*norm(mvec_y) + dot(mvec,mvec_x)*norm(mvec_y) + dot(mvec,mvec_y)*norm(mvec_x) + dot(mvec_x,mvec_y)*norm(mvec);
            Solid(x_i,y_i) = 1/(4*pi*dx*dy)*4*atan2(triple(x_i,y_i),denom);
        end
    end

    % topological charge of previous frame
    Q=0;
    for x_i = 1:N
        for y_i = 1:N
            Q = Q + triple(x_i,y_i)/dx/dy;
        end
    end
    Q
    
    % plot
    quiver(xx,yy,m1,m2)
    drawnow
    
    % reset for new loop
    m1_init=m1;
    m2_init=m2;
    for x_i = 1:N
        for y_i = 1:N
            m3_init(x_i,y_i) = sqrt(1-(m1_init(x_i,y_i))^2-(m2_init(x_i,y_i))^2);
        end
    end
    t=t+dt;
end


