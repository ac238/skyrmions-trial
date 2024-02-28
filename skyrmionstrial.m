

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
lambda=20;

omega = ((xx + yy*1i - z_0)/lambda).^n;
m1_init=4*real(omega)./((abs(omega)).^2+4);
m2_init=4*imag(omega)./((abs(omega)).^2+4);
m3_init=((abs(omega)).^2-4)./((abs(omega)).^2+4);

% plot
quiver(xx,yy,m1_init,m2_init)
axis([xlow xhigh ylow yhigh])

% dynamics
landau = 0;
zeeman = 1;

t=0;
t_final=1;
dt=0.001;
dx=(xhigh-xlow)/(N-1);
dy=(yhigh-ylow)/(N-1);
kappa = 1000;
gmuB = 10;
q = 1;
while t<t_final

    m1 = m1_init;
    m2 = m2_init;
    m3 = m3_init;

    % Zeeman term with RK4
    m1_k1 = gmuB*m2_init;
    m2_k1 = -gmuB*m1_init;
    m1_k2 = gmuB*(m2_init+dt/2*m2_k1);
    m2_k2 = -gmuB*(m1_init+dt/2*m1_k1);
    m1_k3 = gmuB*(m2_init+dt/2*m2_k2);
    m2_k3 = -gmuB*(m1_init+dt/2*m1_k2);
    m1_k4 = gmuB*(m2_init+dt*m2_k3);
    m2_k4 = -gmuB*(m1_init+dt*m1_k3);
    m1 = m1 + zeeman*dt/6*(m1_k1 + 2*m1_k2 + 2*m1_k3 + m1_k4);
    m2 = m2 + zeeman*dt/6*(m2_k1 + 2*m2_k2 + 2*m2_k3 + m2_k4);

    % L-L term with RK4 (only works for dx=dy) (doesn't conserve norm=1) (needs periodic boundaries)
    m1_k1 = kappa/(dx^2) * (m2_init(2:N-1,2:N-1).*(m3_init(1:N-2,2:N-1)+m3_init(3:N,2:N-1)+m3_init(2:N-1,1:N-2)+m3_init(2:N-1,3:N)-4*m3_init(2:N-1,2:N-1)) - m3_init(2:N-1,2:N-1).*(m2_init(1:N-2,2:N-1)+m2_init(3:N,2:N-1)+m2_init(2:N-1,1:N-2)+m2_init(2:N-1,3:N)-4*m2_init(2:N-1,2:N-1)));
    m2_k1 = kappa/(dx^2) * (m3_init(2:N-1,2:N-1).*(m1_init(1:N-2,2:N-1)+m1_init(3:N,2:N-1)+m1_init(2:N-1,1:N-2)+m1_init(2:N-1,3:N)-4*m1_init(2:N-1,2:N-1)) - m1_init(2:N-1,2:N-1).*(m3_init(1:N-2,2:N-1)+m3_init(3:N,2:N-1)+m3_init(2:N-1,1:N-2)+m3_init(2:N-1,3:N)-4*m3_init(2:N-1,2:N-1)));
    m3_k1 = kappa/(dx^2) * (m1_init(2:N-1,2:N-1).*(m2_init(1:N-2,2:N-1)+m2_init(3:N,2:N-1)+m2_init(2:N-1,1:N-2)+m2_init(2:N-1,3:N)-4*m2_init(2:N-1,2:N-1)) - m2_init(2:N-1,2:N-1).*(m1_init(1:N-2,2:N-1)+m1_init(3:N,2:N-1)+m1_init(2:N-1,1:N-2)+m1_init(2:N-1,3:N)-4*m1_init(2:N-1,2:N-1)));
    
    m1_k1_expand = zeros(N,N);
    m1_k1_expand(2:N-1,2:N-1)=m1_k1;
    m2_k1_expand = zeros(N,N);
    m2_k1_expand(2:N-1,2:N-1)=m2_k1;
    m3_k1_expand = zeros(N,N);
    m3_k1_expand(2:N-1,2:N-1)=m3_k1;
    m1_k2arg = m1_init+dt/2*m1_k1_expand;
    m2_k2arg = m2_init+dt/2*m2_k1_expand;
    m3_k2arg = m3_init+dt/2*m3_k1_expand;
    m1_k2 = kappa/(dx^2) * (m2_k2arg(2:N-1,2:N-1).*(m3_k2arg(1:N-2,2:N-1)+m3_k2arg(3:N,2:N-1)+m3_k2arg(2:N-1,1:N-2)+m3_k2arg(2:N-1,3:N)-4*m3_k2arg(2:N-1,2:N-1)) - m3_k2arg(2:N-1,2:N-1).*(m2_k2arg(1:N-2,2:N-1)+m2_k2arg(3:N,2:N-1)+m2_k2arg(2:N-1,1:N-2)+m2_k2arg(2:N-1,3:N)-4*m2_k2arg(2:N-1,2:N-1)));
    m2_k2 = kappa/(dx^2) * (m3_k2arg(2:N-1,2:N-1).*(m1_k2arg(1:N-2,2:N-1)+m1_k2arg(3:N,2:N-1)+m1_k2arg(2:N-1,1:N-2)+m1_k2arg(2:N-1,3:N)-4*m1_k2arg(2:N-1,2:N-1)) - m1_k2arg(2:N-1,2:N-1).*(m3_k2arg(1:N-2,2:N-1)+m3_k2arg(3:N,2:N-1)+m3_k2arg(2:N-1,1:N-2)+m3_k2arg(2:N-1,3:N)-4*m3_k2arg(2:N-1,2:N-1)));
    m3_k2 = kappa/(dx^2) * (m1_k2arg(2:N-1,2:N-1).*(m2_k2arg(1:N-2,2:N-1)+m2_k2arg(3:N,2:N-1)+m2_k2arg(2:N-1,1:N-2)+m2_k2arg(2:N-1,3:N)-4*m2_k2arg(2:N-1,2:N-1)) - m2_k2arg(2:N-1,2:N-1).*(m1_k2arg(1:N-2,2:N-1)+m1_k2arg(3:N,2:N-1)+m1_k2arg(2:N-1,1:N-2)+m1_k2arg(2:N-1,3:N)-4*m1_k2arg(2:N-1,2:N-1)));
    
    m1_k2_expand = zeros(N,N);
    m1_k2_expand(2:N-1,2:N-1)=m1_k2;
    m2_k2_expand = zeros(N,N);
    m2_k2_expand(2:N-1,2:N-1)=m2_k2;
    m3_k2_expand = zeros(N,N);
    m3_k2_expand(2:N-1,2:N-1)=m3_k2;
    m1_k3arg = m1_init+dt/2*m1_k2_expand;
    m2_k3arg = m2_init+dt/2*m2_k2_expand;
    m3_k3arg = m3_init+dt/2*m3_k2_expand;
    m1_k3 = kappa/(dx^2) * (m2_k3arg(2:N-1,2:N-1).*(m3_k3arg(1:N-2,2:N-1)+m3_k3arg(3:N,2:N-1)+m3_k3arg(2:N-1,1:N-2)+m3_k3arg(2:N-1,3:N)-4*m3_k3arg(2:N-1,2:N-1)) - m3_k3arg(2:N-1,2:N-1).*(m2_k3arg(1:N-2,2:N-1)+m2_k3arg(3:N,2:N-1)+m2_k3arg(2:N-1,1:N-2)+m2_k3arg(2:N-1,3:N)-4*m2_k3arg(2:N-1,2:N-1)));
    m2_k3 = kappa/(dx^2) * (m3_k3arg(2:N-1,2:N-1).*(m1_k3arg(1:N-2,2:N-1)+m1_k3arg(3:N,2:N-1)+m1_k3arg(2:N-1,1:N-2)+m1_k3arg(2:N-1,3:N)-4*m1_k3arg(2:N-1,2:N-1)) - m1_k3arg(2:N-1,2:N-1).*(m3_k3arg(1:N-2,2:N-1)+m3_k3arg(3:N,2:N-1)+m3_k3arg(2:N-1,1:N-2)+m3_k3arg(2:N-1,3:N)-4*m3_k3arg(2:N-1,2:N-1)));
    m3_k3 = kappa/(dx^2) * (m1_k3arg(2:N-1,2:N-1).*(m2_k3arg(1:N-2,2:N-1)+m2_k3arg(3:N,2:N-1)+m2_k3arg(2:N-1,1:N-2)+m2_k3arg(2:N-1,3:N)-4*m2_k3arg(2:N-1,2:N-1)) - m2_k3arg(2:N-1,2:N-1).*(m1_k3arg(1:N-2,2:N-1)+m1_k3arg(3:N,2:N-1)+m1_k3arg(2:N-1,1:N-2)+m1_k3arg(2:N-1,3:N)-4*m1_k3arg(2:N-1,2:N-1)));
    
    m1_k3_expand = zeros(N,N);
    m1_k3_expand(2:N-1,2:N-1)=m1_k3;
    m2_k3_expand = zeros(N,N);
    m2_k3_expand(2:N-1,2:N-1)=m2_k3;
    m3_k3_expand = zeros(N,N);
    m3_k3_expand(2:N-1,2:N-1)=m3_k3;
    m1_k4arg = m1_init+dt*m1_k3_expand;
    m2_k4arg = m2_init+dt*m2_k3_expand;
    m3_k4arg = m3_init+dt*m3_k3_expand;
    m1_k4 = kappa/(dx^2) * (m2_k4arg(2:N-1,2:N-1).*(m3_k4arg(1:N-2,2:N-1)+m3_k4arg(3:N,2:N-1)+m3_k4arg(2:N-1,1:N-2)+m3_k4arg(2:N-1,3:N)-4*m3_k4arg(2:N-1,2:N-1)) - m3_k4arg(2:N-1,2:N-1).*(m2_k4arg(1:N-2,2:N-1)+m2_k4arg(3:N,2:N-1)+m2_k4arg(2:N-1,1:N-2)+m2_k4arg(2:N-1,3:N)-4*m2_k4arg(2:N-1,2:N-1)));
    m2_k4 = kappa/(dx^2) * (m3_k4arg(2:N-1,2:N-1).*(m1_k4arg(1:N-2,2:N-1)+m1_k4arg(3:N,2:N-1)+m1_k4arg(2:N-1,1:N-2)+m1_k4arg(2:N-1,3:N)-4*m1_k4arg(2:N-1,2:N-1)) - m1_k4arg(2:N-1,2:N-1).*(m3_k4arg(1:N-2,2:N-1)+m3_k4arg(3:N,2:N-1)+m3_k4arg(2:N-1,1:N-2)+m3_k4arg(2:N-1,3:N)-4*m3_k4arg(2:N-1,2:N-1)));
    m3_k4 = kappa/(dx^2) * (m1_k4arg(2:N-1,2:N-1).*(m2_k4arg(1:N-2,2:N-1)+m2_k4arg(3:N,2:N-1)+m2_k4arg(2:N-1,1:N-2)+m2_k4arg(2:N-1,3:N)-4*m2_k4arg(2:N-1,2:N-1)) - m2_k4arg(2:N-1,2:N-1).*(m1_k4arg(1:N-2,2:N-1)+m1_k4arg(3:N,2:N-1)+m1_k4arg(2:N-1,1:N-2)+m1_k4arg(2:N-1,3:N)-4*m1_k4arg(2:N-1,2:N-1)));
    
    m1(2:N-1,2:N-1) = m1(2:N-1,2:N-1) + landau*dt/6*(m1_k1+2*m1_k2+2*m1_k3+m1_k4);
    m2(2:N-1,2:N-1) = m2(2:N-1,2:N-1) + landau*dt/6*(m2_k1+2*m2_k2+2*m2_k3+m2_k4);
    m3(2:N-1,2:N-1) = m3(2:N-1,2:N-1) + landau*dt/6*(m3_k1+2*m3_k2+2*m3_k3+m3_k4);

    % Pontryagin density
    m1_x=m1(mod(1:N,N)+1,1:N);
    m1_y=m1(1:N,mod(1:N,N)+1);
    m2_x=m2(mod(1:N,N)+1,1:N);
    m2_y=m2(1:N,mod(1:N,N)+1);
    m3_x=m3(mod(1:N,N)+1,1:N);
    m3_y=m3(1:N,mod(1:N,N)+1);
    m_dot_mx = m1.*m1_x + m2.*m2_x + m3.*m3_x;
    m_dot_my = m1.*m1_y + m2.*m2_y + m3.*m3_y;
    mx_dot_my = m1_x.*m1_y + m2_x.*m2_y + m3_x.*m3_y;
    triple = ( m1.*(m2_x.*m3_y) + m2.*(m3_x.*m1_y) + m3.*(m1_x.*m2_y) - m1.*(m3_x.*m2_y) - m3.*(m2_x.*m1_y) - m2.*(m1_x.*m3_y) );
    denom = 1 + m_dot_mx + m_dot_my + mx_dot_my;
    rho = 4*atan2(triple,denom)/(4*pi*dx*dy);

    
    % Coulomb term (not vectorized)
    %coulomb_int_x=zeros(N,N);
    %coulomb_int_y=zeros(N,N);
    %for j_x = 1:N
    %    for j_y = 1:N
    %        x_interact = (xx-xx(j_x,j_y)).*rho(j_x,j_y)./(((xx-xx(j_x,j_y)).^2.+(yy-yy(j_x,j_y)).^2)).^1.5*dx*dy;
    %        y_interact = (yy-yy(j_x,j_y)).*rho(j_x,j_y)./(((xx-xx(j_x,j_y)).^2.+(yy-yy(j_x,j_y)).^2)).^1.5*dx*dy;
    %        x_interact(j_x,j_y)=0;
    %        y_interact(j_x,j_y)=0;
    %        coulomb_int_x = coulomb_int_x + x_interact;
    %        coulomb_int_y = coulomb_int_y + y_interact;
    %    end
    %end


    % check conserved quantities
    Q_top=sum(sum(rho))*dx*dy  % topological charge
    S_z=sum(sum(m3));           % total spin in z-direction

    % check norm is preserved
    m_norm=m1.^2+m2.^2+m3.^2;
    
    % plot
    quiver(xx,yy,m1,m2) % full 2D vector field
    %quiver(xx(N/2,:),zeros(1,N),m1(N/2,:),m3(N/2,:)) % 1D slice
    hold on
    contour(xx,yy,rho,10) % color plot
    hold off
    axis([xlow xhigh ylow yhigh])
    drawnow
    
    % reset for new loop
    m1_init=m1;
    m2_init=m2;
    m3_init=m3;
    t=t+dt;
end


