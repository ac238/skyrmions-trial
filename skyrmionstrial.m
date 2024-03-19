

% create bounds of graph
N=70;
dx=1;
dy=1;
xlow=-(N+1)/2;
ylow=-(N+1)/2;
xhigh=(N-1)/2;
yhigh=(N-1)/2;
[yy,xx]=meshgrid(linspace(xlow,xhigh,N),linspace(ylow,yhigh,N)); %xx and yy are swapped to correspond w/ indices
axis([xlow xhigh ylow yhigh])

% initialize skyrmion according to QHMF 35
n=-2;
z_0=0;
lambda=20;

omega = ((xx + yy*1i - z_0)/lambda).^n;
m_init(:,:,1)=4*real(omega)./((abs(omega)).^2+4);
m_init(:,:,2)=4*imag(omega)./((abs(omega)).^2+4);
m_init(:,:,3)=((abs(omega)).^2-4)./((abs(omega)).^2+4);


% initialize Coulomb distance matrix
for i = 1:N
    for j = 1:N
        for k = 1:3
            dist_x(:,:,k,i,j) = (xx-xx(i,j))./(((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^1.5);
            dist_x(i,j,k,i,j) = 0;
            dist_y(:,:,k,i,j) = (yy-yy(i,j))./(((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^1.5);
            dist_y(i,j,k,i,j) = 0;
        end
    end
end


% constants
B_field=1.1;
E_field=1.1;
q_electron=-3.1;
mass_electron=1.1;
rho_stiff=1.1;
dielectric=1.1;
g_factor=2.002;
casimir=0.5;
nu_level=1.0;

kappa = -rho_stiff*q_electron/(2*pi*casimir*nu_level*B_field);
gmuB = g_factor*B_field*q_electron/(2*mass_electron);
q = q_electron^3/(8*pi^2*casimir*nu_level*B_field*dielectric);
El_amp = -q_electron^2*E_field/(8*pi^2*casimir*nu_level*B_field);


% dynamics
landau = 1; % change to 0/1 to turn off/on effects
zeeman = 1;
coulomb = 1;
electric = 0;

t=0;
t_final=10;
dt=0.001;
t_ind=1;
El_freq = 10;
Q_top = -n;
Q_top_list = []; % store values for final plot
E_B_list = [];
E_LL_list = [];
E_C_list = [];
E_El_list = [];
while t<t_final

    m = m_init;

    % Zeeman term with RK4
    m1_k1 = gmuB*m(:,:,2);
    m2_k1 = -gmuB*m(:,:,1);
    m1_k2 = gmuB*(m(:,:,2)+dt/2*m2_k1);
    m2_k2 = -gmuB*(m(:,:,1)+dt/2*m1_k1);
    m1_k3 = gmuB*(m(:,:,2)+dt/2*m2_k2);
    m2_k3 = -gmuB*(m(:,:,1)+dt/2*m1_k2);
    m1_k4 = gmuB*(m(:,:,2)+dt*m2_k3);
    m2_k4 = -gmuB*(m(:,:,1)+dt*m1_k3);

    m(:,:,1) = m(:,:,1) + zeeman*dt/6*(m1_k1 + 2*m1_k2 + 2*m1_k3 + m1_k4);
    m(:,:,2) = m(:,:,2) + zeeman*dt/6*(m2_k1 + 2*m2_k2 + 2*m2_k3 + m2_k4);

    % Electric field term with RK4
    m_k1 = (deriv_center_y(m))/(2*dy);
    m_k2 = (deriv_center_y(m) + dt/2*(deriv_center_y(m_k1)))/(2*dy);
    m_k3 = (deriv_center_y(m) + dt/2*(deriv_center_y(m_k2)))/(2*dy);
    m_k4 = (deriv_center_y(m) + dt*(deriv_center_y(m_k3)))/(2*dy);

    El = El_amp*cos(El_freq*t);
    m = m + electric*El*dt/6*(m_k1 + 2*m_k2 + 2*m_k3 + m_k4);
    m = m./(sqrt(sum(m.^2,3))); % Renormalize

    % L-L term with RK4 (only works for dx=dy) (doesn't conserve norm=1)
    m_k1 = kappa/(dx^2) * (stiffness(m));
    m_k2arg = m+dt/2*m_k1;
    m_k2 = kappa/(dx^2) * (stiffness(m_k2arg));
    m_k3arg = m+dt/2*m_k2;
    m_k3 = kappa/(dx^2) * (stiffness(m_k3arg));
    m_k4arg = m+dt*m_k3;
    m_k4 = kappa/(dx^2) * (stiffness(m_k4arg));

    m = m + landau*dt/6*(m_k1+2*m_k2+2*m_k3+m_k4);
    m = m./(sqrt(sum(m.^2,3))); % Renormalize

    % Pontryagin density
    m_x=m(mod(1:N,N)+1,1:N,:);
    m_y=m(1:N,mod(1:N,N)+1,:);
    triple = sum( m(:,:,mod(1:3,3)+1).*m_x(:,:,mod(2:4,3)+1).*m_y(:,:,mod(3:5,3)+1) - m(:,:,mod(1:3,3)+1).*m_y(:,:,mod(2:4,3)+1).*m_x(:,:,mod(3:5,3)+1) ,3);
    denom = 1 + sum(m.*m_x,3) + sum(m.*m_y,3) + sum(m_x.*m_y,3);
    rho = 4*atan2(triple,denom)/(4*pi*dx*dy);
    rho_avg = (rho(:,:)+rho(mod(-1:N-2,N)+1,:)+rho(:,mod(-1:N-2,N)+1)+rho(mod(-1:N-2,N)+1,mod(-1:N-2,N)+1))/4;

    % Coulomb term
    if coulomb == 1
        m_k1=zeros(N,N,3);
        m_x=(m(mod(1:N,N)+1,1:N,:)-m)/(dx);
        m_y=(m(1:N,mod(1:N,N)+1,:)-m)/(dy);
        parfor i = 1:N
            for j = 1:N
                m_k1 = m_k1 + (dist_x(:,:,:,i,j).*m_y - dist_y(:,:,:,i,j).*m_x )*rho(i,j)*dx*dy;
            end
        end

        m_k2=zeros(N,N,3);
        m_k2arg = m + dt/2*m_k1;
        m_x=(m_k2arg(mod(1:N,N)+1,1:N,:)-m_k2arg)/(dx);
        m_y=(m_k2arg(1:N,mod(1:N,N)+1,:)-m_k2arg)/(dy);
        parfor i = 1:N
            for j = 1:N
                m_k2 = m_k2 + (dist_x(:,:,:,i,j).*m_y - dist_y(:,:,:,i,j).*m_x )*rho(i,j)*dx*dy;
            end
        end

        m_k3=zeros(N,N,3);
        m_k3arg = m + dt/2*m_k2;
        m_x=(m_k3arg(mod(1:N,N)+1,1:N,:)-m_k3arg)/(dx);
        m_y=(m_k3arg(1:N,mod(1:N,N)+1,:)-m_k3arg)/(dy);
        parfor i = 1:N
            for j = 1:N
                m_k3 = m_k3 + (dist_x(:,:,:,i,j).*m_y - dist_y(:,:,:,i,j).*m_x )*rho(i,j)*dx*dy;
            end
        end

        m_k4=zeros(N,N,3);
        m_k4arg = m + dt*m_k3;
        m_x=(m_k4arg(mod(1:N,N)+1,1:N,:)-m_k4arg)/(dx);
        m_y=(m_k4arg(1:N,mod(1:N,N)+1,:)-m_k4arg)/(dy);
        parfor i = 1:N
            for j = 1:N
                m_k4 = m_k4 + (dist_x(:,:,:,i,j).*m_y - dist_y(:,:,:,i,j).*m_x )*rho(i,j)*dx*dy;
            end
        end

        m_RK4 = (m_k1+2*m_k2+2*m_k3+m_k4)/6;
        m = m + q*coulomb*dt*(m_RK4+m_RK4(mod(-1:N-2,N)+1,:,:)+m_RK4(:,mod(-1:N-2,N)+1,:)+m_RK4(mod(-1:N-2,N)+1,mod(-1:N-2,N)+1,:))/4;
        m = m./(sqrt(sum(m.^2,3))); % Renormalize
    end


    % check conserved quantities
    Q_top_init = Q_top;
    Q_top=sum(sum(rho))*dx*dy;  % topological charge
    S_z=sum(sum(m(:,:,3)));           % total spin in z-direction

    cent_of_mass_x = sum(sum(rho.*xx))/Q_top*dx*dy;
    cent_of_mass_y = sum(sum(rho.*yy))/Q_top*dx*dy;
    variance = sum(sum(rho.*((xx-cent_of_mass_x).^2+(yy-cent_of_mass_y).^2)))/Q_top*dx*dy;
    st_dev = sqrt(variance);

    m_dx=(m(mod(1:N,N)+1,:,:)-m(mod(-1:N-2,N)+1,:,:))/(2*dx); 
    m_dy=(m(:,mod(1:N,N)+1,:)-m(:,mod(-1:N-2,N)+1,:))/(2*dy);

    E_B = gmuB*S_z;    % B energy
    E_LL = sum(sum(-kappa/2 * (sum(m_dx.^2+m_dy.^2))));   % stiffness energy
    E_C = 0;   % Coulomb energy (need to do)
    E_El = El*sum(sum(xx.*rho));   % Applied electric field energy
    Q_top_list(length(Q_top_list)+1)=Q_top;
    E_B_list(length(E_B_list)+1)=E_B;
    E_LL_list(length(E_LL_list)+1)=E_LL;
    E_C_list(length(E_C_list)+1)=E_C;
    E_El_list(length(E_El_list)+1)=E_El;

    % check norm is preserved
    m_norm=sqrt(sum(m.^2,3));
    
    % plot
    quiver(xx,yy,m(:,:,1),m(:,:,2)) % full 2D vector field
    %quiver(xx(:,N/2),zeros(1,N),m(N/2,:,1),m(N/2,:,3))  % 1D slice
    hold on
    %quiver(cent_of_mass_x,cent_of_mass_y,st_dev,0,'r')  %radius vector
    contour(xx-dx/2,yy-dy/2,rho_avg,10) % color plot
    hold off
    axis([xlow xhigh ylow yhigh])
    drawnow
    
    % reset for new loop
    m_init=m;
    t=t+dt;
    m_full(:,:,:,t_ind) = m;
    rho_full(:,:,t_ind) = rho;
    t_ind=t_ind+1;
end

% to plot one of the energies:
% plot(1:length(Q_top_list),Q_top_list)

% draw saved data with no lag, can copypaste to console to do again
for i = 1:t_ind
    quiver(xx,yy,m_full(:,:,1,i),m_full(:,:,2,i))
    hold on
    contour(xx,yy,rho_full(:,:,i),10) % color plot
    hold off
    drawnow
end


