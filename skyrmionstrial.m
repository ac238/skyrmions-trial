

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
n=1;
z_0=0;
lambda=5;

omega = ((xx + yy*1i - z_0)/lambda).^n;
m_init(:,:,1)=4*real(omega)./((abs(omega)).^2+4);
m_init(:,:,2)=4*imag(omega)./((abs(omega)).^2+4);
m_init(:,:,3)=((abs(omega)).^2-4)./((abs(omega)).^2+4);


% optional: set boundary spins to +z
%for i = 1:N
%    m_init(1,i,:) = [0 0 0];
%    m_init(i,1,:) = [0 0 0];
%    m_init(N,i,:) = [0 0 0];
%    m_init(i,N,:) = [0 0 0];
%end


% initialize Coulomb distance matrix
%for i = 1:N
%    for j = 1:N
%        for k = 1:3
%            dist_x(:,:,k,i,j) = (xx-xx(i,j))./(((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^1.5);
%            dist_x(i,j,k,i,j) = 0;
%            dist_y(:,:,k,i,j) = (yy-yy(i,j))./(((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^1.5);
%            dist_y(i,j,k,i,j) = 0;
%        end
%        energy_dist(:,:,i,j) = ((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^-0.5;
%        energy_dist(i,j,i,j) = 0;
%    end
%end


% physical constants
B_field=1.1;
E_field=10.1;
q_electron=-3.1;
mass_electron=1.1;
rho_stiff=1.1;
dielectric=1.1;
g_factor=2.002;
casimir=0.5;
nu_level=1.0;

stiff_val = -rho_stiff*q_electron/(2*pi*casimir*nu_level*B_field);
b_val = g_factor*B_field*q_electron/(2*mass_electron);
alpha_val = q_electron^3/(8*pi^2*casimir*nu_level*B_field*dielectric);
e_val = -q_electron^2*E_field/(8*pi^2*casimir*nu_level*B_field);

%custom parameters
b_val = 0;
stiff_val = 1;
e_val = 6.2;
alpha_val = 0;


t=0;
t_final=100;
dt=0.01;
t_ind=1;
El_freq = 0.05 *2*pi/t_final;
Q_top = -n;
Q_top_list = []; % store values for final plot
E_B_list = [];
E_LL_list = [];
E_C_list = [];
Spin_list = [];

while t<t_final

    m = m_init;

    % Zeeman term with RK4
    m1_k1 = b_val*m(:,:,2);
    m2_k1 = -b_val*m(:,:,1);
    m1_k2 = b_val*(m(:,:,2)+dt/2*m2_k1);
    m2_k2 = -b_val*(m(:,:,1)+dt/2*m1_k1);
    m1_k3 = b_val*(m(:,:,2)+dt/2*m2_k2);
    m2_k3 = -b_val*(m(:,:,1)+dt/2*m1_k2);
    m1_k4 = b_val*(m(:,:,2)+dt*m2_k3);
    m2_k4 = -b_val*(m(:,:,1)+dt*m1_k3);

    m(:,:,1) = m(:,:,1) + dt/6*(m1_k1 + 2*m1_k2 + 2*m1_k3 + m1_k4);
    m(:,:,2) = m(:,:,2) + dt/6*(m2_k1 + 2*m2_k2 + 2*m2_k3 + m2_k4);

    % Electric field term with RK4
    m_k1 = (Elec_change_y(m))/(2*dy);
    m_k2 = (Elec_change_y(m) + dt/2*(Elec_change_y(m_k1)))/(2*dy);
    m_k3 = (Elec_change_y(m) + dt/2*(Elec_change_y(m_k2)))/(2*dy);
    m_k4 = (Elec_change_y(m) + dt*(Elec_change_y(m_k3)))/(2*dy);

    El = -e_val*cos(El_freq*t);
    m = m + El*dt/6*(m_k1 + 2*m_k2 + 2*m_k3 + m_k4);
    m = m./(sqrt(sum(m.^2,3))); % Renormalize

    % Stiffness term with RK4 (only works for dx=dy) (doesn't conserve norm=1)
    m_k1 = stiff_val/(dx^2) * (stiffness(m));
    m_k2arg = m+dt/2*m_k1;
    m_k2 = stiff_val/(dx^2) * (stiffness(m_k2arg));
    m_k3arg = m+dt/2*m_k2;
    m_k3 = stiff_val/(dx^2) * (stiffness(m_k3arg));
    m_k4arg = m+dt*m_k3;
    m_k4 = stiff_val/(dx^2) * (stiffness(m_k4arg));

    m = m + dt/6*(m_k1+2*m_k2+2*m_k3+m_k4);
    m = m./(sqrt(sum(m.^2,3))); % Renormalize

    % Pontryagin density
    m_x=m(mod(1:N,N)+1,1:N,:);
    m_y=m(1:N,mod(1:N,N)+1,:);
    triple = sum( m(:,:,mod(1:3,3)+1).*m_x(:,:,mod(2:4,3)+1).*m_y(:,:,mod(3:5,3)+1) - m(:,:,mod(1:3,3)+1).*m_y(:,:,mod(2:4,3)+1).*m_x(:,:,mod(3:5,3)+1) ,3);
    denom = 1 + sum(m.*m_x,3) + sum(m.*m_y,3) + sum(m_x.*m_y,3);
    rho = 4*atan2(triple,denom)/(4*pi*dx*dy);
    rho_avg = (rho(:,:)+rho(mod(-1:N-2,N)+1,:)+rho(:,mod(-1:N-2,N)+1)+rho(mod(-1:N-2,N)+1,mod(-1:N-2,N)+1))/4;

    % Coulomb term
    if alpha_val ~= 0
        m_k1=coulomb_loop(m,dist_x,dist_y,rho);
        m_k2arg = m + dt/2*m_k1;
        m_k2=coulomb_loop(m_k2arg,dist_x,dist_y,rho);
        m_k3arg = m + dt/2*m_k2;
        m_k3=coulomb_loop(m_k3arg,dist_x,dist_y,rho);
        m_k4arg = m + dt*m_k3;
        m_k4=coulomb_loop(m_k4arg,dist_x,dist_y,rho);

        m_RK4 = (m_k1+2*m_k2+2*m_k3+m_k4)/6;
        m = m - alpha_val*dt*(m_RK4+m_RK4(mod(-1:N-2,N)+1,:,:)+m_RK4(:,mod(-1:N-2,N)+1,:)+m_RK4(mod(-1:N-2,N)+1,mod(-1:N-2,N)+1,:))/4;
        m = m./(sqrt(sum(m.^2,3))); % Renormalize
    end


    % Coulomb energy
    coulomb_energy_field = zeros(N,N);
    parfor i = 1:N
        for j = 1:N
            coulomb_energy_field = coulomb_energy_field + rho.*rho(i,j).*energy_dist(:,:,i,j);
        end
    end
    E_C = alpha_val*2*pi*sum(sum(coulomb_energy_field));



    % check conserved quantities
    Q_top_init = Q_top;
    Q_top=sum(sum(rho))*dx*dy;  % topological charge
    S_z=sum(sum(m(:,:,3)));           % total spin in z-direction
    S_x=sum(sum(m(:,:,1)));      
    S_y=sum(sum(m(:,:,2)));      


    % track skyrmion c.o.m.
    %cent_of_mass_x = sum(sum(rho.*xx))/Q_top*dx*dy;
    %cent_of_mass_y = sum(sum(rho.*yy))/Q_top*dx*dy;
    %variance = sum(sum(rho.*((xx-cent_of_mass_x).^2+(yy-cent_of_mass_y).^2)))/Q_top*dx*dy;
    %st_dev = sqrt(variance);

    m_dx=(m(2:N,1:N-1,:)-m(1:N-1,1:N-1,:))/(dx); 
    m_dy=(m(1:N-1,2:N,:)-m(1:N-1,1:N-1,:))/(dy);

    E_B = b_val*S_z;    % B energy
    E_LL = sum(sum(stiff_val/2 * (sum(m_dx.^2+m_dy.^2))));   % stiffness energy
    Q_top_list(length(Q_top_list)+1)=Q_top;
    E_B_list(length(E_B_list)+1)=E_B;
    E_LL_list(length(E_LL_list)+1)=E_LL;
    E_C_list(length(E_C_list)+1)=E_C;
    Spin_list(length(Spin_list)+1,:)=[S_x S_y S_z];

    % check norm is preserved
    m_norm=sqrt(sum(m.^2,3));
    
    % plot
    quiver(xx,yy,m(:,:,1),m(:,:,2)) % full 2D vector field
    %quiver(xx(:,N/2),zeros(1,N),m(N/2,:,1),m(N/2,:,3))  % 1D slice
    hold on
    %quiver(cent_of_mass_x,cent_of_mass_y,st_dev,0,'r')  %radius vector
    contour(xx(1:N-1,1:N-1)-dx/2,yy(1:N-1,1:N-1)-dy/2,rho(1:N-1,1:N-1),10) % color plot
    hold off
    axis([xlow xhigh ylow yhigh])
    drawnow
    
    % reset for new loop
    m_init=m;
    t=t+dt
    m_full(:,:,:,t_ind) = m;
    rho_full(:,:,t_ind) = rho;
    t_ind=t_ind+1;
end

E_total_list = E_B_list+E_LL_list+E_C_list;

% to plot one of the energies:
% plot((1:length(Q_top_list))*dt,Q_top_list)
% plot((1:length(E_total_list))*dt,E_total_list)
% plot((1:length(E_B_list))*dt,E_B_list,(1:length(E_LL_list))*dt,E_LL_list,(1:length(E_C_list))*dt,E_C_list,(1:length(E_total_list))*dt,E_total_list)
% plot((1:length(Spin_list(:,1)))*dt,Spin_list(:,1),(1:length(Spin_list(:,2)))*dt,Spin_list(:,2),(1:length(Spin_list(:,3)))*dt,Spin_list(:,3))

% draw saved data with no lag, can copypaste to console to do again
for i = 1:t_ind
    i = i*10;
    E_total_list(i)
    quiver(xx,yy,m_full(:,:,1,i),m_full(:,:,2,i))
    hold on
    contour(xx(1:N-1,1:N-1)-dx/2,yy(1:N-1,1:N-1)-dy/2,rho_full(1:N-1,1:N-1,i),10) % color plot
    hold off
    drawnow
end




