
% open file to log outputs
fileID = fopen('skyrmionstrial.out','w');

% begin parallelization
numCores = feature('numcores')
%p = parpool(numCores);

% create bounds of graph
N=100;
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

"initialized"



%custom parameters, all positive!
b_val = 0.13;
stiff_val = 1.0;
e_val = 0;
alpha_val = 70;

% initialize Coulomb distance matrix
if alpha_val ~= 0
    dist_x=zeros(N,N,N,N);
    dist_y=zeros(N,N,N,N);
    for i = 1:N
        for j = 1:N
            dist_x(:,:,i,j) = (xx-xx(i,j))./(((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^1.5);
            dist_x(i,j,i,j) = 0;
            dist_y(:,:,i,j) = (yy-yy(i,j))./(((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^1.5);
            dist_y(i,j,i,j) = 0;
            energy_dist(:,:,i,j) = ((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^-0.5;
            energy_dist(i,j,i,j) = 0;
        end
    end
end


t=0;
t_final=10;
dt=0.01;
t_ind=1;
El_freq = 5 *2*pi/t_final; %num of cycles * 2pi*t_final
Q_top = -n;
Q_top_list = []; % store values for final plot
E_B_list = [];
E_LL_list = [];
E_C_list = [];
Spin_list = [];



while t<t_final
    tic
    m = m_init;
    rho = pontryagin(m);
    rho_avg = (rho(:,:)+rho(mod(-1:N-2,N)+1,:)+rho(:,mod(-1:N-2,N)+1)+rho(mod(-1:N-2,N)+1,mod(-1:N-2,N)+1))/4;

    % RK4
    B_field_thisstep=B_eff(m,b_val,e_val,alpha_val,dist_x,dist_y,rho);
    m_k1 = m(:,:,mod(1:3,3)+1).*B_field_thisstep(:,:,mod(2:4,3)+1)-m(:,:,mod(2:4,3)+1).*B_field_thisstep(:,:,mod(1:3,3)+1);

    m_k2arg=m+dt/2*m_k1;
    B_field=B_eff(m_k2arg,b_val,e_val,alpha_val,dist_x,dist_y,rho);
    m_k2 = m_k2arg(:,:,mod(1:3,3)+1).*B_field(:,:,mod(2:4,3)+1)-m_k2arg(:,:,mod(2:4,3)+1).*B_field(:,:,mod(1:3,3)+1);

    m_k3arg=m+dt/2*m_k2;
    B_field=B_eff(m_k3arg,b_val,e_val,alpha_val,dist_x,dist_y,rho);
    m_k3 = m_k3arg(:,:,mod(1:3,3)+1).*B_field(:,:,mod(2:4,3)+1)-m_k3arg(:,:,mod(2:4,3)+1).*B_field(:,:,mod(1:3,3)+1);

    m_k4arg=m+dt*m_k3;
    B_field=B_eff(m_k4arg,b_val,e_val,alpha_val,dist_x,dist_y,rho);
    m_k4 = m_k4arg(:,:,mod(1:3,3)+1).*B_field(:,:,mod(2:4,3)+1)-m_k4arg(:,:,mod(2:4,3)+1).*B_field(:,:,mod(1:3,3)+1);

    m = m + dt/6*(m_k1+2*m_k2+2*m_k3+m_k4);
    m = m./(sqrt(sum(m.^2,3))); % Renormalize

    % Electric field term with RK4
    m_k1 = (Elec_change_y(m))/(2*dy);
    m_k2 = (Elec_change_y(m) + dt/2*(Elec_change_y(m_k1)))/(2*dy);
    m_k3 = (Elec_change_y(m) + dt/2*(Elec_change_y(m_k2)))/(2*dy);
    m_k4 = (Elec_change_y(m) + dt*(Elec_change_y(m_k3)))/(2*dy);

    El = e_val*cos(El_freq*t);
    m = m - El*dt/6*(m_k1 + 2*m_k2 + 2*m_k3 + m_k4);
    m = m./(sqrt(sum(m.^2,3))); % Renormalize



    % check conserved quantities
    Q_top_init = Q_top;
    Q_top=sum(sum(rho(1:N-1,1:N-1)))*dx*dy;  % topological charge
    S_x=sum(sum(m(:,:,1)));      
    S_y=sum(sum(m(:,:,2)));      
    S_z=sum(sum(m(:,:,3)));           % total spin in z-direction

    m_dx=(m(2:N,1:N-1,:)-m(1:N-1,1:N-1,:))/(dx); 
    m_dy=(m(1:N-1,2:N,:)-m(1:N-1,1:N-1,:))/(dy);

    E_B = -b_val*S_z;    % B energy
    E_LL = sum(sum(stiff_val/2 * (sum(m_dx.^2+m_dy.^2))));   % stiffness energy
    Q_top_list(length(Q_top_list)+1)=Q_top;
    E_B_list(length(E_B_list)+1)=E_B;
    E_LL_list(length(E_LL_list)+1)=E_LL;
    E_C_list(length(E_C_list)+1)=E_C;
    Spin_list(length(Spin_list)+1,:)=[S_x S_y S_z];
    
    % plot
    %quiver(xx,yy,m(:,:,1),m(:,:,2)) % full 2D vector field
    %hold on
    %contour(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),rho(1:N-1,1:N-1),10) % color plot
    %hold off
    %axis([xlow xhigh ylow yhigh])
    %title(t)
    %drawnow
    
    % reset for new loop
    m_init=m;
    t=t+dt;
    m_full(:,:,:,t_ind) = m;
    rho_full(:,:,t_ind) = rho;
    t_ind=t_ind+1;
    toc
end

"completed time evolution"

E_total_list = E_B_list+E_LL_list+E_C_list;

% plotting all the energies
plot((1:length(Q_top_list))*dt,Q_top_list)
drawnow
saveas(gcf,"fig_Q_top.m")
plot((1:length(E_B_list))*dt,E_B_list,(1:length(E_LL_list))*dt,E_LL_list,(1:length(E_C_list))*dt,E_C_list,(1:length(E_total_list))*dt,E_total_list)
legend("Zeeman","Stiffness","Coulomb","Total")
drawnow
saveas(gcf,"fig_Energies.m")
plot((1:length(Spin_list(:,1)))*dt,Spin_list(:,1),(1:length(Spin_list(:,2)))*dt,Spin_list(:,2),(1:length(Spin_list(:,3)))*dt,Spin_list(:,3))
legend("S_x","S_y","S_z")
drawnow
saveas(gcf,"fig_Spin_components")


"created plots"


% draw saved data with no lag, can copypaste to console to do again
for i = 1:(t_ind/10)
    i = i*10;
    quiver(xx,yy,m_full(:,:,1,i),m_full(:,:,2,i))
    drawnow
    saveas(gcf,"zz_quiver_frame"+string(i)+".png")
    pcolor(xx(1:N-1,1:N-1)-dx/2,yy(1:N-1,1:N-1)-dy/2,rho_full(1:N-1,1:N-1,i)) % color plot
    drawnow
    saveas(gcf,"zz_contour_frame"+string(i)+".png")
end

"drew frames of evolution"

fclose(fileID);
%delete(p);


