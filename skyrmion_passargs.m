
function thebigfunction = skyrmion_passargs(jobidp, cpusp, framespacep, b_valp, stiff_valp, e_valp, alpha_valp, damp_valp, damp_falloffp, Np, magnon_ampp, t_finalp)

% args to pass
jobid = jobidp
framespace = framespacep
feature('numcores')
numCores = cpusp
b_val = b_valp
stiff_val = stiff_valp
e_val = e_valp
alpha_val = alpha_valp
damp_val = damp_valp
damp_falloff = damp_falloffp
N=Np
magnon_amp = magnon_ampp
t_final=t_finalp

% make directory to place files
%jobid = string(datetime)
[status,msg] = mkdir(jobid)
cd(jobid)

% open file to log outputs
fileID = fopen('skyrmions.out','w');
diary skyrmions.out

% begin parallelization
%p = parpool(numCores);

% create bounds of graph
dx=1;
dy=1;
xlow=-(N+1)/2;
ylow=-(N+1)/2;
xhigh=(N-1)/2;
yhigh=(N-1)/2;
clow = 0;
chigh = 0.03;
[yy,xx]=meshgrid(linspace(xlow,xhigh,N),linspace(ylow,yhigh,N)); %xx and yy are swapped to correspond w/ indices
axis([xlow xhigh ylow yhigh])

% initialize skyrmion according to QHMF 35
n=1;
z_0=0;
lambda=3;
omega = ((xx + yy*1i - z_0)/lambda).^n;
%omega = (((xx + yy*1i - z_0*ones(N,N))/lambda).^n).*(((xx + yy*1i + z_0*ones(N,N))/lambda).^1);
m_init(:,:,1)=4*real(omega)./((abs(omega)).^2+4);
m_init(:,:,2)=4*imag(omega)./((abs(omega)).^2+4);
m_init(:,:,3)=((abs(omega)).^2-4)./((abs(omega)).^2+4);
% initialize all vertical
%m_init(:,:,1)=zeros(N,N);
%m_init(:,:,2)=zeros(N,N);
%m_init(:,:,3)=ones(N,N);
% initialize all random
%m_init(:,:,1)=2*rand(N)-ones(N,N);
%m_init(:,:,2)=2*rand(N)-ones(N,N);
%m_init(:,:,3)=2*rand(N)-ones(N,N);
% edges vertical
%for i = 1:N
%    m_init(1,i,1)=0;
%    m_init(i,1,1)=0;
%    m_init(N,i,1)=0;
%    m_init(i,N,1)=0;
%    m_init(1,i,2)=0;
%    m_init(i,1,2)=0;
%    m_init(N,i,2)=0;
%    m_init(i,N,2)=0;
%    m_init(1,i,3)=1;
%    m_init(i,1,3)=1;
%    m_init(N,i,3)=1;
%    m_init(i,N,3)=1;
%end
% Initialize annealed
%load("m_final_annealed200.mat");
%load("m_final_annealedDI.mat");
%m_init = m;
%m_init = m(26:75,26:75,:);
clear m;

%m_init = m_init./(sqrt(sum(m_init.^2,3))); % Renormalize

"initialized"



%custom parameters, all positive!
damp_mat=zeros(N,N,3);
for i=1:N
    for j=1:N
        for k=1:3
            %damp_mat(i,N,k)=damp_val;   %one line
            damp_mat(i,j,k)=damp_val*exp((j-N)/damp_falloff);   %exponential
            %damp_mat(i,j,k)=damp_val;   %const everywhere
        end
    end
end
%damp_mat(1,1,1)=damp_val;   %corner only!


% initialize Coulomb distance matrix
dist_x=zeros(N,N,N,N);
dist_y=zeros(N,N,N,N);
energy_dist=zeros(N,N,N,N);
if alpha_val ~= 0
    for i = 1:N
        for j = 1:N
            dist_x(:,:,i,j) = (xx-xx(i,j))./(((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^1.5);
            dist_x(i,j,i,j) = 0;
            dist_y(:,:,i,j) = (yy-yy(i,j))./(((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^1.5);
            dist_y(i,j,i,j) = 0;
            energy_dist(i,j,:,:) = ((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^-0.5;
            energy_dist(i,j,i,j) = 0;
        end
    end
end
center_dist = energy_dist(2:N-1,2:N-1,:,:);
sliced_dist = energy_dist(1:N-1,1:N-1,:,:);


t=0;
dt=0.01;
t_ind=1;
El_freq = 5 *2*pi/t_final; %num of cycles * 2pi*t_final
mag_freq = b_val+2
kappax=2*pi*0.1;
kappay=2*pi*0.11;

Q_top_list = []; % store values for final plot
E_B_list = [];
E_LL_list = [];
E_C_list = [];
E_eff_list = [];
E_loss1_list = [];
E_loss1=0;
E_loss2_list = [];
E_loss2=0;
Spin_list = [];
mean_x_list = [];
mean_y_list = [];
stdev_list = [];



while t<t_final
    tic
    m = m_init;
    rho = pontryagin(m);
    %rho = DI_Pontryagin(0,m)/(8*pi);
    %rho_avg = (rho(:,:)+rho(mod(-1:N-2,N)+1,:)+rho(:,mod(-1:N-2,N)+1)+rho(mod(-1:N-2,N)+1,mod(-1:N-2,N)+1))/4;  %oops still periodic

    % set electric field
    B_Zeeman = zeros(N,N,3);
    B_Zeeman(:,:,3)=b_val*ones(N,N);
    El_x = e_val*cos(El_freq*t);
    El_y = e_val*cos(0.8*El_freq*t);

    % Gilbert damping
    %for posx = 1:N
    %    for posy = N-4:N
    %        m(posx,posy,:) = m(posx,posy,:) - dt*damp_val*(sum(m(posx,posy,:).*B_field_thisstep(posx,posy,:))*m(posx,posy,:)-B_field_thisstep(posx,posy,:));
    %    end
    %end
    %m(1:N,d_rows,:) = m(1:N,d_rows,:) - dt*damp_val*((m(1:N,d_rows,1).*B_field_thisstep(1:N,d_rows,1) + m(1:N,d_rows,2).*B_field_thisstep(1:N,d_rows,2) + m(1:N,d_rows,3).*B_field_thisstep(1:N,d_rows,3)).*m(1:N,d_rows,:)-B_field_thisstep(1:N,d_rows,:));

    % RK4
    B_field_thisstep=B_eff_solid(m,b_val,stiff_val,El_x,El_y,alpha_val,energy_dist,center_dist,sliced_dist);
    %B_field_thisstep=DI_Beffective2( N, alpha_val, stiff_val, B_Zeeman, zeros(N,N), m);
    B_mod_thisstep=B_field_thisstep;
    m_k1 = m(:,:,mod(1:3,3)+1).*B_field_thisstep(:,:,mod(2:4,3)+1)-m(:,:,mod(2:4,3)+1).*B_field_thisstep(:,:,mod(1:3,3)+1) - damp_mat.*((m(:,:,1).*B_mod_thisstep(:,:,1) + m(:,:,2).*B_mod_thisstep(:,:,2) + m(:,:,3).*B_mod_thisstep(:,:,3)).*m-B_mod_thisstep);

    m_k2arg=m+dt/2*m_k1;
    B_field=B_eff_solid(m_k2arg,b_val,stiff_val,El_x,El_y,alpha_val,energy_dist,center_dist,sliced_dist);
    %B_field=DI_Beffective2( N, alpha_val, stiff_val, B_Zeeman, zeros(N,N), m_k2arg);
    B_mod=B_field;
    m_k2 = m_k2arg(:,:,mod(1:3,3)+1).*B_field(:,:,mod(2:4,3)+1)-m_k2arg(:,:,mod(2:4,3)+1).*B_field(:,:,mod(1:3,3)+1) - damp_mat.*((m_k2arg(:,:,1).*B_mod(:,:,1) + m_k2arg(:,:,2).*B_mod(:,:,2) + m_k2arg(:,:,3).*B_mod(:,:,3)).*m_k2arg-B_mod);

    m_k3arg=m+dt/2*m_k2;
    B_field=B_eff_solid(m_k3arg,b_val,stiff_val,El_x,El_y,alpha_val,energy_dist,center_dist,sliced_dist);
    %B_field=DI_Beffective2( N, alpha_val, stiff_val, B_Zeeman, zeros(N,N), m_k3arg);
    B_mod=B_field;
    m_k3 = m_k3arg(:,:,mod(1:3,3)+1).*B_field(:,:,mod(2:4,3)+1)-m_k3arg(:,:,mod(2:4,3)+1).*B_field(:,:,mod(1:3,3)+1) - damp_mat.*((m_k3arg(:,:,1).*B_mod(:,:,1) + m_k3arg(:,:,2).*B_mod(:,:,2) + m_k3arg(:,:,3).*B_mod(:,:,3)).*m_k3arg-B_mod);

    m_k4arg=m+dt*m_k3;
    B_field=B_eff_solid(m_k4arg,b_val,stiff_val,El_x,El_y,alpha_val,energy_dist,center_dist,sliced_dist);
    %B_field=DI_Beffective2( N, alpha_val, stiff_val, B_Zeeman, zeros(N,N), m_k4arg);
    B_mod=B_field;
    m_k4 = m_k4arg(:,:,mod(1:3,3)+1).*B_field(:,:,mod(2:4,3)+1)-m_k4arg(:,:,mod(2:4,3)+1).*B_field(:,:,mod(1:3,3)+1) - damp_mat.*((m_k4arg(:,:,1).*B_mod(:,:,1) + m_k4arg(:,:,2).*B_mod(:,:,2) + m_k4arg(:,:,3).*B_mod(:,:,3)).*m_k4arg-B_mod);

    dmdt = (m_k1+2*m_k2+2*m_k3+m_k4)/6;
    m = m + dmdt*dt;

    % set edges to 0
    %m(:,1,1:2) = zeros(N,1,2);
    %m(1,:,1:2) = zeros(1,N,2);
    %m(:,N,1:2) = zeros(N,1,2);
    %m(N,:,1:2) = zeros(1,N,2);

    % drive magnons
    if magnon_amp ~= 0
        for i=1:N
            m(i,1,1) = magnon_amp*sin(mag_freq*t);
            m(i,1,2) = magnon_amp*cos(mag_freq*t);
            %diagonal:
            %m(i,1,1) = magnon_amp*sin(mag_freq*t+i*0.5*pi);
            %m(i,1,2) = magnon_amp*cos(mag_freq*t+i*0.5*pi);
            %charged:
            %m(i,1,1) = cos(kappay.*yy(i,1)-mag_freq*t)*sin(kappax.*xx(i,1));
            %m(i,1,2) = cos(kappay.*yy(i,1)-mag_freq*t)*cos(kappax.*xx(i,1));
            %m(i,1,3) = sin(kappay.*yy(i,1)-mag_freq*t);
        end
    end

    % Renormalize
    m = m./(sqrt(sum(m.^2,3))); 


    % check conserved quantities
    Q_top=sum(sum(rho(1:N-1,1:N-1)))*dx*dy;  % topological charge
    S_x=sum(sum(m(:,:,1)));      
    S_y=sum(sum(m(:,:,2)));      
    S_z=sum(sum(m(:,:,3)));           % total spin in z-direction

    m_dx=(m(2:N,1:N-1,:)-m(1:N-1,1:N-1,:))/(dx); 
    m_dy=(m(1:N-1,2:N,:)-m(1:N-1,1:N-1,:))/(dy);

    % Coulomb energy
    if alpha_val ~= 0
        coulomb_energy_field = zeros(N,N);
        for i = 1:N
            for j = 1:N
                coulomb_energy_field = coulomb_energy_field + rho.*rho(i,j).*energy_dist(:,:,i,j);
            end
        end
        %coulomb_energy_field = distributed(tensorprod(rho,energy_dist,[1 2])).*rho;
        E_C = alpha_val*2*pi*sum(sum(coulomb_energy_field));
    else
        E_C = 0;
    end

    E_B = -b_val*(S_z-N*N);    % B energy
    E_LL = sum(sum(stiff_val/2 * (sum(m_dx.^2+m_dy.^2))));   % stiffness energy
    E_eff = -sum(sum(sum(B_field_thisstep.*m_init))) + b_val*N*N;
    %general energy loss %warning: backward derivative
    E_loss1 = E_loss1 + sum(sum(sum(B_field_thisstep.*dmdt)))*dt; 
    %This one works!!
    E_loss2 = E_loss2 + sum(-damp_mat(:,:,1).*(dot(B_field_thisstep,m_init,3).^2-dot(B_field_thisstep,B_field_thisstep,3)),[1 2])*dt;
    %Torque dot omega
    %E_loss1 = E_loss1 + sum(-damp_mat(:,:,1).*dot(m_init,B_field_thisstep,3).*dot(m_init,cross(B_field_thisstep,dmdt),3),[1 2])*dt;

    Q_top_list(length(Q_top_list)+1)=Q_top;
    E_B_list(length(E_B_list)+1)=E_B;
    E_LL_list(length(E_LL_list)+1)=E_LL;
    E_C_list(length(E_C_list)+1)=E_C;
    E_eff_list(length(E_eff_list)+1)=E_eff;
    E_loss1_list(length(E_loss1_list)+1)=E_loss1;
    E_loss2_list(length(E_loss2_list)+1)=E_loss2;
    Spin_list(length(Spin_list)+1,:)=[S_x S_y S_z];

    
    % mean and stdev
    mean_x = sum(sum(rho.*xx))/Q_top*dx*dy;
    mean_y = sum(sum(rho.*yy))/Q_top*dx*dy;
    variance = sum(sum(rho.*((xx-mean_x).^2+(yy-mean_y).^2)))/Q_top*dx*dy;
    st_dev = sqrt(variance);
    mean_x_list(length(mean_x_list)+1)=mean_x;
    mean_y_list(length(mean_y_list)+1)=mean_y;
    stdev_list(length(stdev_list)+1)=st_dev;

    
    %plot
    if t_ind==1
        contour(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),rho(1:N-1,1:N-1),10)
        climsave=clim;
        %climsave=[-0.000000000001 0.000000000001]; %for charged magnons
        contour(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),m(1:N-1,1:N-1,3),10)
        climsave_z=[-1 1];
    end

    %rho_deriv = pontryagin_deriv(m,1,1) - pontryagin_deriv(m,-1,1) - pontryagin_deriv(m,1,-1) + pontryagin_deriv(m,-1,-1);
    
    if mod(t_ind,framespace)==0

        pc=pcolor(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),rho(1:N-1,1:N-1)); % color plot
        %pc=pcolor(xx,yy,rho);
        pc.EdgeColor='none';
        clim(climsave);
        colorbar
        title("Topological charge at t="+num2str(t_ind*dt));
        drawnow
        %saveas(gcf,"zz_charge_frame"+string(t_ind)+".png")
        %movefile("zz_charge_frame"+string(t_ind)+".png",jobid);
        exportgraphics(gcf,"charge.gif",'Append',true)


        pc=pcolor(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),m(1:N-1,1:N-1,3)); % color plot
        pc.EdgeColor='none';
        clim(climsave_z);
        thebar=colorbar;
        thebar.Label.String = 'z component';
        title("Magnetization field components at t="+num2str(t_ind*dt));
        hold on
        quiver(xx,yy,m(:,:,1),m(:,:,2),'r')
        quiver(mean_x,mean_y,st_dev/sqrt(2),st_dev/sqrt(2),'g')
        hold off
        drawnow
        %saveas(gcf,"zz_quiver_frame"+string(t_ind)+".png")
        %movefile("zz_quiver_frame"+string(t_ind)+".png",jobid);
        exportgraphics(gcf,"quiver.gif",'Append',true)

        m_full(:,:,:,t_ind/framespace) = m;
        rho_full(:,:,t_ind/framespace) = rho;
    end

    % reset for new loop
    m_init=m;
    t=t+dt
    t_ind=t_ind+1;
    toc

    %if t>50
    %    magnon_amp=0;
    %end
end

"completed time evolution"

E_total_list = E_B_list+E_LL_list+E_C_list;

% plotting all the energies
plot((1:length(Q_top_list))*dt,Q_top_list)
drawnow
saveas(gcf,"fig_Q_top.m")
plot((1:length(E_B_list))*dt,E_B_list,(1:length(E_LL_list))*dt,E_LL_list,(1:length(E_C_list))*dt,E_C_list,(1:length(E_total_list))*dt,E_total_list,(1:length(E_loss1_list))*dt,E_loss1_list,(1:length(E_loss2_list))*dt,E_loss2_list)
legend("Zeeman","Stiffness","Coulomb","Total","Total loss?","Damping?")
drawnow
saveas(gcf,"fig_Energies.m")
plot((1:length(Spin_list(:,1)))*dt,Spin_list(:,1),(1:length(Spin_list(:,2)))*dt,Spin_list(:,2),(1:length(Spin_list(:,3)))*dt,Spin_list(:,3))
legend("S_x","S_y","S_z")
drawnow
saveas(gcf,"fig_Spin_components")
plot((1:length(mean_x_list))*dt,mean_x_list,(1:length(mean_y_list))*dt,mean_y_list,(1:length(stdev_list))*dt,stdev_list)
legend("mean_x","mean_y","stdev")
drawnow
saveas(gcf,"fig_Pos_radius")


pc=pcolor(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),rho(1:N-1,1:N-1)); % color plot
%pc=pcolor(xx,yy,rho);
pc.EdgeColor='none';
%clim(climsave);
colorbar
drawnow
saveas(gcf,"fig_charge_frame"+string(t_ind))
saveas(gcf,"fig_charge_frame"+string(t_ind)+".png")

pc=pcolor(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),m(1:N-1,1:N-1,3)); % color plot
pc.EdgeColor='none';
%clim(climsave_z);
colorbar
hold on
quiver(xx,yy,m(:,:,1),m(:,:,2),'r')
quiver(mean_x,mean_y,st_dev/sqrt(2),st_dev/sqrt(2),'g')
hold off
drawnow
%saveas(gcf,"fig_quiver_frame"+string(t_ind))
%saveas(gcf,"fig_quiver_frame"+string(t_ind)+".png")
exportgraphics(gcf,"quiver.gif",'Append',true)


"created plots"

% copypaste this to view animation
%for i = 1:(t_ind/10)
%    i = i*10;
%    quiver(xx,yy,m_full(:,:,1,i),m_full(:,:,2,i))
%    hold on
%    contour(xx(1:N-1,1:N-1)-dx/2,yy(1:N-1,1:N-1)-dy/2,rho_full(1:N-1,1:N-1,i),10) % color plot
%    hold off
%    drawnow
%end

"drew frames of evolution"

fclose(fileID);
%delete(p);
cd ..

end



function B_field = B_eff_solid(m_arg,b_val,stiff_val,El_x,El_y,alpha_val,energy_dist,center_dist,sliced_dist)
    N=length(m_arg(:,1,1));
    N_1=N-1;
    dx=1;
    dy=1;
    rho = pontryagin(m_arg);

    %Zeeman
    B_Zeeman = zeros(N,N,3);
    B_Zeeman(:,:,3)=b_val*ones(N,N);

    %Stiffness
    B_stiffness=zeros(N,N,3);
    %Centered stiffness:
    B_stiffness(2:N-1,2:N-1,:) = (m_arg(1:N-2,2:N-1,:)+m_arg(3:N,2:N-1,:)+m_arg(2:N-1,1:N-2,:)+m_arg(2:N-1,3:N,:)-4*m_arg(2:N-1,2:N-1,:))/(dx^2)*stiff_val;
    %Wrapped stiffness:
    %B_stiffness(:,2:N-1,:) = (m_arg(mod(-1:N-2,N)+1,2:N-1,:)+m_arg(mod(1:N,N)+1,2:N-1,:)+m_arg(:,1:N-2,:)+m_arg(:,3:N,:)-4*m_arg(:,2:N-1,:))/(dx^2)*stiff_val;

    %Coulomb
    centered_int = zeros(N,N);
    sliced_int = zeros(N,N);
    full_int=zeros(N,N,3);
    full_intx=zeros(N,N,3);
    full_inty=zeros(N,N,3);
    full_intxy=zeros(N,N,3);

    %could cause slowdown...
    %center_dist = energy_dist(2:N-1,2:N-1,:,:);
    %sliced_dist = energy_dist(1:N-1,1:N-1,:,:);
    center_rho = rho(2:N_1,2:N_1);
    sliced_rho = rho(1:N_1,1:N_1);

    if alpha_val ~=0
        %centered_intx = zeros(N,N);
        %centered_inty = zeros(N,N);
        %centered_intxy = zeros(N,N);
        for i = 1:N
            for j = 1:N
                centered_int(i,j) = sum(sum(center_dist(:,:,i,j).*center_rho));
                %sliced_int(i,j) = sum(sum(center_dist(i,j).*center_rho));
            end
        end
        %centered_intx = distributed(tensorprod(dist_x(2:N-1,2:N-1,2:N-1,2:N-1),centered_rho(1:N-2,1:N-2),[1 2]));
        %centered_inty = distributed(tensorprod(dist_y(2:N-1,2:N-1,2:N-1,2:N-1),centered_rho(1:N-2,1:N-2),[1 2]));
    end
   
    for k=1:3
        full_int(:,:,k)=centered_int;
        full_intx(2:N,:,k)=centered_int(1:N-1,:);
        full_inty(:,2:N,k)=centered_int(:,1:N-1);
        full_intxy(2:N,2:N,k)=centered_int(1:N-1,1:N-1);
        %sliced version
        %full_int(1:N-1,1:N-1,k)=centered_int(1:N-1,1:N-1);
        %full_intx(2:N,1:N-1,k)=centered_int(1:N-1,1:N-1);
        %full_inty(1:N-1,2:N,k)=centered_int(1:N-1,1:N-1);
        %full_intxy(2:N,2:N,k)=centered_int(1:N-1,1:N-1);
    end

    %Deepak's averaging line
    %centered_int = 0.25*(centered_int + circshift(centered_int,[-1 0]) + circshift(centered_int,[0 -1])+ circshift(centered_int,[-1 -1]));

    %B_coulomb=-4*pi*alpha_val*(centered_int.*pontryagin_deriv(m_arg,1,1) - circshift(centered_int,[1 0]).*pontryagin_deriv(m_arg,-1,1) - circshift(centered_int,[0 1]).*pontryagin_deriv(m_arg,1,-1) + circshift(centered_int,[1 1]).*pontryagin_deriv(m_arg,-1,-1));
    B_coulomb=-4*pi*alpha_val*(full_int.*pontryagin_deriv2(m_arg,1,1) - full_intx.*pontryagin_deriv2(m_arg,-1,1) - full_inty.*pontryagin_deriv2(m_arg,1,-1) + full_intxy.*pontryagin_deriv2(m_arg,-1,-1));

    B_electric = (El_x+El_y)*pontryagin_deriv2(m_arg,1,1) + (-El_x+El_y)*pontryagin_deriv2(m_arg,-1,1) + (El_x-El_y)*pontryagin_deriv2(m_arg,1,-1) + (-El_x-El_y)*pontryagin_deriv2(m_arg,-1,-1);

    B_field = B_Zeeman+B_stiffness+B_coulomb+B_electric;

end


function rho = pontryagin(m)
    N=length(m);
    dx=1;
    dy=1;
    m_x=m(mod(1:N,N)+1,1:N,:);
    m_y=m(1:N,mod(1:N,N)+1,:);
    m_xy=m(mod(1:N,N)+1,mod(1:N,N)+1,:);
    %triple1 = sum( m(:,:,mod(1:3,3)+1).*m_x(:,:,mod(2:4,3)+1).*m_xy(:,:,mod(3:5,3)+1) - m(:,:,mod(1:3,3)+1).*m_xy(:,:,mod(2:4,3)+1).*m_x(:,:,mod(3:5,3)+1) ,3);
    triple1 = dot(m,cross(m_x,m_xy,3),3);
    denom1 = 1 + sum(m.*m_x,3) + sum(m.*m_xy,3) + sum(m_x.*m_xy,3);
    %triple2 = sum( m(:,:,mod(1:3,3)+1).*m_xy(:,:,mod(2:4,3)+1).*m_y(:,:,mod(3:5,3)+1) - m(:,:,mod(1:3,3)+1).*m_y(:,:,mod(2:4,3)+1).*m_xy(:,:,mod(3:5,3)+1) ,3);
    triple2 = dot(m,cross(m_xy,m_y,3),3);
    denom2 = 1 + sum(m.*m_xy,3) + sum(m.*m_y,3) + sum(m_xy.*m_y,3);
    rho = 4*(atan2(triple1,denom1)+atan2(triple2,denom2))/(8*pi*dx*dy);
end




function rho_prime = pontryagin_deriv2(m,x_dir,y_dir)
    N=length(m);
    dx=1;
    dy=1;
    m_x=m(2+x_dir:N-1+x_dir,2:N-1,:);
    m_y=m(2:N-1,2+y_dir:N-1+y_dir,:);
    m_xy=m(2+x_dir:N-1+x_dir,2+y_dir:N-1+y_dir,:);
    m_c=m(2:N-1,2:N-1,:);
    triple = sum( m_c(:,:,mod(1:3,3)+1).*m_x(:,:,mod(2:4,3)+1).*m_xy(:,:,mod(3:5,3)+1) - m_c(:,:,mod(1:3,3)+1).*m_xy(:,:,mod(2:4,3)+1).*m_x(:,:,mod(3:5,3)+1) ,3);
    denom = 1 + sum(m_c.*m_x,3) + sum(m_c.*m_xy,3) + sum(m_x.*m_xy,3);
    triple2 = sum( m_c(:,:,mod(1:3,3)+1).*m_xy(:,:,mod(2:4,3)+1).*m_y(:,:,mod(3:5,3)+1) - m_c(:,:,mod(1:3,3)+1).*m_y(:,:,mod(2:4,3)+1).*m_xy(:,:,mod(3:5,3)+1) ,3);
    denom2 = 1 + sum(m_c.*m_xy,3) + sum(m_c.*m_y,3) + sum(m_xy.*m_y,3);
    rho_prime1 = zeros(N,N,3);
    rho_prime1(2:N-1,2:N-1,:) = 1/(pi*dx*dy) * ((m_x(:,:,mod(1:3,3)+1).*m_xy(:,:,mod(2:4,3)+1) - m_xy(:,:,mod(1:3,3)+1).*m_x(:,:,mod(2:4,3)+1)).*denom - (m_x+m_xy).*triple ) ./ (triple.^2+denom.^2);
    rho_prime2 = zeros(N,N,3);
    rho_prime2(2:N-1,2:N-1,:) = 1/(pi*dx*dy) * ((m_xy(:,:,mod(1:3,3)+1).*m_y(:,:,mod(2:4,3)+1) - m_y(:,:,mod(1:3,3)+1).*m_xy(:,:,mod(2:4,3)+1)).*denom2 - (m_xy+m_y).*triple2 ) ./ (triple2.^2+denom2.^2);
    rho_prime=(rho_prime1+rho_prime2)*0.5;
end