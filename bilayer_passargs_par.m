
function thebigfunction = bilayer_passargs(jobid, framespace, b_val, stiff_val, e_val, alpha1, alpha2, alpha12, time_to_reach, dt, damp_val, damp_falloff, N, lambda, magnon_amp_real, mag_freq, t_final, Lspace, z_1, z_2)

[status,msg] = mkdir(jobid)
cd(jobid)

% open file to log outputs
fileID = fopen('skyrmions.out','w');
diary skyrmions.out


% begin parallelization
cpus = feature('numcores')
p = parpool(cpus);

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

% initialize all vertical
m_init(:,:,1,1)=zeros(N,N);
m_init(:,:,2,1)=zeros(N,N);
m_init(:,:,3,1)=ones(N,N);
m_init(:,:,1,2)=zeros(N,N);
m_init(:,:,2,2)=zeros(N,N);
m_init(:,:,3,2)=ones(N,N);
% initialize skyrmion according to QHMF 35
n=1;
z_0=0;
omega = ((xx + yy*1i - z_1)/lambda);
omega2 = ((xx + yy*1i - z_2)/lambda);
%omega2 = ((yy + xx*1i - z_0)/lambda^2).^(-1); %flips charge
%omega = (((xx + yy*1i - z_0*ones(N,N))/lambda).^n).*(((xx + yy*1i + z_0*ones(N,N))/lambda).^1);
m_init(:,:,1,1)=4*real(omega)./((abs(omega)).^2+4);
m_init(:,:,2,1)=4*imag(omega)./((abs(omega)).^2+4);
m_init(:,:,3,1)=((abs(omega)).^2-4)./((abs(omega)).^2+4);
m_init(:,:,1,2)=4*real(omega2)./((abs(omega2)).^2+4);
m_init(:,:,2,2)=4*imag(omega2)./((abs(omega2)).^2+4);
m_init(:,:,3,2)=((abs(omega2)).^2-4)./((abs(omega2)).^2+4);
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
clear m;

%m_init = m_init./(sqrt(sum(m_init.^2,3))); % Renormalize

"initialized"



%custom parameters, all positive!
damp_mat_real=zeros(N,N,3,2);
for i=1:N
    for j=1:N
        for k=1:3
            %damp_mat(i,N,k,l)=damp_val;   %one line
            damp_mat_real(i,j,k,1)=damp_val*exp((j-N)/damp_falloff);  %exponential L1
            %damp_mat(i,j,k,2)=damp_val*exp((j-N)/damp_falloff);  %exponential L2
            %damp_mat(i,j,k,l)=damp_val;   %const everywhere
        end
    end
end


% initialize Coulomb distance matrix
dist_x=zeros(N,N,N,N);
dist_y=zeros(N,N,N,N);
energy_dist=zeros(N,N,N,N);
energy_distL=zeros(N,N,N,N);
for i = 1:N
    for j = 1:N
        energy_dist(i,j,:,:) = ((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^-0.5;
        energy_dist(i,j,i,j) = 0;
        energy_distL(i,j,:,:) = ((xx-xx(i,j)).^2+(yy-yy(i,j)).^2 + Lspace^2).^-0.5;
        %energy_distL(i,j,i,j) = 0;
    end
end
center_dist = energy_dist(2:N-1,2:N-1,:,:);
center_distL = energy_distL(2:N-1,2:N-1,:,:);


t=0;
t_ind=1;
El_freq = 5 *2*pi/t_final; %num of cycles * 2pi*t_final
El_rad = 10;
V_ext = zeros(N,N,2);
%V_ext(:,:,1) = e_val*exp(-(xx.^2+yy.^2)/(2*(El_rad^2))); %Gaussian
%V_ext(:,:,2) = e_val*exp(-(xx.^2+yy.^2)/(2*(El_rad^2)));
rho = pontryagin_bilayer(m_init);
rho_init=rho;
alphadouble=zeros(1,1,2);
alphadouble(1,1,1)=alpha1;
alphadouble(1,1,2)=alpha2;
center_rho = rho(2:N-1,2:N-1,:);
for i = 1:N  %exact match
    for j = 1:N
        V_ext(i,j,:) = -e_val*sum(sum(alphadouble.*center_dist(:,:,i,j).*center_rho + alpha12*center_distL(:,:,i,j).*center_rho(:,:,[2;1])));
    end
end


"Curve fit potential"

gaussEqn = 'a*exp(-((x-b)/c)^2)+d'
thefit=fit(xx(:,N/2,1),V_ext(:,N/2,1),gaussEqn)
plot(thefit,xx(:,N/2,1),V_ext(:,N/2,1))
legend("Coulomb Potential","Gaussian fit")
drawnow
saveas(gcf,"Vfit")
fit_vals = coeffvalues(thefit)
V_ext(:,:,1) = e_val*fit_vals(1)*exp(-((xx-fit_vals(2)).^2+(yy-fit_vals(2)).^2)/(fit_vals(3))); %Gaussian
V_ext(:,:,2) = e_val*fit_vals(1)*exp(-((xx-fit_vals(2)).^2+(yy-fit_vals(2)).^2)/(fit_vals(3))); 


kappax=2*pi*0.1;
kappay=2*pi*0.11;

Q_top_list = []; % store values for final plot
E_B_list = [];
E_LL_list = [];
E_C_list = [];
E_C12_list = [];
E_El_list = [];
E_eff_list = [];
E_loss1_list = [];
E_loss1=0;
E_loss2_list = [];
E_loss2=0;
Spin_list = [];
mean_x_list1 = [];
mean_y_list1 = [];
stdev_list1 = [];
mean_x_list2 = [];
mean_y_list2 = [];
stdev_list2 = [];




while t<t_final
    tic
    m = m_init;
    rho = pontryagin_bilayer(m);
    %rho_avg = (rho(:,:)+rho(mod(-1:N-2,N)+1,:)+rho(:,mod(-1:N-2,N)+1)+rho(mod(-1:N-2,N)+1,mod(-1:N-2,N)+1))/4;  %oops still periodic

    % set electric field
    % oscillating
    %El_x = e_val*cos(El_freq*t);
    %El_y = e_val*cos(0.8*El_freq*t);
    % Gaussian
    %El_x = e_val/(El_rad^2)*yy.*exp(-(xx.^2+yy.^2)/(2*(El_rad^2)));
    %El_y = e_val/(El_rad^2)*xx.*exp(-(xx.^2+yy.^2)/(2*(El_rad^2)));

    % turn on stuff once time passes
    if t<time_to_reach
        magnon_amp=0;
        damp_mat=zeros(N,N,3,2);
    else
        magnon_amp=magnon_amp_real;
        damp_mat=damp_mat_real;
    end

    % check conserved quantities
    Q_top=sum(sum(rho(1:N-1,1:N-1,:),2),1)*dx*dy;  % topological charge
    Q_top=permute(Q_top,[4,3,2,1]);
    S_x=sum(sum(m(:,:,1,:),2),1);      
    S_y=sum(sum(m(:,:,2,:),2),1);      
    S_z=sum(sum(m(:,:,3,:),2),1); 
    S_x=permute(S_x,[4,3,2,1]);
    S_y=permute(S_y,[4,3,2,1]);
    S_z=permute(S_z,[4,3,2,1]);

    m_dx=(m(2:N,1:N-1,:,:)-m(1:N-1,1:N-1,:,:))/(dx); 
    m_dy=(m(1:N-1,2:N,:,:)-m(1:N-1,1:N-1,:,:))/(dy);

    B_field_thisstep=B_eff_bilayer(m,b_val,stiff_val,V_ext,alpha1,alpha2,alpha12,center_dist,center_distL);

    

    
    % mean and stdev
    mean_x1 = sum(sum(rho(:,:,1).*xx))/Q_top(1)*dx*dy;
    mean_y1 = sum(sum(rho(:,:,1).*yy))/Q_top(1)*dx*dy;
    variance1 = sum(sum(rho(:,:,1).*((xx-mean_x1).^2+(yy-mean_y1).^2)))/Q_top(1)*dx*dy;
    st_dev1 = sqrt(variance1);
    mean_x_list1(length(mean_x_list1)+1)=mean_x1;
    mean_y_list1(length(mean_y_list1)+1)=mean_y1;
    stdev_list1(length(stdev_list1)+1)=st_dev1;

    mean_x2 = sum(sum(rho(:,:,2).*xx))/Q_top(2)*dx*dy;
    mean_y2 = sum(sum(rho(:,:,2).*yy))/Q_top(2)*dx*dy;
    variance2 = sum(sum(rho(:,:,2).*((xx-mean_x2).^2+(yy-mean_y2).^2)))/Q_top(2)*dx*dy;
    st_dev2 = sqrt(variance2);
    mean_x_list2(length(mean_x_list2)+1)=mean_x2;
    mean_y_list2(length(mean_y_list2)+1)=mean_y2;
    stdev_list2(length(stdev_list2)+1)=st_dev2;

    
    %plot
    if t_ind==1
        climsave_z = [-1 1];
        contour(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),rho(1:N-1,1:N-1,1),10)
        climsave1=clim;
        %climsave1=[-0.0001 0.0001];
        contour(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),rho(1:N-1,1:N-1,2),10)
        climsave2=clim;
    end

    if mod(t_ind,framespace)==1

        pc=pcolor(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),rho(1:N-1,1:N-1,1)); % color plot
        %pc=pcolor(xx,yy,rho);
        pc.EdgeColor='none';
        clim(climsave1);
        colorbar
        title("L1 Topological charge at t="+num2str(t_ind*dt));
        drawnow
        %saveas(gcf,"zz1_charge_frame"+string(t_ind)+".png")
        exportgraphics(gcf,"L1_top.gif",'Append',true)
        %movefile("zz1_charge_frame"+string(t_ind)+".png",jobid);
        "L1charge"


        pc=pcolor(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),m(1:N-1,1:N-1,3,1)); % color plot
        pc.EdgeColor='none';
        clim(climsave_z);
        thebar=colorbar;
        thebar.Label.String = 'z component';
        title("L1 Magnetization field components at t="+num2str(t_ind*dt));
        hold on
        quiver(xx,yy,m(:,:,1,1),m(:,:,2,1),'r')
        quiver(mean_x1,mean_y1,st_dev1/sqrt(2),st_dev1/sqrt(2),'g')
        hold off
        drawnow
        %saveas(gcf,"zz1_quiver_frame"+string(t_ind)+".png")
        exportgraphics(gcf,"L1_mag.gif",'Append',true)
        %movefile("zz1_quiver_frame"+string(t_ind)+".png",jobid);
        "L1quiver"



        pc=pcolor(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),rho(1:N-1,1:N-1,2)); % color plot
        %pc=pcolor(xx,yy,rho);
        pc.EdgeColor='none';
        clim(climsave2);
        colorbar
        title("L2 Topological charge at t="+num2str(t_ind*dt));
        drawnow
        %saveas(gcf,"zz2_charge_frame"+string(t_ind)+".png")
        exportgraphics(gcf,"L2_top.gif",'Append',true)
        %movefile("zz2_charge_frame"+string(t_ind)+".png",jobid);
        "L2charge"

        pc=pcolor(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),rho(1:N-1,1:N-1,2)-rho_init(1:N-1,1:N-1,2)); % color plot
        %pc=pcolor(xx,yy,rho);
        pc.EdgeColor='none';
        clim([-0.001 0.001]);
        colorbar
        title("L2 Topological charge diff at t="+num2str(t_ind*dt));
        drawnow
        %saveas(gcf,"zz2_charge_frame"+string(t_ind)+".png")
        exportgraphics(gcf,"L2_top_diff.gif",'Append',true)
        %movefile("zz2_charge_frame"+string(t_ind)+".png",jobid);
        "L2charge_diff"


        pc=pcolor(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),m(1:N-1,1:N-1,3,2)); % color plot
        pc.EdgeColor='none';
        clim(climsave_z);
        thebar=colorbar;
        thebar.Label.String = 'z component';
        title("L2 Magnetization field components at t="+num2str(t_ind*dt));
        hold on
        quiver(xx,yy,m(:,:,1,2),m(:,:,2,2),'r')
        quiver(mean_x2,mean_y2,st_dev2/sqrt(2),st_dev2/sqrt(2),'g')
        hold off
        drawnow
        %saveas(gcf,"zz2_quiver_frame"+string(t_ind)+".png")
        exportgraphics(gcf,"L2_mag.gif",'Append',true)
        %movefile("zz2_quiver_frame"+string(t_ind)+".png",jobid);
        "L2quiver"


        plot(yy(N/2,:,1),V_ext(N/2,:,1),yy(N/2,:,1),rho(N/2,:,1)*thefit(1)/rho_init(N/2,N/2,1))
        legend("V_ext","rho")
        drawnow
        exportgraphics(gcf,"Vfit_frame.gif",'Append',true)
        "Vfit_frame"

    end
    
    % dynamics!

    % RK4
    %B_field_thisstep=B_eff_bilayer(m,b_val,stiff_val,El_x,El_y,alpha1,alpha2,alpha12,center_dist,center_distL);
    B_mod_thisstep=B_field_thisstep;
    m_k1 = m(:,:,mod(1:3,3)+1,:).*B_field_thisstep(:,:,mod(2:4,3)+1,:)-m(:,:,mod(2:4,3)+1,:).*B_field_thisstep(:,:,mod(1:3,3)+1,:) - damp_mat.*((m(:,:,1,:).*B_mod_thisstep(:,:,1,:) + m(:,:,2,:).*B_mod_thisstep(:,:,2,:) + m(:,:,3,:).*B_mod_thisstep(:,:,3,:)).*m-B_mod_thisstep);

    m_k2arg=m+dt/2*m_k1;
    B_field=B_eff_bilayer(m_k2arg,b_val,stiff_val,V_ext,alpha1,alpha2,alpha12,center_dist,center_distL);
    B_mod=B_field;
    m_k2 = m_k2arg(:,:,mod(1:3,3)+1,:).*B_field(:,:,mod(2:4,3)+1,:)-m_k2arg(:,:,mod(2:4,3)+1,:).*B_field(:,:,mod(1:3,3)+1,:) - damp_mat.*((m_k2arg(:,:,1,:).*B_mod(:,:,1,:) + m_k2arg(:,:,2,:).*B_mod(:,:,2,:) + m_k2arg(:,:,3,:).*B_mod(:,:,3,:)).*m_k2arg-B_mod);

    m_k3arg=m+dt/2*m_k2;
    B_field=B_eff_bilayer(m_k3arg,b_val,stiff_val,V_ext,alpha1,alpha2,alpha12,center_dist,center_distL);
    B_mod=B_field;
    m_k3 = m_k3arg(:,:,mod(1:3,3)+1,:).*B_field(:,:,mod(2:4,3)+1,:)-m_k3arg(:,:,mod(2:4,3)+1,:).*B_field(:,:,mod(1:3,3)+1,:) - damp_mat.*((m_k3arg(:,:,1,:).*B_mod(:,:,1,:) + m_k3arg(:,:,2,:).*B_mod(:,:,2,:) + m_k3arg(:,:,3,:).*B_mod(:,:,3,:)).*m_k3arg-B_mod);

    m_k4arg=m+dt*m_k3;
    B_field=B_eff_bilayer(m_k4arg,b_val,stiff_val,V_ext,alpha1,alpha2,alpha12,center_dist,center_distL);
    B_mod=B_field;
    m_k4 = m_k4arg(:,:,mod(1:3,3)+1,:).*B_field(:,:,mod(2:4,3)+1,:)-m_k4arg(:,:,mod(2:4,3)+1,:).*B_field(:,:,mod(1:3,3)+1,:) - damp_mat.*((m_k4arg(:,:,1,:).*B_mod(:,:,1,:) + m_k4arg(:,:,2,:).*B_mod(:,:,2,:) + m_k4arg(:,:,3,:).*B_mod(:,:,3,:)).*m_k4arg-B_mod);

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
            m(i,1,1,1) = magnon_amp*sin(mag_freq*t);
            m(i,1,2,1) = magnon_amp*cos(mag_freq*t);
            %diagonal:
            %m(i,1,1,1) = magnon_amp*sin(mag_freq*t+i*0.5*pi);
            %m(i,1,2,1) = magnon_amp*cos(mag_freq*t+i*0.5*pi);
            %charged:
            %m(i,1,1,1) = cos(kappay.*yy(i,1)-mag_freq*t)*sin(kappax.*xx(i,1));
            %m(i,1,2,1) = cos(kappay.*yy(i,1)-mag_freq*t)*cos(kappax.*xx(i,1));
            %m(i,1,3,1) = sin(kappay.*yy(i,1)-mag_freq*t);
        end
    end

    % Renormalize
    m = m./(sqrt(sum(m.^2,3))); 

    
    %save energies 
    % Coulomb energy
    %tic
    if alpha1 ~= 0
        coulomb_energy_field = zeros(N,N,2);
        for i = 1:N
            for j = 1:N
                coulomb_energy_field = coulomb_energy_field + rho.*rho(i,j,:).*energy_dist(:,:,i,j);
            end
        end
        %coulomb_energy_field = distributed(tensorprod(rho,energy_dist,[1 2])).*rho;
        sum_coulombfield = sum(sum(coulomb_energy_field));
        sum_coulombfield=permute(sum_coulombfield,[3,2,1]);
        E_C = [alpha1;alpha2]*2*pi.*sum_coulombfield;
    else
        E_C = zeros(2,1);
    end
    %toc

    % interlayer Coulomb energy
    %tic
    if alpha12 ~= 0
        coulomb_energy_field = zeros(N,N);
        for i = 1:N
            for j = 1:N
                coulomb_energy_field = coulomb_energy_field + rho(:,:,1).*rho(i,j,2).*energy_dist(:,:,i,j);
            end
        end
        %coulomb_energy_field = distributed(tensorprod(rho,energy_dist,[1 2])).*rho;
        E_C12 = alpha12*2*pi*sum(sum(coulomb_energy_field));
    else
        E_C12 = 0;
    end
    %toc

    E_B = -b_val*(S_z-N*N);    % B energy
    E_LL = sum(sum(stiff_val/2 * (sum(m_dx.^2+m_dy.^2))));   % stiffness energy
    E_LL=permute(E_LL,[4,3,2,1]);
    E_El = 4*pi*sum(sum(rho.*V_ext)); % electric field potential energy (Gaussian)
    E_El=permute(E_El,[4,3,2,1]);
    %general energy loss %warning: backward derivative
    E_loss1 = E_loss1 + sum(sum(sum(B_field_thisstep.*dmdt)))*dt; 
    %friction due to damping
    E_loss2 = E_loss2 + permute(sum(-damp_mat(:,:,1,:).*(dot(B_field_thisstep,m_init,3).^2-dot(B_field_thisstep,B_field_thisstep,3)),[1 2]),[4 3 2 1])*dt;

    Q_top_list(length(Q_top_list)+1,:)=Q_top;
    E_B_list(length(E_B_list)+1,:)=E_B;
    E_LL_list(length(E_LL_list)+1,:)=E_LL;
    E_C_list(length(E_C_list)+1,:)=E_C;
    E_C12_list(length(E_C12_list)+1)=E_C12;
    E_El_list(length(E_El_list)+1,:)=E_El;
    E_loss1_list(length(E_loss1_list)+1,:)=E_loss1;
    E_loss2_list(length(E_loss2_list)+1,:)=E_loss2;
    Spin_list(length(Spin_list)+1,:,:)=[S_x S_y S_z];
    

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


E_C12_list=permute(E_C12_list,[2,1]);
E_total_list = E_B_list+E_LL_list+E_C_list+E_El_list;
E_both_list = (E_total_list(1:t_final/dt,1)+E_total_list(1:t_final/dt,2)+E_C12_list);

% plotting all the energies
plot((1:length(Q_top_list))*dt,Q_top_list(:,1))
drawnow
saveas(gcf,"fig_Q_top1")
saveas(gcf,"fig_Q_top1.png")
plot((1:length(Q_top_list))*dt,Q_top_list(:,2))
drawnow
saveas(gcf,"fig_Q_top2")
saveas(gcf,"fig_Q_top2.png")

plot((1:length(E_B_list(:,1)))*dt,E_B_list(:,1),(1:length(E_LL_list(:,1)))*dt,E_LL_list(:,1),(1:length(E_C_list(:,1)))*dt,E_C_list(:,1),(1:length(E_C12_list))*dt,E_C12_list,(1:length(E_El_list(:,1)))*dt,E_El_list(:,1),(1:length(E_total_list(:,1)))*dt,E_total_list(:,1),(1:length(E_loss1_list(:,1)))*dt,E_loss1_list(:,1),(1:length(E_loss2_list(:,1)))*dt,E_loss2_list(:,1))
legend("Zeeman","Stiffness","Coulomb","Interlayer Coulomb","Electric","Total","Change","Friction")
title("L1 Energies");
drawnow
saveas(gcf,"fig_Energies1")
saveas(gcf,"fig_Energies1.png")
plot((1:length(E_B_list(:,2)))*dt,E_B_list(:,2),(1:length(E_LL_list(:,2)))*dt,E_LL_list(:,2),(1:length(E_C_list(:,2)))*dt,E_C_list(:,2),(1:length(E_C12_list))*dt,E_C12_list,(1:length(E_El_list(:,2)))*dt,E_El_list(:,2),(1:length(E_total_list(:,2)))*dt,E_total_list(:,2),(1:length(E_loss1_list(:,2)))*dt,E_loss1_list(:,2),(1:length(E_loss2_list(:,2)))*dt,E_loss2_list(:,2))
legend("Zeeman","Stiffness","Coulomb","Interlayer Coulomb","Electric","Total","Change","Friction")
title("L2 Energies");
drawnow
saveas(gcf,"fig_Energies2")
saveas(gcf,"fig_Energies2.png")
plot((1:length(E_total_list(:,1)))*dt,E_total_list(:,1),(1:length(E_total_list(:,2)))*dt,E_total_list(:,2),(1:length(E_C12_list))*dt,E_C12_list,(1:length(E_C12_list))*dt,E_both_list)
legend("L1","L2","Interlayer","Both")
title("Energies for full system");
drawnow
saveas(gcf,"fig_Energies_both")
saveas(gcf,"fig_Energies_both.png")

plot((1:length(Spin_list(:,1,1)))*dt,Spin_list(:,1,1),(1:length(Spin_list(:,1,2)))*dt,Spin_list(:,1,2),(1:length(Spin_list(:,1,3)))*dt,Spin_list(:,1,3))
legend("S_x","S_y","S_z")
drawnow
saveas(gcf,"fig_Spin_components1")
saveas(gcf,"fig_Spin_components1.png")
plot((1:length(Spin_list(:,2,1)))*dt,Spin_list(:,2,1),(1:length(Spin_list(:,2,2)))*dt,Spin_list(:,2,2),(1:length(Spin_list(:,2,3)))*dt,Spin_list(:,2,3))
legend("S_x","S_y","S_z")
drawnow
saveas(gcf,"fig_Spin_components2")
saveas(gcf,"fig_Spin_components2.png")

plot((1:length(mean_x_list1))*dt,mean_x_list1,(1:length(mean_y_list1))*dt,mean_y_list1,(1:length(stdev_list1))*dt,stdev_list1)
legend("mean_x","mean_y","stdev")
drawnow
saveas(gcf,"fig_Pos_radius1")
saveas(gcf,"fig_Pos_radius1.png")
plot((1:length(mean_x_list2))*dt,mean_x_list2,(1:length(mean_y_list2))*dt,mean_y_list2,(1:length(stdev_list2))*dt,stdev_list2)
legend("mean_x","mean_y","stdev")
drawnow
saveas(gcf,"fig_Pos_radius2")
saveas(gcf,"fig_Pos_radius2.png")


%pc=pcolor(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),rho(1:N-1,1:N-1,1)); % color plot
%pc.EdgeColor='none';
%clim(climsave1);
%colorbar
%drawnow
%saveas(gcf,"fig_charge_frame1"+string(t_ind))
%saveas(gcf,"fig_charge_frame1"+string(t_ind)+".png")
%pc=pcolor(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),rho(1:N-1,1:N-1,2)); % color plot
%pc.EdgeColor='none';
%clim(climsave2);
%colorbar
%drawnow
%saveas(gcf,"fig_charge_frame2"+string(t_ind))
%saveas(gcf,"fig_charge_frame2"+string(t_ind)+".png")


%pc=pcolor(xx(1:N-1,1:N-1),yy(1:N-1,1:N-1),m(1:N-1,1:N-1,3,1)); % color plot
%pc.EdgeColor='none';
%clim(climsave_z2);
%colorbar
%hold on
%quiver(xx,yy,m(:,:,1),m(:,:,2),'r')
%quiver(mean_x,mean_y,st_dev/sqrt(2),st_dev/sqrt(2),'g')
%hold off
%drawnow
%saveas(gcf,"fig_quiver_frame"+string(t_ind))
%saveas(gcf,"fig_quiver_frame"+string(t_ind)+".png")

save("energyvars.mat","E_total_list","E_B_list","E_LL_list","E_C_list","E_C12_list","E_El_list")
save("m_final.mat","m")


"drew frames of evolution"



fclose(fileID);
delete(p);
cd ..

end



function B_field = B_eff_bilayer(m_arg,b_val,stiff_val,V_ext,alpha1,alpha2,alpha12,center_dist,other_dist)
    N=length(m_arg);
    N_1=N-1;
    dx=1;
    dy=1;
    rho = pontryagin_bilayer(m_arg);
    alphadouble=zeros(1,1,2);
    alphadouble(1,1,1)=alpha1;
    alphadouble(1,1,2)=alpha2;

    %Zeeman
    B_Zeeman = zeros(N,N,3,2);
    B_Zeeman(:,:,3,:)=b_val*ones(N,N,1,2);

    %Stiffness
    B_stiffness=zeros(N,N,3,2);
    %Centered stiffness:
    B_stiffness(2:N-1,2:N-1,:,:) = (m_arg(1:N-2,2:N-1,:,:)+m_arg(3:N,2:N-1,:,:)+m_arg(2:N-1,1:N-2,:,:)+m_arg(2:N-1,3:N,:,:)-4*m_arg(2:N-1,2:N-1,:,:))/(dx^2)*stiff_val;
    %Wrapped stiffness:
    %B_stiffness(:,2:N-1,:,1) = (m_arg(mod(-1:N-2,N)+1,2:N-1,:,1)+m_arg(mod(1:N,N)+1,2:N-1,:,1)+m_arg(:,1:N-2,:,1)+m_arg(:,3:N,:,1)-4*m_arg(:,2:N-1,:,1))/(dx^2)*stiff_val;
    %B_stiffness(2:N-1,2:N-1,:,2) = (m_arg(1:N-2,2:N-1,:,2)+m_arg(3:N,2:N-1,:,2)+m_arg(2:N-1,1:N-2,:,2)+m_arg(2:N-1,3:N,:,2)-4*m_arg(2:N-1,2:N-1,:,2))/(dx^2)*stiff_val;
    
    %Coulomb
    B_coulomb = zeros(N,N,3,2);
    centered_int = zeros(N,N,2);
    full_int=zeros(N,N,2,3);
    full_intx=zeros(N,N,2,3);
    full_inty=zeros(N,N,2,3);
    full_intxy=zeros(N,N,2,3);

    %could cause slowdown...
    %center_dist = energy_dist(2:N-1,2:N-1,:,:);
    %sliced_dist = energy_dist(1:N-1,1:N-1,:,:);
    center_rho = rho(2:N_1,2:N_1,:);
    %sliced_rho = rho(1:N_1,1:N_1,:);

    %tic
    if alpha12 ~=0
        for i = 1:N
            for j = 1:N
                %centered_int(i,j,1) = sum(sum(alpha1*center_dist(:,:,i,j).*center_rho(:,:,1) + alpha12*other_dist(:,:,i,j).*center_rho(:,:,2)));
                %centered_int(i,j,2) = sum(sum(alpha2*center_dist(:,:,i,j).*center_rho(:,:,2) + alpha12*other_dist(:,:,i,j).*center_rho(:,:,1)));
                centered_int(i,j,:) = sum(sum(alphadouble.*center_dist(:,:,i,j).*center_rho + alpha12*other_dist(:,:,i,j).*center_rho(:,:,[2;1])));
            end
        end
    end
    %toc

    centered_int = centered_int + V_ext;
   
    for k=1:3
        full_int(:,:,:,k)=centered_int;
        full_intx(2:N,:,:,k)=centered_int(1:N-1,:,:);
        full_inty(:,2:N,:,k)=centered_int(:,1:N-1,:);
        full_intxy(2:N,2:N,:,k)=centered_int(1:N-1,1:N-1,:);
    end

    %switch xyz to position 3, layer to position 4
    full_int=permute(full_int,[1,2,4,3]);
    full_intx=permute(full_intx,[1,2,4,3]);
    full_inty=permute(full_inty,[1,2,4,3]);
    full_intxy=permute(full_intxy,[1,2,4,3]);


    B_coulomb=-4*pi*(full_int.*pontryagin_deriv_bi(m_arg,1,1) - full_intx.*pontryagin_deriv_bi(m_arg,-1,1) - full_inty.*pontryagin_deriv_bi(m_arg,1,-1) + full_intxy.*pontryagin_deriv_bi(m_arg,-1,-1));

    %B_electric = (El_x+El_y).*pontryagin_deriv_bi(m_arg,1,1) + (-El_x+El_y).*pontryagin_deriv_bi(m_arg,-1,1) + (El_x-El_y).*pontryagin_deriv_bi(m_arg,1,-1) + (-El_x-El_y).*pontryagin_deriv_bi(m_arg,-1,-1);

    %B_field = B_Zeeman+B_stiffness+B_coulomb+B_electric;
    B_field = B_Zeeman+B_stiffness+B_coulomb;

end

function rho_old = pontryagin_bilayer(m)
    N=length(m);
    dx=1;
    dy=1;
    tripleL = zeros(N,N,3,2);
    denomL = zeros(N,N,3,2);
    tripleR = zeros(N,N,3,2);
    denomR = zeros(N,N,3,2);
    m_x=m(mod(1:N,N)+1,1:N,:,:);
    m_y=m(1:N,mod(1:N,N)+1,:,:);
    m_xy=m(mod(1:N,N)+1,mod(1:N,N)+1,:,:);
    tripleL = dot(m,cross(m_x,m_xy,3),3);
    denomL = 1 + sum(m.*m_x,3) + sum(m.*m_xy,3) + sum(m_x.*m_xy,3);
    tripleR = dot(m,cross(m_xy,m_y,3),3);
    denomR = 1 + sum(m.*m_xy,3) + sum(m.*m_y,3) + sum(m_xy.*m_y,3);
    %puts layer in index 3
    rho_old(:,:,1) = 4*(atan2(tripleL(:,:,1,1),denomL(:,:,1,1))+atan2(tripleR(:,:,1,1),denomR(:,:,1,1)))/(8*pi*dx*dy);
    rho_old(:,:,2) = 4*(atan2(tripleL(:,:,1,2),denomL(:,:,1,2))+atan2(tripleR(:,:,1,2),denomR(:,:,1,2)))/(8*pi*dx*dy);
end


function rho_prime = pontryagin_deriv_bi(m,x_dir,y_dir)
    N=length(m);
    dx=1;
    dy=1;
    m_x=m(2+x_dir:N-1+x_dir,2:N-1,:,:);
    m_y=m(2:N-1,2+y_dir:N-1+y_dir,:,:);
    m_xy=m(2+x_dir:N-1+x_dir,2+y_dir:N-1+y_dir,:,:);
    m_c=m(2:N-1,2:N-1,:,:);
    triple = dot(m_c,cross(m_x,m_xy,3),3);
    denom = 1 + sum(m_c.*m_x,3) + sum(m_c.*m_xy,3) + sum(m_x.*m_xy,3);
    triple2 = dot(m_c,cross(m_xy,m_y,3),3);
    denom2 = 1 + sum(m_c.*m_xy,3) + sum(m_c.*m_y,3) + sum(m_xy.*m_y,3);
    rho_prime1 = zeros(N,N,3,2);
    rho_prime1(2:N-1,2:N-1,:,:) = 1/(pi*dx*dy) * (cross(m_x,m_xy,3).*denom - (m_x+m_xy).*triple ) ./ (triple.^2+denom.^2);
    rho_prime2 = zeros(N,N,3,2);
    rho_prime2(2:N-1,2:N-1,:,:) = 1/(pi*dx*dy) * (cross(m_xy,m_y,3).*denom2 - (m_xy+m_y).*triple2 ) ./ (triple2.^2+denom2.^2);
    rho_prime=(rho_prime1+rho_prime2)*0.5;
    %rho_prime(:,:,1)=(rho_prime1(:,:,1,1)+rho_prime2(:,:,1,1))*0.5;
    %rho_prime(:,:,2)=(rho_prime1(:,:,1,2)+rho_prime2(:,:,1,2))*0.5;
end