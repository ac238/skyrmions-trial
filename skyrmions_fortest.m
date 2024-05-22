
function time_for = skyrmions_fortest(N)
    % create bounds of graph
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
    
    
    %custom parameters, all positive!
    b_val = 0.13;
    stiff_val = 0;
    e_val = 0;
    alpha_val = 0.7;
    
    % initialize Coulomb distance matrix
    if alpha_val ~= 0
        dist_x=zeros(N,N,3,N,N);
        dist_y=zeros(N,N,3,N,N);
        for i = 1:N
            for j = 1:N
                for k = 1:3
                    dist_x(:,:,k,i,j) = (xx-xx(i,j))./(((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^1.5);
                    dist_x(i,j,k,i,j) = 0;
                    dist_y(:,:,k,i,j) = (yy-yy(i,j))./(((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^1.5);
                    dist_y(i,j,k,i,j) = 0;
                end
                energy_dist(:,:,i,j) = ((xx-xx(i,j)).^2+(yy-yy(i,j)).^2).^-0.5;
                energy_dist(i,j,i,j) = 0;
            end
        end
    end
    
    
    
    
    t=0;
    t_final=1;
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
    
        El = e_val*cos(El_freq*t);
        m = m - El*dt/6*(m_k1 + 2*m_k2 + 2*m_k3 + m_k4);
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
            m_k1=coulomb_forloop(m,dist_x,dist_y,rho);
            m_k2arg = m + dt/2*m_k1;
            m_k2=coulomb_forloop(m_k2arg,dist_x,dist_y,rho);
            m_k3arg = m + dt/2*m_k2;
            m_k3=coulomb_forloop(m_k3arg,dist_x,dist_y,rho);
            m_k4arg = m + dt*m_k3;
            m_k4=coulomb_forloop(m_k4arg,dist_x,dist_y,rho);
    
            m_RK4 = (m_k1+2*m_k2+2*m_k3+m_k4)/6;
            m = m - alpha_val*dt*(m_RK4+m_RK4(mod(-1:N-2,N)+1,:,:)+m_RK4(:,mod(-1:N-2,N)+1,:)+m_RK4(mod(-1:N-2,N)+1,mod(-1:N-2,N)+1,:))/4;
            m = m./(sqrt(sum(m.^2,3))); % Renormalize
    
    
            % Coulomb energy
            coulomb_energy_field = zeros(N,N);
            for i = 1:N
                for j = 1:N
                    coulomb_energy_field = coulomb_energy_field + rho.*rho(i,j).*energy_dist(:,:,i,j);
                end
            end
            E_C = alpha_val*2*pi*sum(sum(coulomb_energy_field));
        else
            E_C = 0;
        end
    
    
    
        % check conserved quantities
        Q_top_init = Q_top;
        Q_top=sum(sum(rho(1:N-1,1:N-1)))*dx*dy;  % topological charge
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
    
        E_B = -b_val*S_z;    % B energy
        E_LL = sum(sum(stiff_val/2 * (sum(m_dx.^2+m_dy.^2))));   % stiffness energy
        Q_top_list(length(Q_top_list)+1)=Q_top;
        E_B_list(length(E_B_list)+1)=E_B;
        E_LL_list(length(E_LL_list)+1)=E_LL;
        E_C_list(length(E_C_list)+1)=E_C;
        Spin_list(length(Spin_list)+1,:)=[S_x S_y S_z];
    
        % check norm is preserved
        m_norm=sqrt(sum(m.^2,3));
        
        % plot
        %quiver(xx,yy,m(:,:,1),m(:,:,2)) % full 2D vector field
        %quiver(xx(:,N/2),zeros(1,N),m(N/2,:,1),m(N/2,:,3))  % 1D slice
        %hold on
        %quiver(cent_of_mass_x,cent_of_mass_y,st_dev,0,'r')  %radius vector
        %contour(xx(1:N-1,1:N-1)-dx/2,yy(1:N-1,1:N-1)-dy/2,rho(1:N-1,1:N-1),10) % color plot
        %hold off
        %axis([xlow xhigh ylow yhigh])
        %title(t)
        %drawnow
        
        % reset for new loop
        m_init=m;
        t=t+dt;
        m_full(:,:,:,t_ind) = m;
        rho_full(:,:,t_ind) = rho;
    
        time_thisrun(t_ind)=toc;
        t_ind=t_ind+1;
    end
    time_for=min(time_thisrun);
end



