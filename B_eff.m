
function B_field = B_eff(m_arg,b_val,stiff_val,El_x,El_y,alpha_val,dist_x,dist_y,rho)
    N=length(m_arg(:,1,1));
    dx=1;
    dy=1;

    %Zeeman
    B_Zeeman = zeros(N,N,3);
    B_Zeeman(:,:,3)=b_val*ones(N,N);

    %Stiffness
    B_stiffness=zeros(N,N,3);
    B_stiffness(2:N-1,2:N-1,:) = (m_arg(1:N-2,2:N-1,:)+m_arg(3:N,2:N-1,:)+m_arg(2:N-1,1:N-2,:)+m_arg(2:N-1,3:N,:)-4*m_arg(2:N-1,2:N-1,:))/(dx^2)*stiff_val;

    %Coulomb
    centered_m_x=(m_arg(3:N,:,:)-m_arg(1:N-2,:,:))/(2);  %N-2xN-2x3
    centered_m_y=(m_arg(:,3:N,:)-m_arg(:,1:N-2,:))/(2);
    far_edge_m_x = (m_arg(N,:,:)-m_arg(N-1,:,:));
    near_edge_m_x = (m_arg(2,:,:)-m_arg(1,:,:));
    far_edge_m_y = (m_arg(:,N,:)-m_arg(:,N-1,:));
    near_edge_m_y = (m_arg(:,2,:)-m_arg(:,1,:));

    m_x = zeros(N,N,3);
    m_x(2:N-1,:,:)=centered_m_x;
    m_x(N,:,:)=far_edge_m_x;
    m_x(1,:,:)=near_edge_m_x;
    m_y = zeros(N,N,3);
    m_y(:,2:N-1,:)=centered_m_y;
    m_y(:,N,:)=far_edge_m_y;
    m_y(:,1,:)=near_edge_m_y;

    full_intx=zeros(N,N,3);
    full_inty=zeros(N,N,3);
    if alpha_val ~=0
        centered_rho = (rho(2:N-1,2:N-1)+rho(1:N-2,2:N-1)+rho(2:N-1,1:N-2)+rho(1:N-2,1:N-2))/4;
        centered_intx = zeros(N-2,N-2);
        centered_inty = zeros(N-2,N-2);
        for i = 2:N-1
            for j = 2:N-1
                %full_intx(2:N-1,2:N-1,:) = full_intx(2:N-1,2:N-1,:) + dist_x(2:N-1,2:N-1,i,j)*centered_rho(i-1,j-1);
                %full_inty(2:N-1,2:N-1,:) = full_inty(2:N-1,2:N-1,:) + dist_y(2:N-1,2:N-1,i,j)*centered_rho(i-1,j-1);
                centered_intx = centered_intx + dist_x(2:N-1,2:N-1,i,j).*centered_rho(i-1,j-1);
                centered_inty = centered_inty + dist_y(2:N-1,2:N-1,i,j).*centered_rho(i-1,j-1);
            end
        end
        %centered_intx = distributed(tensorprod(dist_x(2:N-1,2:N-1,2:N-1,2:N-1),centered_rho(1:N-2,1:N-2),[1 2]));
        %centered_inty = distributed(tensorprod(dist_y(2:N-1,2:N-1,2:N-1,2:N-1),centered_rho(1:N-2,1:N-2),[1 2]));
        for k = 1:3
            full_intx(2:N-1,2:N-1,k)=centered_intx;
            full_inty(2:N-1,2:N-1,k)=centered_inty;
        end
    end
    B_coulomb=((m_arg(:,:,mod(1:3,3)+1).*m_y(:,:,mod(2:4,3)+1)-m_arg(:,:,mod(2:4,3)+1).*m_y(:,:,mod(1:3,3)+1)).*(El_x*ones(N,N,3)+alpha_val*full_intx) - (m_arg(:,:,mod(1:3,3)+1).*m_x(:,:,mod(2:4,3)+1)-m_arg(:,:,mod(2:4,3)+1).*m_x(:,:,mod(1:3,3)+1)).*(El_y*ones(N,N,3)+alpha_val*full_inty));

    B_field = B_Zeeman+B_stiffness+B_coulomb;

end