
function B_field = B_eff(m_arg,b_val,e_val,alpha_val,dist_x,dist_y,rho)
    N=length(m_arg(:,1,1));
    dx=1;
    dy=1;

    %Zeeman
    B_Zeeman = zeros(N,N,3);
    B_Zeeman(:,:,3)=b_val*ones(N,N);

    %Stiffness
    B_stiffness=zeros(N,N,3);
    B_stiffness(2:N-1,2:N-1,:) = (m_arg(1:N-2,2:N-1,:)+m_arg(3:N,2:N-1,:)+m_arg(2:N-1,1:N-2,:)+m_arg(2:N-1,3:N,:)-4*m_arg(2:N-1,2:N-1,:))/(dx^2);

    %Coulomb
    m_x=(m_arg(3:N,2:N-1,:)-m_arg(1:N-2,2:N-1,:))/(2);  %N-2xN-2x3
    m_y=(m_arg(2:N-1,3:N,:)-m_arg(2:N-1,1:N-2,:))/(2);
    centered_intx = zeros(N-2,N-2);
    centered_inty = zeros(N-2,N-2);
    centered_rho = (rho(2:N-1,2:N-1)+rho(1:N-2,2:N-1)+rho(2:N-1,1:N-2)+rho(1:N-2,1:N-2))/4;
    for i = 2:N-1
        for j = 2:N-1
            centered_intx = centered_intx + dist_x(2:N-1,2:N-1,i,j).*centered_rho(i-1,j-1);
            centered_inty = centered_inty + dist_y(2:N-1,2:N-1,i,j).*centered_rho(i-1,j-1);
        end
    end
    full_intx=zeros(N-2,N-2,3);
    full_inty=zeros(N-2,N-2,3);
    for k = 1:3
        full_intx(:,:,k)=centered_intx;
        full_inty(:,:,k)=centered_inty;
    end
    B_coulomb=zeros(N,N,3);
    B_coulomb(2:N-1,2:N-1,:)=alpha_val*((m_arg(2:N-1,2:N-1,mod(1:3,3)+1).*m_y(:,:,mod(2:4,3)+1)-m_arg(2:N-1,2:N-1,mod(2:4,3)+1).*m_y(:,:,mod(1:3,3)+1)).*full_intx - (m_arg(2:N-1,2:N-1,mod(1:3,3)+1).*m_x(:,:,mod(2:4,3)+1)-m_arg(2:N-1,2:N-1,mod(2:4,3)+1).*m_x(:,:,mod(1:3,3)+1)).*full_inty);

    %need to add e field

    B_field = B_Zeeman+B_stiffness+B_coulomb;

end