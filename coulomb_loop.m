
function coulomb = coulomb_loop(m_arg,dist_x,dist_y,rho)
    N=length(m_arg(:,1,1));
    m_x=(m_arg(3:N,2:N-1,:)-m_arg(1:N-2,2:N-1,:))/(2);
    m_y=(m_arg(2:N-1,3:N,:)-m_arg(2:N-1,1:N-2,:))/(2);
    centered_m = zeros(N-2,N-2,3);
    centered_rho = (rho(2:N-1,2:N-1)+rho(1:N-2,2:N-1)+rho(2:N-1,1:N-2)+rho(1:N-2,1:N-2))/4;
    parfor i = 1:N-2
        for j = 1:N-2
            centered_m = centered_m + (dist_x(2:N-1,2:N-1,:,i,j).*m_y - dist_y(2:N-1,2:N-1,:,i,j).*m_x )*centered_rho(i,j);
        end
    end
    coulomb = zeros(N,N,3);
    coulomb(2:N-1,2:N-1,:)=centered_m;
end