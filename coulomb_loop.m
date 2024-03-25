
function coulomb = coulomb_loop(m_arg,dist_x,dist_y,rho)
    N=length(m_arg(:,1,1));
    m_x=(m_arg(3:N,2:N-1,:)-m_arg(2:N-1,2:N-1,:));
    m_y=(m_arg(2:N-1,3:N,:)-m_arg(2:N-1,2:N-1,:));
    centered_m = zeros(N-2,N-2,3);
        parfor i = 2:N-1
            for j = 2:N-1
                centered_m = centered_m + (dist_x(2:N-1,2:N-1,:,i-1,j-1).*m_y - dist_y(2:N-1,2:N-1,:,i-1,j-1).*m_x )*rho(i,j);
            end
        end
    coulomb = zeros(N,N,3);
    coulomb(2:N-1,2:N-1,:)=centered_m;
end