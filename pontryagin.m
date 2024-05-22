
function rho = pontryagin(m)
    N=length(m);
    dx=1;
    dy=1;
    m_x=m(mod(1:N,N)+1,1:N,:);
    m_y=m(1:N,mod(1:N,N)+1,:);
    triple = sum( m(:,:,mod(1:3,3)+1).*m_x(:,:,mod(2:4,3)+1).*m_y(:,:,mod(3:5,3)+1) - m(:,:,mod(1:3,3)+1).*m_y(:,:,mod(2:4,3)+1).*m_x(:,:,mod(3:5,3)+1) ,3);
    denom = 1 + sum(m.*m_x,3) + sum(m.*m_y,3) + sum(m_x.*m_y,3);
    rho = 4*atan2(triple,denom)/(4*pi*dx*dy);
end
