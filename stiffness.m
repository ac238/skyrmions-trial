
function stiff = stiffness(m)
    N=length(m(:,1,1));
    stiffed_m = m(1:N,1:N,mod(1:3,3)+1).*(m(mod(-1:N-2,N)+1,1:N,mod(2:4,3)+1)+m(mod(1:N,N)+1,1:N,mod(2:4,3)+1)+m(1:N,mod(-1:N-2,N)+1,mod(2:4,3)+1)+m(1:N,mod(1:N,N)+1,mod(2:4,3)+1)-4*m(1:N,1:N,mod(2:4,3)+1)) - m(1:N,1:N,mod(2:4,3)+1).*(m(mod(-1:N-2,N)+1,1:N,mod(1:3,3)+1)+m(mod(1:N,N)+1,1:N,mod(1:3,3)+1)+m(1:N,mod(-1:N-2,N)+1,mod(1:3,3)+1)+m(1:N,mod(1:N,N)+1,mod(1:3,3)+1)-4*m(1:N,1:N,mod(1:3,3)+1));
    centered_m = m(2:N-1,2:N-1,mod(1:3,3)+1).*(m(mod(0:N-3,N)+1,2:N-1,mod(2:4,3)+1)+m(mod(2:N-1,N)+1,2:N-1,mod(2:4,3)+1)+m(2:N-1,mod(0:N-3,N)+1,mod(2:4,3)+1)+m(2:N-1,mod(2:N-1,N)+1,mod(2:4,3)+1)-4*m(2:N-1,2:N-1,mod(2:4,3)+1)) - m(2:N-1,2:N-1,mod(2:4,3)+1).*(m(mod(0:N-3,N)+1,2:N-1,mod(1:3,3)+1)+m(mod(2:N-1,N)+1,2:N-1,mod(1:3,3)+1)+m(2:N-1,mod(0:N-3,N)+1,mod(1:3,3)+1)+m(2:N-1,mod(2:N-1,N)+1,mod(1:3,3)+1)-4*m(2:N-1,2:N-1,mod(1:3,3)+1));
    stiff = zeros(N,N,3);
    stiff(2:N-1,2:N-1,:)=centered_m;
end
