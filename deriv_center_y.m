
function E_field = deriv_center_y(m)
    N=length(m(:,1,1));
    centered_m = m(:,mod(0:N-3,N)+1,:)-m(:,mod(2:N-1,N)+1,:);
    E_field = zeros(N,N,3);
    E_field(:,2:N-1,:)=centered_m;
end
