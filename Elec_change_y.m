
function E_field = Elec_change_y(m)
    N=length(m(:,1,1));

    centered_m = (m(:,mod(2:N-1,N)+1,:)-m(:,mod(0:N-3,N)+1,:))/2;
    far_edge_m = (m(:,N,:)-m(:,N-1,:));
    near_edge_m = (m(:,2,:)-m(:,1,:));

    E_field = zeros(N,N,3);
    E_field(:,2:N-1,:)=centered_m;
    E_field(:,N,:)=far_edge_m;
    E_field(:,1,:)=near_edge_m;
end
