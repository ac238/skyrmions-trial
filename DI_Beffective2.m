function Beff = DI_Beffective2( NL, wm, Jm, B, U, S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% function Beffective is used to calculate the effective magnetic
%%%% field arisig from all terms in the Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% INPUT: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% NL : lattice size NL*NL
%%%%
%%%% S : spin configuration of the lattice at time t
%%%%%%%%%% constants:
%%%% wm: Coulomb coupling strength
%%%% Jm: stiffness coupling strength
%%%% B: external field
%%%% U: electric potential?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Fx(nx,ny) gives the force in x direction on spin (nx,ny)
%%%% Fy(nx,ny) gives the force in y direction on spin (nx,ny)
%%%% Fz(nx,ny) gives the force in z direction on spin (nx,ny)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% preallocate

  %cutoff = 5;
BCoul = zeros(NL,NL,3);

BZeeStiff =  B +  Jm * (circshift(S,[0 1 0])+circshift(S,[1 0 0]) + ...
		       circshift(S,[-1 0 0]) + circshift(S,[0 -1 0])-4*S);

    function f = fN(S1,S2,S3)
    s1s2 = dot(S1,S2,3);
    s2s3 = dot(S2,S3,3);
    s3s1 = dot(S3,S1,3);
    cs2s3 = cross(S2,S3,3);
    f = (cs2s3.*(1+s1s2 + s2s3 + s3s1)...
        - (S2 + S3).*dot(S1,cs2s3,3))./ ...
        ((1+s1s2).*(1+s2s3).*(1+s3s1));
    end

if wm==0
    BCoul = 0;
    Beff = BZeeStiff;
else
    bigU = zeros(NL,NL);
    rhoP = DI_Pontryagin(NL,S); % Potryagin density matrix
    
    % generate generic matrix for 1/|m-n| with n=(1,1). from here all other
    % matrices can be generated through appropriate cyclic shifts that take
    % into account n and the offset for neighboring site.
    
    dist =-(NL-1):NL-1; 
    distMatBig = 1./sqrt(dist.^2 + dist'.^2); % matrix with n=(1,1)
    %distMatBig(distMatBig < 1/cutoff) = 0;
    distMatBig(NL,NL) = 0; % set self-force to 0.
     
    
    
    % shifted matrices for nearest neighbor spins
    
    SPx = circshift(S,[-1 0 0]); % n+x
    SMx = circshift(S,[1 0 0]); % n-x
    SPy = circshift(S,[0 -1 0]); % n+y
    SMy = circshift(S,[0 1 0]); % n-y
    
    
    
    fPxPy = fN(S,SPx,SPy);
    fPyMx = fN(S,SPy,SMx);
    fMxMy = fN(S,SMx,SMy);
    fMyPx = fN(S,SMy,SPx);
    
    % Ux = circshift(U,[1 0]);
    % Uy = circshift(U,[0 1]);
    % Uxy = circshift(U,[1 1]);
    
    
    % Am = EEnvelope(type, t, T_e );
    
    % using a loop below to control memory requirement (other option
    % is to use a NL^4 dim matrix, which gets too large too quickly)
    % explore using convolution conv2
    
    for nx=1:NL
        for ny=1:NL
            distMat = distMatBig(1+NL-nx:2*NL-nx,1+NL-ny:2*NL-ny);
            bigU(nx,ny) = U(nx,ny) + wm*sum(sum(rhoP .* distMat,2),1);
        end
    end

    bigU = 0.25 * (bigU + circshift(bigU,[-1 0]) + circshift(bigU,[0 -1]) ...
    + circshift(bigU,[-1 -1]));

    BCoul = bigU.*fPxPy + circshift(bigU,[1 0]).*fPyMx + ...
        circshift(bigU,[1 1]).*fMxMy + circshift(bigU,[0 1]).*fMyPx;
    
Beff = -2 * BCoul + BZeeStiff;

end

end
