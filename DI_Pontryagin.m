function [rhoP,Q,Qdiag,Qquad,Qtri] = DI_Pontryagin(dummy,S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% function pontryagin is used to calculate the Pontryagin density
%%%% at each point from a given spin configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% INPUT: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% dummy: this used to be the size, which is now computed from S to accommodate boundary conditions. 
%%%% S : spin configuration of the lattice at time t 
%%%%%%%%%% constants:
%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% rhoP is an array for the Pontryagin density at each point
%%%% Q is the charge = sum over the density at all points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% RELEVANT FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% NearestNeighbour : find each site's nearest neighbour points
%%%% the neighbour points are stored in the order : n+y, n-y, n+x, n-x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOTE on indices:
%%% in this code, and other related code in this project, the first
%%% index (row index) is x, and the second (column) index is y. 


  NL = length(S);
%%%% preallocate

  rhoP = zeros(NL, NL);

% shifted matrices: +x, +y, +x+y in that order

SPx = circshift(S,[-1 0 0]);
SPy = circshift(S,[0 -1 0]);
SPxy = circshift(S,[-1 -1 0]);

omega1 = 2 * atan2(dot(S,cross(SPx,SPxy,3),3),...
    (1 + dot(S,SPx,3) + dot(SPx,SPxy,3) + dot(SPxy,S,3)));

omega2 = 2 * atan2(dot(S,cross(SPxy,SPy,3),3),...
    (1 + dot(S,SPxy,3) + dot(SPxy,SPy,3) + dot(SPy,S,3)));

rhoP = 2*(omega1 + omega2);

% adjust for the derivative taking only forward difference
% average over 4 plaquettes - check with and wihout this. 
rhoP = 0.25 * (rhoP + circshift(rhoP,[1 0]) + circshift(rhoP,[0 1]) ...
    + circshift(rhoP,[1 1]));

Q = sum(sum(rhoP))/(8*pi);

Qquad = 1/(8*pi)*[sum(sum(rhoP(1:NL/2,1:NL/2))), sum(sum(rhoP(1:NL/2,NL/2+1:NL))),sum(sum(rhoP(NL/2+1:NL,1:NL/2))),sum(sum(rhoP(NL/2+1:NL,NL/2+1:NL)))];

Qtri = 1/(8*pi)*[sum(sum(triu(rhoP,1))),sum(sum(tril(rhoP,-1))) ];
Qdiag = 1/(8*pi)*sum(sum(diag(rhoP)));
end
