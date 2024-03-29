# skyrmions-trial
A trial project in which I simulate the dynamics of skyrmions using MATLAB.



% create bounds of graph
Sets up the set of coordinate points upon which the magnetization vectors will be measured.
N is the number of points in each direction, and the following four values are the bounds of the graph in units of position.
[yy,xx] are two NxN arrays which give the x- and y-positions at each index, respectively; they are swapped to correspond with the positions of the indices.

% initialize skyrmion according to QHMF 35
Creates a single skyrmion with charge -n at the point (Re(z_0),Im(z_0)).
The magnetization field is defined in this block as an NxNx3 matrix which give the magnetization components m(:,:,1)_init, m(:,:,2)_init, m(:,:,3)_init, at each index.
First, the x- and y-positions are redefined as the real and imaginary axis of z. Then, omega is found from z, and m_init is found from omega.

% dynamics
Does numerical calculations to the components of m_init and saves them in m. Plots m(:,:,1) and m(:,:,2), then saves them all in m_init to repeat.
Each term of the equations of motion are added separately so that they may be turned off by setting the values "zeeman" and "landau" to zero.

% Zeeman term
For the Zeeman term, B was chosen to be uniform in the z-direction, causing the magnetization field vectors to rotate in the xy-plane.
The strength coefficient is gmuB=g*mu*B, where g is the g-factor, mu is the magneton, and B is the magnetic field strength.
RK4 is used, which preserves the norm of m as 1 and conserves topological charge.

% Electric field term
For the electric field term, E was chosen to be uniform in the x-direction, causing the skyrmion to move to the right. RK4 is used.
The strength coefficient is q=e*E/(4*pi*s*n), where e is the charge, E is the applied electric field strength, s is the spin Casimir, and n is the particle number density.

% Landau-Lifshitz term
This term causes the skyrmion to spread out.
The strength coefficient is kappa=(rho_s)/(s*n), where rho_s is the spin stiffness, s is the spin Casimir, and n is the particle number density.
RK4 is used, but it does not conserve |m|=1; this is not physical and needs to be addressed.

% Pontryagin density
I used the solid-angle method with the intent to exactly calculate the numerical Pontryagin density.
The numerical Pontryagin density is equal to the solid angle between m and its nearest neighbors divided by 4*pi.
The solid-angle formula is complicated, so I calculated it in pieces. For any parts that used the magnitude, I simply substituted 1.

% Coulomb term
I added the Coulomb term with RK4.
The unchanging parts of the calculation within the for loop were condensed into two larger non-dynamic matrices, dist_x and dist_y, for efficiency.

% Topological charge
I numerically integrated rho over all space to find Q_top, the topological charge.
I used periodic indexing so that Q_top should always be exactly an integer. If this isn't the case when a dynamical term is turned on, then there is an issue.
The sign of Q_top is always positive; this needs to be fixed.

% plot
For each index, quiver creates a field of vectors whose components are proportional to m(:,:,1) and m(:,:,1) at positions xx and yy.
The contour plot is of rho, the Pontryagin density, but it can be changed to check other unseen values such as m(:,:,3) or norm_m (to ensure it is preserved).
A line of code is present which can be uncommented to replace the graphs with a single "slice" of the skyrmion along the x-axis.

% check conserved quantities
Calculates conserved quantities such as S_z and various types of energy.
Only the Zeeman, A, and stiffness terms have been implemented; Coulomb and electric field are coming soon.
Plots them over time after all time is elapsed.
E_El doesn't seem to be conserved, but when it oscillates, it does seem to depend on the skyrmion's position, much like a charge in a potential.
Q_top is not conserved by the electric or Coulomb terms.
E_C still needs to be implemented.

% draw saved data
Records each time-frame of m and rho in a large matrix so that the dynamics can be replayed in constant time.
Copy-paste this into console to replay it without generating the data again.



