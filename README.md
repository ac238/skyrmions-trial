# skyrmions-trial
A trial project in which I simulate the dynamics of skyrmions using MATLAB.



% create bounds of graph
Sets up the set of coordinate points upon which the magnetization vectors will be measured.
N is the number of points in each direction, and the following four values are the bounds of the graph in units of position.
[xx,yy] are two NxN arrays which give the x- and y-positions at each index, respectively.

% initialize skyrmion according to QHMF 35
Creates a single skyrmion with charge -n at the point (Re(z_0),Im(z_0)).
The magnetization field is defined in this block as a pair of NxN matrices which give the magnetization components m1_init, m2_init, m3_init at each index.
First, the x- and y-positions are redefined as the real and imaginary axis of z. Then, omega is found from z, and m1_init and m2_init are found from omega.

% dynamics
Does numerical calculations to m1_init, m2_init, m3_init and saves them in m1, m2, m3. Plots m1 and m2, then saves them all in m1_init, m2_init, m3_init to repeat.
Each term of the equations of motion are added separately so that they may be turned off by setting the values "zeeman" and "landau" to zero.

% Zeeman term
For the Zeeman term, B was chosen to be uniform in the z-direction, causing the magnetization field vectors to rotate in the xy-plane.
The strength coefficient is gmuB=g*mu*B, where g is the g-factor, mu is the magneton, and B is the magnetic field strength.
RK4 is used, which preserves the norm of m as 1 and conserves topological charge.

% Electric field term
For the electric field term, E was chosen to be uniform in the x-direction, causing the skyrmion to move to the right.
RK4 is used, but the skyrmion's shape is not preserved. This needs to be fixed.
El = e*E/(4*pi*s*n)

% Landau-Lifshitz term
This term causes the skyrmion to spread out.
The strength coefficient is kappa=(rho_s)/(s*n), where rho_s is the spin stiffness, s is the spin Casimir, and n is the particle number density.
RK4 is used, but it does not conserve abs(m)=1; this is not physical and needs to be addressed.

% Pontryagin density
I used the solid-angle method with the intent to exactly calculate the numerical Pontryagin density.
The numerical Pontryagin density is equal to the solid angle between m and its nearest neighbors divided by 4*pi.
The solid-angle formula is complicated, so I calculated it in pieces. For any parts that used the magnitude, I simply substituted 1.

% Topological charge
I numerically integrated rho over all space to find Q_top, the topological charge.
I used periodic indexing so that Q_top should always be exactly an integer. If this isn't the case when a dynamical term is turned on, then there is an issue.
The sign of Q_top is always positive; this needs to be fixed.

% plot
For each index, quiver creates a field of vectors whose components are proportional to m1 and m2 at positions xx and yy.
The contour plot is of rho, the Pontryagin density, but it can be changed to check other unseen values such as m3 or norm_m (to ensure it is preserved).
A line of code is present which can be uncommented to replace the graphs with a single "slice" of the skyrmion along the x-axis.

% check conserved quantities
Calculates conserved quantities such as S_z and various types of energy.
Only the Zeeman, A, and stiffness terms have been implemented; Coulomb and electric field are coming soon.
Plots them over time after all time is elapsed.




