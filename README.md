# skyrmions-trial
A trial project in which I simulate the dynamics of skyrmions using MATLAB.



% create bounds of graph
Sets up the set of coordinate points upon which the magnetization vectors will be measured.
N is the number of points in each direction, and the following four values are the bounds of the graph in units of position.
[xx,yy] are two NxN arrays which give the x- and y-positions at each index, respectively.

% initialize skyrmion according to QHMF 35
Creates a single skyrmion with charge -n at the point (Re(z_0),Im(z_0)).
The magnetization field is defined in this block as a pair of NxN matrices which give the physical components m1_init and m2_init at each index.
First, the x- and y-positions are redefined as the real and imaginary axis of z. Then, omega is found from z, and m1_init and m2_init are found from omega.

% plot
For each index, quiver creates a field of vectors whose components are proportional to m1 and m2 at positions xx and yy.

% dynamics
Does numerical calculations to m1_init and m2_init and saves them in m1 and m2. Plots m1 and m2, then saves them in m1_init and m2_init to repeat.
Each term of the equations of motion are added separately so that they may be turned off by setting the values "zeeman" and "landau" to zero.

% Landau-Lifshitz term (broken)
The Landau-Lifshitz term is broken. It does not preserve the norm of m as 1. Periodic boundary conditions are implemented using the mod function.

% Zeeman term
For the Zeeman term, B was chosen to be uniform in the z-direction, causing the skyrmion to rotate.
RK4 is used, which preserves the norm of m as 1 and conserves topological charge.




