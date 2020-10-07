%==========================================================================
%
% Example case using the FVMLab framework 4 students
%
% Purpose: Provides an example for setting up a case and calling a solver.
%          This involves creating a mesh, defining materials, defining
%          boundary conditions, defining iteration parameters, and finally
%          calling the solver.
%
% by Frederik Rogiers
%
%==========================================================================

% TIP: use "clear variables" instead of "clear all" to clear variables
%      use "clear classes" when the interface of a class has changes
%      use "close all" to close figures
%      use "clc" to clear the command window
% 
% TIP: pressing CTRL+D while the cursor is on a function opens that function
%      in the m-editor. This is the most convenient way of browsing through
%      your source code.

clear variables
clc
close all

% Create a mesh
% seedI = LineSeed.lineSeedOneWayBias([0 0],[0.75 0],10,1.00,'o');
% seedJ = LineSeed.lineSeedOneWayBias([0 0],[0.5 0.5],5,1.00,'o');
Lx=1; Nx=10;
Ly=1; Ny=10;
seedI = LineSeed.lineSeedOneWayBias([0 0],[Lx 0],Nx,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 Ly],Ny,1.00,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');


% Set up initial fields
T = Field(casedef.dom.allCells,0);     % Temperature [K] (scalar); empty field
% reset(T,0);                          % Reset with all zeros
randomdata = rand(T.elsize,T.elcountzone)-0.5;
set(T,randomdata);                     % Set with random numbers

U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[10*ones(1,U.elcountzone);zeros(1,U.elcountzone)]);
%reset(U,[1;0.2]);


% Define material properties
casedef.material.k = 16;  % Thermal conductivity [W/(m K)]
casedef.vars.U = U;


% Define boundary conditions
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 10;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;


% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.TTol     = 1e-6;


% Call solver
result = examplesolver(casedef);

% Plot result
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 0; quiver=1;
fvmplotfield(result.T,scale,lw);
Uoost = restrictto(faceInterpolate(casedef.dom,U),getzone(casedef.dom,'OOSTRAND'));
fvmplotvectorfield(Uoost,quiver);
% fvmplotmesh(casedef.dom,lw);
% fvmplotcellnumbers(casedef.dom,8);
% fvmplotfacenumbers(casedef.dom,8);
% fvmplotvertexnumbers(casedef.dom,8);
% fvmplotnormals(casedef.dom,lw);


% Cross section

t = result.T.data;

xmesh = linspace(0,Lx,Nx);
t_plot = t(1:Ny:casedef.dom.nPc);

figure
plot(xmesh,t_plot,'linewidth',2)
set(gca,'Fontsize',12);xlabel('x','Fontsize',14);ylabel('\phi','Fontsize',14)
title('Cross section of species value')
grid on
