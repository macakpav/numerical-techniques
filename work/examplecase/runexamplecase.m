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

%% Case parameters
Lx=1; Nx=10;
Ly=1; Ny=10;
uniformU=[0,0]; %initial velocity field

casedef.boundarynames = { 'WEST', 'EAST', 'SOUTH', 'NORTH' };
BCtype = { 'Neumann', 'Neumann', 'Dirichlet', 'Dirichlet'  };
BCval = [ 0 0; 0 0; 0 0; 0 0 ];

%% Create a mesh
seedI = LineSeed.lineSeedOneWayBias([0 0],[Lx 0],Nx,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 Ly],Ny,1.00,'o');

mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);

%% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');


%% Set up initial fields
%T = Field(casedef.dom.allCells,0);     % Temperature [K] (scalar); empty field
%reset(T,0);                            % Reset with all zeros

U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[uniformU(1)*ones(1,U.elcountzone);uniformU(2)*ones(1,U.elcountzone)]);


%% Define material properties
casedef.material.nu = 0.1;
casedef.material.rho = 1;

casedef.vars.Uinit = U;
casedef.vars.pIn = 100;
casedef.vars.pOut = 0;
casedef.vars.L = Lx;

%% Define boundary conditions
for i = 1:length(casedef.boundarynames)
  casedef.BC{i}.zoneID = casedef.boundarynames(i);
  casedef.BC{i}.kind   = BCtype(i);
  casedef.BC{i}.data.bcval = BCval(i,:);    
end

%% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.TTol     = 1e-6;
casedef.iteration.dt       = 50;

%% Call solver
result = examplesolver(casedef);

%% Plot result
% figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 0; quiver=1;
% fvmplotfield(result.T,scale,lw);
% Uoost = restrictto(faceInterpolate(casedef.dom,U),getzone(casedef.dom,'EAST'));
fvmplotvectorfield(U,quiver);
% fvmplotmesh(casedef.dom,lw);
% fvmplotcellnumbers(casedef.dom,8);
% fvmplotfacenumbers(casedef.dom,8);
% fvmplotvertexnumbers(casedef.dom,8);
% fvmplotnormals(casedef.dom,lw);


%% Cross section

% xmesh = linspace(0,Lx,Nx);
% t_plot = result.T.data(1:Ny:casedef.dom.nPc);
% 
% figure
% plot(xmesh,t_plot,'linewidth',2)
% set(gca,'Fontsize',12);xlabel('x','Fontsize',14);ylabel('\phi','Fontsize',14)
% title('Cross section of species value')
% grid on
% hold on
% plot(xmesh, analyticalConvectionDiffusion(xmesh, Lx, uniformU(1), 1, casedef.material.k, BCval(1), BCval(2)),'*')

