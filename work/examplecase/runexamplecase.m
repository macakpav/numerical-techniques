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
% clc
close all

%% Case parameters
Lx=1; Nx=10;
Ly=1; Ny=10;
uniformU=[0,0]; %initial velocity field

rotation_vector=[1,0]/sqrt(1);

northUx = 5;
dp = -100;
gradP = dp * rotation_vector;

casedef.boundarynames = { 'WEST', 'EAST', 'SOUTH', 'NORTH' };
BCtype = { 'Neumann', 'Neumann', 'Dirichlet', 'Dirichlet'  };
BCval = [ 0 0; 0 0; 0 0; northUx*rotation_vector ];

% convection diffusion tests
% variant = "default";
% variant = "hotPlate";
% switch lower(variant)
%   case "hotplate"
%     %hot plate example
%     T1 = 100;
%     T2 = 10;
%     BCval = [T1, T1, T1, T2];
%     BCtype = {'Dirichlet', 'Dirichlet', 'Dirichlet', 'Dirichlet'};
%   otherwise
%     % left to right diffusiom
%     BCval = [200, 300, 0, 0];
%     BCtype = {'Dirichlet', 'Dirichlet', 'Neumann', 'Neumann'};
% end

%% Create a mesh
seedI = LineSeed.lineSeedOneWayBias([0, 0], [Lx, 0], Nx, 1.00, 'o');
seedJ = LineSeed.lineSeedOneWayBias([0, 0], [0, Ly], Ny, 1.00, 'o');

mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);

%% Create domain from mesh
casedef.dom = newdomain(mesh, 'MyDomain');

%% Set up initial fields
%T = Field(casedef.dom.allCells,0);     % Temperature [K] (scalar); empty field
%reset(T,0);                            % Reset with all zeros

U = Field(casedef.dom.allCells, 1); % Velocity [m/s] (vector);
set(U, [uniformU(1) * ones(1, U.elcountzone); uniformU(2) * ones(1, U.elcountzone)]);

%% Define material properties
casedef.material.nu = 0.15;
casedef.material.rho = 10;
casedef.material.k = 16; % Thermal conductivity [W/(m K)]

casedef.vars.Uinit = U;
casedef.vars.gradP = gradP;

%% Define boundary conditions
for i = 1:length(casedef.boundarynames)
  casedef.BC{i}.zoneID = casedef.boundarynames(i);
  casedef.BC{i}.kind   = BCtype(i);
  casedef.BC{i}.data.bcval = BCval(i,:);
end

%% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.TTol     = 1e-9;
casedef.iteration.dt       = 50;

%% Call solver
result = examplesolver(casedef);

u = Field(casedef.dom.allCells,0);     % Velocity [m/s] (vector);
set(u,result.U.data(1,:));

%% Plot result
figure; hold on; axis off; axis equal; colormap(jet(50));
title("Values of u_x");
scale = 'lin'; lw = 0; quiver=1;
fvmplotfield(u,scale,lw);
% Uoost = restrictto(faceInterpolate(casedef.dom,U),getzone(casedef.dom,'EAST'));
fvmplotvectorfield(faceInterpolate(casedef.dom,U),quiver);
fvmplotmesh(casedef.dom,lw);
% fvmplotcellnumbers(casedef.dom,8);
% fvmplotfacenumbers(casedef.dom,8);
% fvmplotvertexnumbers(casedef.dom,8);
% fvmplotnormals(casedef.dom,lw);

%% Cross section

ymesh = linspace(Ly/Ny/2,Ly-Ly/Ny/2,Ny);
u_sim = u.data(3*Ny+1:4*Ny);

u_anal = analyticalCanalFlow(ymesh,-norm(gradP),Ly,casedef.material.nu*casedef.material.rho, northUx);

rel_error = norm(u_sim-u_anal)/norm(u_anal);
fprintf("Relative error of numerical solution: %10.4e\n", rel_error);

figure()
plot(u_sim,ymesh,'linewidth',2)
set(gca,'Fontsize',12);ylabel('y','Fontsize',14);xlabel('U_x','Fontsize',14)
title('Cross section of u_x')
grid on
hold on
plot(u_anal,ymesh, '*')
legend("numerical solution", "analytical solution",'Location', 'northwest' )

%% Convection diffusion
% switch lower(variant)
%   case "hotplate"
%     t_anal = Field(casedef.dom.allCells,0);
%     set(t_anal, analyticalHotPlateDiffusion(casedef.dom, T1, T2, Ly, Lx));
%     %%
%     figure;
%     hold on;
%     axis off;
%     axis equal;
%     title("Analytical solution");
%     colormap(jet(50));
%     fvmplotfield(t_anal, scale, lw);
%    %% 
%     diff = Field(casedef.dom.allCells,0);
%     set(diff, abs(t_anal.data-result.T.data) ./ abs(t_anal.data))
%     figure;
%     hold on;
%     axis off;
%     axis equal;
%     title("Difference in analytical and numerical solution");
%     colormap(jet(50));
%     fvmplotfield(diff, scale, lw);
%     
% %     rel_error = norm(t_anal.data-result.T.data)./norm(t_anal.data);
%     fprintf("Relative error, max: %10.4e, mean: %10.4e\n", max(diff.data(1:casedef.dom.nPc)), mean(diff.data(1:casedef.dom.nPc)));
%     
%   otherwise
%     xmesh = linspace(Lx/Nx/2, Lx-Lx/Nx/2, Nx);
%     t_sim = result.T.data(1:Ny:casedef.dom.nPc);
%     t_anal = analyticalConvectionDiffusion(xmesh, Lx, uniformU(1), 1, casedef.material.k, BCval(1), BCval(2));
% 
% 
%     rel_error = norm(t_sim-t_anal) / norm(t_anal);
%     fprintf("Relative error: %10.4e\n", rel_error);
% 
%     figure
%     plot(xmesh, t_sim, 'linewidth', 2)
%     set(gca, 'Fontsize', 12);
%     xlabel('x', 'Fontsize', 14);
%     ylabel('\phi', 'Fontsize', 14)
%     title('Cross section of values')
%     grid on
%     hold on
%     plot(xmesh, t_anal, '*')
%     legend('Numerical solution', 'Analytical solution')
% end
