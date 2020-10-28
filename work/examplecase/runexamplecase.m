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
Lx = 1;
Nx = 10;
Ly = 1;
Ny = 10;
MMS_func = 'MMS'; % MMS linear onlyX default
uniformU = [0, 0]; %initial velocity field

casedef.boundarynames = {'WEST', 'EAST', 'SOUTH', 'NORTH'};
% BCtype = {'Neumann', 'Neumann', 'Neumann', 'Neumann'};
BCtype = {'Dirichlet', 'Dirichlet', 'Dirichlet', 'Dirichlet'};
% BCtype = {'Neumann', 'Neumann', 'Dirichlet', 'Dirichlet'};
BCval = [0, 0; 0, 0; 0, 0; 0, 0];

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
casedef.material.nu = 1;
casedef.material.rho = 1;

casedef.vars.Uinit = U;
casedef.vars.pIn = 100;
casedef.vars.pOut = 0;
casedef.vars.L = Lx;

%% Define boundary conditions
for i = 1:length(casedef.boundarynames)
  casedef.BC{i}.zoneID = casedef.boundarynames(i);
  casedef.BC{i}.kind = BCtype(i);
  casedef.BC{i}.data.bcval = BCval(i, :);
end

%% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.TTol = 1e-6;
casedef.iteration.dt = 50;

%% MMS function pointers

switch MMS_func
  case 'MMS'
    % solution
    sol_x = @(x, y) x.^2;
    sol_y = @(x, y) -y.^2;
    % source
    src_x = @(x, y) 4 .* x.^3 - 2 .* x.^2 .* y - 2; %* casedef.material.nu;
    src_y = @(x, y) -2 .* x .* y.^2 + 4 .* y.^3 + 2; %* casedef.material.nu;
  case 'linear'
    % solution
    sol_x = @(x, y) x;
    sol_y = @(x, y) y;
    % source
    src_x = @(x, y) 2*x+x;
    src_y = @(x, y) 2*y+y;
  case 'onlyX'
    % solution
    sol_x = @(x, y) x.^2;
    sol_y = @(x, y) 0*y;
    % source
    src_x = @(x, y) 4*x.^3-2; %*casedef.material.nu;
    src_y = @(x, y) 0*y;
  case 'default'
    % solution
    sol_x = @(x, y) 0*x;
    sol_y = @(x, y) 0*y;
    % source
    src_x = @(x, y) -ones(1,length(x))*(casedef.vars.pOut-casedef.vars.pIn)/casedef.material.rho/casedef.vars.L;
    src_y = @(x, y) 0*y;
  otherwise
    
end

%set MMS funcitons
Src = Field(casedef.dom.allCells,1);
set(Src,[src_x(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:));
         src_y(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:))])
casedef.S.Src = Src;
casedef.S.src_x = src_x;
casedef.S.src_y = src_y;
Sol = Field(casedef.dom.allCells,1);
set(Sol, [sol_x(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:));
          sol_y(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:))])
casedef.S.Sol = Sol;
casedef.S.sol_x = sol_x;
casedef.S.sol_y = sol_y;

%% Call solver
result = examplesolver(casedef);

%% 
u = Field(casedef.dom.allCells, 0); % Velocity [m/s] (vector);
set(u, result.U.data(1, :));
v = Field(casedef.dom.allCells, 0); % Velocity [m/s] (vector);
set(v, result.U.data(2, :));
uS = Field(casedef.dom.allCells, 0); % Velocity [m/s] (vector);
set(uS, casedef.S.Sol.data(1, :));
vS = Field(casedef.dom.allCells, 0); % Velocity [m/s] (vector);
set(vS, casedef.S.Sol.data(2, :));

%% Plot result
figure;
hold on;
axis off;
axis equal;
colormap(jet(50));
scale = 'lin';
lw = 0;
quiver = 1;

fvmplotfield(uS, scale, lw);
fvmplotvectorfield(faceInterpolate(casedef.dom,casedef.S.Sol),quiver)

figure;
hold on;
axis off;
axis equal;
colormap(jet(50));
scale = 'lin';
lw = 0;
quiver = 1;
fvmplotfield(u, scale, lw);
fvmplotvectorfield(faceInterpolate(casedef.dom,U),quiver);

% Uoost = restrictto(faceInterpolate(casedef.dom,U),getzone(casedef.dom,'EAST'));
% fvmplotmesh(casedef.dom, lw);
% fvmplotcellnumbers(casedef.dom,8);
% fvmplotfacenumbers(casedef.dom,8);
% fvmplotvertexnumbers(casedef.dom,8);
% fvmplotnormals(casedef.dom,lw);

%% Cross section
figure()
dif = abs(result.U.data(1,:) - casedef.S.Sol.data(1,:));
%dif = sum(dif,1);
diffield = Field(casedef.dom.allCells, 0); % Velocity [m/s] (vector);
set(diffield, dif);
fvmplotfield(diffield, scale, lw);
title('Difference in X-component')

figure()
dif = abs(result.U.data(2,:) - casedef.S.Sol.data(2,:));
%dif = sum(dif,1);
diffield = Field(casedef.dom.allCells, 0); % Velocity [m/s] (vector);
set(diffield, dif);
fvmplotfield(diffield, scale, lw);
title('Difference in Y-component')


% 
% ymesh = linspace(0, Ly, Ny);
% u_plot = u.data(3*Ny+1:4*Ny);
% 
% figure()
% plot(u_plot, ymesh, 'linewidth', 2)
% set(gca, 'Fontsize', 12);
% ylabel('y', 'Fontsize', 14);
% xlabel('U_x', 'Fontsize', 14)
% title('Cross section of U_x')
% grid on
% hold on
% dp = (casedef.vars.pOut - casedef.vars.pIn) / casedef.vars.L;
% plot(analyticalCanalFlow(ymesh, dp, Ly, casedef.material.nu * casedef.material.rho), ymesh, '*') %analyticalConvectionDiffusion(xmesh, Lx, uniformU(1), 1, casedef.material.k, BCval(1), BCval(2)),'*')


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
