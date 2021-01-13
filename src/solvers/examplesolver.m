%==========================================================================
%
% Example solver using the FVMLab framework 4 students
%
% Purpose: Provides code structure for solving a scalar conservation
%          equation using data structures provided by the framework.
%
% by Frederik Rogiers
%
%==========================================================================
function result = examplesolver(casedef)

dom = casedef.dom;

% Create field objects
U = casedef.vars.Uinit;

% Create all other required data structures
nu = casedef.material.nu;
rho = casedef.material.rho;
dt = casedef.iteration.dt;
L = casedef.vars.L;
gradP = casedef.vars.gradP;

% Create an equation object for holding a scalar conservation equation
eqn_u = ScalarFvEqn2(dom);
eqn_v = ScalarFvEqn2(dom);

iterate = true;
niter = 0;
Uold = U;
while iterate

  niter = niter + 1;

  % Set all terms to zero
  reset(eqn_u);
  reset(eqn_v);

  % Create sparse matrix parts
  adiag_u = zeros(dom.nC, 1);
  bdata_u = zeros(dom.nC, 1);
  anb_internal_u = zeros(2*dom.nIf, 1);
  anb_boundary_u = zeros(2*dom.nBf, 1);

  adiag_v = zeros(dom.nC, 1);
  bdata_v = zeros(dom.nC, 1);
  anb_internal_v = zeros(2*dom.nIf, 1);
  anb_boundary_v = zeros(2*dom.nBf, 1);

  % Calculate U_(n,f) from old iteration
  Uface = faceInterpolate(dom, U);
  Unorm = scalarProduct(Uface, dom.fNormal);

  % Compute coefficients for physical cell eqns and add them to eqn object
  for i = 1:dom.nF % loop over internal faces
    PC = dom.fNbC(2*i-1); % physical cell index
    NBC = dom.fNbC(2*i); % neighbor cell
    xiMag = dom.fXiMag(i); % length of Xi vector (distance between centers)
    fArea = dom.fArea(i); % area of the face
    fUnorm = Unorm(i); % scalar product of faceNormal and Uface(interpolated U on the face)
    cVols = dom.cVol([PC,NBC]);

    anb_conv = fArea * fUnorm / 2; % convective term
    anb_diff = nu * fArea / xiMag; % diffusive term

    if i <= dom.nIf
      adiag_u(PC) = adiag_u(PC) - anb_diff - anb_conv + (cVols(1)/dt)/4; % diagonal element
      adiag_u(NBC) = adiag_u(NBC) - anb_diff + anb_conv + (cVols(2)/dt)/4;
      
      adiag_v(PC) = adiag_v(PC) - anb_diff - anb_conv + (cVols(1)/dt)/4; % diagonal element
      adiag_v(NBC) = adiag_v(NBC) - anb_diff + anb_conv + (cVols(2)/dt)/4;
      
      anb_internal_u(2*i-1) = anb_diff - anb_conv; % offdiagonal element
      anb_internal_u(2*i) = anb_diff + anb_conv;
      
      anb_internal_v(2*i-1) = anb_diff - anb_conv; % offdiagonal element
      anb_internal_v(2*i) = anb_diff + anb_conv;
      
      
      bdata_u(PC) = bdata_u(PC) + (Uold.data(1,PC) * cVols(1) / dt - gradP(1)/rho * cVols(1)) /4;
      bdata_v(PC) = bdata_v(PC) + (Uold.data(2,PC) * cVols(1) / dt - gradP(2)/rho * cVols(1)) /4;
      bdata_u(NBC) = bdata_u(NBC) + (Uold.data(1,NBC) * cVols(2) / dt - gradP(1)/rho * cVols(2)) /4;
      bdata_v(NBC) = bdata_v(NBC) + (Uold.data(2,NBC) * cVols(2) / dt - gradP(2)/rho * cVols(2)) /4;
    else
      adiag_u(PC) = adiag_u(PC) - anb_diff - anb_conv + (cVols(1)/dt)/4;
      anb_boundary_u(2*(i - dom.nIf)-1) = anb_diff - anb_conv;
      adiag_v(PC) = adiag_v(PC) - anb_diff - anb_conv + (cVols(1)/dt)/4;
      anb_boundary_v(2*(i - dom.nIf)-1) = anb_diff - anb_conv;
      bdata_u(PC) = bdata_u(PC) + (Uold.data(1,PC) * cVols(1) / dt - gradP(1)/rho * cVols(1)) /4;
      bdata_v(PC) = bdata_v(PC) + (Uold.data(2,PC) * cVols(1) / dt - gradP(2)/rho * cVols(1)) /4;
    end
  end

  % Compute coefficients for ghost cell eqns and add them to eqn object
  for j = 1:length(casedef.BC) % loop over all boundary zones
    Bzone = dom.getzone(casedef.BC{j}.zoneID);
    Brange = Bzone.range(1):Bzone.range(2); %range of boundary face indicies
    Bkind = casedef.BC{j}.kind;
    Bcoef = casedef.BC{j}.data.bcval;

    for i = Brange
      GC = dom.fNbC(2*i); % ghost cell index
      bdata_u(GC) = Bcoef(1);
      bdata_v(GC) = Bcoef(2);

      if (Bkind == "Dirichlet")
        adiag_u(GC) = 0.5;
        adiag_v(GC) = 0.5;
        anb_boundary_u(2*(i - dom.nIf)) = 0.5;
        anb_boundary_v(2*(i - dom.nIf)) = 0.5;
      elseif (Bkind == "Neumann")
        adiag_u(GC) = 1 / xiMag;
        adiag_v(GC) = 1 / xiMag;
        anb_boundary_u(2*(i - dom.nIf)) = -1 / xiMag;
        anb_boundary_v(2*(i - dom.nIf)) = -1 / xiMag;
      else
        printf("Unknown kind of boundary condition for zone %s", Bzone.id)
      end
    end
  end
  
  % False time stepping and volume source terms
%   for i = 1:dom.nC
%     PC = i;
%     cVol = dom.cVol(PC);
%     
%     adiag_u(PC) = adiag_u(PC) + cVol/dt; % diagonal element
%     adiag_v(PC) = adiag_v(PC) + cVol/dt;
%     bdata_u(PC) = bdata_u(PC) + Uold.data(1,PC) * cVol / dt - dp/rho * cVol;
%     bdata_v(PC) = bdata_v(PC) + Uold.data(2,PC) * cVol / dt;%- dp/rho * cVol;
%   end

  eqn_u.adata = [adiag_u; anb_internal_u; anb_boundary_u];
  eqn_u.bdata = bdata_u;
  
  eqn_v.adata = [adiag_v; anb_internal_v; anb_boundary_v];
  eqn_v.bdata = bdata_v;

  % Create a matlab sparse linear system from the eqn object
  [Au, bu] = to_msparse(eqn_u);
  [Av, bv] = to_msparse(eqn_v);
  %spy(Au);
  u = U.data(1,:)';
  v = U.data(2,:)';

  % Check tolerance and iteration count
  uRes = bu - Au * u;
  vRes = bv - Av * v;
  UResnorm = max( [norm(uRes) norm(vRes)] );
  if UResnorm < casedef.iteration.TTol
    Uconverged = true;
    iterate = false;
  elseif niter > casedef.iteration.maxniter
    Uconverged = false;
    iterate = false;
    %    elseif checkstoprequest(stopmon)
    %       Tconverged = false;
    %       iterate = false;
  else
    u = Au \ bu;
    v = Av \ bv;
    % Alternatives: gmres, bicgstabb, ...
    set(U, [u';v']); % Put algebraic solution in the Field
  end


end % iterate

result.endtime = now; % call datestr(now) for displaying this time
result.Uconverged = Uconverged;
result.niter = niter;
result.UResnorm = UResnorm;
result.uRes = Field(dom.allCells, 0); set(result.uRes, uRes');
result.vRes = Field(dom.allCells, 0); set(result.vRes, vRes');
result.U = U;


end
