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
T = Field(dom.allCells,0);      % Temperature [K] (scalar); empty field
reset(T,0);                     % Reset with all zeros

% Create all other required data structures
k=casedef.material.k;
U=casedef.vars.U;
Uface=faceInterpolate(dom,U);
Unorm=scalarProduct(Uface,dom.fNormal);

% Create an equation object for holding a scalar conservation equation
eqn = ScalarFvEqn2(dom);

iterate = true;
niter = 0;
while iterate   
   
   niter = niter+1;
   
   % Set all terms to zero
   reset(eqn); 
   
   adiag = zeros(dom.nC,1);
   bdata = zeros(dom.nC,1);
   anb_internal = zeros(2*dom.nIf,1);
   anb_boundary = zeros(2*dom.nBf,1);
   
% Compute coefficients for physical cell eqns and add them to eqn object
   for i = 1:dom.nIf % loop over internal faces
      PC=dom.fNbC(2*i-1); %physical cell index
      NBC=dom.fNbC(2*i); %neighbor cell
      xiMag=dom.fXiMag(i); % length of Xi vector (distance between centers)
      fArea=dom.fArea(i); % area of the face
      fUnorm=Unorm(i); % scalar product of faceNormal and Uface(interpolated U on the face)
      
      anb_conv=fArea*fUnorm; %convective term
      anb_diff=k*fArea/xiMag; %diffusive term
      
      adiag(PC) = adiag(PC)-anb_diff-anb_conv; %diagonal elements
      adiag(NBC) = adiag(NBC)-anb_diff+anb_conv;
      anb_internal(2*i-1) = anb_diff-anb_conv; %offdiagonal elements
      anb_internal(2*i) = anb_diff+anb_conv;
   end

% Compute coefficients for ghost cell eqns and add them to eqn object
    for i = 1:length(casedef.BC) %loop over all boundary zones
      Bzone = dom.getzone(casedef.BC{i}.zoneID);
      Brange = Bzone.range(1):Bzone.range(2);
      Bkind = casedef.BC{i}.kind;
      Bcoef = casedef.BC{i}.data.bcval;

      for i = Brange
        Cind=dom.fNbC(2*(i-1)+1:2*i); % indicies of neighbour cells
        PC=Cind(1);
        GC=Cind(2);
        fUnorm=Unorm(i);

        xiMag=dom.fXiMag(i); % length of Xi vector (distance between centers)
        fArea=dom.fArea(i); % area of the face

        anb_conv=fArea*fUnorm;
        anb_diff=k*fArea/xiMag;
        adiag(PC) = adiag(PC)-anb_diff-anb_conv;
        anb_boundary(2*(i-dom.nIf)-1) = anb_diff-anb_conv;
        bdata(GC) = Bcoef;

        if (Bkind == "Dirichlet")
          adiag(GC) = 0.5;
          anb_boundary(2*(i-dom.nIf)) = 0.5;
        elseif (Bkind == "Neumann")
          adiag(GC) = 1/xiMag;
          anb_boundary(2*(i-dom.nIf)) = -1/xiMag;
        else
          printf("Unknown kind of boundary condition for zone %s", Bzone.id)
        end
     end
     
     
   end

   eqn.adata = [adiag; anb_internal; anb_boundary]; % just a meaningless example adata
   eqn.bdata = bdata;   % just a meaningless example bdata
      
   % Create a matlab sparse linear system from the eqn object
   [A,b] = to_msparse(eqn);
   %spy(A);
   x = get(T);
   x = x';
   
   % Check tolerance and iteration count
   TRes = b-A*x;
   TResnorm = norm(TRes);         
   if TResnorm < casedef.iteration.TTol
      Tconverged = true;
      iterate = false;
   elseif niter > casedef.iteration.maxniter
      Tconverged = false;
      iterate = false;
%    elseif checkstoprequest(stopmon)
%       Tconverged = false;
%       iterate = false;
   else
      x = A\b; % Direct sparse solver.
               % Alternatives: gmres, bicgstabb, ...
      set(T,x'); % Put algebraic solution in the Field
   end
      
   
end % iterate

result.endtime = now; % call datestr(now) for displaying this time 
result.Tconverged = Tconverged;
result.niter = niter;
result.TResnorm = TResnorm;
result.TRes = Field(dom.allCells,0);
   set(result.TRes,TRes');
result.T = T;


end



