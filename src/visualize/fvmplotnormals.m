function fvmplotnormals(dom,lw)

fCoord = dom.fCoord;
fNormal = dom.fNormal;

  gh = quiver(gca,fCoord(1,:),fCoord(2,:),fNormal(1,:),fNormal(2,:)); 
  set(gh,'Color','k',...
          'LineWidth',lw,...
          'MaxHeadSize',0.08,...
          'AutoScaleFactor',0.75);

end