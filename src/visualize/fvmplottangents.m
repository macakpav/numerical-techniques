function fvmplottangents(dom,lw)

fCoord = dom.fCoord;
fTangent = dom.fTangent;

  gh = quiver(gca,fCoord(1,:),fCoord(2,:),fTangent(1,:),fTangent(2,:)); 
  set(gh,'Color','k',...
          'LineWidth',lw,...
          'MaxHeadSize',0.08,...
          'AutoScaleFactor',0.75);

end