function fvmplotcellnumbers(dom,fs)

cCoord = dom.cCoord;

for jC = 1:dom.nC
      
   text(cCoord(1,jC),cCoord(2,jC),num2str(jC),...
        'HorizontalAlignment','center',...
        'HitTest','off',...
        'FontAngle','normal',...
        'FontWeight','bold',...
        'FontUnits','points',...
        'FontSize',fs);
           
end


end