function DirichletBoundaryValue = u_d(x,t)
%Code from https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf
% if t<.5
% DirichletBoundaryValue = (.0001)*ones(size(x,1),1);
% 
% elseif t>=.5 && t<1
%     
%     DirichletBoundaryValue = (.0001)*ones(size(x,1),1);
%     
% elseif t>=1 && t<2
%     
%     DirichletBoundaryValue = (.0001)*ones(size(x,1),1);
%     
% else
%    
% DirichletBoundaryValue = (.0001+.0000001*t)*ones(size(x,1),1);
% end
DirichletBoundaryValue = 0.00085*ones(size(x,1),1);
end