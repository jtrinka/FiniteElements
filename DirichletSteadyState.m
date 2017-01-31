function [Dnodes] = DirichletSteadyState(nodes,xmin,xmax,ymin,ymax,xy)
%Code from https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf
Dnodes=[];
for j=1:max(nodes(:))
 
        
        %just x values
   if xy(j,1)==xmin
       
       Dnodes=[Dnodes;j];
   elseif xy(j,1)==xmax
       Dnodes=[Dnodes;j];
       
   elseif xy(j,2)==ymin
       Dnodes=[Dnodes;j];
       
   elseif xy(j,2)==ymax
        Dnodes=[Dnodes;j];
   end
        
        
        Dnodes=unique(Dnodes);
    
    
end