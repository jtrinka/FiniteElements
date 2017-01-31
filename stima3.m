function [M] = stima3(vertices)
% Code from https://www.math.hu-berlin.de/~cc/cc_homepage/download/1999-AJ_CC_FS-50_Lines_of_Matlab.pdf
d=size(vertices,2);
G=[ones(1,d+1);vertices']\[zeros(1,d);eye(d)];
M=det([ones(1,d+1);vertices'])*G*G'/prod(1:d);
end