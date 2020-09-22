function [dd mdd] = focusProfile2(fcen,coords,cfl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 22, 2019
% focal profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dd=zeros(size(coords,1),1);
for i=1:length(fcen)
  dd = dd + (coords(:,i)-fcen(i)).^2;
end
dd=sqrt(dd);
mdd=dd/cfl;
dd=mdd;
dd=dd-min(dd);

