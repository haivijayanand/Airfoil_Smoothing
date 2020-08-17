function [dist]=distance(yori,ynew,np)
%
% -------------------------------------------------------------------------
% PURPOSE - Compute distance at each x between original and
%           smoothed airfoil. Distance=0.0 means they are identical.
% -------------------------------------------------------------------------
% Coded by: Massimiliano Di Muzio - November 2008
dist=0.0;
for i=1:np
   dist=dist+sqrt((yori(i)-ynew(i))^2.);
end
end
  


