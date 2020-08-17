function [first_derivative]=dS(xold,K,r,n,N1,N2)
%
% -------------------------------------------------------------------------
% PURPOSE - Calculate 1st derivative of z(xold) function
% -------------------------------------------------------------------------
% Coded by: Massimiliano Di Muzio - November 2008
% -------------------------------------------------------------------------

first_derivative=(1-xold)^(N2 + n - r)* K*(-xold^(N1 - 1 + r)*N1+xold^(N1 + r)*N1+xold^(N1 + r)*N2-xold^(N1 - 1 + r)*r+xold^(N1 + r)*n)/(-1+xold);