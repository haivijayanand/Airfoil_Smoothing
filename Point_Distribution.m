function [ptd]=Point_Distribution(xstart,xend,npt,SpacingLaw)
%
% -------------------------------------------------------------------------
% PURPOSE - Compute a set of points between LE and TE. The intermediate 
%           points are computed according to various spacing rules.
% -------------------------------------------------------------------------
% Coded by: Massimiliano Di Muzio - November 2008
% -------------------------------------------------------------------------
% Basic description:
% SpacingLaw    1 = uniform
%               2 = full cosine, dense near both ends
%               3 = half cosine, dense near start                                          
%               4 = half sine, dense near end      
%               5 = read from file                                                                                                                   
% -------------------------------------------------------------------------

ptd(npt)=xend;
ptd(1)=xstart;
  
for i=1:npt-2
    tmp(i)=i;
    tmp(i)=tmp(i)/(npt-1);
end

if (SpacingLaw<1) || (SpacingLaw>5)
    SpacingLaw=1;
end

for i=2:npt-1
    switch SpacingLaw
        case 1            
            dx=tmp(i-1);                                % linear
            ptd(i)=xstart + (xend-xstart)*dx;
        case 2            
            dx=0.5*(1.0-cos(pi*tmp(i-1)));              % full cosine, dense near both ends
            ptd(i)=xstart + (xend-xstart)*dx;
        case 3            
            dx=(1.0-cos(pi/2.*tmp(i-1)));               % half cosine, dense near start
            ptd(i)=xstart + (xend-xstart)*dx;
        case 4
            dx=sin(pi/2.*tmp(i-1));                     % half sine, dense near end
            ptd(i)=xstart + (xend-xstart)*dx;
    end
end

if (SpacingLaw == 5)    
    ptd=load('xnew.dat');
    ptd=ptd';
end
end
