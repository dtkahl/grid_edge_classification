function [cx cy crow ccol]=poly2rastermask_change_density(xpoly,ypoly,xll,yll,nrows,ncols,cellsize)
%POLY2RASTERMASK Computes the rows and columns of a raster grid that are
%intersected by a polyline
%   Given a list of points with x-y coordinates stored in xpoly and ypoly,
%   and information about a raster grid (xll,yll,cellsize,nrows,ncols),
%   this code will return a list of row/column pairs corresponding to
%   pixels intersected by the polyline. The code relies on the function
%   "improfile" (part of the image processing toolbox) which creates a set
%   of points along each polyline and then finds the cell contained by the
%   point. The density of points can be made higher by changing the
%   parameter "pointdensity".

pointdensity=5;

%test for NaN
pointlist=~isnan(xpoly) | ~isnan(ypoly);
xpoly=xpoly(pointlist);
ypoly=ypoly(pointlist);
npoly=length(xpoly);

%set up grid of cell-center coordinates
xc=xll-cellsize/2+cellsize*[1:ncols];
yc=yll-cellsize/2+cellsize*[nrows:-1:1];

%create image for input into improfile
A=zeros(nrows,ncols);

%find distance to determine the number of required points for profile
dist=0;
for p=1:npoly-1
    dx=xpoly(p+1)-xpoly(p);
    dy=ypoly(p+1)-ypoly(p);
    dist=dist+sqrt(dx^2+dy^2);
end
np=pointdensity*floor(dist/cellsize);
if np <=  1
    np =2;
end
%np
[cx,cy,c]=improfile(xc,yc,A,xpoly,ypoly,np);

for p=1:np
    ccol(p)=floor((cx(p)-xll)/cellsize)+1;
    crow(p)=nrows-floor((cy(p)-yll)/cellsize);   
end

end



