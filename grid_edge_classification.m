%% Grid Edge Classification
% Authors:Daniel Kahl &  Brett Sanders
%         University of California, Irvine
% Source: https://doi.org/10.1016/j.advwatres.2022.104287

% Contact: bsanders@uci.edu, dkahl1@uci.edu

%The purpose of this program is to create a list of grid edges that are
%classified to represent the effects of levees (or walls) in a raster grid
%flood inundation model. We assume that every grid edge in a flood model
%can be assigned to one of the following three classes:
%    -1 : levee/wall located within one cell to the left/bottom of the edge
%    1 : levee/wall located within one cell to the right/top of the edge
%    0 : no levee/wall within one cell to the left or right

%The input to this program includes two files:
%    1) a levee polyline file (here supplied as a shape file, *.shp) 
%    2) and a DEM file (here supplied as a *.tif file).

%Within the code, several additional inputs must be specified by the user
%    1) the location and input filenames for the DEM and levee shapefile.
%    2) the upscale factor, which represents the ratio of the coarse grid
%       size to the fine grid size. 
%    3) the row and column indices of grid cell "plugs" that are needed to
%       fully encircle interior channel areas with cells. That is, before
%       running the "imfill.m" command, it is important to make sure there
%       are not any openings/gaps in loop of cells around each channel.
%    4) the row and column indices of cells that fall within "holes"
%       See Fig. 2b in Kahl et al. (2022). This information supports the
%       "imfill.m" routine (hole filling operation)
%    5) the row and column indices of grid edges that are incorrectly
%       classified as levees and, if implemented, would create a
%       non-physical blockage effect within a channel. We term these edges
%       as "weirs" since they would present a barrier similar to a weir
%       if not removed. See Fig. 2c in Kahl et al. (2022).
%       Note that there are separate lists for horizontal and vertical
%       "weir" edges that neee to be removed.

%The output of this program is an ascii file (*.levee) containing two 
%lists of edges:
%    1) A list of vertical edges defined by the row and cell index
%    2) A list of horizontal edges defined by the row and cell index

%These lists contain only those edges classified as -1 or 1.

%The method that is implemented in this code is described in a research
%paper by Kahl et al. (2022) appearing in Advances in Water Resources
% https://doi.org/10.1016/j.advwatres.2022.104287

%The code was developed to improve the representation of levees in
%dual-grid flood inundation models - which resolve topographic data on a
%fine grid and perform solution updates on a relatively coarse grid. In
%particular, the authors developed the code in support of the PRIMo model
%described by Sanders and Schubert (2019)
% https://doi.org/10.1016/j.advwatres.2019.02.007

%However, we anticipate that the method (and code) could support any raster
%grid flood inundation model. That is, the code could be applied to come up
%with a list of grid edges that, when connected, can approximately
%represent a levee that would otherwise fall between neighboring grid
%edges.

%It is left to the user to reformat the output for supporting other
%applications.

%This code relies on commands from the following matlab toolboxes
%       Image Processing Toolbox
%       Mapping Toolbox

%======================================================================

%This code is set up to demonstrate the method in a testcase. 
%First, it is configured to load the following files:
%    1) ex_levees.shp (levee polyline data)
%    2) ex_DTM.tif  (DEM file)

%Second, it is configured to use an upscale factor of 10
%    usf=10

%Third, it is configured to add zero "plugs".

%fourth, it is configured to use two holes
%    hole 1: row=224, col=1
%    hole 2: row=1, col=172

%Fifth, it is configured to remove four vertical weirs and two horizontal
%weirs with the following grid edge indices
%    Vweir 1: row=160, col=173
%    Vweir 2: row=161, col=173
%    Vweir 3: row=96, col=47
%    Vweir 4: row=96, col=48

%    Hweir 1: row=97, col=46
%    Hweir 2: row=96, col=47

%We note that the rows and columns required of plugs, holes and weirs are
%dependent on the upscale factor. 

%It is left as an exercise to the user to configure plugs, holes and weirs
%for other upscale factors such as 5, 20, and 50

clear;close all;clc
%% Controls
usf=10;
% Specify directory and file names for DTM/DEM, levee, ends and plugs
% shapefiles.
dir1 = 'data/';
% Specify DEM/DTM
fname_DEM='ex_DTM';
% Levee .shp File
fname_levee= 'ex_levees.shp';
% Save output .levee file
makeleveefile=1; %0=off, 1=on
%% Manual Input of plugs, holes and weirs
% Plugs are used to fully encircle channel areas (create closed loops)
% within the logical mask, where "true" means the cell contains a levee
plugs = []; %default

% Holes are used to guide the "imfill.m" command, telling it where filling
% begins. 

holes=[];  %default
% Settings for usf=10 
holes=[224 1; 1 172];

% Weirs are edges incorrectly assigned to be a levee, and should be removed
% from the list of grid edges classified as levees. There are both
% horizontal and vertical weirs

Vweir=[]; %default
Hweir=[]; %default

% Settings for usf=10
Vweir=[160 173;
   161 173;
   96 47
   96 48];
Hweir=[97 46;
   96 47];

%% Grid Edge Classification Method

%Step 1 of the 5-step method described in Kahl et al. (2022)
%===========================================================
% read DEM (fine grid) and shapefile (polyline) data 
% edit polyline data if needed
fname_poly = strcat(dir1, fname_levee);
fname_DEM = strcat(dir1, fname_DEM);

[dem, R] = geotiffread(fname_DEM);
[levees,A] = shaperead(fname_poly, 'UseGeoCoords', true);

xll=R.LongitudeLimits(1);
yll=R.LatitudeLimits(1);
nrows=R.RasterSize(1);
ncols=R.RasterSize(2);
dx=(R.LongitudeLimits(2)-R.LongitudeLimits(1))/ncols;
dy=(R.LatitudeLimits(2)-R.LatitudeLimits(1))/nrows;
cellsize=dx;
nline=length(levees);

%Create upscale (coarse) grid
nrowu=floor(nrows/usf);
ncolu=floor(ncols/usf);
cellsizeu=cellsize*usf;
%Check: nrows and ncols should both be a multiple of usf;
if (nrows-nrowu*usf > 0 || ncols-ncolu*usf > 0)
    %will trim raster and reset xll,yll
    xll_0=xll; %save original xll
    yll_0=yll; %save original yll
    xll=xll_0; %new xll
    yll=yll_0+cellsize*(nrows-nrowu*usf); %new yll
end
z=dem(1:nrowu*usf,1:ncolu*usf); 
%If DEM is trimmed, data in upper left corner is retained.

%Specify cell center and cell vertex coordinates
xuc=xll+cellsizeu*[1:ncolu]-0.5*cellsizeu; %upscale grid
yuc=yll+cellsizeu*[nrowu:-1:1]-0.5*cellsizeu; %upscale grid
xuv=xll+cellsizeu*[0:ncolu]; %upscale grid
yuv=yll+cellsizeu*[nrowu:-1:0]; %upscale grid
xc=xll+cellsize*[1:ncolu*usf]-0.5*cellsize; %trimmed fine-res grid
yc=yll+cellsize*[nrowu*usf:-1:1]-0.5*cellsize; %trimmed fine-res grid

%Make a list of edges to plot grid
for col=1:ncolu+1
    VertEdgex(1,col)=xuv(col);
    VertEdgex(2,col)=xuv(col);
    VertEdgey(1,col)=yuv(1);
    VertEdgey(2,col)=yuv(nrowu+1);
end
for row=1:nrowu+1
    HorEdgex(1,row)=xuv(1);
    HorEdgex(2,row)=xuv(ncolu+1);
    HorEdgey(1,row)=yuv(row);
    HorEdgey(2,row)=yuv(row);
end

%Create a more granular DEM for visualization purposes
for row=1:nrowu
    for col=1:ncolu
        sumz=0;
        for r=1:usf
            for c=1:usf
                ri=usf*(row-1)+r;
                ci=usf*(col-1)+c;
                sumz=sumz+z(ri,ci);
            end
        end
        zu(row,col)=sumz/(usf*usf);
    end
end

figure(1)
imagesc(xc,yc,z)
axis square
set(gca,'YDir','normal')
title('Step 1: DEM and Levee Polylines')
% axis(XYlims)
hold on
plot(levees(1).Lon,levees(1).Lat,'r-','LineWidth',2)
for i=2:nline
    plot(levees(i).Lon,levees(i).Lat,'r-','LineWidth',2)
end
hold off

%Step 2 of the 5-step method described in Kahl et al. (2022)
%===========================================================
%Intersect polylines with raster grid and create a mask of affected cells
nline=length(levees);
leveemask1=zeros(nrowu,ncolu);
for l=1:nline
    xpoly=levees(l).Lon;
    ypoly=levees(l).Lat;
    [px py prow pcol]=poly2rastermask(xpoly,ypoly,xll,yll,nrowu,ncolu,cellsizeu);
    np=length(prow);
    for p=1:np
        if (prow(p) >= 1 && prow(p) <= nrowu && pcol(p) >= 1 && pcol(p) <= ncolu)
           leveemask1(prow(p),pcol(p))=1;
       end
    end
%     end %
end

%Step 1A: Add plugs to fully encircle channel areas by the levee mask
for p=1:size(plugs,1)
    leveemask1(plugs(p,1),plugs(p,2))=1;
end
LM1=logical(leveemask1);

figure(2)
imagesc(LM1)
title('Step 2: Mask of cells intersected by levees')

%Step 3 of the 5-step method described in Kahl et al. (2022)
%===========================================================
%Hole-filling step
LM2=imfill(LM1,'holes'); %first fill without holes specified
%Tip: stop the code here when developing a new application, and visualize
%the mask "LM2" to see where holes will be needed to fill all channel
%areas. Can use "imagesc(LM2)" to visualize mask
nholes=size(holes,1);
for i=1:nholes
    LM2=imfill(LM2,holes(i,:));
end

figure(3)
imagesc(LM2)
title('Step 3: Mask of areas covering flood channels')

%Step 4 and 5 of the 5-step method described in Kahl et al. (2022)
%===========================================================
% identify edges that fall along the perimeter of the LM2 region

%First sweep over vertical edges
notVweir=ones(nrowu,ncolu+1); %this will be used to exclude channel ends
for e=1:size(Vweir,1)
    notVweir(Vweir(e,1),Vweir(e,2))=0;
end
notVweir=logical(notVweir);
Vedge=zeros(nrowu,ncolu+1);
e=0;
for row=1:nrowu
    for col=2:ncolu
        if xor(LM2(row,col-1),LM2(row,col)) && notVweir(row,col) %edge here
            if LM2(row,col)
                Vedge(row,col)=1; %levee on the right
            else
                Vedge(row,col)=-1; %levee on the left
            end
            e=e+1;
            VedgeX(1,e)=xuv(col);
            VedgeX(2,e)=xuv(col);
            VedgeY(1,e)=yuv(row);
            VedgeY(2,e)=yuv(row+1);
        end
    end
end
nVedge=e;

%Second sweep over horizontal edges
notHweir=ones(nrowu+1,ncolu); %this will be used to exclude channel ends
for e=1:size(Hweir,1)
    notHweir(Hweir(e,1),Hweir(e,2))=0;
end
notHweir=logical(notHweir);
Hedge=zeros(nrowu+1,ncolu);
e=0;
for col=1:ncolu
    for row=2:nrowu
        if xor(LM2(row-1,col),LM2(row,col)) && notHweir(row,col) %edge here === notHweir
            if LM2(row,col)
                Hedge(row,col)=-1; %levee below edge
            else
                Hedge(row,col)=1; %levee above edge
            end
            e=e+1;
            HedgeX(1,e)=xuv(col);
            HedgeX(2,e)=xuv(col+1);
            HedgeY(1,e)=yuv(row);
            HedgeY(2,e)=yuv(row);
        end
    end
end
nHedge=e;

figure(4)
imagesc(xc,yc,z)
axis square
set(gca,'YDir','normal')
title('Step 4: Grid edges classified as levees (in white)')
% axis(XYlims)
hold on
plot(levees(1).Lon,levees(1).Lat,'r-','LineWidth',2)
for i=2:nline
    plot(levees(i).Lon,levees(i).Lat,'r-','LineWidth',2)
end
line(HorEdgex,HorEdgey,'Color','k')
line(VertEdgex,VertEdgey,'Color','k')
plot(VedgeX,VedgeY,'w-','LineWidth',2)
plot(HedgeX,HedgeY,'w-','LineWidth',2)
hold off


%% now write a .levee file for raster grid flood inundation model
fname='leveetest.levee';
%fname=[fileroot '/' fileroot '.levee'];
fileID=fopen(fname,'w');
fprintf(fileID,'usf \n');
fprintf(fileID,'%d  \n',usf);
fprintf(fileID,'nVedge \n');
fprintf(fileID,'%d  \n',nVedge);
fprintf(fileID,'row, col, leveeloc (-1=left, 1=right) \n');
for row=1:nrowu
    for col=1:ncolu+1
        if Vedge(row,col) ~= 0
                A=[row col Vedge(row,col)];
                fprintf(fileID,'%d %d %d\n',A);
        end 
    end
end
fprintf(fileID,'nHedge \n');
fprintf(fileID,'%d  \n',nHedge);
fprintf(fileID,'row, col, leveeloc (-1=below, 1=above) \n');
for col=1:ncolu
    for row=1:nrowu+1
        if Hedge(row,col) ~= 0
                A=[row col Hedge(row,col)];
                fprintf(fileID,'%d %d %d\n',A);
        end 
    end
end
fclose(fileID);
