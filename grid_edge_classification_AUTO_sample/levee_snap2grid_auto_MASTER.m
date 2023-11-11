%% Grid edge classification: Auto
%This Matlab program is designed to input two files:
%   (1) a shape file with polyline data
%   (2) a DEM file with topographic heights
%
%and then output a next file (leveetest.levee) that lists out edges of a
%primo computational grid that can be used to approximate the position of
%the levee/wall. The file "leveetest.levee" is designed to be loaded into
%PRIMO.
%
%The user of this code needs to do several things:
%
%   (1) input a DEM and shape file with levee polylines
%   (2) specify the upscale factor
%   (3) manually specify the row/col where holes exist for hole-filling
%   purposes using the "imfill" function in Matlab
%   (4) maniually specify the row/col at the end of channels that would
%   otherwise "plug" the channel and not allow water to flow properly.
%
%This code is set up to run a test case called "leveetest" using an upscale
%factor of either 10, 20 or 50. The holes and plugs need to be manually
%edited to successfully build a .levee file for each case.
clear; close all; clc

%% Controls
usf=5;
dir_gis = []%'data\';
makeleveefile=1; %0=off, 1=on

fname_DEM= []; %'Holy_Fire_lidar_0918_PRIMo_Leach_ALLculverts_project_ND9999';
fname_poly= []%'Levee_Leach_combine.shp';
fname_ends=[]%'Leach_ends.shp';
fname_plugs=[]%'Leach_plugs.shp';

%% Levee file processing
%holes and plugs (optional)
holes=[];
plugs = [];

%ends of channels need to be opened manually
Vends=[];
Hends=[];

%read raster and shapefile (polyline) data 
[dem, R] = geotiffread([dir_gis,fname_DEM]);
[levees,A] = shaperead([dir_gis,fname_poly], 'UseGeoCoords', true);
[levees_e,A_e] = shaperead([dir_gis,fname_ends], 'UseGeoCoords', true);
[levees_p,A_p] = shaperead([dir_gis,fname_plugs], 'UseGeoCoords', true);

% Specify loaded DEM attributes
% xll=R.LongitudeLimits(1);
% yll=R.LatitudeLimits(1);
xll=R.XWorldLimits(1);
yll=R.YWorldLimits(1);
nrows=R.RasterSize(1);
ncols=R.RasterSize(2);
% dx=(R.LongitudeLimits(2)-R.LongitudeLimits(1))/ncols;
% dy=(R.LatitudeLimits(2)-R.LatitudeLimits(1))/nrows;
dx=(R.XWorldLimits(2)-R.XWorldLimits(1))/ncols;
dy=(R.YWorldLimits(2)-R.YWorldLimits(1))/nrows;
cellsize=dx;

%plot both raster and shapefile
% figure,
% hold on;
% % geoshow(dem,R,'DisplayType','texturemap');
% geoshow(levees,'Color', 'red','LineWidth',2);
% hold off

%Create upscale grid
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

%Create a more granular DEM for visualization
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

% ====================================
% Creating mask of Ends Locations
nline2=length(levees_e);
leveemask1_e=zeros(nrowu,ncolu);
for l=1:nline2
    xpoly2=levees_e(l).Lon;
    ypoly2=levees_e(l).Lat;
    [px2 py2 prow2 pcol2]=poly2rastermask_change_density(xpoly2,ypoly2,xll,yll,nrowu,ncolu,cellsizeu);
    np2=length(prow2);
    for p=1:np2
        if (prow2(p) >= 1 && prow2(p) <= nrowu && pcol2(p) >= 1 && pcol2(p) <= ncolu)
           leveemask1_e(prow2(p),pcol2(p))=1;
       end
    end
%     end %
end


leveemask1_e=imfill(leveemask1_e,'holes');

figure
imagesc(leveemask1_e)
axis square
% === Invert mask here, because want 0 in location to not want edge
leveemask1_ee = not(leveemask1_e);

% creating Vends, Hends 
[ends_row, ends_col] = find(leveemask1_ee==0);
Vends = [ends_row, ends_col+1;ends_row, ends_col];
Hends = [ends_row+1, ends_col;ends_row, ends_col];

% figure
% imagesc(leveemask1_ee)
%filling holes...

% ==========================
% === plugs AUTO
% ==========================
% PLUGS Creating mask of Plugs Locations
nline3=length(levees_p);
leveemask1_p=zeros(nrowu,ncolu);
for l=1:nline3
    xpoly3=levees_p(l).Lon;
    ypoly3=levees_p(l).Lat;
    [px3 py3 prow3 pcol3]=poly2rastermask(xpoly3,ypoly3,xll,yll,nrowu,ncolu,cellsizeu);
    np3=length(prow3);
    for p=1:np3
        if (prow3(p) >= 1 && prow3(p) <= nrowu && pcol3(p) >= 1 && pcol3(p) <= ncolu)
           leveemask1_p(prow3(p),pcol3(p))=1;
       end
    end
%     end %
end

for ii = 1:nrowu
    for ww = 1:ncolu
        if leveemask1_p(ii,ww) ==1;
            leveemask1(ii,ww)=1;
        end
    end
end
% ==========================



%Add plugs
for p=1:size(plugs,1)
    leveemask1(plugs(p,1),plugs(p,2))=1;
end

%Fill in space bounded by levees to create a levee mask
LM1=logical(leveemask1);
LM2=imfill(LM1,'holes');
nholes=size(holes,1);
for i=1:nholes
    LM2=imfill(LM2,holes(i,:));
end

imagesc(LM2)
%Compute perimeter of the levee mask (LM2)

%First sweep over vertical edges
notVends=ones(nrowu,ncolu+1); %this will be used to exclude channel ends
for e=1:size(Vends,1)
    notVends(Vends(e,1),Vends(e,2))=0;
end
notVends=logical(notVends);
Vedge=zeros(nrowu,ncolu+1);
e=0;
for row=1:nrowu
    for col=2:ncolu
        if xor(LM2(row,col-1),LM2(row,col)) && notVends(row,col) %edge here
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
notHends=ones(nrowu+1,ncolu); %this will be used to exclude channel ends
for e=1:size(Hends,1)
    notHends(Hends(e,1),Hends(e,2))=0;
end
notHends=logical(notHends);
Hedge=zeros(nrowu+1,ncolu);
e=0;
for col=1:ncolu
    for row=2:nrowu
        if xor(LM2(row-1,col),LM2(row,col)) && notHends(row,col) %edge here === notHends
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

%%
figure
imagesc(xc,yc,z)
axis square
set(gca,'YDir','normal')
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

figure,
imagesc(LM2)
axis square

%now write a .levee file for primo

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
    


if leveemask1_ee(1,1)
    disp('hey')
else
    disp('not')
end


