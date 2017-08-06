function [area, thickmap] = glacier_cross_sect_area( xx, yy, thickmap )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
addpath /Users/trevorhillebrand/Documents/MATLAB/Toolboxes/AntarcticMappingTools_v5.00/AntarcticMappingTools
% Read your data
[Data,RefMat]= arcgridread('/Users/trevorhillebrand/Documents/Antarctica/Darwin-Hatherton/Data/QGIS Data/Bed and Ice thickness data (from M. Riger-Kusk)/icethic.asc');

[nrows,ncols,~]=size(Data);
[row,col]=ndgrid(1:nrows,1:ncols);
[ygrid,xgrid]=pix2latlon(RefMat,row,col);
x = xgrid(1,:);
y = ygrid(:,1)';
xres = gradient(x);

x1_index = find(abs(x - xx(1)) ==min(min(abs(x-xx(1)))));
x2_index = find(abs(x - xx(2)) ==min(min(abs(x-xx(2)))));
y1_index = find(abs(y - yy(1)) ==min(min(abs(y-yy(1)))));
y2_index = find(abs(y - yy(2)) ==min(min(abs(y-yy(2)))));

if thickmap == 0
figure; 
thickmap = imagesc(x,y, Data);
set (gca, 'Ydir', 'normal')
end

hold on; line ([x(x1_index) x(x2_index)],[y(y1_index) y(y2_index)],'linewidth', 2)

[cx, cy, cross_section] = improfile(Data, [x1_index x2_index], [y1_index y2_index]);

distance = sqrt((x(x1_index)-x(x2_index)).^2 + (y(y1_index)-y(y2_index)).^2)./length(cross_section);

area = sum(cross_section.*distance);

end

