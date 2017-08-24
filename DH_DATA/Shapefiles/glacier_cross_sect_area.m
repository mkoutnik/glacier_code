function [area, thickmap] = glacier_cross_sect_area( x, y, xx, yy, Data, thickmap)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% find the x and y indexes of the endpoints of the cross-section
x1_index = find(abs(x - xx(1)) ==min(min(abs(x-xx(1)))));
x2_index = find(abs(x - xx(2)) ==min(min(abs(x-xx(2)))));
y1_index = find(abs(y - yy(1)) ==min(min(abs(y-yy(1)))));
y2_index = find(abs(y - yy(2)) ==min(min(abs(y-yy(2)))));

%plot the ice thickness map only once
if thickmap == 0
figure; 
thickmap = imagesc(x,y, Data);
set (gca, 'Ydir', 'normal')
colorbar
end

% plot the location of all cross sections
hold on; line ([x(x1_index) x(x2_index)],[y(y1_index) y(y2_index)],'linewidth', 2)

%extract the ice thickness values along the profile. cx and cy are the
%spatial coordinates of the points. They are not used here, but might be
%useful at some point
[cx, cy, cross_section] = improfile(Data, [x1_index x2_index], [y1_index y2_index]);

%for the Riemann-type integration, calculate the width of each box
distance = sqrt((x(x1_index)-x(x2_index)).^2 + (y(y1_index)-y(y2_index)).^2)./length(cross_section);

%finally, calculate the area of the cross-section
area = sum(cross_section.*distance);

end

