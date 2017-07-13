function [width, distance_along_centerline] = glacier_widths(flowline_path, buffer_path, nodes)
%glacier_widths automatically measures widths of a polygon (glacier outline
%shapefile) perpendicular to a flowline (also a shapefile)
%  
%flowline = Polar stereographic pairs of flowline coordinates, starting at
%the glacier mouth. Use Qgis to make this, then use interp1 to densify.
%Read into matlab with shaperead
%buffer = the polygon of the glacier outline, made in qgis. use shaperead
%nodes = number of evenly-spaced points at which to measure width.
%% get shapefiles into usable format. Buffer is fine the way it is.


buffer = shaperead(buffer_path);

flowline_struct = shaperead(flowline_path);
flowline_struct.X = flowline_struct.X(1:end-1);
flowline_struct.Y = flowline_struct.Y(1:end-1);

flowx = linspace(flowline_struct.X(1), flowline_struct.X(end), nodes);
flowy = interp1(flowline_struct.X, flowline_struct.Y, flowx);
flowline = [flowx; flowy]';


%% Normal lies in the null space of the matrix A - B, where A and B are each
%(x,y) pairs
%calculate normals to left and right of start and end of flowline (looking
%upstream)
temp_halfwidth = 5e4; %this just needs to be wider than any point along the flowline. 50 km is a good halfwidth

Lnormal_start = flowline(1,:) + temp_halfwidth.*(null(flowline(2,:) - flowline(1,:)))';
Rnormal_start = flowline(1,:) - temp_halfwidth.*(null(flowline(2,:) - flowline(1,:)))';
Lnormal_end = flowline(end,:) + temp_halfwidth.*(null(flowline(end,:) - flowline(end-1,:)))';
Rnormal_end = flowline(end,:) - temp_halfwidth.*(null(flowline(end,:) - flowline(end-1,:)))';

Lnormal = Lnormal_start;
Rnormal = Rnormal_start;

for jj = 2:size(flowline,1)-1
Lnormal_temp = flowline(jj, :) + temp_halfwidth.*(null(flowline(jj+1,:)-flowline(jj-1,:)))';
Rnormal_temp = flowline(jj, :) - temp_halfwidth.*(null(flowline(jj+1,:)-flowline(jj-1,:)))';

Lnormal = [Lnormal; Lnormal_temp];
Rnormal = [Rnormal; Rnormal_temp];
end

Lnormal = [Lnormal; Lnormal_end];
Rnormal = [Rnormal; Rnormal_end];

%% Now find the points where these normal vectors intersect the buffer
count = 1;
for jj = 1:size(flowline,1)
  
    normalx = [Rnormal(jj,1), Lnormal(jj,1)];
    normaly = [Rnormal(jj,2), Lnormal(jj,2)];
    
    % This gives the x and y coordinates of the points where widthlines
    % intersect the buffer
    [xx,yy] = polyxpoly(normalx, normaly, buffer.X, buffer.Y); 
    
    plot(xx,yy)
    hold on
    
    %The width is then the euclidian distance between these two points
   
    if isempty(xx) == 0 && isempty(yy) == 0 
        width(count) = pdist2([xx(1) yy(1)],[xx(2) yy(2)]);
    else
        width(count) = NaN;
    end
    count = count+1;
end

plot(buffer.X, buffer.Y)
plot(flowline(:,1), flowline(:,2), 'k', 'linewidth', 2)
daspect([1 1 1])

distance_vec = [0];
for jj = 2:size(flowline,1)
distance_temp= [pdist2([flowline(jj-1,1) flowline(jj-1,2)], [flowline(jj,1) flowline(jj,2)])'];

distance_vec  = [distance_vec; distance_temp];
end

distance_along_centerline = cumsum(distance_vec);
end






