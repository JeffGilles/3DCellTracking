function [cx, cy] = centroid(polygon)

% return the centroid of the polygon

    shape = polyshape(polygon(:,1),polygon(:,2));
    [cx, cy] = shape.centroid;
end