function reshaped = reshapePolygon(polygon1, polygon2)
    
    [cx1, cy1] = Functions.centroid(polygon1);
    [cx2, cy2] = Functions.centroid(polygon2);
    
    drift = [cx2-cx1, cy2-cy1];
    driftPoly2 = polygon2 - drift;
    
    % get the biggest corresponding
    if length(polygon1) >= length(polygon2)
        ind = dsearchn(driftPoly2, polygon1);
        reshaped = polygon2(ind, :);
    else
        ind = dsearchn(polygon1,driftPoly2);
        reshaped = polygon1(ind, :);
    end
    
end