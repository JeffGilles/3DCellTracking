function [bins, f] = pdf_distance_n1(txy)

% --- Compute distances

D = [];
for t = unique(txy(:,1))'
    
    I = txy(:,1)==t;
    x = txy(I,2);
    y = txy(I,3);
    
    % Tesselation
    T = delaunay(x, y);
    
    % Links
    tmp = [T(:,1) T(:,2) ; T(:,1) T(:,3) ; T(:,2) T(:,3)];
    L = unique([min(tmp, [], 2) max(tmp, [], 2)], 'rows');
    
    % Distances
    D = [D ; sqrt((x(L(:,1)) - x(L(:,2))).^2 + (y(L(:,1)) - y(L(:,2))).^2)];

end

% --- Output

[f, bins] = ksdensity(D, linspace(min(D), max(D), 500));

% --- Display

if ~nargout
    
    clf
    hold on
    
    plot(bins, f, '-');
    
    xlabel('$d_{n1}$', 'Interpreter', 'LaTeX');
    ylabel('pdf', 'Interpreter', 'LaTeX');
    
    box on
    set(gca, 'YScale', 'log')
    
end