function [dx, dy, Cxy] = piv2Stack(stack, gridSize, extSize, min_Int, min_corr)

% PIV on a stack between slices z

    % --- Init ------------------------------------------------------------
    
    dim = size(stack);
    xg = 1:gridSize:dim(2);
    yg = 1:gridSize:dim(1);
    nL = dim(3);
    
    [X_, Y_] = meshgrid(xg, yg);
    Gxy = [X_(:) Y_(:)];

    Nxy = size(Gxy,1);
    
    % --- XY displacements --------------------------------------------

    % --- All layers

    fprintf('Computing xy displacements ...');
    tic

    dx = zeros(Nxy, nL);
    dy = zeros(Nxy, nL);
    Cxy = zeros(Nxy, nL);

    for l = 1:nL-1

        for k = 1:Nxy

            i1 = max(Gxy(k,2)-extSize, 1);
            i2 = min(Gxy(k,2)+extSize, dim(1));
            j1 = max(Gxy(k,1)-extSize, 1);
            j2 = min(Gxy(k,1)+extSize, dim(2));

            Sub1 = stack(i1:i2, j1:j2, l);
            Sub2 = stack(i1:i2, j1:j2, l+1);

            if mean(Sub1(:))<min_Int || mean(Sub2(:))<min_Int
                continue
            end

            [dx(k,l), dy(k,l), Cxy(k,l)] = fcorr(Sub1, Sub2);

        end

    end

    Cxy(Cxy<min_corr) = 0;

    fprintf(' %.02f sec\n', toc);
end