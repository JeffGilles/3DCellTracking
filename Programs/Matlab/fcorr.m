function [dx, dy, c] = fcorr(Img1, Img2, varargin)
%fcorr FFT-based image correlation
%   FCORR(IMG1, IMG2) computes the positions of the best image correlation 
%   between IMG1 and IMG2.
%

% === Input variables =====================================================

in = inputParser;
in.addParameter('method', 'parabolic', @ischar);
in.parse(varargin{:});

% =========================================================================

% --- Perform FFT
res = ifftshift(ifft2(fft2(Img1).*conj(fft2(Img2))));

% --- Get maximum
[c, I] = max(abs(res(:)));
[y0, x0] = ind2sub(size(res), I);

if nargout==3
    r1 = ifftshift(ifft2(fft2(Img1).*conj(fft2(Img1))));
    r2 = ifftshift(ifft2(fft2(Img2).*conj(fft2(Img2))));
    c = c/sqrt(max(abs(r1(:))))/sqrt(max(abs(r2(:))));
end

switch in.Results.method
    
    case 'rough'
        
        % --- Rough maximum estimation
        
        dx = size(res,2)/2 + 1.5 - x0;
        dy = size(res,1)/2 + 1.5 - y0;

    case 'parabolic'
        
        % --- Parabolic estimation
        
        u = x0 + [-1 0 1];
        if u(1)<1 || u(3)>size(res,2)
            dx = size(res,2)/2 + 1.5 - x0;
        else
        v = sqrt(res(y0, u));
        dx = size(res,2)/2 + 1.5 + real((u(3)^2*(v(1)-v(2)) + u(2)^2*(v(3)-v(1)) + u(1)^2*(v(2)-v(3)))/...
                                 (u(3)*(v(2)-v(1)) + u(2)*(v(1)-v(3)) + u(1)*(v(3)-v(2)))/2);
        end
                             
        u = y0 + [-1 0 1];
        if u(1)<1 || u(3)>size(res,1)
            dy = size(res,1)/2 + 1.5 - y0;
        else
            v = sqrt(res(u, x0));
            dy = size(res,1)/2 + 1.5 + real((u(3)^2*(v(1)-v(2)) + u(2)^2*(v(3)-v(1)) + u(1)^2*(v(2)-v(3)))/...
                (u(3)*(v(2)-v(1)) + u(2)*(v(1)-v(3)) + u(1)*(v(3)-v(2)))/2);
        end

        % --- Refine correlation value

% % %         if nargout==3   
% % %             [X, Y] = meshgrid(1:size(res,2), 1:size(res,1));
% % %             c = interp2(X,Y,V,x0)
% % %         end
        
end

