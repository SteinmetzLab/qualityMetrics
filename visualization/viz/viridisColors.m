function cmap = viridisColors(n)
% viridisColors  Return an [n x 3] viridis colormap (purple -> yellow).
%
%   cmap = viridisColors(n)
%
% Dependency-free reproduction of matplotlib's "viridis" so plots match the
% reference figure regardless of MATLAB version. Row 1 is dark purple (low),
% row n is yellow (high).

if nargin < 1 || isempty(n), n = 256; end

% 11 evenly spaced anchor points from matplotlib viridis.
anchors = [ ...
    0.267004 0.004874 0.329415;
    0.282623 0.140926 0.457517;
    0.253935 0.265254 0.529983;
    0.206756 0.371758 0.553117;
    0.163625 0.471133 0.558148;
    0.127568 0.566949 0.550556;
    0.134692 0.658636 0.517649;
    0.266941 0.748751 0.440573;
    0.477504 0.821444 0.318195;
    0.741388 0.873449 0.149561;
    0.993248 0.906157 0.143936];

xa = linspace(0, 1, size(anchors,1));
if n == 1
    xq = 0.5;
else
    xq = linspace(0, 1, n);
end
cmap = [interp1(xa, anchors(:,1), xq, 'linear')', ...
        interp1(xa, anchors(:,2), xq, 'linear')', ...
        interp1(xa, anchors(:,3), xq, 'linear')'];
cmap = min(max(cmap, 0), 1);
end
