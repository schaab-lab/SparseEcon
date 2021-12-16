
% FIGURE FORMATTING
fontname = 'Times';

set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultLineLineWidth', 3.5);

set(groot, 'defaultAxesXGrid', 'on')
set(groot, 'defaultAxesYGrid', 'on')
set(groot, 'defaultAxesFontSize', 12)
set(groot, 'defaultAxesFontName', fontname);
set(groot, 'defaultTextFontName', fontname);

colorPalette = [000,076,153; ...
                255,188,000; ...
                039,169,65; ...
                000,205,108; ...
                064,173,090; ...
                202,091,035] ./ 255;

set(groot, 'defaultAxesColorOrder', colorPalette);
