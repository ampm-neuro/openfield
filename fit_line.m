function [R, p, poly] = fit_line(X,Y, varargin)
%fits a line to a scatterplot of X,Y

if nargin ==3
    degree = varargin{1};
else
    degree = 1;
end

%force column vectors
X = X(:);
Y = Y(:);

%common points index
nidx = ~isnan(X) & ~isnan(Y);

%polynomial
poly=polyfit(X(nidx),Y(nidx), degree);

%fit line coordinates
fit_x = (min(X(nidx))-(range(X(nidx))/15)):(range(X(nidx))/100):...
    (max(X(nidx))+(range(X(nidx))/15));
fit_y=polyval(poly, fit_x);

%figure;
hold on; plot(X,Y,'o', 'color', [.7 .7 .7]);plot(fit_x, fit_y, 'k-', 'LineWidth', 2)
xlim([min(fit_x) max(fit_x)])

if degree == 1
    [R, p] = corr(X(nidx),Y(nidx));
elseif degree == 2
    [~,gof,~] = fit(X(nidx), Y(nidx), 'poly2', 'normalize', 'on');
    R = gof;
    p = nan;
else
    R = nan;
    p = nan;
end

end