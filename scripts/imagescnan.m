
function imagescnan(IM, col)
% function imagescnan(IM, col)
% -- to display NaNs in imagesc as white/black

% white
if ~isempty(col)
    nanjet = [ 1,1,1; col  ];
else
    nanjet = [ 1,1,1; jet  ];
end
nanjetLen = length(nanjet); 
pctDataSlotStart = 2/nanjetLen;
pctDataSlotEnd   = 1;
pctCmRange = pctDataSlotEnd - pctDataSlotStart;

dmin = nanmin(IM(:));
dmax = nanmax(IM(:));
dRange = dmax - dmin;   % data range, excluding NaN

cLimRange = dRange / pctCmRange;
cmin = dmin - (pctDataSlotStart * cLimRange);
cmax = dmax;
imagesc(IM);
set(gcf,'colormap',nanjet);
caxis([cmin cmax]);

end