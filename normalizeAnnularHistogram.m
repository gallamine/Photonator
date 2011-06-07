%%

% Determine correct normalization for a 2d histogram of points on the x/y
% plane. By taking the histogram of the radius, the data is skewed by the
% area of the annular rings. This needs to be corrected. - WCC 6/7/11

% figure(4);
% cnt = hist3([rec_loc_final(:,1) rec_loc_final(:,2)],'Edges',{-2:1e-3:2 -2:1e-3:2});
% finalPdf = cnt./max(max(cnt));
% 
% imshow(finalPdf);
% colormap(jet);


figure(4);
cnt = hist3([recPosY' recPosX'],'Edges',{-0.3:1e-2:0.3 -0.3:1e-2:0.3});
% finalPdf = cnt./max(max(cnt));
imshow(finalPdf,[0 1])
colormap(jet);

[cnt,bins] = hist(distances,100);
binDelta = [bins] - [0 bins(1:end-1)];
cntNorm = cnt./(pi.*(2.*bins.*binDelta - binDelta.^2));
bar(bins,cntNorm);