function [hdDir] = compute_hd(x,y)

xy(:,1) = x;
xy(:,2) = y;
% calculate the angle:
myDiff = diff(xy);
hdDir = mod(atan2(myDiff(:,1), myDiff(:,2))*180/pi, 360);

end