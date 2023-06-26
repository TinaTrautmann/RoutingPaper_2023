function [lats,lons]=PixelCoord2LatLon(x,y,ul_lat,ul_lon,res)
% [lats,lons]=PixelCoord2LatLon(x,y,ul_lat,ul_lon,res)
% ul_lat=90;
% ul_lon=-180;
% res=0.0833333333;

lats=ul_lat-y.*res+res/2;
lons=ul_lon+x.*res-res/2;

end