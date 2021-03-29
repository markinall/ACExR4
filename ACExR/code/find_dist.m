function distance=find_dist(lat1,lat2,long1,long2,units)

% Calculates the distance in metres between two points on earth
% Assumes a spherical earth
% Default Units are metres, specify km or nm if required
% Radius of the Earth is 6378.160 km
% lkm =  1.853nm

faktor = pi/180;

R=1000*12756.32/2;
c1=(cos(lat1*faktor)*cos(long1*faktor)*cos(lat2*faktor)*cos(long2*faktor));
c2=(cos(lat1*faktor)*sin(long1*faktor)*cos(lat2*faktor)*sin(long2*faktor));
c3=(sin(lat1*faktor)*sin(lat2*faktor));
d = R*acos(c1+c2+c3); 

if not(exist('units','var'))
   units='m';
end

switch units
case{'km'}
   distance=d/1000;
case{'nm'}
   distance=(d/1000)/1.853;
otherwise
   distance=d;
end

%This problem can be most easily solved by using spherical coordinates on the
%earth.  Have you dealt with those before?  Here's the transformation from
%spherical coordinates to normal rectangular coordinates, where a=latitude
%and b=longitude, and r is the radius of the earth:x = r Cos[a] Cos[b]
%y = r Cos[a] Sin[b]z = r Sin[a]
%Then we'll use the following property of the dot product (notated [p,q]):
%[p,q] = Length[p] * Length[q] * Cos[angle between p & q]
%Now, any vector that points to a point on the surface of the earth will have
%length r.  So the right side we have r^2 * Cos[angle between p & q].  On the
%left side, we can pull the r's out of the dot product, and cancel them with
%the r's on the right side.  Let t represent the angle between p and q.  Then
%if the latitude and longitude of our two cities, p and q, are (a1,b1) and
%(a2,b2), we have 
%Cos[a1] Cos[b1] Cos[a2] Cos[b2] + Cos[a1] Sin[b1] Cos[a2] Sin[b2] 
%    + Sin[a1] Sin[a2]  =  Cos[t]
%So you can compute the angle t as a function of a1, b1, a2, b2, which are
%the latitudes and longitudes of our cities p and q.  Then visualize what
%you've got:  draw a great circle through the points p and q.  This is just a
%plain old Joe-Schmoe circle of radius r, and the angle t is the angle of the
%arc that subtends p and q.  The problem from here on out is just figuring
%out what the arc length between p and q is.  The relevant formula is 
%Arc length = t/360 * 2Pi* r.  So that's your formula.  By substitution, wehave 
%Arccos[Cos[a1] Cos[b1] Cos[a2] Cos[b2] + Cos[a1] Sin[b1] Cos[a2] Sin[b2]
%    + Sin[a1] Sin[a2]]/360 * 2Pi * r
%Oh, by the way, West longitude means negative values of b, and South
%latitude means negative values of a.%define the (surface) distance we seek as D. 

