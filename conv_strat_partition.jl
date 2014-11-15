
# conv_strat_partition.jl
# to partition pixels in radar observations of hurricanes into convective, stratiform and weak echo.
# Based on Appendix A of Didlake and Houze (MWR 2009), referred to here as DH; initially leaving tuned parameters to their 
# values, though they used airborne ELDORA observations, may not be appropriate for ground based NEXRAD.
# 
#   Bonnie Brown, November 2014
#    University of Hawaii, Manoa
#

using NetCDF

# haversine function from http://rosettacode.org/wiki/Haversine_formula#Julia to calculate distance between points
haversine(lat1,lon1,lat2,lon2) = 2 .* 6372.8 .* asin(sqrt(sind((lat2-lat1)./2).^2 + cosd(lat1).* cosd(lat2).* sind((lon2 - lon1)./2).^2))

nc = "/Users/brbrown/julia_scripts/test.nc"
#nc = "/Users/brbrown/data/arthur2014/radar/KLTX/output/20140703/20140703/ncf_20140703_210457.nc"
x = ncread(nc,"x0"); #km
y = ncread(nc,"y0"); #km
z = ncread(nc,"z0"); #km
REF = ncread(nc,"REF"); #dBZ
lat = ncread(nc,"lat0");
lon = ncread(nc,"lon0");

ncclose(nc)

level = 2 # low-level altitude to conduct partitioning at (km)
a = 9     # arbitrary parameter
b = 45    # arbitrary parameter
Zti = 42  # convective threshold intensity (dBZ)
Zwe = 20  # weak echo threshold (dBZ)


# Zbg (background reflectivity; dBZ) is average of nonnegative and nonzero reflectivity 
# within a radius of 11km around the grid point
refl = REF[:,:,find(z.==level),:];
refl = squeeze(refl,[3 4]);
Zbg = NaN*ones(refl);
s = size(refl);
for n = 1:s[1]
    for m = 1:s[2]
        dist = haversine(lat[n,m],lon[n,m],lat,lon);  # find great circle distance from each point
        tmp = refl[find(dist.<=11.0)];                # find only points within 11km
        Zbg[n,m] = mean(tmp[find(tmp.>0)]);           # take mean of points within radius only if reflectivity is nonnegative and nonzero
    end
end

print("background reflectivity calculation complete")

# Now define the convective center criterion delta Zcc first introduced by Steiner et al 1995. The cosine function used
# by DH was introduced in Yuter and Houze (1997). If Z exceeds Zbg by delta Zcc, it is a convective center.

dZcc = a*cos( (1/b) * (pi.*Zbg/2) );
delZ = refl - Zbg;
cc = find(delZ.>=dZcc);

print("convective center calculation complete")

# define the convective radius R - this is the radius of points around a convective center which are also classified as convective
R = zeros(Zbg);
R[find(Zbg.<20)] = 0.5;
R[find(Zbg.>=20)] = 0.5 + 3.5 * (Zbg[find(Zbg.>20)] - 20)./15;
R[find(Zbg.>=35)] = 4.0;      # this overwrites any values that were defined immediately above but where Zbg was over 35 (cannot put two logical statements in find function)

print("convective radius calculation complete\n")

# Classify all convective centers, points within a convective radius, and points exceeding the convective intensity threshold 
# as convective points (1). Classify points beneath the Zwe threshold as weak echoes (2), and everything else as stratiform (0). 
# Missing data should remain missing (-9999)

csmask = zeros(refl);
csmask[cc] = 1;
ri = R[cc];
for r = 1:length(ri)
    dist = haversine(lat[cc[r]],lon[cc[r]],lat,lon);  # find great circle distance from each point
    pts = find(dist.<=ri[r]);
    csmask[pts] = 1;
end
csmask[find(refl.>=Zti)] = 1;
csmask[find(refl.<Zwe)] = 2;
csmask[find(refl.<-900)] = refl[find(refl.<-900)]; #keep missing missing

# the remaining points stay at a value of 0 (zero) indicating stratiform

print("classification into convective, weak echo, stratiform, and missing complete\n")

# write the partition mask (csmask) to a new variable in the gridded radar or WRF output file
nccreate(nc,"CSMASK","x0",x,"y0",y)
ncwrite(csmask,nc,"CSMASK")

print("convective-stratiform partition written to new variable in netcdf file")

