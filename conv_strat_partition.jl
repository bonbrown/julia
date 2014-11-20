# conv_strat_partition.jl
# to partition pixels in radar observations of hurricanes into convective, stratiform and weak echo.
# Based on Appendix A of Didlake and Houze (MWR 2009), referred to here as DH; initially leaving tuned parameters to their 
# values, though they used airborne ELDORA observations, may not be appropriate for ground based NEXRAD.
# 
#   Bonnie Brown, November 2014
#    University of Hawaii, Manoa
#

using NetCDF

# define functions
# haversine function from http://rosettacode.org/wiki/Haversine_formula#Julia to calculate distance between points
haversine(lat1,lon1,lat2,lon2) = 2 .* 6372.8 .* asin(sqrt(sind((lat2-lat1)./2).^2 + cosd(lat1).* cosd(lat2).* sind((lon2 - lon1)./2).^2))
# get unstaggered height field from wrfout variable geopotential (base + perturbation) IN METERS
height(PHI) = ( PHI[:,:,2:end,:]+ PHI[:,:,1:end-1,:] ) ./ (2*9.81) ;

#------------------- INPUT OPTIONS ------------------------------------#
level = 2 # low-level altitude to conduct partitioning at (km)
a = 9     # arbitrary parameter
b = 45    # arbitrary parameter
Zti = 42  # convective threshold intensity (dBZ)
Zwe = 20  # weak echo threshold (dBZ)
km = 1000;

# path and name of file
nc = "/Users/brbrown/wrf/arthur2014/wsm6/wrfout_d04_2014-07-03_12:00:00"

# doing gridded radar or wrf?
wrf = bool(1)

#----------------------- WRF OUTPUT FIELDS ---------------------------#
if wrf
    print("doing WRF\n")
    REF = ncread(nc,"REFL_10CM");
    lat  = ncread(nc,"XLAT");
    lon = ncread(nc,"XLONG");
    PHB = ncread(nc,"PHB");
    PH = ncread(nc,"PH") + PHB;

    z = height(PH)./km;
    s = size(REF);
    ti = s[4];
    
    reflec = zeros(lat);
    
    print("Interpolating WRF reflectivity to height level\n")
    #interpolate reflectivity to level - perhaps put this in a separate function in the future
    for t = 1:ti
        data = squeeze(REF[:,:,:,t],4);
        for ii = 1:s[1] 
            for jj = 1:s[2]
                if ((z[ii,jj,1]< level) || (z[ii,jj,s[3]] > level))
                    klev = 0;
                    # find the nearest height levels
                    for kk = 2:s[3]
                        if z[ii,jj,kk] > level
                            klev = kk-1;
                            break
                        end
                    end

                    # linearly interpolate
                    m = (data[ii,jj,klev+1] - data[ii,jj,klev]) ./ (z[ii,jj,klev+1] - z[ii,jj,klev]); 
                    reflec[jj,ii,t] = m .* (level - z[ii,jj,klev]) + data[ii,jj,klev];

                else # if height is above or below model, set to NaN
                    reflec[ii,jj,t] = NaN;
                end
            end
        end
    end
#-------------------- GRIDDED RADAR FIELDS -------------------------#
else
    print("doing gridded radar\n")
    x = ncread(nc,"x0"); #km
    y = ncread(nc,"y0"); #km
    z = ncread(nc,"z0"); #km
    REF = ncread(nc,"REF"); #dBZ
    lat = ncread(nc,"lat0");
    lon = ncread(nc,"lon0");
    
    reflec = squeeze(REF[:,:,find(z.==level),:],3);
    s = size(reflec);
    ti = s[3]
end
#--------------------------------------------------------------------#
ncclose(nc)

#preallocate final product
csmask_write = zeros(reflec);
size(reflec)


for t = 1:ti
#---------------------------- TIME LOOP ----------------------------------------------------#
	# Zbg (background reflectivity; dBZ) is average of nonnegative and nonzero reflectivity 
	# within a radius of 11km around the grid point
	print("Beginning background reflectivity calculation\n")
	refl = squeeze(reflec[:,:,t],3);
	Zbg = NaN*ones(refl);
	s = size(refl);
	for n = 1:s[1]
	    for m = 1:s[2]
    	    dist = haversine(lat[n,m,t],lon[n,m,t],lat[:,:,t],lon[:,:,t]);  # find great circle distance from each point
        	tmp = refl[find(dist.<=11.0)];                # find only points within 11km
        	Zbg[n,m] = mean(tmp[find(tmp.>0)]);           # take mean of points within radius only if reflectivity is nonnegative and nonzero
    	end
	end

	print("Background reflectivity calculation complete\n")

	# Now define the convective center criterion delta Zcc first introduced by Steiner et al 1995. The cosine function used
	# by DH was introduced in Yuter and Houze (1997). If Z exceeds Zbg by delta Zcc, it is a convective center.

	dZcc = a*cos( (1/b) * (pi.*Zbg/2) );
	delZ = refl - Zbg;
	cc = find(delZ.>=dZcc);

	print("Convective center calculation complete\n")

	# define the convective radius R - this is the radius of points around a convective center which are also classified as convective
	R = zeros(Zbg);
	R[find(Zbg.<20)] = 0.5;
	R[find(Zbg.>=20)] = 0.5 + 3.5 * (Zbg[find(Zbg.>=20)] - 20)./15;
	R[find(Zbg.>=35)] = 4.0;      # this overwrites any values that were defined immediately above but where Zbg was over 35 (cannot put two logical statements in find function)
	
	print("Convective radius calculation complete\n")
	
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
	
	csmask_write[:,:,t] = csmask;
	# the remaining points stay at a value of 0 (zero) indicating stratiform

	print("classification into convective, weak echo, stratiform, and missing complete\n")
	print("TIME ",string(t)," COMPLETE\n\n")
#--------------------------- END TIME LOOP ---------------------------------------------------------#
end

# write the partition mask (csmask) to a new variable in the gridded radar or WRF output file
if wrf
	nccreate(nc,"CSMASK","west-east",0:s[1]-1,"south-north",0:s[2]-1,"Time",0:ti-1);
	ncwrite(csmask_write,nc,"CSMASK")
else	
	nccreate(nc,"CSMASK","x0",x,"y0",y)
	ncwrite(csmask,nc,"CSMASK")
end

print("Convective-stratiform partition written to new variable CSMASK in netcdf file\n")

