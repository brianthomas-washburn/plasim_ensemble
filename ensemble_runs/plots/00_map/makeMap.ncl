begin

    ; --- Historical Orography NetCDF File ---
    
    ; Silurian period (425 Mya) map
    ;fi = addfile("../../historical_geography/Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc_v2/Map74_PALEOMAP_1deg_Late_Silurian_425Ma.nc","r")
    
    ; Present day (0 Mya) map
    fi = addfile("../../historical_geography/Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc_v2/Map01_PALEOMAP_1deg_Holocene_0Ma.nc","r")
    
    fidat = fi->z
    
    dsizes = dimsizes(fidat)
    nlat = dsizes(0)
    nlon = dsizes(1)
    
    dat_lsm = new((/nlat,nlon/),"double")
    copy_VarMeta(fidat,dat_lsm)
    
    do it = 0,nlat-1
        do it2 = 0,nlon-1
            if (fidat(it,it2) .gt. 0)
                dat_lsm(it,it2) = 1.0; land mass
            else
                dat_lsm(it,it2) = 0.0; sea position
            end if
        end do
    end do
    
    printVarSummary(dat_lsm)
    oro = dat_lsm
   
    ;; --- Single Run File
    ;fi = addfile("~/plasimfiles/data/silurian/Silurian450_ensemble/ensembleRun.01/PostProcessing/section.1.nc","r")
    ;oro = fi->lsm
    
    ; --- Output Land Map File for Contour Plots ---
    orolat= fi->lat
    orolon= fi->lon
    
    basename = "./customMap"
    landmask = 1
    paleo_outline (oro(:,:), orolat, orolon, landmask, basename)
    
end
