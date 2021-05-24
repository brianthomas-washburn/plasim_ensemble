;*** This NCL code is for analyzing netCDF files produced by PlaSim postprocessor
;    from PlaSim output files.
;    Started by BCThomas, Oct 2019
;     for project with Cody Ratterman to look at climate changes associated with O3 variation
;
; Use needs to specify:
;  1. Where code should look for input files (indir).
;  2. Files to read (in_file1 and in_file2)
;  3. Where to put any output (e.g. plots)
;  4. What to plot (what data, what times, what level, etc - see plotting section)
;     Also note you can specify the range of values for plots
;  5. Whether to do a running time average and what timescale to use.


; You can make a series of plots at different time values,
;  which can then be combined into an 
;  animation using "convert" (ImageMagick); eg: convert -delay 50 *.png anim.gif

begin
;--- Read in file/data:
    fi = addfile("../Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc_v2/Map74_PALEOMAP_1deg_Late_Silurian_425Ma.nc","r") ;silurian - 425.2 Ma
    ;fi = addfile("../Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc_v2/Map01_PALEOMAP_1deg_Holocene_0Ma.nc","r") ;modern day - 0 Ma
    
    printVarSummary(fi)
    print(str_join(getfilevarnames(fi),", ")) ;
    
    fidat = fi->z
    printVarSummary(fidat) 
    
;--- Coordinate Transformation Using Area Averaging 
    dat = g2fsh_Wrap (fidat, (/32,64/)) ; T21 Resolution
    ;dat = g2fsh_Wrap (fidat, (/64,128/)) ; T42 Resolution
    ;printVarSummary(dat)
    ;print(dat(0,0))
    ;print(dat(31,0))
    
    dat = dat(::-1,:)  ; Flip the latitude dimension
    
    printVarSummary(dat)
    ;print(dat(0,0))
    ;print(dat(31,0))
    
; --- Surface Geopotential
    standard_g = 9.8 ;standar gravity at mean sea level
    dat_sg = dat(:,:)*standard_g
    copy_VarMeta(dat,dat_sg)
    printVarSummary(dat_sg)
    ;print(min(dat_sg))
    ;print(max(dat_sg))
    
; --- Land Sea Mask
    dat_lsm = new((/32,64/),"double")
    copy_VarMeta(dat,dat_lsm)
    
    do it = 0,31
        do it2 = 0,63
            if (dat(it,it2) .gt. 0)
                dat_lsm(it,it2) = 1.0; land mass
            else
                dat_lsm(it,it2) = 0.0; sea position
            end if
        end do
    end do
    
    printVarSummary(dat_lsm)
    ;print(min(dat_lsm))
    ;print(max(dat_lsm))
    
; --- Sea Ice Cover
    dat_sic = new((/32,64/),"double")
    copy_VarMeta(dat,dat_sic)
    dat_sic(:,:) = 0.0
    printVarSummary(dat_sic)
    
; --- Glaciers
    dat_glac = new((/32,64/),"double")
    copy_VarMeta(dat,dat_glac)
    dat_glac(:,:) = 0.0
    printVarSummary(dat_glac)
    
; --- Vegetation
    dat_vegc = new((/32,64/),"double")
    copy_VarMeta(dat,dat_vegc)
    dat_vegc(:,:) = 0.0
    printVarSummary(dat_vegc)
    
    dat_vegf = new((/32,64/),"double")
    copy_VarMeta(dat,dat_vegf)
    dat_vegf(:,:) = 0.0
    printVarSummary(dat_vegf)
    
; --- Maximum Soil Water Holding (Field) Capacity
    dat_mrfc = new((/32,64/),"double")
    copy_VarMeta(dat,dat_mrfc)
    dat_mrfc(:,:) = 0.2 ;<---- set the capacity here (m) 
    printVarSummary(dat_mrfc)
    
; --- Surface Albedo
    as_sea = 0.06
    as_land = 0.24
    dat_as = new((/32,64/),"double")
    copy_VarMeta(dat,dat_as)
    dat_as(:,:) = dat_lsm(:,:)*(as_land-as_sea)+as_sea ; lsm is 0 or 1, all values are set to as_sea or as_land
    printVarSummary(dat_as)
    
    as0 = 0.0
    dat_as0 = new((/32,64/),"double")
    copy_VarMeta(dat,dat_as0)
    dat_as0(:,:) = 0.0
    printVarSummary(dat_as0)
    
; --- Make Contour Plots    
    res = True
    res@cnFillOn = True     ;filled contours?
    res@cnLinesOn = False   ;contour lines?
    ;res@lbOrientation = "Vertical"

    res@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
    res@cnMinLevelValF       = min(fidat)          ; set min contour level
    res@cnMaxLevelValF       = max(fidat)          ; set max contour level
    res@cnLevelSpacingF      = (max(fidat)-min(fidat))/20         ; set contour spacing
;      
    ;-- File Data
    res@tiMainString = "File Data" ;title to write on plot
    wks = gsn_open_wks("png","fidat") ;sends to X11 window or file
    plot = gsn_csm_contour_map(wks,lonFlip(fidat(1:180,1:360)),res) ;make the plot
;    
    ;-- Coordinate Transformation Data
    res@tiMainString = "New Grid Data" ;title to write on plot
    wks = gsn_open_wks("png","dat") ;sends to X11 window or file
    plot = gsn_csm_contour_map(wks,dat,res) ;make the plot
    
    ;-- Land Sea Mask Data
    res@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
    res@cnMinLevelValF       = 0          ; set min contour level
    res@cnMaxLevelValF       = 1          ; set max contour level
    res@cnLevelSpacingF      = 0.5         ; set contour spacing
    
    res@tiMainString = "Land Sea Mask Data" ;title to write on plot
    wks = gsn_open_wks("png","dat_lsm") ;sends to X11 window or file
    plot = gsn_csm_contour_map(wks,dat_lsm,res) ;make the plot
;

; --- Output to Surface Files
   dsizes = dimsizes(dat) ;use either dat1 or dat2 here, if they are same dimensions
   nlat = dsizes(0)
   nlon = dsizes(1)
   
   lines = new(257,"string")
   
   ; --- surface geopotential
   dat_sg_string = new((/257,8/),"string")
   
   dat_sg_string(0,0) = "       129"
   dat_sg_string(0,1) = "         0"
   dat_sg_string(0,2) = "  20070101"
   dat_sg_string(0,3) = "        64"
   dat_sg_string(0,4) = "         0"
   dat_sg_string(0,5) = "        32"
   dat_sg_string(0,6) = "         0"
   dat_sg_string(0,7) = "         0"
   
   do i = 0, 31
      do j = 1, 8
         dat_sg_string(8*i+j,0:7) = tostring_with_format(dat_sg(i, 8*(j-1):(8*j)-1), "%10.3f")
      end do
      ;dat_sg_string(1,0:7) = tostring_with_format(dat_sg(0, 0: 7), "%10.3f")
      ;dat_sg_string(2,0:7) = tostring_with_format(dat_sg(0, 8:15), "%10.3f")
      ;dat_sg_string(3,0:7) = tostring_with_format(dat_sg(0,16:23), "%10.3f")
      ;dat_sg_string(4,0:7) = tostring_with_format(dat_sg(0,24:31), "%10.3f")
      ;dat_sg_string(5,0:7) = tostring_with_format(dat_sg(0,32:39), "%10.3f")
      ;dat_sg_string(6,0:7) = tostring_with_format(dat_sg(0,40:47), "%10.3f")
      ;dat_sg_string(7,0:7) = tostring_with_format(dat_sg(0,48:55), "%10.3f")
      ;dat_sg_string(8,0:7) = tostring_with_format(dat_sg(0,56:63), "%10.3f")
   end do
   
   do i = 0, 256
      lines(i) = str_join(dat_sg_string(i,:), "")
   end do
   
   asciiwrite ("N032_surf_0129.sra", lines)
   
   ; --- land sea mask
   dat_lsm_string = new((/257,8/),"string")
   lines(:) = ""
   
   dat_lsm_string(0,0) = "       172"
   dat_lsm_string(0,1) = "         0"
   dat_lsm_string(0,2) = "  20090101"
   dat_lsm_string(0,3) = "        64"
   dat_lsm_string(0,4) = "         0"
   dat_lsm_string(0,5) = "        32"
   dat_lsm_string(0,6) = "         0"
   dat_lsm_string(0,7) = "         0"
   
   do i = 0, 31
      do j = 1, 8
         dat_lsm_string(8*i+j,0:7) = tostring_with_format(dat_lsm(i, 8*(j-1):(8*j)-1), "%10.3f")
      end do
   end do
   
   do i = 0, 256
      lines(i) = str_join(dat_lsm_string(i,:), "")
   end do
   
   asciiwrite ("N032_surf_0172.sra", lines)
   
   ; --- sea ice cover
   dat_sic_string = new((/257,8/),"string")
   lines(:) = ""
   
   dat_sic_string(0,0) = "       210"
   dat_sic_string(0,1) = "         0"
   dat_sic_string(0,2) = "  20090000"
   dat_sic_string(0,3) = "        -1"
   dat_sic_string(0,4) = "        64"
   dat_sic_string(0,5) = "        32"
   dat_sic_string(0,6) = "         0"
   dat_sic_string(0,7) = "         0"
   
   do i = 0, 31
      do j = 1, 8
         dat_sic_string(8*i+j,0:7) = tostring_with_format(dat_sic(i, 8*(j-1):(8*j)-1), "%7.3f")
      end do
   end do
   
   do i = 0, 256
      lines(i) = str_join(dat_sic_string(i,:), "")
   end do
   
   asciiwrite ("N032_surf_0210.sra", lines)
   
   ; --- glaciers
   delete(lines) ; resize lines array because there are now 4 data points per-line instead of 8
   lines = new(513,"string")
   dat_glac_string = new((/513,4/),"string")
   lines(:) = ""
   
   dat_glac_string(0,0) = "       232         0"
   dat_glac_string(0,1) = "    990100        -1"
   dat_glac_string(0,2) = "        64        32"
   dat_glac_string(0,3) = "         0         0"
   
   do i = 0, 31
      do j = 1, 16
        dat_glac_string(16*i+j,0:3) = tostring_with_format(dat_glac(i, 4*(j-1):(4*j)-1), "%10.3f")
      end do
   end do
   
   do i = 0, 512
      lines(i) = str_join(dat_glac_string(i,:), "")
   end do
   
   asciiwrite ("N032_surf_0232.sra", lines)
   
   ; --- vegetation (c & f)                                 ERROR HERE!!!
   dat_vegc_string = new((/513,4/),"string")
   lines(:) = ""
   
   dat_vegc_string(0,0) = "       199         0"
   dat_vegc_string(0,1) = "  20090001        -1"
   dat_vegc_string(0,2) = "        64        32"
   dat_vegc_string(0,3) = "         0         0"
   
   do i = 0, 31
      do j = 1, 16
        dat_vegc_string(16*i+j,0:3) = tostring_with_format(dat_vegc(i, 4*(j-1):(4*j)-1), "%16.6f")
      end do
   end do
   
   do i = 0, 512
      lines(i) = str_join(dat_vegc_string(i,:), "")
   end do
   
   asciiwrite ("N032_surf_0199.sra", lines)
   ; ---
   dat_vegf_string = new((/513,4/),"string")
   lines(:) = ""
   
   dat_vegf_string(0,0) = "       212         0"
   dat_vegf_string(0,1) = "    990100        -1"
   dat_vegf_string(0,2) = "        64        32"
   dat_vegf_string(0,3) = "         0         0"
   
    do i = 0, 31
      do j = 1, 16
        dat_vegf_string(16*i+j,0:3) = tostring_with_format(dat_vegf(i, 4*(j-1):(4*j)-1), "%16.6f")
      end do
   end do
   
   do i = 0, 512
      lines(i) = str_join(dat_vegf_string(i,:), "")
   end do
   
   asciiwrite ("N032_surf_0212.sra", lines)
   
   ; --- Surface Albedo
   dat_as_string = new((/513,4/),"string")
   lines(:) = ""
   
   dat_as_string(0,0) = "       174         0"
   dat_as_string(0,1) = "  20090001         0"
   dat_as_string(0,2) = "        64        32"
   dat_as_string(0,3) = "         0         0"
   
   do i = 0, 31
      do j = 1, 16
        dat_as_string(16*i+j,0:3) = tostring_with_format(dat_as(i, 4*(j-1):(4*j)-1), "%10.3f")
      end do
   end do
   
   do i = 0, 512
      lines(i) = str_join(dat_as_string(i,:), "")
   end do
   
   asciiwrite ("N032_surf_0174.sra", lines)
   
   ; ---
   dat_as0_string = new((/513,4/),"string")
   lines(:) = ""
   
   dat_as0_string(0,0) = "      1740         0"
   dat_as0_string(0,1) = "  20090101        -1"
   dat_as0_string(0,2) = "        64        32"
   dat_as0_string(0,3) = "         0         0"
   
   do i = 0, 31
      do j = 1, 16
        dat_as0_string(16*i+j,0:3) = tostring_with_format(dat_as0(i, 4*(j-1):(4*j)-1), "%10.3f")
      end do
   end do
   
   do i = 0, 512
      lines(i) = str_join(dat_as0_string(i,:), "")
   end do
   
   asciiwrite ("N032_surf_1740.sra", lines)
   
   ; --- Contour Map Outline
    oro = dat_lsm
    orolat = fi->lat
    orolon = fi->lon
    
    basename = "./silurianMap"
    landmask = 1
    paleo_outline (oro(:,:), orolat, orolon, landmask, basename)
   
end
