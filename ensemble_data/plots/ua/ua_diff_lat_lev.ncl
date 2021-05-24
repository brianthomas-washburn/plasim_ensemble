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

;--- Set input data file attributes:
; directory:
  indir = "../../"
  runName = "ensembleRun"
  nRuns = tointeger(systemfunc("ls -d "+indir+runName+".* | wc -l"))
  
  control_equilibriumStart = 0
  supernova_equilibriumStart = 1
  nYears = 1
  
  var = "ua"
  
  do i = 1,nRuns
    if (i .lt. 10) then
        dir_base = ("ensembleRun.0"+i)
    else 
        dir_base = ("ensembleRun."+i)
    end if
    files = addfiles( systemfunc("ls "+indir+dir_base+"/PostProcessing/*.nc") , "r")
    ListSetType (files, "cat")
    filedat = files[:]->$var$
    ;printVarSummary(filedat)
  
    if (i .eq. 1) then
        dsizes = dimsizes(filedat) ;use either dat1 or dat2 here, if they are same dimensions
        dat = new((/nRuns,dsizes(0),dsizes(1),dsizes(2),dsizes(3)/), typeof(filedat))
        delete(dsizes)
      ;read in arrays for latitude, longitude, level (pressure) and time:
        time = files[:]->time
        lat = files[0]->lat
        lon = files[0]->lon
        lev = files[0]->lev
    end if
    
    dat(i-1,:,:,:,:) = filedat(:,:,:,:)
    
  end do
  dat!0 = "run"
  
  printVarSummary(dat)
  print(min(dat))
  print(max(dat))
  
  runAv = dim_avg_n_Wrap(dat,0)
  printVarSummary(runAv)
  
  dat1 = runAv(12*control_equilibriumStart:12*control_equilibriumStart+12*nYears-1,:,:,:)
  dat2 = runAv(12*supernova_equilibriumStart:12*supernova_equilibriumStart+12*nYears-1,:,:,:)
  
  diff = dat2-dat1 ;(dat2-dat1)/dat1 ;dat2-dat1
  copy_VarMeta(runAv(12*control_equilibriumStart:12*control_equilibriumStart+12*nYears-1,:,:,:), diff(:,:,:,:))
  printVarSummary(diff)
  print(min(diff))
  print(max(diff))
  
  dsizes = dimsizes(diff)
  
  diff_jan = new((/dsizes(0)/12,dsizes(1),dsizes(2)/), typeof(diff))
  diff_april = new((/dsizes(0)/12,dsizes(1),dsizes(2)/), typeof(diff))
  diff_july = new((/dsizes(0)/12,dsizes(1),dsizes(2)/), typeof(diff))
  diff_oct = new((/dsizes(0)/12,dsizes(1),dsizes(2)/), typeof(diff))
  
  do i = 0,dsizes(0)/12-1
    diff_jan(i,:,:) = dim_avg_n_Wrap(diff(12*i,:,:,:), 2)
    diff_april(i,:,:) = dim_avg_n_Wrap(diff(12*i+3,:,:,:), 2)
    diff_july(i,:,:) = dim_avg_n_Wrap(diff(12*i+6,:,:,:), 2)
    diff_oct(i,:,:) = dim_avg_n_Wrap(diff(12*i+9,:,:,:), 2)
    i = i + 12
  end do
  
  printVarSummary(diff_jan)


;----------------------------------------------------------------------
undef("gsnColorRange")
function gsnColorRange(lower:numeric, upper:numeric, step:numeric, \
                       center:numeric, color_end:integer, center_color:integer)
local nboxes_left, nboxes_right, ncolors_left, ncolors_right, output, \
      color_start
begin

 color_start = 2 ; All of the color maps begin with 2.
;
; Calculate the number of color bar boxes to the left and right of
; the "center" value.
;
 nboxes_left  = (center - (lower - step))/step
 nboxes_right = ((upper + step) - center)/step

;
; Calculate the number of colors in the map on the left and right hand sides.
;
 ncolors_left = (center_color - color_start + 1)
 ncolors_right = (color_end - center_color)

 output = True

;
;
; Either the lower or upper extent of the color map will be adjusted. If
; the magnitude of the lower limit is less than the magnitude of the
; upper limit, then the lower limit has to be "moved in" towards the
; center color. Oppositely, the upper limit will be moved. If both the
; lower and upper numerical values are the same, then pass back 2 as the
; lower extent and color_end (number of colors) as the upper extent (use
; the whole map in other words).
;
 if(abs(lower) .lt. abs(upper))then

  output@ColorStart = round(center_color - (ncolors_right/nboxes_right)*nboxes_left, 3)
  output@ColorEnd   = color_end

 else

;---no "else if" in NCL :(
  if(abs(lower) .gt. abs(upper))then
   output@ColorStart = 2
   output@ColorEnd   = round(center_color + (ncolors_left/nboxes_left)*nboxes_right, 3)
  else
   output@ColorStart = 2
   output@ColorEnd   = color_end
  end if

 end if

 return(output)

end
;

  
; --Contour plot of averages

; ---January
   ;What variable to plot:    
    plotvar = dim_avg_n_Wrap(diff_jan , 0)
    ;
    printVarSummary(plotvar)
    print(min(plotvar))
    print(max(plotvar))
    ;set the plot output file name:
    ;uses the input file name and variable's long name attribute
    pfile = "ua_diff_lat_lev_contour"
    pout = "png" ;"x11" ;where to make plot (x window or file type)
    wks = gsn_open_wks(pout,pfile) ;sends to X11 window or file
    ;
    ;set plot attributes:
    res = True
    res@cnFillOn = True     ;filled contours?
    res@cnLinesOn = False   ;contour lines?
    ;
    contMin = floor(min(diff)); 100
    contMax = ceil(max(diff)); -100
    contSpac =  0.5 ; 10
    contCenter = 0
    ;
    res@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
    res@cnMinLevelValF       = contMin          ; set min contour level
    res@cnMaxLevelValF       = contMax          ; set max contour level
    res@cnLevelSpacingF      = contSpac         ; set contour spacing
    ;
    gsn_range = gsnColorRange(contMin, contMax, contSpac, contCenter, 97, 49) 
    cmap      = read_colormap_file("BlRe")
    res@cnFillPalette = cmap(gsn_range@ColorStart-2:gsn_range@ColorEnd-2,:)
    ;
    res@trYTensionF = 5.0
    ;
    res@tiMainString = "January" ;title to write on plot
    plotvar@long_name = "Eastward Wind Difference"
    plotvar@units = "m/s"
    ;      
    res@tiXAxisString = "Latitude"
    res@tiYAxisString = "Pressure"
    res@trYReverse = True ; reverse y-axis  
    ;      
    plot = gsn_csm_pres_hgt(wks, plotvar, res) ;make the plot
; 

; ---April 
    plotvar = dim_avg_n_Wrap(diff_april , 0)
    print(min(plotvar))
    print(max(plotvar))
    res@tiMainString = "April" ;title to write on plot
    plotvar@long_name = "Eastward Wind Difference"
    plotvar@units = "m/s"
    plot = gsn_csm_pres_hgt(wks, plotvar, res) ;make the plot
    
; ---July 
    plotvar = dim_avg_n_Wrap(diff_july , 0)
    print(min(plotvar))
    print(max(plotvar))
    res@tiMainString = "July" ;title to write on plot
    plotvar@long_name = "Eastward Wind Difference"
    plotvar@units = "m/s"
    plot = gsn_csm_pres_hgt(wks, plotvar, res) ;make the plot
    
; ---October
    plotvar = dim_avg_n_Wrap(diff_oct , 0)
    print(min(plotvar))
    print(max(plotvar))
    res@tiMainString = "October" ;title to write on plot
    plotvar@long_name = "Eastward Wind Difference"
    plotvar@units = "m/s"
    plot = gsn_csm_pres_hgt(wks, plotvar, res) ;make the plot

;end contour plots

end
