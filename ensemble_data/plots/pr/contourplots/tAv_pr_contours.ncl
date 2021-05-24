;*** This NCL code is for analyzing netCDF files produced by PlaSim postprocessor
;    from PlaSim output files.
;    Started by BCThomas, Oct 2019
;    Edited by CLRatterman, Spring 2021
;
; Use needs to specify:
;  1. Where code should look for input files (indir).
;  2. Runs read in (runName: ensembleRun would read ensembleRun.01, ensembleRun.02 , etc.)
;  3. Number of completed runs to analyze (nRuns)
;  4. Start and length of equilibrium periods (control_equilibriumStart, supernova_equilibriumStart, nYears)

begin

;--- Set input data file attributes:
; directory:
  indir = "../../../"
  runName = "ensembleRun"
  nRuns = tointeger(systemfunc("ls -d "+indir+runName+".* | wc -l"))
  
  control_equilibriumStart = 0
  supernova_equilibriumStart = 1
  nYears = 1
  
  var = "pr"
  
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
        dat = new((/nRuns,dsizes(0),dsizes(1),dsizes(2)/), typeof(filedat))
        delete(dsizes)
      ;read in arrays for latitude, longitude, level (pressure) and time:
        time = files[:]->time
        lat = files[0]->lat
        lon = files[0]->lon
        lev = files[0]->lev
    end if
    
    dat(i-1,:,:,:) = filedat(:,:,:)
    
  end do
  dat!0 = "run"
  
  printVarSummary(dat)
  print(min(dat))
  print(max(dat))
  
  runAv = dim_avg_n_Wrap(dat,0)
  printVarSummary(runAv)
  
  dat1 = dim_avg_n_Wrap(runAv(12*control_equilibriumStart:12*control_equilibriumStart+12*nYears-1,:,:), 0)
  dat2 = dim_avg_n_Wrap(runAv(12*supernova_equilibriumStart:12*supernova_equilibriumStart+12*nYears-1,:,:), 0)
  
  diff = dat2-dat1
  copy_VarMeta(runAv(0,:,:), diff(:,:))
  printVarSummary(diff)
  print(min(diff))
  print(max(diff))
; --- 
;*** Plotting:
;
;----------------------------------------------------------------------
; This function sets the gsnSpreadXXX resources necessary to correctly
; span a two-color colortable.
;
; lower  = cnMinLevelValF
; upper  = cnMaxLevelValF
; step   = cnLevelSpacingF
; center = The numerical value the colormap is centered on. For
;          anomalies or trends, it's common to use 0.0, so blue means
;          cold or cooling and red means warm or warming.
; color_end = The number of colors in colormap (ex. 97 for BlRe, 253 
;             for BlueRed)
; center_color = Color value on the left of the "center" value
; (see above).
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
; --- Diff Contour Plot ---
  plotvar = diff                                                  
      plotvar@long_name = "                                         Precipitation Difference (m/s)"
      plotvar@units = ""
      
  pfile = "tav_diff_contour"
      pout = "png" ;"x11" ;where to make plot (x window or file type)
      wks = gsn_open_wks(pout,pfile) ;sends to X11 window or file
      
  res = True
      res@cnFillOn = True     ;filled contours?
      res@cnLinesOn = False   ;contour lines?
      
  contMin = min(plotvar)
  contMax = max(plotvar)
  contSpac = 1e-8
  contCenter = 0.0
      
  gsn_range = gsnColorRange(contMin, contMax, contSpac, contCenter, 97, 49) 
  cmap      = read_colormap_file("BlRe")
  res@cnFillPalette         = cmap(gsn_range@ColorStart-2:gsn_range@ColorEnd-2,:)
      
  res@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
        res@cnMinLevelValF       = contMin          ; set min contour level
        res@cnMaxLevelValF       = contMax          ; set max contour level
        res@cnLevelSpacingF      = contSpac         ; set contour spacing
        
  res@tiMainString = "SN and Custom Control Sections";title to write on plot
  
  res@mpDataBaseVersion = "Ncarg4_1"
  res@mpDataSetName = "../../00_map/customMap"
  res@mpGeophysicalLineThicknessF = 5.0
  
  plot = gsn_csm_contour_map_ce(wks,plotvar,res) ;make the plot

; -- Control Contour Plot
plotvar = dat1                                                  
      plotvar@long_name = "                                              Precipitation (m/s)"
      plotvar@units = ""
      
  pfile = "tav_control_contour"
      pout = "png" ;"x11" ;where to make plot (x window or file type)
      wks = gsn_open_wks(pout,pfile) ;sends to X11 window or file
      
  res = True
      res@cnFillOn = True     ;filled contours?
      res@cnLinesOn = False   ;contour lines?
      
   if (min(dat1) .lt. min(dat2)) then
    contMin = min(dat1)
  else
    contMin = min(dat2)
  end if
  
  if (max(dat1) .gt. max(dat2)) then
    contMax = max(dat1)
  else
    contMax = max(dat2)
  end if
  
  contSpac = 1e-8
   
  ;res@cnLevelSelectionMode = "AutomaticLevels"
  res@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
        res@cnMinLevelValF       = contMin          ; set min contour level
        res@cnMaxLevelValF       = contMax          ; set max contour level
        res@cnLevelSpacingF      = contSpac         ; set contour spacing
        
  res@tiMainString = "Custom Control Case";title to write on plot
  
  res@mpDataBaseVersion = "Ncarg4_1"
  res@mpDataSetName = "../../00_map/customMap"
  res@mpGeophysicalLineThicknessF = 5.0
  
  plot = gsn_csm_contour_map_ce(wks,plotvar,res) ;make the plot
  
; -- Supernova Contour Plot
plotvar = dat2                                                  
      plotvar@long_name = "                                              Precipitation (m/s)"
      plotvar@units = ""
      
  pfile = "tav_supernova_contour"
      pout = "png" ;"x11" ;where to make plot (x window or file type)
      wks = gsn_open_wks(pout,pfile) ;sends to X11 window or file
      
  res = True
      res@cnFillOn = True     ;filled contours?
      res@cnLinesOn = False   ;contour lines?

  contMin = contMin
  contMax = contMax
  contSpace = contSpac
   
  ;res@cnLevelSelectionMode = "AutomaticLevels"
  res@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
        res@cnMinLevelValF       = contMin          ; set min contour level
        res@cnMaxLevelValF       = contMax          ; set max contour level
        res@cnLevelSpacingF      = contSpac         ; set contour spacing
        
  res@tiMainString = "Supernova Case";title to write on plot
  
  res@mpDataBaseVersion = "Ncarg4_1"
  res@mpDataSetName = "../../00_map/customMap"
  res@mpGeophysicalLineThicknessF = 5.0
  
  plot = gsn_csm_contour_map_ce(wks,plotvar,res) ;make the plot
  
end
