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
  indir1 = "../ensembleRun.01/PostProcessing/"
  in_file1 = "section.1"
  print("in_file1 = "+indir1+in_file1)
  
  indir2 = "../ensembleRun.01/PostProcessing/"
  in_file2 = "section.2"
  print("in directory2 = "+indir2+in_file2)
  
  indir3 = "../ensembleRun.01/PostProcessing/"
  in_file3 = "section.2"
  print("in directory3 = "+indir3+in_file3)
;-- end
;
;Specify what time period to average over.  Setup to average over the *last* tav_yrn years
  tav_yrn = 30 ;e.g. 30 = last 30 years
  tav_tblock = 12 ; 360 for daily data, 12 for monthly data
;-- 
; 
;--- Plotting options: 
; make multiple plots at various times? Currently this makes contour plots
  pmt = "yes";
; **** (Note! Some other plot options must be modified below in plotting section.)
; where to put plot files:
  pfiledir = ""
  system("if ! -d" + pfiledir +" ; then mkdir " + pfiledir +" ; fi"); make pfiledir if one does not exist
; 
;Currently setup to plot from some chosen time to the end of the run
  ;make line plot x-axis go from 0 until yrn years from end:
  yrn = 100
  ;define how many time blocks per year (monthly = 12; daily = 365
  tblock = tav_tblock
  ;deine timestep of plots (e.g. every 30 days, every 1 month, etc.)
  itstep = 1
;  
;-- end
;
;--- Do a running average on time?
  runav = "no" ;
; if yes, what time period to use?
  nave = 1 ;to average over a year, if data is monthly
        ;31 ;to average over ~ month
        ;365 ;to average over year, if data is daily
;-- end
;
;--- Read in data:
;setup files for readin:
  f1 = addfile(indir1+in_file1+".nc","r")
  f2 = addfile(indir2+in_file2+".nc","r")
  f3 = addfile(indir3+in_file3+".nc","r")
;
;read in data of interest (to dat1 variable):
  dat1 = f1->var237 ;customized for ozone levels
  printVarSummary(dat1)
  print("min, max values for "+in_file1+":")
  print(min(dat1))
  print(max(dat1))
  ;print(dat1(0,0,0,:))
  
  dat2 = f2->var237 ;customized for ozone levels
  printVarSummary(dat2)
  print("min, max values for "+in_file2+":")
  print(min(dat2))
  print(max(dat2))
  
  dat3 = f3->var237 ;customized for ozone levels
  printVarSummary(dat3)
  print("min, max values for "+in_file3+":")
  print(min(dat3))
  print(max(dat3))
;
 
;also read in arrays for latitude, longitude, level (pressure) and time:
; use either f1 or f2 here, if they are same dimensions
  lat = f1->lat
;  printVarSummary(lat)
  lon = f1->lon
;  printVarSummary(lon)
  lev = f1->lev
;  printVarSummary(lev)
  time = f1->time
;  printVarSummary(time)
;
;Get the dimensions of the array
  dsizes = dimsizes(dat1) ;use either dat1 or dat2 here, if they are same dimensions
;set the maximum number of time, level, lat, lon points
  tmax   = dsizes(0)
  levmax = dsizes(1)
  latmax = dsizes(2)
  lonmax = dsizes(3)
;
;-- end read in block

;***Average over longitude
  longitudeAv1 = dim_avg_n_Wrap(dat1, 3)
  ;printVarSummary(longitudeAv1)
  longitudeAv2 = dim_avg_n_Wrap(dat2, 3)
  ;printVarSummary(longitudeAv2)
  longitudeAv3 = dim_avg_n_Wrap(dat3, 3)
  ;printVarSummary(longitudeAv3)
;-- end averaging over longitude

; Full ratio of data sets:
  ratio = longitudeAv1/longitudeAv2
  copy_VarMeta(longitudeAv1,ratio)
  ratio@long_name = "ratio_"+in_file1+"-"+in_file2
  ratio@units = "Ratio "+in_file1+"/"+in_file2 
  printVarSummary(ratio)
  print("min, max values for full ratio:")
  print(min(ratio)) 
  print(max(ratio))    
;


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

;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
;-- Contour plots:

;You can specify continour limits (max, min and spacing)
  ; If you want values chosen automatically, use -9999.99 as the limit values
  contMin = 4.0*10^(-7) ;-9999.99
  contMax = 6.5*10^(-6) ;-9999.99
  contSpac = 0.5*10^(-6) ; -9999.99
 
        ;--- tav plot 1 ---
        print("tav plot 1")
        
        plotvar = dim_avg_n_Wrap(longitudeAv1, 0) ;data from file 1
        ;printVarSummary(plotvar)
        
        plotvar@long_name = "";"O~B~3~N~ Concentration Mixing Ratio"
        plotvar@units = "";"kg/kg"
        
        pfiledir = ""
        pfile = "control_custom_"+dat1@long_name+pfiledir+"_tav"
        
        pout = "png" ;"x11" ;where to make plot (x window or file type)
        wks = gsn_open_wks(pout,pfile) ;sends to X11 window or file
        
        res = True
        res@cnFillOn = True     ;filled contours?
        res@cnLinesOn = False   ;contour lines?
        
        if(contMin.ne.-9999.99) then
            res@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
            res@cnMinLevelValF       = contMin          ; set minControl contour level
            res@cnMaxLevelValF       = contMax          ; set max contour level
            res@cnLevelSpacingF      = contSpac         ; set contour spacing
        end if
        res@trYTensionF = 5.0
        res@tiMainString = "Custom Control" ;title to write on plot
        res@tiXAxisString = "Latitude"
        res@tiYAxisString = "Pressure (hPa)"
        
        ;plot = gsn_csm_contour(wks,plotvar,res) ;make the plot
        plot = gsn_csm_pres_hgt(wks, plotvar, res)
        
        ;--- tav plot 2 ---
        print("tav plot 2")
        
        plotvar = dim_avg_n_Wrap(longitudeAv2, 0) ;data from file 2
        ;printVarSummary(plotvar)
        
        plotvar@long_name = "";"O~B~3~N~ Concentration Mixing Ratio"
        plotvar@units = "";"kg/kg"
        
        pfiledir = ""
        pfile = "supernova_custom_"+dat1@long_name+pfiledir+"_tav"
        
        pout = "png" ;"x11" ;where to make plot (x window or file type)
        wks = gsn_open_wks(pout,pfile) ;sends to X11 window or file
        

        res@tiMainString = "Custom Supernova" ;title to write on plot
        
        ;plot = gsn_csm_contour(wks,plotvar,res) ;make the plot
        plot = gsn_csm_pres_hgt(wks, plotvar, res)
 
        ;--- tav plot 3 ---
        print("tav plot 3")
        
        plotvar = dim_avg_n_Wrap(longitudeAv3, 0) ;data from file 3
        ;printVarSummary(plotvar)
        
        plotvar@long_name = "";"O~B~3~N~ Concentration Mixing Ratio"
        plotvar@units = "";"kg/kg"
        
        pfiledir = ""
        pfile = "default_approximation_"+dat1@long_name+pfiledir+"_tav"
        
        pout = "png" ;"x11" ;where to make plot (x window or file type)
        wks = gsn_open_wks(pout,pfile) ;sends to X11 window or file
        
        res = True
        res@cnFillOn = True     ;filled contours?
        res@cnLinesOn = False   ;contour lines?
        if(contMin.ne.-9999.99) then
            res@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
            res@cnMinLevelValF       = contMin          ; set min contour level
            res@cnMaxLevelValF       = contMax          ; set max contour level
            res@cnLevelSpacingF      = contSpac         ; set contour spacing
        end if
        res@trYTensionF = 5.0
        res@tiMainString = "Default Approximation" ;title to write on plot
        
        ;plot = gsn_csm_contour(wks,plotvar,res) ;make the plot
        plot = gsn_csm_pres_hgt(wks, plotvar, res)
    
        ;;--- percent difference custom(1) & default(3) ---
        ;print("percent difference custom(1) & default(3)")
        
        ;plotvar = 100 * dim_avg_n_Wrap((longitudeAv1-longitudeAv3)/longitudeAv3, 0)
        ;plotvar!0 = "lev"
        ;plotvar!1 = "lat"
        ;plotvar&lev = lev
        ;plotvar&lat = lat
        
        ;plotvar@long_name = "";"O~B~3~N~ Concentration Mixing Ratio"
        ;plotvar@units = "";"kg/kg"
        ;;printVarSummary(plotvar)
        
        ;;You can specify continour limits (max, min and spacing)
        ;; If you want values chosen automatically, use -9999.99 as the limit values
        ;contMin = -94 ; min(plotvar)
        ;contMax = 27 ;max(plotvar)
        ;contSpac = 10 ;(contMax-contMin)/10
        ;contCenter = 0.0
        
        ;pfiledir = ""
        ;pfile = "percent_difference_control_custom_default"+dat1@long_name+pfiledir+"_tav"
        
        ;pout = "png" ;"x11" ;where to make plot (x window or file type)
        ;wks = gsn_open_wks(pout,pfile) ;sends to X11 window or file
        
        ;;
        ;; For the "BlRe" colormap, we have 97 colors and the color 
        ;; transitions from blue to red at 49. See :
        ;;
        ;;  http://www.ncl.ucar.edu/Document/Graphics/ColorTables/BlRe.shtml
        ;;
        ;; For the "BlueRed" or "GreenYellow" color maps, we have 253 colors and
        ;; a transition at color 97.
        ;;
        ;gsn_range = gsnColorRange(contMin, contMax, contSpac, contCenter, 97, 49) 
        ;cmap      = read_colormap_file("BlRe")
        
        ;res = True
        ;res@cnFillOn = True     ;filled contours?
        ;;Note: subtract 2 from color start and end since read_colormap_file doesn't return the first two colors
        ;res@cnFillPalette         = cmap(gsn_range@ColorStart-2:gsn_range@ColorEnd-2,:)
        
        ;res@cnLinesOn = False   ;contour lines?
        
        ;if(contMin.ne.-9999.99) then
            ;res@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
            ;res@cnMinLevelValF       = contMin          ; set min contour level
            ;res@cnMaxLevelValF       = contMax          ; set max contour level
            ;res@cnLevelSpacingF      = contSpac         ; set contour spacing
        ;end if
        ;res@trYTensionF = 5.0
        ;res@tiMainString = "Percent Difference Control Custom & Default" ;title to write on plot
       
        
        ;;plot = gsn_csm_contour(wks,plotvar,res) ;make the plot
        ;plot = gsn_csm_pres_hgt(wks, plotvar, res)
        
        ;--- percent difference supernova(2) & custom(1) ---
        print("percent difference supernova(2) & custom(1)")
        
        plotvar = 100 * dim_avg_n_Wrap((longitudeAv2-longitudeAv1)/longitudeAv2, 0) ;file 2 vs file 1
        plotvar!0 = "lev"
        plotvar!1 = "lat"
        plotvar&lev = lev
        plotvar&lat = lat 
        
        plotvar@long_name = "";"O~B~3~N~ Concentration Mixing Ratio"
        plotvar@units = "";"kg/kg"
        ;printVarSummary(plotvar)
        
        ;You can specify continour limits (max, min and spacing)
        ; If you want values chosen automatically, use -9999.99 as the limit values
        contMin = min(plotvar) ;-146 ;-52829 
        contMax = max(plotvar) ;16 
        contSpac = (contMax-contMin)/10
        contCenter = 0.0

        print(contMin) 
        print(contMax)
        print(contSpac)
        print(contCenter)
        
        pfiledir = ""
        pfile = "percent_difference_supernova_control_custom"+dat1@long_name+pfiledir+"_tav"
        
        pout = "png" ;"x11" ;where to make plot (x window or file type)
        wks = gsn_open_wks(pout,pfile) ;sends to X11 window or file
        
        ;
        ; For the "BlRe" colormap, we have 97 colors and the color 
        ; transitions from blue to red at 49. See :
        ;
        ;  http://www.ncl.ucar.edu/Document/Graphics/ColorTables/BlRe.shtml
        ;
        ; For the "BlueRed" or "GreenYellow" color maps, we have 253 colors and
        ; a transition at color 97.
        ;
        gsn_range = gsnColorRange(contMin, contMax, contSpac, contCenter, 97, 49) 
        cmap      = read_colormap_file("BlRe")
        
        res = True
        res@cnFillOn = True     ;filled contours?
        ;Note: subtract 2 from color start and end since read_colormap_file doesn't return the first two colors
        res@cnFillPalette = cmap(gsn_range@ColorStart-2:gsn_range@ColorEnd-2,:)
        
        res@cnLinesOn = False   ;contour lines?
        
        if(contMin.ne.-9999.99) then
            res@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
            res@cnMinLevelValF       = contMin          ; set min contour level
            res@cnMaxLevelValF       = contMax          ; set max contour level
            res@cnLevelSpacingF      = contSpac         ; set contour spacing
        end if
        res@trYTensionF = 5.0
        res@tiMainString = "Percent Difference Supernova & Control Custom" ;title to write on plot
       
        
        ;plot = gsn_csm_contour(wks,plotvar,res) ;make the plot
        plot = gsn_csm_pres_hgt(wks, plotvar, res)
;end contour plots
end
