;*** This NCL code is for analyzing netCDF files produced by PlaSim postprocessor
;    from PlaSim output files.
;    Started by BCThomas, Oct 2019
;    Edited by CLRatterman, Spring 2021
;
; Use needs to specify:
;  1. Where code should look for input files (indir).
;  2. Runs read in (runName: ensembleRun would read ensembleRun.01, ensembleRun.02 , etc.)
;  3. Number of completed runs to analyze (nRuns)

; You can make a series of plots at different time values,
;  which can then be combined into an 
;  animation using "convert" (ImageMagick); eg: convert -delay 50 *.png anim.gif

begin

;--- Set input data file attributes:
; directory:
  indir = "../../../"
  runName = "ensembleRun"
  nRuns = tointeger(systemfunc("ls -d "+indir+runName+".* | wc -l"))
  var = "ts"
  
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

;Get the dimensions of the array
  dsizes = dimsizes(dat)
;set the maximum number of time, level, lat, lon points
  runmax = dsizes(0)
  tmax   = dsizes(1)
  latmax = dsizes(2)
  lonmax = dsizes(3)
  
; For area averaging (from https://www.ncl.ucar.edu/Document/Functions/Built-in/wgt_areaave.shtml):
  jlat  = dimsizes(lat) 
  rad = 4.0*atan(1.0)/180.0 ;pi/180
  re = 6371220.0 ;radius of the Earth in meters
  rr = re*rad
  dlon = abs(lon(2)-lon(1))*rr
  dx = dlon*cos(lat*rad)
; lat can have variable spacing
  dy = new(jlat,typeof(dx))
  dy(0) = abs(lat(2)-lat(1))*rr
  dy(1:jlat-2) = abs(lat(2:jlat-1)-lat(0:jlat-3))*rr*0.5   
  dy(jlat-1) = abs(lat(jlat-1)-lat(jlat-2))*rr
  area = dx*dy   ; cell area function of latitude only

  printVarSummary(area)
  
; Area average (that is, average over latitude and longitude, so a global average)
  areaAv_dat = wgt_areaave_Wrap(dat,area,1.0,1)
  printVarSummary(areaAv_dat)
  
  ; --Line plot of area difference in area averages over time
        nMonths = tmax
        
        plotvar = new((/runmax+1,nMonths/),float)
        
        do i = 0,runmax-1
            plotvar(i,:) =  areaAv_dat(i,:) 
        end do
        
        ;lastly, plot plot run averages for dat1 and dat2 area averages
        plotvar(runmax,:) = dim_avg_n_Wrap(areaAv_dat(:,:) ,0)
        
        printVarSummary(plotvar)
        
        ;convert from months to years
        temp = new((/runmax+1,nMonths/12/),float)
        do i = 0,runmax
            temp(i,:) = month_to_annual(plotvar(i,:), 1)
        end do
        delete(plotvar)
        plotvar = temp
        delete(temp)
        
        printVarSummary(plotvar)
        ;print(plotvar)
        
        years = ispan (1,nMonths/12,1) ; create a years array (/1, 2, ..., 100/)
        printVarSummary(years)
        
        ;Name to save plot as...
        plotvarName = "areaAv_dat"
;;
        ;uses the input file name and variable's long name attribute
        pfile = "./areaAv_overtime" ;custom plot file name
;;
        pout = "png" ;"x11" ;where to make plot (x window or file type)
        wks = gsn_open_wks(pout,pfile) ;open workspace for plot
;;    
        ;set other plot parameters:
        res = True
        
        res@trXMaxF = max(years)
        res@trXMinF = min(years)
        res@trYMaxF = max(plotvar)
        res@trYMinF = min(plotvar)
        
        res@xyLineColors = new(runmax+1,string)  ; change line color
        res@xyLineThicknesses = new(runmax+1,float) ;change line thickness
        res@xyDashPattern = new(runmax+1,integer)   ;change line dash pattern (solid = 0)
        
        do i = 0,runmax-1
            res@xyLineColors(i)  = "gray16"  ; change line color
            res@xyLineThicknesses(i) = 1.0
            res@xyDashPattern(i) = 0
        end do
        res@xyLineColors(runmax)  = "black"  ; change line color
        ;print(res@xyLineColors)
        res@xyLineThicknesses(runmax) = 5.0
        res@xyDashPattern(runmax) = 0
        
        res@tiMainString  = "Annual Global Average Surface Temperatures" ;write a title on the plot
        
        
        res@tiXAxisString = "Year"                              ; xaxis string
        res@tiYAxisString = "Surface Temperature (K)"  ; yaxis string
        
        plot = gsn_csm_xy(wks,years,plotvar,res) ;*** make the plot
        
        gres = True
        gres@YPosPercent = 20    ; expressed as %, 0->100, sets position of top border of legend 
                            ;;  when gres@Position is set to its default setting of "Top" (Default = 95.)
        gres@XPosPercent = 60      ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

        
        lineres = True
        lineres@lgLineColors = (/res@xyLineColors(0),res@xyLineColors(runmax)/) ; line colors
        lineres@lgLineThicknesses = (/res@xyLineThicknesses(0),res@xyLineThicknesses(runmax)/)                       ; line thicknesses
        lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

        textres = True
        textres@lgLabels = (/"individual runs","run average"/)  ; legend labels (required)

        plot = simple_legend(wks,plot,gres,lineres,textres) ;*** make the legend
        draw(plot)
        frame(wks)
end
