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
  
  ; --Histogram of monthly global average surface temperatures
        plotvar = new((/2,12*runmax*nYears/),float)
        ; convert from nRuns and time dimensioned to two separate one dimensional arrays
        plotvar(0,:) =  ndtooned(areaAv_dat(:,12*control_equilibriumStart:12*control_equilibriumStart+12*nYears-1))
        plotvar(1,:) =  ndtooned(areaAv_dat(:,12*supernova_equilibriumStart:12*supernova_equilibriumStart+12*nYears-1))
        
        printVarSummary(plotvar)
        
        ;uses the input file name and variable's long name attribute
        pfile = "./areaAv_hist_sections" ;custom plot file name

        pout = "png" ;"x11" ;where to make plot (x window or file type)
        wks = gsn_open_wks(pout,pfile) ;open workspace for plot
   
        ;set other plot parameters:
        res = True
        
        ; Set histogram bin and label parameters here
        xMin = floor(min(plotvar)) ; xaxis minimum
        xMax = ceil(max(plotvar)) ; xaxis maximum
        binWidth = 0.1 ; width of bins
        
        nBins = tointeger((xMax-xMin)*(1/binWidth)) ; number of bins
        ;print(nBins)
        labelStride = 1/binWidth ; how many ticks to skip between x-axis labels
        
        res@gsnHistogramBinIntervals = fspan(xMin, xMax, nBins+1) 
        ;print(res@gsnHistogramBinIntervals)
        res@tmXBLabelStride = labelStride
        res@tmXBMajorLengthF = 0.01
        
        res@tiMainString  = "Monthly Global Average Surface Temperatures" ; write a title on the plot
        
        res@tiXAxisString = "Surface Temperature (K)" ; xaxis string
        res@tiYAxisString = "Frequency" ; yaxis string
        
        res@gsnHistogramBarColors = "blue" ; control bar color
        res@gsFillOpacityF = 0.8 ; control bar opacity (see-throughness, 0-1)
        plot1 = gsn_histogram(wks,plotvar(0,:),res) ;*** make the plot
        
        res@gsnHistogramBarColors = "red" ; supernova bar color
        res@gsFillOpacityF = 0.4 ; supernova bar opacity (see-throughness, 0-1)
        plot2 = gsn_histogram(wks,plotvar(1,:),res) ;*** make the plot
        
        overlay(plot1,plot2) ; plot both histograms together
        draw(plot1)
        frame(wks)
        
    ; -- T-test for significant difference between control and supernova cases
        siglvl  = 0.01 ; significance level test
        aveC = dim_avg(plotvar(0,:))
        aveS = dim_avg(plotvar(1,:))
        print("averages = " + aveC + ", " + aveS)
        varC = variance(plotvar(0,:))
        varS = variance(plotvar(1,:))
        print("variances = " + varC + ", " + varS)
        print("standard deviations = " + sqrt(varC) + ", " + sqrt(varS))
        sC = dimsizes(plotvar(0,:))
        sS = dimsizes(plotvar(1,:))
        print("sample sizes = " + sC + ", " + sS)
        
        iflag = False ; false when we can assume the two populations have the same variance
        prob = ttest(aveC,varC,sC, aveS,varS,sS, iflag, False)
        
        print("null hypothesis prob (0 to 1) = " + sprintf("%.10f", prob))
        
        print("significance level = " + siglvl)
        if (prob .lt. siglvl) then
            print("There is a significant difference between the two sample averages.")
        else
            print("There is likely no significant difference between the two sample averages.")
        end if
end
