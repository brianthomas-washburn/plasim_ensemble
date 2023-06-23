;This code reads in a specified variable from a PlaSim
; post-processed file (generated using "burn7")

;This version assumes data is lat-lon-level-time dependent

;Code can find how many ensemble members there are
; and add files from each section within each ensemble member run


;-------------------------------- User specified settings: -----------------------;
;Name to use for this case:
cname ="coupledControl-20pcSNCR_5member_25July2022"
;"test_ControlOnly_March2022"

;Location of data to read in:
indir = "../../runsOutput/"+cname+"/"
runName = "ensembleRun"

;Variable to look at (e.g. ta = temperature); 
var = "ta"

;Plot settings:

;What kind of plot?
; 1 Global average at a specific level, vs time (line plot)
; 2 lat-lon at a specific level and time points; 
plotType = 1
plotLev = 0 ;0 is top level

;For single-time plots, specify which time range (to plot multiple times)
plotTimeStart = 0
plotTimeStop = 0

;Annual averaging?
annAv = "yes"

;Ensemble average?
ensemAv = "yes"


;----------------------------- END User specified settings: -----------------------;


;------------------- Read-in of data ----------------------------------------------;
;determine how many ensemble members:
nRuns = tointeger(systemfunc("ls -d "+indir+runName+".* | wc -l"))
;loop over ensemble members:
do i = 1,nRuns
  if (i .lt. 10) then
    dir_base = ("ensembleRun.0"+i)
  else 
    dir_base = ("ensembleRun."+i)
  end if
  ;add all post-processed files for this member:
  files = addfiles( systemfunc("ls "+indir+dir_base+"/PostProcessing/*.nc") , "r")
  ListSetType (files, "cat")
  filedat = files[:]->$var$
  ;printVarSummary(filedat)

  ;setup input data array, with dimensions number of runs, time, level, lat, lon
  if (i .eq. 1) then
    dsizes = dimsizes(filedat) ;use either dat1 or dat2 here, if they are same dimensions
    dat = new((/nRuns,dsizes(0),dsizes(1),dsizes(2),dsizes(3)/), typeof(filedat))
    delete(dsizes)
    ;read in arrays for time, latitude, longitude, level (pressure):
    time = files[:]->time
    lat = files[0]->lat
    lon = files[0]->lon
    lev = files[0]->lev
  end if
  ;fill run slot with this run's data:
  dat(i-1,:,:,:,:) = filedat(:,:,:,:)
  
end do
dat!0 = "run" ;label first dimension as run number
dat@_FillValue = -9999.99 ;change fill value
printVarSummary(dat)
printVarSummary(time)
;------------------- END Read-in of data -------------------------------------------;

if(ensemAv .eq. "yes") then
;Create ensemble average:
  runAv = dim_avg_n_Wrap(dat,0)
  printVarSummary(runAv)
end if 


if(plotType .eq. 1) then
  ;Average over lat-lon:

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

  datp = wgt_areaave_Wrap(runAv,area,1.0,1)
  datp@_FillValue = -9999.99
  printVarSummary(datp)

  times = time
  tlabel = "Months"
  if(annAv .eq. "yes") then
    datAnnAv = month_to_annual(datp,1)
    datp := datAnnAv
    times :=  ispan(0,-1+dimsizes(time)/12,1)
    tlabel = "Years"
  end if
;  printVarSummary(datp)
;  printVarSummary(times)
  

  ;Plotting:
  axisLabel=var+" lev "+plotLev+" globally averaged "+datp@units
  plotFile := var+"_lev"+plotLev+"_globalAvg-"+cname
  wks   := gsn_open_wks ("png",plotFile )  ; send graphics to PNG file                            
  res                  = True                ; plot mods desired                                  
  res@xyLineThicknessF = 4.0                                                                      
  res@tiXAxisString = "Months"         
  res@tiYAxisString = axisLabel
  res@tiMainString  = " " ;plotFile                    
  res@tiMainFontHeightF = 0.020 ;changes font size       
  ;res@trXMinF = 1200-24    
  ;res@trXMaxF = 1200+24
  
  res@tiMainString     = cname 
  
  plot = gsn_csm_xy(wks,times,datp(:,plotLev),res)
  
end if

