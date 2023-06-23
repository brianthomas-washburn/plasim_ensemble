;This code to handle plotting pretty generic input from PlaSim
; Assumes any data input is *time dependent* - then asks about
;  space dimentsionality

;Specify how many ensemble members to include:
;ensemN = 5

;How many sections for each ensemble?
;secN = 2





;Set up input array:
; First dimension is for files (data from file 1, 2, etc)
; Next dimensions are for data from each file
;latDimIn = 

;What dimensionality are we looking at (besides time): lat, lon, alt?
;dimIn = 3 ;3 for lat+lon+alt; 2 for just lat+lon (ie. single alt data)

;dimIn = dimIn+1 ;add in time dimension



;dataIn_global =


  indir = "../../runsOutput/coupledControl-20pcSNCR_5member_25July2022/"
  runName = "ensembleRun"
  nRuns = tointeger(systemfunc("ls -d "+indir+runName+".* | wc -l"))
  
  control_equilibriumStart = 0
  supernova_equilibriumStart = 1
  nYears = 1
  
  var = "ta"
  
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