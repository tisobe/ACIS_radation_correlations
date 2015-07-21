;+++++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION NO_AXIS_LABELS, axis, index, value
;+++++++++++++++++++++++++++++++++++++++++++++++++++
;
return, string(" ")
end
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++
FUNCTION MTA_RAD, startpar, endpar, PLOTX=plotx
;+++++++++++++++++++++++++++++++++++++++++++++++++++
;
;---  can supply startpar and endpar in format yyyy:doy
;---  to run subset of data

;--- this is updated version of eph_cti.pro
;--- Ver 2.0 Feb 2002 BDS
;
;---   updated w/ mrdfits,rdfloat, better timing plots 
;---   fixed eph limits (they were originally speced wrong
;
;--- Ver 2.1 Feb 2002 BDS
;---   use ephemeris file instead of guessing when perigee is
;---         Mar 2002 BDS
;---   add crm region to config plot
;
;--- May 2010
;---  adjust for new ephin channels E1300 and E150
;
;--- Mar 2015
;--- automated config file array / modified ylog option / added skipping plot for 0 data cases
;---- modified by t. isobe (tisobe@cfa.harvard.edu)
;
; ********** param start *************************
; ************ here's a few things that are easily changed ******
;--- set operational limits

e1300_lim = 1000.0/65.6/1.17
e150_lim  = 800000.0/65.6/1.17
p_4gm_lim = 300
p41gm_lim = 8.47
;
;--- config_files are split into years, so only read the ones needed
;
jtime = systime(/julian, /utc)
caldat, jtime, cyear, m, d, ht, mt, st

tyear         = cyear - 1999
config_files  = strarr(tyear+1)

for i = 0, tyear do begin
    pyear = 1999 + i
    syear = string(pyear)
    syear = strtrim(syear, 1)

    config_files(i) = '/data/mta_www/mta_temp/mta_states/MJ/' + syear + '/comprehensive_data_summary' + syear

;    print, config_files(i)                 ;--- debug
endfor
;
; ********** param end *************************
;
;
; ***** data files *******************
;
ephfile     = '/data/mta4/www/DAILY/mta_rad/Data/mk_eph_avg_2010.out'
ephscafile  = '/data/mta4/www/DAILY/mta_rad/Data/ephsca.fits'
ctifile     = '/data/mta4/www/DAILY/mta_rad/cti_data.txt'
goesfile    = '/data/mta4/www/DAILY/mta_rad/goes_data.txt'
acefile     = '/data/mta4/www/DAILY/mta_rad/ace_data.txt'
cnt_file    = strarr(3)
cnt_file[0] = '/data/mta/www/mta_dose_count/*/ccd5'
cnt_file[1] = '/data/mta/www/mta_dose_count/*/ccd6'
cnt_file[2] = '/data/mta/www/mta_dose_count/*/ccd7' 
;/*/    ;---- this is here just for convenience for editing.
ephemfile   = '/data/mta/DataSeeker/data/repository/dephem.rdb'
;xmmfile    = '/data/mta4/space_weather/XMM/xmm.archive'

; ********** param end ***************************

exit_status = [1,1]
;
;--- set starting and ending dates of interest
;
if (n_params() eq 0) then begin
  startpar = string('1999:203')
  endpar   = string(strmid(systime(), 20, 4) + ':366')
endif

if (n_params() eq 1) then $
  endpar = string(strmid(systime(), 20, 4) + ':366')

ttmp      = fltarr(2)
ttmp      = strsplit(startpar,':', /extract)
startpars = float(((ttmp(0) - 1998) * 31536000) + ((ttmp(1)-1) * 86400))
leap      = 2000

while (ttmp(0) gt leap) do begin
  startpars = startpars + 86400
  leap      = leap + 4
endwhile

ttmp    = fltarr(2)
ttmp    = strsplit(endpar,':', /extract)
endpars = float(((ttmp(0) - 1998) * 31536000) + ((ttmp(1)-1) * 86400))
leap    = 2000

while (ttmp(0) gt leap) do begin
  endpars = endpars + 86400
  leap    = leap + 4
endwhile
;   
; ********** ACE  **************
;
rdfloat, acefile, yy,mm,dd,hm,jd,sec, $
         x,x,x,x,ace_ch1,ace_ch2,ace_ch3,ace_ch4,ace_ch5,x

ace_time = float((jd - 50814) * 86400) + sec

b = where(ace_time ge startpars and ace_time le endpars, j)

if (j gt 0) then begin
  ace_time = ace_time(b)
  ace_ch1  = ace_ch1(b)
  ace_ch2  = ace_ch2(b)
  ace_ch3  = ace_ch3(b)
  ace_ch4  = ace_ch4(b)
  ace_ch5  = ace_ch5(b)
endif else begin
  ace_time = [0, 0]
  ace_ch1  = [0, 0]
  ace_ch2  = [0, 0]
  ace_ch3  = [0, 0]
  ace_ch4  = [0, 0]
  ace_ch5  = [0, 0]
endelse
;
; ********** GOES **************
;
rdfloat, goesfile, yy,mm,dd,hm,jd,sec, $
         goes_p1,goes_p2,x,x,goes_p5,x,x,x,x,x,x

goes_time = float((jd - 50814) * 86400) + sec

b = where(goes_time ge startpars and goes_time le endpars, j)

if (j gt 0) then begin
  goes_time = goes_time(b)
  goes_p1   = goes_p1(b)
  goes_p2   = goes_p2(b)
  goes_p5   = goes_p5(b)
endif else begin
  goes_time = [0, 0]
  goes_p1   = [0, 0]
  goes_p2   = [0, 0]
  goes_p5   = [0, 0]
endelse
;
; ********** EPHIN **************
;
rdfloat, ephfile, time, e1300, e150, p_4gm, p41gm, scint,skipline=2

b = where(time ge startpars and time le endpars, done)

print, startpars, endpars
print, min(time)
print, max(time)
print, startpars
print, endpars

if (done eq 0) then begin
;
;--- no files to close now, but add them here if that changes
;---stop, 'No observations found between ', startpar, ' and ', endpar
;
  print, 'No observations found between ', startpar, ' and ', endpar
endif

if (done gt 0) then begin       ;---- this if goes all the way to the end of the script
time  = time(b)

e1300 = e1300(b)
e150  = e150(b)
;
;--- sca00 is not updated any more
;
ephsca = mrdfits(ephscafile, 1)

print,"Hey you ",max(ephsca.time)

b = where(ephsca.time ge startpars and ephsca.time le endpars,j)

if (j gt 0) then begin
  ephscatime = ephsca(b).time
  ephsca00   = ephsca(b).sca00_avg
endif else begin
  ephscatime = [0,0]
  ephsca00   = [0,0]
endelse
;
; ********** CTI *****************
;
readcol, ctifile, cti_start, x,x,x,x,cti_obsid,cti_end,x,x, $
         format='A,A,A,A,A,L,A,A,A'

cti_start = cxtime(cti_start,'cal','sec')
cti_end   = cxtime(cti_end,'cal','sec')

b = where(cti_start ge startpars and cti_end le endpars, j)

if (j gt 0) then begin
;
;--- extract non-zero data
;
  cti_start = cti_start(b)
  cti_end   = cti_end(b)
  cti_obsid = cti_obsid(b)
endif else begin
  cti_start = [0]
  cti_end   = [1]
  cti_obsid = [0]
endelse
;
; ********** XMM *****************
;
readcol,"/data/mta4/space_weather/XMM/xmm.archive", $
  time_xmm,le0,le1,le2, $
  hes0,hes1,hes2,hec, skipline=2, $
  format='F,F,F,F,F,F,F,F'
;
; ************************** read counts data *************************
;
for i = 0, n_elements(cnt_file)-1 do begin
  name    = strmid(cnt_file[i], rstrpos(cnt_file[i], '/')+1, strlen(cnt_file[i])-1)
  command = 'cat '+cnt_file[i]+' > '+name+'.dat'
  spawn, command, xnum 
endfor

rdfloat, 'ccd5.dat', ccd5time, ccd5cnts
rdfloat, 'ccd6.dat', ccd6time, ccd6cnts
rdfloat, 'ccd7.dat', ccd7time, ccd7cnts

; convert from old DOM
ccd5time = (ccd5time+201+365)*86400
ccd6time = (ccd6time+201+365)*86400
ccd7time = (ccd7time+201+365)*86400

s        = sort(ccd5time)
ccd5time = ccd5time(s)
ccd5cnts = ccd5cnts(s)

s        = sort(ccd6time)
ccd6time = ccd6time(s)
ccd6cnts = ccd6cnts(s)

s        = sort(ccd7time)
ccd7time = ccd7time(s)
ccd7cnts = ccd7cnts(s)
; 
; **************** plots of all data ****************************
;
xwidth  = 750
yheight = 700
win     = 0
;
;--- choose device
;
if (keyword_set(plotx)) then begin
  set_plot, 'X'
  window, win, xsize=xwidth, ysize=yheight, retain=2
  win=win+1
endif else begin
  set_plot, 'Z'
  device, set_resolution = [xwidth, yheight]
endelse

loadct,   39
white  = 255
bgrd   =   0
red    = 230 
yellow = 190 
blue   = 100 
green  = 150 
purple =  50  
orange = 215 
;
;--- for historic reasons, all x-axes based on ephin data range
;
xmin   = time(0)
xmax   = time(n_elements(time)-1) 
xdiff  = xmax - xmin
;
;---- if this is a long term plots, set plotting range slightly differently
;
if xdiff gt 31622400 then xmax     = xmax + 0.05 * xdiff

nticks = 6              ;---- effective only for the short term (< 3 yrs)

doyticks = pick_doy_ticks(xmin,xmax-43200,num=nticks)
;
;--- for the case there are more ticks than the originally set
;
ntest    = n_elements(doyticks)
if ntest gt nticks then nticks = ntest
;
;--- start plots
;
!P.MULTI = [0, 1, 5, 0, 0]
;
;------------ goes rates ----------------
;
;
;--- find y range
;
if (n_elements(goes_p1) gt 2) then begin
  min1 = min(goes_p1(where(goes_p1 gt 0)))
  min2 = min(goes_p2(where(goes_p2 gt 0)))
  min3 = min(goes_p5(where(goes_p5 gt 0)))
  ymin = min([min1, min2, min3])
  ymax = max([goes_p1, goes_p2, goes_p5])
  no_data = 0
endif else begin
  ymin = 0
  ymax = 1
  no_data = 1
endelse

if (ymax - ymin gt 100) then setlog = 1
if (ymax - ymin le 100) then setlog = 0

plot, goes_time, replicate(0, n_elements(time)), psym = 3, $
      ytitle = 'particles/cm^2 Sr s MeV', $
      yrange = [0.001, ymax], ystyle = 1, xstyle = 1, ylog=setlog, $
      xrange = [xmin, xmax], xmargin= [8, 15], ymargin= [0, 2], $
      charsize = 1.5, charthick =1, color = white, background = bgrd, $
      xticks=nticks-1, xtickv = doyticks, $
      xtickformat='no_axis_labels', xminor=10

if (no_data eq 0) then begin
  oplot, goes_time, goes_p1, color = purple
  oplot, goes_time, goes_p2, color = green
  oplot, goes_time, goes_p5, color = blue
endif

if (no_data eq 1) then begin
  xyouts, (xmin+xmax)/2, ymax/2, 'No valid data', alignment = 0.5, /data
endif

xyouts, 1.0, 0.95, '0.8-4', alignment = 1.0, color = purple, /normal
xyouts, 1.0, 0.93, '4-9', alignment = 1.0, color = green, /normal
xyouts, 1.0, 0.91, '40-80', alignment = 1.0, color = blue, /normal
xyouts, 1.0, 0.89, '(keV)', alignment = 1.0, color = white, /normal
;
;--- switch to GOES10 at 2003:098:15:00:00
;--- now Goes 13
;
xyouts, xmax+((xmax-xmin)*0.025), ymax, 'GOES13 RATES', alignment=1.0, $
        orientation = 90, charsize=1.1
;
;--- write calender day for reference
;
xyouts, xmin, ymax*1.3, cxtime(xmin,'sec','cal'), alignment = 0.0, /data
xyouts, xmax, ymax*1.3, cxtime(xmax,'sec','cal'), alignment = 1.0, /data
;
;----------- ace rates -------------
;
;
;--- find y range
;
if (n_elements(ace_ch1) gt 2) then begin
  min1 = max([min(ace_ch1),0.01])
  min2 = max([min(ace_ch2),0.01])
  min3 = max([min(ace_ch3),0.01])
  min4 = max([min(ace_ch4),0.01])
  min5 = max([min(ace_ch5),0.01])
  ymin = min([min1, min2, min3, min4, min5])
  ymax = max([ace_ch1, ace_ch2, ace_ch3, ace_ch4, ace_ch5])

  no_data = 0
endif else begin
  ymin    = 0
  ymax    = 1
  no_data = 1
endelse


if (ymax - ymin gt 100) then setlog = 1
if (ymax - ymin le 100) then setlog = 0

plot, ace_time, replicate(0, n_elements(time)), psym = 3, $
      ;xtitle = 'time (DOM)', $
      ytitle = 'particles/cm^2 Sr s MeV', $
      yrange = [0.001, ymax], ystyle = 1, xstyle = 1, ylog=setlog, $
      xrange = [xmin, xmax], xmargin= [8, 15], ymargin= [0, 0], $
      charsize = 1.5, charthick =1, color = white, background = bgrd, $
      xticks=nticks-1, xtickv = doyticks, $
      xtickformat='no_axis_labels', xminor=10

if (no_data eq 0) then begin
  oplot, ace_time, ace_ch1, color = green
  oplot, ace_time, ace_ch2, color = yellow
  oplot, ace_time, ace_ch3, color = orange
  oplot, ace_time, ace_ch4, color = white
  oplot, ace_time, ace_ch5, color = red
endif

if (no_data eq 1) then begin
  xyouts, (xmin+xmax)/2, ymax/2, 'No valid data', alignment = 0.5, /data
endif

xyouts, 1.0, 0.77, '47-65',     alignment = 1.0, color = green,  /normal
xyouts, 1.0, 0.75, '112-187',   alignment = 1.0, color = yellow, /normal
xyouts, 1.0, 0.73, '310-580',   alignment = 1.0, color = orange, /normal
xyouts, 1.0, 0.71, '761-1220',  alignment = 1.0, color = white,  /normal
xyouts, 1.0, 0.69, '1060-1910', alignment = 1.0, color = red,    /normal
xyouts, 1.0, 0.67, '(keV)',     alignment = 1.0, color = white,  /normal
xyouts, xmax+((xmax-xmin)*0.025), ymax, 'ACE RATES', alignment=1.0, $
        orientation = 90, charsize=1.1
;
;--------------- xmm rates -----------------------
;
le1_color  = purple
le2_color  = blue
hes1_color = red
hes2_color = orange
hec_color  = yellow

ymin       = min(le1(where(le1 gt 0)))
ymax       = max([le1,le2,hes1,hes2,hec])
ymax       = 2.0e5

plot,time_xmm,le2, $
  xstyle   = 1,ystyle=1,xrange=[xmin,xmax],yrange=[ymin,ymax],/nodata, $
  xtitle   = '', ytitle = 'counts/sec', title = '', $
  xmargin  = [8, 15], ymargin= [0, 0], /ylog, $
  charsize = 1.5, charthick =1, color = white, background = bgrd, $
  xticks   = nticks-1, xtickv = doyticks, $
  xtickformat='no_axis_labels', xminor=10

oplot,time,le1, color=le1_color
oplot,time,le2, color=le2_color
oplot,time,hes1,color=hes1_color
oplot,time,hes2,color=hes2_color
oplot,time,hec, color=hec_color

xyouts, 1.0, 0.56,"LE1",color=le1_color,align=1.0, /norm
xyouts, 1.0, 0.54,"LE2",color=le2_color,align=1.0, /norm
xyouts, 1.0, 0.52,"HES1",color=hes1_color,align=1.0, /norm
xyouts, 1.0, 0.50,"HES2",color=hes2_color,align=1.0, /norm
xyouts, 1.0, 0.48,"HEC",color=hec_color,align=1.0, /norm
xyouts, xmax+((xmax-xmin)*0.025), 1e5, 'XMM RATES', alignment=1.0, $
        orientation = 90, charsize=1.1
;
;---------------  ephin rates -----------------------
;
;; note to self, recall rich logan conversation, 11/13/02:
;; here is what is REALLY happening.  The level 1 lightcurve data 
;; that goes into the plots is actually the sum of several counts. 
;; That is, P4 = P4GM+P4GR+P4S, likewise for P41, while E1300 is just E1300.  
;; Therefore, what I have labeled as P41GM (now changed the SCP41)
;; is really the sum of the three P41 counts.
;;   Now, the geometry factor can be up to 3.8 for P41GR and P41S in large 
;; format, which gives a minimum fluence of 0.004.  But, what I am finally 
;; plotting are 5 minute averages, so some points can be even lower than
;; 0.004 (eg. if all counts are 0, save 1 minimum reading in the 5 minute
;; sample).
;
;---- find y range
;
tmp1 = min(e1300(where(e1300 gt 0)))
tmp4 = min(e1300(where(e1300 gt 0)))
x    = where(ephsca00 gt 0, xn)
if (xn gt 0) then tmp4 = min(ephsca00(x))

ymin   = min([tmp1, tmp4])
tmp1   = max(e1300)
tmp4   = max(ephsca00)
ymax   = max([tmp1, tmp4])

if (ymax - ymin gt 100) then setlog = 1
if (ymax - ymin le 100) then setlog = 0

plot, time, replicate(0, n_elements(time)), psym = 3, $
      ytitle = 'particles/cm^2 Sr s', $
      yrange = [0.1, ymax], ystyle = 1, xstyle = 1, /ylog, $
      xrange = [xmin, xmax], xmargin= [8, 15], ymargin= [1, 0], $
      charsize = 1.5, charthick =1, color = white, background = bgrd, $
      xticks=nticks-1, xtickv = doyticks, $
      xtickformat='no_axis_labels', xminor=10
;
;--- overplot ephin rates
;
oplot, time, e1300, color = red
xyouts, 1.0, 0.36, 'E1300', alignment = 1.0, color = red, /normal
oplot, time, e150, color = yellow
xyouts, 1.0, 0.34, 'E150', alignment = 1.0, color = yellow, /normal
;
;--- plot operation limits
;
oplot, [xmin,xmax], [e1300_lim, e1300_lim], color = red
oplot, [xmin,xmax], [e150_lim, e150_lim], color = yellow

xyouts, xmax+((xmax-xmin)*0.025), ymax, 'EPHIN RATES', alignment=1.0, $
        orientation = 90, charsize=1.1
;
;---------------- acis rates ---------
;
plot, ccd5time, replicate(0, n_elements(ccd5time)), psym=3, $
      xtitle = 'time (DOY)', ytitle = 'counts/sec', title = '', $
      yrange = [0, 50], ystyle = 1, xstyle = 1, $
      xrange = [xmin, xmax], xmargin= [8, 15], ymargin= [4, -1], $
      charsize = 1.5, charthick =1, color = white, background = bgrd, $
      xticks=nticks-1, xtickv = doyticks, xminor=10, $
      xtickformat='s2doy_axis_labels'

oplot, ccd5time, ccd5cnts/300, psym=3, color = red
oplot, ccd6time, ccd6cnts/300, psym=3, color = yellow
oplot, ccd7time, ccd7cnts/300, psym=3, color = blue

xyouts, 1.0, 0.15, 'ccd5', alignment = 1.0, color = red, /normal
xyouts, 1.0, 0.13, 'ccd6', alignment = 1.0, color = yellow, /normal
xyouts, 1.0, 0.11, 'ccd7', alignment = 1.0, color = blue, /normal
xyouts, xmax+((xmax-xmin)*0.025), 50, 'ACIS RATES', alignment=1.0, $
        orientation = 90, charsize=1.1

ycon = convert_coord([0],[0], /device, /to_data)

label_year, xmin, xmax, ycon(1), csize=1.2
write_gif, 'rad_cnts.gif', tvrd()

;
; ************** Config plots ****************
;
xwidth  = 750
yheight = 420
;
;--- choose device
;
if (keyword_set(plotx)) then begin
  set_plot, 'X'
  window, win, xsize=xwidth, ysize=yheight,retain=2
  win = win + 1
endif else begin
  set_plot, 'Z'
  device, set_resolution = [xwidth, yheight]
endelse

!P.MULTI = [0, 1, 6, 0, 0]
xyouts, 0.1, 0.5, 'Configuration', alignment = 0.5, color = white, $
         orientation = 90, charsize=1.5, /normal
;
;--- plot s/c config
;
hetg_start   = fltarr(1)
hetg_end     = fltarr(1)
letg_start   = fltarr(1)
letg_end     = fltarr(1)
acisi_start  = fltarr(1)
acisi_end    = fltarr(1)
aciss_start  = fltarr(1)
aciss_end    = fltarr(1)
hrci_start   = fltarr(1)
hrci_end     = fltarr(1)
hrcs_start   = fltarr(1)
hrcs_end     = fltarr(1)
radmon_start = fltarr(1)
radmon_end   = fltarr(1)
fmt_start    = fltarr(1)
fmt_list     = strarr(1)
hetg_in      = 0
letg_in      = 0
acisi_in     = 0
aciss_in     = 0
hrci_in      = 0
hrcs_in      = 0
radmon_on    = 0
fmt_on       = 'FMT0'
;
;--- config_files are split into years, so only read the ones needed
;
secs_per_year = 31536000
start_year    = 63158400
needed        = 0

while (startpars gt start_year) do begin
  start_year = start_year + secs_per_year
  needed = needed + 1
endwhile

b    = intarr(1)
b(0) = needed

while (endpars gt start_year) do begin
  needed     = needed + 1
  b          = [b, needed]
  start_year = start_year + secs_per_year
endwhile

if ((b(0)) eq n_elements(config_files)-1) then b=b(0)

b = indgen(n_elements(config_files))

for k = 0, n_elements(b) - 1 do begin
;
;---- figure num lines input
;
  xnum    = strarr(1)
  command = 'wc -l '+config_files(b(k))

  spawn, command, xnum 

  xxnum  = fltarr(2)
  xxnum  = strsplit(xnum(0),' ', /extract)
  numobs = long(xxnum(0))
  
  get_lun, iunit
  openr,   iunit, config_files(b(k))
  
  array = strarr(numobs)
  print, "Reading ", config_files(b(k))         ;---  debug
  readf, iunit, array
    
  free_lun, iunit
   
  for i = 45L, numobs-1 do begin        ;--- don't start at beginning, 
                                        ;--- there is bad data there
    tmp  = strarr(16)
    tmp  = strsplit(array(i), '	', /extract)
    ttmp = fltarr(5)
    ttmp = strsplit(tmp(0),':', /extract)
;
;--- time in seconds since 1998:00:00:00
;
    ctime = float(((ttmp(0) - 1998) * 31536000) $
                    + (ttmp(1) * 86400) + (double(ttmp(2)) * 3600) $
                    + (ttmp(3) * 60) + ttmp(4) - 86400)

    leap = 2000
    while (ttmp(0) gt leap) do begin
      ctime = ctime + 86400             ;---- add leap day
      leap = leap + 4 
    endwhile
  
    last_time = ctime
    if ((ctime ge startpars) and (ctime le endpars)) then begin
      if (float(tmp[1]) gt 89000 and acisi_in eq 0) then begin
        acisi_start = [acisi_start, ctime]
        acisi_in = 1
      endif

      if (float(tmp[1]) lt 89000 and acisi_in eq 1) then begin
        acisi_end = [acisi_end, ctime]
        acisi_in = 0 
      endif

      if (float(tmp[1]) lt 80000 and float(tmp[1]) gt 71000 and aciss_in eq 0) then begin
        aciss_start = [aciss_start, ctime]
        aciss_in = 1
      endif

      if ((float(tmp[1]) gt 80000 or float(tmp[1]) lt 71000) and aciss_in eq 1) then begin
        aciss_end = [aciss_end, ctime]
        aciss_in = 0
      endif

      if (float(tmp[1]) lt -45000 and float(tmp[1]) gt -55000 and hrci_in eq 0) then begin
        hrci_start = [hrci_start, ctime]
        hrci_in = 1
      endif

      if ((float(tmp[1]) gt -45000 or float(tmp[1]) lt -55000) and hrci_in eq 1) then begin
        hrci_end = [hrci_end, ctime]
        hrci_in = 0
      endif

      if (float(tmp[1]) lt -90000 and hrcs_in eq 0) then begin
        hrcs_start = [hrcs_start, ctime]
        hrcs_in = 1
      endif

      if (float(tmp[1]) gt -90000 and hrcs_in eq 1) then begin
        hrcs_end = [hrcs_end, ctime]
        hrcs_in = 0
      endif

      if (tmp[9] lt 20 and hetg_in eq 0) then begin
        hetg_start = [hetg_start, ctime]
        hetg_in = 1
      endif

      if (tmp[9] gt 60 and hetg_in eq 1) then begin
        hetg_end = [hetg_end, ctime]
        hetg_in = 0
      endif

      if (tmp[11] lt 20 and letg_in eq 0) then begin
        letg_start = [letg_start, ctime]
        letg_in = 1
      endif

      if (tmp[11] gt 60 and letg_in eq 1) then begin
        letg_end = [letg_end, ctime]
        letg_in = 0
      endif

      if (strtrim(tmp(5),2) eq 'DISA' and radmon_on eq 0) then begin
        radmon_start = [radmon_start, ctime]
        radmon_on = 1
      endif

      if (strtrim(tmp(5),2) eq 'ENAB' and radmon_on eq 1) then begin
        radmon_end = [radmon_end, ctime]
        radmon_on = 0
      endif

      fmt_test = strtrim(tmp(7),2) 
      if (fmt_test ne fmt_on) then begin
        fmt_start = [fmt_start, ctime]
        fmt_list = [fmt_list, fmt_test]
        fmt_on = fmt_test
      endif
    endif
  endfor
endfor

hetg_end   = [hetg_end,   last_time]
letg_end   = [letg_end,   last_time]
acisi_end  = [acisi_end,  last_time]
aciss_end  = [aciss_end,  last_time]
hrci_end   = [hrci_end,   last_time]
hrcs_end   = [hrcs_end,   last_time]
radmon_end = [radmon_end, last_time]
fmt_start  = [fmt_start,  last_time]
;
;------------ plot grating intervals ------------------
;
plot, time, replicate(0, n_elements(time)), psym = 3, $
        yticks = 1, yticklen = 0, $
        yrange = [ymin, ymax], ystyle = 1, xstyle = 1, ylog=setlog, $
        xrange = [xmin, xmax], xmargin= [8, 15], ymargin= [0, 2], $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks=nticks-1, xtickv = doyticks, xminor=10, $
        ytickformat='no_axis_labels', $
        xtickformat='no_axis_labels'
  
  for i = 1, n_elements(hetg_start)-1 do begin
    tstart = max([hetg_start(i), xmin])
    tend   = min([hetg_end(i), xmax])
    if (tend gt tstart) then $
    polyfill, [tstart, tstart, tend, tend], [ymin, ymax, ymax, ymin], $
                      color = yellow, /data
  endfor

  xyouts, 1.0, 0.93, 'HETG', alignment = 1.0, $
                    color = yellow, /normal
  for i = 1, n_elements(letg_start)-1 do begin
    tstart = max([letg_start(i), xmin])
    tend   = min([letg_end(i), xmax])
    if (tend gt tstart) then $
    polyfill, [tstart, tstart, tend, tend], [ymin, ymax, ymax, ymin], $
                      color = blue, /data
  endfor

  xyouts, 1.0, 0.90, 'LETG', alignment = 1.0, $
                    color = blue, /normal
;
;------------ plot detector intervals  -------------------
;
plot, time, replicate(0, n_elements(time)), psym = 3, $
        yticks = 1, yticklen = 0, $
        yrange = [ymin, ymax], ystyle = 1, xstyle = 1, ylog=setlog, $
        xrange = [xmin, xmax], xmargin= [8, 15], ymargin= [0, 0], $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks=nticks-1, xtickv = doyticks, xminor=10, $
        ytickformat='no_axis_labels', $
        xtickformat='no_axis_labels'
  
  for i = 1, n_elements(acisi_start)-1 do begin
    tstart = max([acisi_start(i), xmin])
    tend   = min([acisi_end(i), xmax])
    if (tend gt tstart) then $
    polyfill, [tstart, tstart, tend, tend], [ymin, ymax, ymax, ymin], $
                      color = purple, /data
  endfor

  xyouts, 1.0, 0.78, 'ACIS-I', alignment = 1.0, $
                    color = purple, /normal

  for i = 1, n_elements(aciss_start)-1 do begin
    tstart = max([aciss_start(i), xmin])
    tend   = min([aciss_end(i), xmax])
    if (tend gt tstart) then $
    polyfill, [tstart, tstart, tend, tend], [ymin, ymax, ymax, ymin], $
                      color = yellow, /data
  endfor
  xyouts, 1.0, 0.75, 'ACIS-S', alignment = 1.0, $
                    color = yellow, /normal

  for i = 1, n_elements(hrci_start)-1 do begin
    tstart = max([hrci_start(i), xmin])
    tend   = min([hrci_end(i), xmax])
    if (tend gt tstart) then $
    polyfill, [tstart, tstart, tend, tend], [ymin, ymax, ymax, ymin], $
                      color = blue, /data
  endfor
  xyouts, 1.0, 0.72, 'HRC-I', alignment = 1.0, $
                    color = blue, /normal

  for i = 1, n_elements(hrcs_start)-1 do begin
    tstart = max([hrcs_start(i), xmin])
    tend   = min([hrcs_end(i), xmax])
    if (tend gt tstart) then $
    polyfill, [tstart, tstart, tend, tend], [ymin, ymax, ymax, ymin], $
                      color = green, /data
  endfor
  xyouts, 1.0, 0.69, 'HRC-S', alignment = 1.0, $
                    color = green, /normal
;
;------ find RAD/CTI intervals  --------------------
;
t_max     = 500000                  ;---- sort of works to find missing data (gaps)
t_min     = 600                     ;---- minimum time in rad zone, to avoid false detections
num       = n_elements(time)
j         = 0                       ;---- gti indexer
rad_start = fltarr(1)
rad_end   = fltarr(1)

for i = 0L, num - 1 do begin
  rad_in = 0
  tstart = time(i)  ; initialize 
  tend   = time(i)  ; initialize 

  if (i lt num - 2) then i = i + 1

  while (((e1300(i) gt e1300_lim) or $
         (e150(i) gt e150_lim) or $
         (p_4gm(i) gt p_4gm_lim) or $
         (p41gm(i) gt p41gm_lim)) and (i lt num - 2)) do begin
    if (rad_in eq 0) then begin
      tstart = time(i)
      rad_in = 1
    endif

    i = i + 1
  endwhile

  if (rad_in eq 1) then begin
    tend = time(i-1)                ;--- we've been in the rad zone, what time did we exit
  endif                             ;--- rad_in is reset at top of loop

  if ((tend - tstart lt t_max) and (tend-tstart gt t_min)) then begin
    rad_start = [rad_start, tstart]
    rad_end = [rad_end, tend]
    j = j + 1
  endif

endfor

j = 1
if (j gt 0) then begin
  
  plot, time, replicate(0, n_elements(time)), psym = 3, $
        yticks = 1, yticklen = 0, $
        ;xtitle = 'time (DOM)', $
        ;ytitle = '', title = 'Usable time intervals', $
        yrange = [ymin, ymax], ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin= [8, 15], ymargin= [0, 0], $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks=nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='no_axis_labels', $
        ytickformat='no_axis_labels'
;  
;--- plot cti measurement intervals
;
  for i = 0, n_elements(cti_obsid)-1 do begin
    tstart = cti_start(i)
    tend   = cti_end(i)
    polyfill, [tstart, tstart, tend, tend], [ymin, ymax, ymax, ymin], $
                      color = yellow, /data
  endfor

  xyouts, 1.0, 0.58, 'CTI', alignment = 1.0, $
                    color = yellow, /normal
;
;--- oplot altitude
;
  rdfloat, ephemfile, $
           alt_time,alt,alt_x,alt_y,alt_z,mag_x,mag_y,mag_z,crm_reg, skipline=2
  b = where(alt_time ge startpars and alt_time le endpars, cnt)
  if (cnt gt 0) then begin
    alt_time=alt_time(b)
    alt=alt(b)
    alt_x=alt_x(b)
    alt_y=alt_y(b)
    alt_z=alt_z(b)
    mag_x=mag_x(b)
    mag_y=mag_y(b)
    mag_z=mag_z(b)
    crm_reg=crm_reg(b)
  endif 
  alt_s = ((alt-min(alt))*(ymax-ymin)/(max(alt)-min(alt))) + ymin ; rescale
  oplot, alt_time, alt_s, color=blue
  xyouts, 1.0, 0.55, 'ALTITUDE', alignment = 1.0, $
                    color = blue, /normal
;
;--- find crm regions
;
  wind_start=fltarr(1)
  wind_stop=fltarr(1)
  shth_start=fltarr(1)
  shth_stop=fltarr(1)
  sphr_start=fltarr(1)
  sphr_stop=fltarr(1)
  r = 0L
  while (r lt cnt-2) do begin
    reg=crm_reg(r)

    case reg of
      1: wind_start=[wind_start,alt_time(r)]
      2: shth_start=[shth_start,alt_time(r)]
      3: sphr_start=[sphr_start,alt_time(r)]
      else: print, "Found CRM_REG ", reg
    endcase

    while (crm_reg(r) eq reg and r lt cnt-2) do begin
      r = r + 1
    endwhile

    case reg of
      1: wind_stop=[wind_stop,alt_time(r)]
      2: shth_stop=[shth_stop,alt_time(r)]
      3: sphr_stop=[sphr_stop,alt_time(r)]
      else: print, "end CRM_REG ", reg
    endcase

  endwhile  
;
;--- plot crm regions
;
  plot, time, replicate(0, n_elements(time)), psym = 3, $
        yticks = 1, yticklen = 0, $
        yrange = [ymin, ymax], ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin= [8, 15], ymargin= [0, 0], $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks=nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='no_axis_labels', $
        ytickformat='no_axis_labels'
  
  for i = 1, n_elements(wind_start)-1 do begin
    tstart = max([wind_start(i),xmin])
    tend   = min([wind_stop(i),xmax])
    if (tend gt tstart) then $
      polyfill, [tstart, tstart, tend, tend], [ymin, ymax, ymax, ymin], $
                      color = green, /data
  endfor
  xyouts, 1.0, 0.48, 'Solar Wind', alignment = 1.0, color = green, /normal
  for i = 1, n_elements(shth_start)-1 do begin
    tstart = max([shth_start(i),xmin])
    tend   = min([shth_stop(i),xmax])
    if (tend gt tstart) then $
      polyfill, [tstart, tstart, tend, tend], [ymin, ymax, ymax, ymin], $
                      color = blue, /data
  endfor
  xyouts, 1.0, 0.45, 'Magnetosheath', alignment = 1.0, color = blue, /normal
  for i = 1, n_elements(sphr_start)-1 do begin
    tstart = max([sphr_start(i),xmin])
    tend   = min([sphr_stop(i),xmax])
    if (tend gt tstart) then $
      polyfill, [tstart, tstart, tend, tend], [ymin, ymax, ymax, ymin], $
                      color = yellow, /data
  endfor
  xyouts, 1.0, 0.42, 'Magnetosphere', alignment = 1.0, color = yellow, /normal
endif       ;---- j gt 0
;
;----------------- plot radmon ------------------------
;
plot, time, replicate(0, n_elements(time)), psym = 3, $
        yticks = 1, yticklen = 0, $
        ;xtitle = 'time (DOM)', $
        ;ytitle = '', title = 'Usable time intervals', $
        yrange = [ymin, ymax], ystyle = 1, xstyle = 1, ylog=setlog, $
        xrange = [xmin, xmax], xmargin= [8, 15], ymargin= [0, 0], $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks=nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='no_axis_labels', $
        ytickformat='no_axis_labels'
;  
;------------- plot rad ------------------------
;
  for i = 0, n_elements(radmon_start)-1 do begin
    tstart = max([radmon_start(i), xmin])
    tend   = min([radmon_end(i), xmax])

    if (tend gt tstart) then $
    polyfill, [tstart, tstart, tend, tend], [ymin, ymax, ymax, ymin], $
                      color = orange, /data
  endfor

  xyouts, 1.0, 0.33, 'RADMON', alignment = 1.0, $
                    color = orange, /normal
  xyouts, 1.0, 0.30, 'DISABLED', alignment = 1.0, $
                    color = orange, /normal
;
;------------ find sca00 saturation intervals  -------------------------------
;
t_max     = 500000      ;---- sort of works to find missing data (gaps)
t_min     = 600         ;---- minimum time in rad zone, to avoid false detections
num       = n_elements(ephsca.time)
j         = 0 ; gti indexer
sca_start = fltarr(1)
sca_end   = fltarr(1)

for i = 0L, num - 1 do begin
  rad_in = 0
  tstart = ephsca(i).time  ; initialize 
  tend   = ephsca(i).time  ; initialize 

  if (i lt num - 2) then i = i + 1

  while ((ephsca(i).sca00_avg eq 8372224) and (i lt num - 2)) do begin
    if (rad_in eq 0) then begin
      tstart = ephsca(i).time
      rad_in = 1
    endif
    i = i + 1
  endwhile

;
;---  we've been in the rad zone, what time did we exit
;---  rad_in is reset at top of loop
;
  if (rad_in eq 1) then begin
    tend = ephsca(i-1).time 
  endif              

  if ((tend - tstart lt t_max) and (tend-tstart gt t_min)) then begin
    sca_start = [sca_start, tstart]
    sca_end   = [sca_end, tend]
    j         = j + 1
  endif

endfor

if (j gt 0) then begin
  sca_start = sca_start[1:*]            ;--- delete first dummy element
  sca_end   = sca_end[1:*]              ;--- delete first dummy element
;
;--- plot sca
;
  for i = 0, n_elements(sca_start)-1 do begin
    tstart = max([sca_start(i), xmin])
    tend   = min([sca_end(i), xmax])
    if (tend gt tstart) then $
    polyfill, [tstart, tstart, tend, tend], [ymin, ymax, ymax, ymin], $
                      color = red, /data
  endfor

  xyouts, 1.0, 0.27, 'SCA00', alignment = 1.0, $
                    color = red, /normal
  xyouts, 1.0, 0.24, 'SAT', alignment = 1.0, $
                    color = red, /normal
endif ; j gt 0
; 
;----------  plot format--------------------
;
plot, time, replicate(0, n_elements(time)), psym = 3, $
        yticks   = 1, yticklen = 0, $
        xtitle   = 'time (DOY)', $
        ;ytitle  = '', title = 'Usable time intervals', $
        yrange   = [ymin, ymax], ystyle = 1, xstyle = 1, ylog=setlog, $
        xrange   = [xmin, xmax], xmargin= [8, 15], ymargin= [4, 0], $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticklen = -0.15, $
        xticks   = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat = 's2doy_axis_labels', $
        ytickformat = 'no_axis_labels'
;  
;--- plot fmt
;
  colors = [blue, green, orange, yellow, red]

  for i = 1, 5 do begin
    format = 'FMT'+strtrim(string(i), 2)
    b      = where(fmt_list eq format, num)

    if (num gt 0) then begin
      f_start = fmt_start(b)
      f_end   = fmt_start(b+1)

        for j = 0, n_elements(b)-1 do begin
          tstart = max([f_start(j), xmin])
          tend   = min([f_end(j), xmax])

          if (tend gt tstart) then $
          polyfill, [tstart, tstart, tend, tend], [ymin, ymax, ymax, ymin], $
                      color = colors(i-1), /data
        endfor
    endif
    xyouts, 1.0, 0.20 - (.03*i), format, alignment = 1.0, $
               color = colors(i-1), /normal
  endfor

  xyouts, 1.0, 0.0, 'Updated: '+systime(), alignment = 1, color = white, /normal

;  ycon = convert_coord([0],[0.90], /norm, /to_data)
  ycon = convert_coord([0],[0], /norm, /to_data)

  label_year, xmin, xmax, ycon(1), csize=1.2

  write_gif, 'rad_use.gif', tvrd()
;
;---- summarize rad zone timing
;
ref_err = 50000             ;---- max radstart and end offset from perigee

a = {mv, peri:      0.0, $
         mon_disa:  0.0, n_mon_disa:  0, mon_enab:    0.0, n_mon_enab:    0, $
         eph_viol:  0.0, n_eph_viol:  0, eph_safe:    0.0, n_eph_safe:    0, $
         sca00_sat: 0.0, n_sca00_sat: 0, sca00_unsat: 0.0, n_sca00_unsat: 0, $
         x: 'N'}                ;----- x is marker for later filtering
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;exit
;;;--- the following line generate "Floating divide by 0" 03/12/15
deriv, alt_time, alt, dv_time,dv
mn = local_min(dv)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;exit

print, "peri cnt ", n_elements(mn)                      ;---- debug
;print, alt_time(mn)                                    ;---- debug

if (mn(0) gt -1) then begin
radtab = replicate({mv}, n_elements(mn))

radmon_start = radmon_start(where(radmon_start gt 0))
radmon_end   = radmon_end(where(radmon_end gt 0))

for k = 0, n_elements(mn)-1 do begin
  peri_time = alt_time(mn(k))
  radtab(k).peri = peri_time
  b = where(rad_start gt peri_time-ref_err and $
            rad_start lt peri_time+ref_err, num)

  if (num gt 0) then $
    radtab(k).eph_viol = rad_start(b(0))

  radtab(k).n_eph_viol = num
  b = where(rad_end gt peri_time-ref_err and $
            rad_end lt peri_time+ref_err, num)

  if (num gt 0) then $
    radtab(k).eph_safe = rad_end(b(num-1))

  radtab(k).n_eph_safe = num
  b = where(sca_start gt peri_time-ref_err and $
            sca_start lt peri_time+ref_err, num)

  if (num gt 0) then $
    radtab(k).sca00_sat = sca_start(b(0))

  radtab(k).n_sca00_sat = num
  b = where(sca_end gt peri_time-ref_err and $
            sca_end lt peri_time+ref_err, num)

  if (num gt 0) then $
    radtab(k).sca00_unsat = sca_end(b(num-1))

  radtab(k).n_sca00_unsat = num
  b = where(radmon_start gt peri_time-ref_err and $
            radmon_start lt peri_time, num)

  if (num gt 0) then $
    radtab(k).mon_disa = radmon_start(b(0))

  radtab(k).n_mon_disa = num
  b = where(radmon_end gt peri_time and $
            radmon_end lt peri_time+ref_err, num)

  if (num gt 0) then $
    radtab(k).mon_enab = radmon_end(b(num-1))

  radtab(k).n_mon_enab = num
;   print, radtab(k)                    ;---- debug
  radtab(k).x = 'X' ; mark as seen
endfor ;for k = 0, n_elements(mn)-1 do begin

get_lun, sunit                          ;---- debug
openw, sunit, 'radtab.out'              ;---- debug

for i = 0, k-1 do begin  ; debug
  printf, sunit, radtab(i)  ; debug
endfor  ; debug
close, sunit                            ;---- debug
;
;---------  timing plots ---------------------
;
xwidth  = 750
yheight = 700
;
;---- choose device
;
if (keyword_set(plotx)) then begin
  set_plot, 'X'
  window, win, xsize=xwidth, ysize=yheight,retain=2
  win=win+1
endif else begin
  set_plot, 'Z'
  device, set_resolution = [xwidth, yheight]
endelse

!P.MULTI = [0, 1, 4, 0, 0]
loadct, 39
;
;----- time interval stats
;
start_diff = [radtab.eph_viol,0] - [0,radtab.eph_viol]
end_diff   = [radtab.eph_safe,0] - [0,radtab.eph_safe]
orbit_diff = [radtab.eph_viol,0] - [0,radtab.eph_safe]
rad_diff   = radtab.eph_safe-radtab.eph_viol
num        = n_elements(radtab)
start_diff = start_diff[1:num-1]
end_diff   = end_diff[1:num-1]
orbit_diff = orbit_diff[1:num-1]
rad_diff   = rad_diff[1:num-1]

b = where(start_diff lt 90.0*3600.0 and start_diff gt 0,bnum)
if (bnum gt 2) then begin
  yplot = start_diff(b)/3600.0
  ymin  = min(yplot)
  ymax  = max(yplot)
  ywide = ymax-ymin
  ymin  = ymin - 0.1 * ywide
  ymax  = ymax + 0.1 * ywide
  plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'start_diff (hr)', $
        title  = 'Time between consecutive RADZONE entries', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

b = where(end_diff lt 90.0*3600.0 and end_diff gt 0,bnum)
if (bnum gt 2) then begin
  yplot = end_diff(b)/3600.0
  ymin  = min(yplot)
  ymax  = max(yplot)
  ywide = ymax-ymin
  ymin  = ymin-0.1*ywide
  ymax  = ymax+0.1*ywide
  plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'end_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time between consecutive RADZONE exits', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

b = where(rad_diff lt 30.0*3600.0 and rad_diff gt 0,bnum)
if (bnum gt 2) then begin
  yplot = rad_diff(b)/3600.0
  ymin  = min(yplot)
  ymax  = max(yplot)
  ywide = ymax-ymin
  ymin  = ymin-0.1*ywide
  ymax  = ymax+0.1*ywide
  plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'rad_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time in RADZONE', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks=nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

b = where(orbit_diff lt 90.0*3600.0 and orbit_diff gt 0,bnum)
if (bnum gt 2) then begin
  yplot = orbit_diff(b)/3600.0
  ymin  = min(yplot)
  ymax  = max(yplot)
  ywide = ymax-ymin
  ymin  = ymin-0.1*ywide
  ymax  = ymax+0.1*ywide
  plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'orbit_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time out of RADZONE', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

ycon = convert_coord([0],[0], /device, /to_data)
label_year, xmin, xmax, ycon(1), csize=1.2

write_gif, 'eph_diff.gif', tvrd()
;
;------------ radmon diff -------------------
;
xwidth  = 750
yheight = 1050
;
;----- choose device
;
if (keyword_set(plotx)) then begin
  set_plot, 'X'
  window, win, xsize=xwidth, ysize=yheight,retain=2
  win = win+1
endif else begin
  set_plot, 'Z'
  device, set_resolution = [xwidth, yheight]
endelse

!P.MULTI = [0, 1, 6, 0, 0]
loadct, 39
start_diff = [radtab.mon_disa,0] - [0,radtab.mon_disa]
end_diff   = [radtab.mon_enab,0] - [0,radtab.mon_enab]
orbit_diff = [radtab.mon_disa,0] - [0,radtab.mon_enab]
rad_diff   = radtab.mon_enab-radtab.mon_disa
entr_diff  = radtab.eph_viol-radtab.mon_disa
exit_diff  = radtab.mon_enab-radtab.eph_safe
num        = n_elements(radtab)
start_diff = start_diff[1:num-1]
end_diff   = end_diff[1:num-1]
orbit_diff = orbit_diff[1:num-1]
rad_diff   = rad_diff[1:num-1]

b = where(start_diff lt 70.0*3600.0 and start_diff gt 0,bnum)
if (bnum gt 2) then begin
  yplot = start_diff(b)/3600.0
  ymin  = min(yplot)
  ymax  = max(yplot)
  ywide = ymax-ymin
  ymin  = ymin-0.1*ywide
  ymax  = ymax+0.1*ywide
  plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'start_diff (hr)', $
        title  = 'Time between consecutive RADMON disables', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

b = where(end_diff lt 70.0*3600.0 and end_diff gt 0,bnum)
if (bnum gt 2) then begin
  yplot = end_diff(b)/3600.0
  ymin  = min(yplot)
  ymax  = max(yplot)
  ywide = ymax-ymin
  ymin  = ymin-0.1*ywide
  ymax  = ymax+0.1*ywide
  plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'end_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time between consecutive RADMON enables', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif
  
b = where(rad_diff lt 30.0*3600.0 and rad_diff gt 0,bnum)
if (bnum gt 2) then begin
    yplot = rad_diff/3600.0
    ymin  = min(yplot)
    ymax  = max(yplot)
    ywide = ymax-ymin
    ymin  = ymin-0.1*ywide
    ymax  = ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'rad_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time in RADMON disable', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

b = where(orbit_diff lt 90.0*3600.0 and orbit_diff gt 0, bmum)
if (bnum gt 2) then begin
    yplot = orbit_diff(b)/3600.0
    ymin  = min(yplot)
    ymax  = max(yplot)
    ywide = ymax-ymin
    ymin  = ymin-0.1*ywide
    ymax  = ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'orbit_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time out of RADMON disable', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

b = where(entr_diff lt 15.0*3600.0 and entr_diff gt 0, bnum)
if (bnum gt 2) then begin
    yplot = entr_diff(b)/3600.0
    ymin  = min(yplot)
    ymax  = max(yplot)
    ywide = ymax-ymin
    ymin  = ymin-0.1*ywide
    ymax  = ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'entr_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time between RADMON disable and RADZONE entry', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

b = where(exit_diff lt 15.0*3600.0 and exit_diff gt 0, bnum)
if (bnum gt 2) then begin
    yplot = exit_diff(b)/3600.0
    ymin  = min(yplot)
    ymax  = max(yplot)
    ywide = ymax-ymin
    ymin  = ymin-0.1*ywide
    ymax  = ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'exit_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time between RADZONE exit and RADMON enable', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

ycon = convert_coord([0],[0], /device, /to_data)
label_year, xmin, xmax, ycon(1), csize=1.2

write_gif, 'mon_diff.gif', tvrd()
;
;---- sca00 diff ---------
;---- sca00 is not ploted any more -----
xwidth  = 750
yheight = 1400
;
;---- choose device
;
if (keyword_set(plotx)) then begin
  set_plot, 'X'
  window, win, xsize=xwidth, ysize=yheight,retain=2
  win=win+1
endif else begin
  set_plot, 'Z'
  device, set_resolution = [xwidth, yheight]
endelse

!P.MULTI = [0, 1, 8, 0, 0]
loadct, 39
start_diff = [radtab.sca00_sat,0] - [0,radtab.sca00_sat]
end_diff   = [radtab.sca00_unsat,0] - [0,radtab.sca00_unsat]
orbit_diff = [radtab.sca00_sat,0] - [0,radtab.sca00_unsat]
rad_diff   = radtab.sca00_unsat-radtab.sca00_sat
entr_diff  = radtab.sca00_sat-radtab.mon_disa
exit_diff  = radtab.mon_enab-radtab.sca00_unsat
ephs_diff  = radtab.sca00_sat-radtab.eph_viol
ephe_diff  = radtab.sca00_unsat-radtab.eph_safe
num        = n_elements(radtab)
start_diff = start_diff[1:num-1]
end_diff   = end_diff[1:num-1]
orbit_diff = orbit_diff[1:num-1]
rad_diff   = rad_diff[1:num-1]
;
;-------- perigee diff --------------------
;
xwidth  = 750
yheight = 1225
;
;--- choose device
;
if (keyword_set(plotx)) then begin
  set_plot, 'X'
  window, win, xsize=xwidth, ysize=yheight,retain=2
  win = win + 1
endif else begin
  set_plot, 'Z'
  device, set_resolution = [xwidth, yheight]
endelse

!P.MULTI = [0, 1, 7, 0, 0]
loadct, 39
orbit_diff = [radtab.peri,0] - [0,radtab.peri]
entr_diff  = radtab.peri-radtab.mon_disa
exit_diff  = radtab.mon_enab-radtab.peri
ephs_diff  = radtab.peri-radtab.eph_viol
ephe_diff  = radtab.eph_safe-radtab.peri
scas_diff  = radtab.peri-radtab.sca00_sat
scae_diff  = radtab.sca00_unsat-radtab.peri
num        = n_elements(radtab)
orbit_diff = orbit_diff[1:num-1]

b = where(orbit_diff lt 90.0*3600.0 and orbit_diff gt 0, bnum)
if (bnum gt 2) then begin
    yplot = orbit_diff(b)/3600.0
    ymin  = min(yplot)
    ymax  = max(yplot)
    ywide = ymax-ymin
    ymin  = ymin-0.1*ywide
    ymax  = ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'orbit_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time between consecutive Perigees', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

b = where(entr_diff gt 0, bnum)
if (bnum gt 2) then begin
    yplot = entr_diff(b)/3600.0
    ymin  = min(yplot)
    ymax  = max(yplot)
    ywide = ymax-ymin
    ymin  = ymin-0.1*ywide
    ymax  = ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'entr_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time between RADMON disable and Perigee', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

b = where(exit_diff gt 0, bnum)
if (bnum gt 2) then begin
    yplot = exit_diff(b)/3600.0
    ymin  = min(yplot)
    ymax  = max(yplot)
    ywide = ymax-ymin
    ymin  = ymin-0.1*ywide
    ymax  = ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'exit_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time between Perigee and RADMON enable', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

b = where(abs(ephs_diff) lt 15.0*3600.0,bnum)
if (bnum gt 2) then begin
    yplot = ephs_diff(b)/3600.0
    ymin  = min(yplot)
    ymax  = max(yplot)
    ywide = ymax-ymin
    ymin  = ymin-0.1*ywide
    ymax  = ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'ephs_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time between RADZONE entry and Perigee', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

b = where(abs(ephe_diff) lt 15.0*3600.0, bnum)
if (bnum gt 2) then begin
    yplot = ephe_diff(b)/3600.0
    ymin  = min(yplot)
    ymax  = max(yplot)
    ywide = ymax-ymin
    ymin  = ymin-0.1*ywide
    ymax  = ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'ephe_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time between Perigee and RADZONE exit', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'

;b = where(abs(scas_diff) lt 15.0*3600.0)
    ymin  = min(yplot)
    ymax  = max(yplot)
    ywide = ymax-ymin
    ymin  = ymin-0.1*ywide
    ymax  = ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'scas_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time between SCA00 sat and Perigee', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels',/nodata

;b = where(abs(scae_diff) lt 15.0*3600.0)
;yplot=scae_diff(b)/3600.0
    ymin=min(yplot)
    ymax=max(yplot)
    ywide=ymax-ymin
    ymin=ymin-0.1*ywide
    ymax=ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'scae_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        ;max_value = mean(0)+(2*sqrt(mean(1))), $
        ;min_value = mean(0)-(2*sqrt(mean(1))), $
        title = 'Time between Perigee and SCA00 unsat', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks=nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels',/nodata
endif

ycon = convert_coord([0],[0], /device, /to_data)
label_year, xmin, xmax, ycon(1), csize=1.2

write_gif, 'per_diff.gif', tvrd()
;
;--- make separate plot for monthly report, xv does not crop-----------
;---  long gif well
;
xwidth  = 750
yheight = 700
;
;---- choose device
;
if (keyword_set(plotx)) then begin
  set_plot, 'X'
  window, win, xsize=xwidth, ysize=yheight,retain=2
  win=win+1
endif else begin
  set_plot, 'Z'
  device, set_resolution = [xwidth, yheight]
endelse
!P.MULTI = [0, 1, 4, 0, 0]

rad_diff = radtab.eph_safe-radtab.eph_viol
b = where(rad_diff lt 30.0*3600.0 and rad_diff gt 0, bnum)
;yplot=rad_diff(b)/3600.0
ymin=min(yplot)
ymax=max(yplot)
ywide=ymax-ymin
ymin=ymin-0.1*ywide
ymax=ymax+0.1*ywide
if (bnum gt 2) then begin
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'rad_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        ;max_value = mean(0)+(2*sqrt(mean(1))), $
        ;min_value = mean(0)-(2*sqrt(mean(1))), $
        title = 'Time in RADZONE', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks=nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

rad_diff = radtab.mon_enab-radtab.mon_disa

b = where(rad_diff gt 0, bnum)
if (bnum gt 2) then begin
    yplot = rad_diff(b)/3600.0
    ymin  = min(yplot)
    ymax  = max(yplot)
    ywide = ymax-ymin
    ymin  = ymin-0.1*ywide
    ymax  = ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'rad_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time in RADMON disable', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

entr_diff = radtab.eph_viol-radtab.mon_disa
exit_diff = radtab.mon_enab-radtab.eph_safe

b = where(entr_diff gt 0, bnum)
if (bnum gt 2) then begin
    yplot = entr_diff(b)/3600.0
    ymin  = min(yplot)
    ymax  = max(yplot)
    ywide = ymax-ymin
    ymin  = ymin-0.1*ywide
    ymax  = ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'entr_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time between RADMON disable and RADZONE entry', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif

b = where(exit_diff gt 0, bnum)
if (bnum gt 2) then begin
    yplot = exit_diff(b)/3600.0
    ymin  = min(yplot)
    ymax  = max(yplot)
    ywide = ymax-ymin
    ymin  = ymin-0.1*ywide
    ymax  = ymax+0.1*ywide
    plot, radtab(b).peri, yplot, psym=1, $
        xtitle = 'time (DOY)', ytitle = 'exit_diff (hr)', $
        ystyle = 1, xstyle = 1, $
        xrange = [xmin, xmax], xmargin = [12,2], $
        yrange = [ymin, ymax], $
        title  = 'Time between RADZONE exit and RADMON enable', $
        charsize = 1.8, charthick =1, color = white, background = bgrd, $
        xticks = nticks-1, xtickv = doyticks, xminor=10, $
        xtickformat='s2doy_axis_labels'
endif 
ycon = convert_coord([0],[0], /device, /to_data)
label_year, xmin, xmax, ycon(1), csize=1.2

write_gif, 'xper_diff.gif', tvrd()

endif else begin  ;if (mn gt -1) then begin
  print, "Can't do timing plots, no rad zones in this time interval"
endelse
;
;---- end 
;
endif ; done = 0
end
