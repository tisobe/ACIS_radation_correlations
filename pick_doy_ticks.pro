FUNCTION PICK_DOY_TICKS, min, max, NUM=num
;
;--- give start and end times in secs
;--- number of tickmarks is optional (default provided)
;--- return tick values to be used to give round DOY values
;---  see also s2doy_axis_labels
;--- 24. Jan 2002 BDS
;--- 16.Mar 2004 BDS -try to get easy times, ie. Jan01, Feb01 ...
;
;--- this one is updated by t. isobe (tisobe@cfa.harvard.edu) Mar 10, 2015
;---    * it finds the beginning year
;---    * a tick starts from that year and then creates ticks every 2 years
;---    note, since it assumes that year is 366 days, there will be
;---          an error in a very long run.
;

if (NOT keyword_set(num)) then num=5
;
;--- use a conversion to met to get a round value
;--- but first calibrate met conversion, in case definition has changed
;
met_calib = cxtime('1999-01-01T00:00:00', 'cal', 'met')
rnd_day   = met_calib - fix(met_calib)

min_met   = (1.0*fix(cxtime(min, 'sec', 'met')))+rnd_day+1.0
max_met   = (1.0*fix(cxtime(max, 'sec', 'met')))+rnd_day+1.0

range     = max_met - min_met
;
;--- if the range is longer than three years, do the next
;
if (range gt (366*3)) then begin
    sperdy  = 86400.0
    speryr  = sperdy * 366.0
    ybefore = long((min - 31536000) / speryr)
    yend    = long((max - 31536000) / speryr)
    range   = long((yend - ybefore +1)/2)
    start   = 1999 + ybefore
    array   = strarr(range+1)
    for i = 0, range do begin
        lyear = strtrim(string(start + 2 * i)) + '-01-01T00:00:00'
        array[i] = lyear
    endfor
    ticks = cxtime(array, 'cal', 'sec')
        
  return, ticks
endif 
;
;--- if the range is shorter than three years, do the follow
;
interval = 0
min_int  = 1.0
i = 0
;
;--- find best scale
;
while (interval lt min_int) do begin
  interval = long((1.0*fix(range*10^i/(num-1)))/10.0^i)
  min_int  = min_int/10.0
  i = i + 1
endwhile

ticks = cxtime(min_met+(interval*indgen(num)),'met','sec')
return, ticks
end
