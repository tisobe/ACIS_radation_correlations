x = mta_rad('2015:040')
spawn, 'mv rad_cnts.gif rad_cnts_curr.gif'
spawn, 'mv rad_use.gif rad_use_curr.gif'
print, 'Done current'
retall
x = mta_rad('2015:061.0','2015:090')
spawn, 'mv rad_cnts.gif rad_cnts_mar15.gif'
spawn, 'mv rad_use.gif rad_use_mar15.gif'
spawn, 'mv eph_diff.gif eph_diff_mar15.gif'
spawn, 'mv mon_diff.gif mon_diff_mar15.gif'
spawn, 'mv per_diff.gif per_diff_mar15.gif'
spawn, 'mv xper_diff.gif mon_per_diff_mar15.gif'
print, 'Done mar15'
x = mta_rad('2014:001','2015:001')
spawn, 'mv rad_cnts.gif rad_cnts_14.gif'
spawn, 'mv rad_use.gif rad_use_14.gif'
spawn, 'mv eph_diff.gif eph_diff_14.gif'
spawn, 'mv mon_diff.gif mon_diff_14.gif'
spawn, 'mv per_diff.gif per_diff_14.gif'
spawn, 'mv xper_diff.gif mon_per_diff_14.gif'
print, 'Done 14'
x = mta_rad('2014:070','2015:070')
spawn, 'mv rad_cnts.gif rad_cnts_last_one_year.gif'
spawn, 'mv rad_use.gif rad_use_last_one_year.gif'
spawn, 'mv eph_diff.gif eph_diff_last_one_year.gif'
spawn, 'mv mon_diff.gif mon_diff_last_one_year.gif'
spawn, 'mv per_diff.gif per_diff_last_one_year.gif'
spawn, 'mv xper_diff.gif mon_per_diff_last_one_year.gif'
print, 'Done Last One Year
x = mta_rad()
spawn, 'mv rad_cnts.gif rad_cnts_all.gif'
spawn, 'mv rad_use.gif rad_use_all.gif'
spawn, 'mv eph_diff.gif eph_diff_all.gif'
spawn, 'mv mon_diff.gif mon_diff_all.gif'
spawn, 'mv per_diff.gif per_diff_all.gif'
print, 'Done all'
retall
exit
