#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################################
#                                                                                                           #
#       update_html_page.py: update radiation related html page                                             #
#                                                                                                           #
#           author: t. isobe    (tiosbe@cfa.harvard.edu)                                                    #
#                                                                                                           #
#           last update: Mar 02, 2015                                                                       #
#                                                                                                           #
#############################################################################################################

import os
import sys
import re
import subprocess

mta_dir   = '/data/mta/Script/Python_script2.7/'
sys.path.append(mta_dir)

import convertTimeFormat as tcnv    #---- contains MTA time conversion routines

mon_list1 = ['031', '060', '091', '121', '152', '182', '213', '244', '274', '305', '335', '366']
mon_list2 = ['031', '060', '090', '120', '151', '181', '212', '243', '273', '304', '334', '365']

fmon_list = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']


#----------------------------------------------------------------------------------------------------------
#--  print_html: create and/or update radiation related html page                                       ---
#----------------------------------------------------------------------------------------------------------

def print_html(year, mon):
    """
    create and/or update radiation related html page
    """
#
#--- find today's date
#
    if year == '':
        [year, mon, day, hours, min, sec, weekday, yday, dst] = tcnv.currentTime()
#
#--- find out the last month
#
#        cyear = year
#        lmon      = mon -1
#        if lmon < 1:
#            lmon  = 12
#            cyear = year -1

#
#--- for the case year and mon are given
#
    cyear = year
    lmon  = mon
#
#--- choose a correct month list depending on whether this is the leap year
#
    if tcnv.isLeapYear(cyear) == 1:
        mon_list = mon_list1
    else:
        mon_list = mon_list2

    last_day  = mon_list[lmon-1]
#
#--- convert the month from a numeric to letter
#
    umon      = tcnv.changeMonthFormat(lmon)
    smon      = umon.lower()

    lmon_year =  str(cyear)
    syear     =  lmon_year[2] + lmon_year[3]
    last_year =  str(year -1)
    syear2    =  last_year[2] + last_year[3]
    monyear   =  smon + syear
#
#--- set output html page names
#
    year_html = 'all' + syear + '.html'
    mon_html  = monyear + '.html'
    rad_html  = 'rad_time_' + monyear + '.html'
#
#--- read yearly html page template
#
    data = open('./Template/yearly_template', 'r').read()

    data = data.replace('$#FYEAR#$', str(year))
    data = data.replace('$#SYEAR#$', syear)

    fo   = open(year_html, 'w')
    fo.write(data)
    fo.close()
#
#--- read monthly html page template
#
    data = open('./Template/monthly_template', 'r').read()

    data = data.replace('$#FYEAR#$', str(year))
    data = data.replace('$#UMON#$',umon )
    data = data.replace('$#MONYEAR#$', monyear)

    fo   = open(mon_html, 'w')
    fo.write(data)
    fo.close()
#
#--- read rad_time html page template
#
    data = open('./Template/rad_time_template', 'r').read()

    data = data.replace('$#LMONTH#$', fmon_list[mon-2])
    data = data.replace('$#FYEAR#$', str(year))
    data = data.replace('$#MONYEAR#$', monyear)

    fo   = open(rad_html, 'w')
    fo.write(data)
    fo.close()


#----------------------------------------------------------------------------------------------------------
#--  print_index:html: update index.html page                                                           ---
#----------------------------------------------------------------------------------------------------------

def print_index_html():
    """
    create and/or update radiation related html page
    """
#
#--- find today's date
#
    [year, mon, day, hours, min, sec, weekday, yday, dst] = tcnv.currentTime()
#
#--- read header part and the first part of the index page
#
    line = open('./Template/index_top', 'r').read()
#
#--- start creating a link table
#
    line = line + '<table border=1 cellpadding=10 cellspacing=2>\n'
    line = line + '<tr>\n'
    line = line + '<td colspan=13>\n'
    line = line + '<table border=1 width=100%>\n'
    line = line + '<tr><th> <a href="all.html" style="font-size:120%">Mission since JAN2010</a></th></tr>\n'
    line = line + '</table>\n'
    line = line + '</td>\n'
    line = line + '</tr>\n'

    line = line + '<tr>'
    line = line + '<th>Year</th><th>Jan</th><th>Feb</th><th>Mar</th><th>Apr</th><th>May</th><th>Jun</th>\n'
    line = line + '<th>Jul</th><th>Aug</th><th>Sep</th><th>Oct</th><th>Nov</th><th>Dec</th>\n'
    line = line + '</tr>\n'
    for eyear in range(year, 1999, -1):
        line = line + '<tr>\n' 
        line = line + '<th>' + str(eyear) + '</th>\n'

        lyear    = str(eyear)
        syear    = lyear[2] + lyear[3]

        for emon in range(1, 13):
            lmon     = tcnv.changeMonthFormat(emon)
            monyear  = lmon.lower() + syear
            cmonyear = monyear.upper()

            if eyear == year and emon > mon:
                line = line +  '<td>' + cmonyear + '</td>\n'
            else:
                line = line + '<td><a href="./' + monyear + '.html">'+ cmonyear + '</a></td>\n'

        line = line + '</tr>\n\n'

    line  = line + '</table>\n'
#
#--- table finished. add a closing part
#
    line  = line + '<div style="padding-top: 15px"></div>\n'
    line  = line + '<hr />'
    line  = line + '<div style="padding-top: 15px"></div>\n'
    line  = line + '<p>This page is maintained by B. Spitzbart (<a href="bspitzbart@cfa.harvard.edu">bspitzbart@cfa.harvard.edu</a>).\n'
    line  = line + '</body</html>\n'
#
#--- now write out the page
#
    fo  = open('/data/mta4/www/DAILY/mta_rad/index.html', 'w')
    fo.write(line)
    fo.close()

#----------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
#
#--- if you provide year and month (in format of 2015 3), it will create the html pages for that month. 
#
    if len(sys.argv) == 2:
        year = argv[1]
        mon  = argv[2]
    else:
        year = ''
        mon  = ''
    print_html(year, mon)
#
#--- index page is always written up to this month of this year
#
    print_index_html()
