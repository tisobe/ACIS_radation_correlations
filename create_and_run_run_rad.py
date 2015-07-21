#!/usr/bin/env /proj/sot/ska/bin/python

#############################################################################################################
#                                                                                                           #
#       create_and_run_run_rad.py: create and run idl script today's "run_rad.pro"A                         #
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


idl_path = '+/data/mta/Script/Ephem/Scripts/:+/usr/local/rsi/user_contrib/astron_Oct09/pro:+/home/mta/IDL:/home/nadams/pros:+/data/swolk/idl_libs:/home/mta/IDL/tara:widget_tools:utilities:event_browser'


#----------------------------------------------------------------------------------------------------------
#-- print_run_rad: create today's run_rad.pro idl script                                                ---
#----------------------------------------------------------------------------------------------------------

def print_run_rad():
    """
    create today's run_rad.pro idl script
    """
#
#--- find today's date
#
    [year, mon, day, hours, min, sec, weekday, yday, dst] = tcnv.currentTime()
#
#---if this is the first of the month, compute the last month
#
    if day < 2:
        subtract = 1
    else:
        subtract = 0
    cyear = year
    lmon      = mon - subtract
    if lmon < 1:
        lmon  = 12
        cyear = year -1
#
#--- choose a correct month list depending on whether this is the leap year
#
    if tcnv.isLeapYear(cyear) == 1:
        mon_list = mon_list1
    else:
        mon_list = mon_list2

    if lmon == 1:
        start_day = '001'
    else:
        sday = float(mon_list[lmon-2]) + 1
        start_day = str(sday)
        if sday < 10:
            start_day = '00' + start_day
        elif sday < 100:
            start_day = '0' + start_day

    last_day  = mon_list[lmon-1]
#
#--- convert the month from a numeric to letter
#
    smon      = tcnv.changeMonthFormat(lmon)
    smon      = smon.lower()
#
#--- set date format to the appropriate ones
#
    today     =  str(yday)
    if yday < 100:
        today = '0' + today
    elif today < 10:
        today = '00' + today

    day30     = yday - 30
    if day30 < 0:
        day30 = 366 + day30
        tday_year = str(year -1)
    else:
        tday_year =  str(year)

    if day30 < 10:
        day30 = '00' + str(day30)
    elif day30 < 100:
        day30 = '0' + str(day30)
    else:
        day30 = str(day30)


    lmon_year =  str(cyear)
    syear     =  lmon_year[2] + lmon_year[3]
    last_year =  str(year -1)
    syear2    =  last_year[2] + last_year[3]
    monyear   =  smon + syear
#
#--- now print out the run_rad.pro
#
    line = "x = mta_rad('" + tday_year + ":" + str(day30) +"')\n"
    line = line + "spawn, 'mv rad_cnts.gif rad_cnts_curr.gif'\n"
    line = line + "spawn, 'mv rad_use.gif rad_use_curr.gif'\n"
    line = line + "print, 'Done current'\n"
    line = line + "retall\n"
    
    line = line + "x = mta_rad('" + lmon_year + ":" + start_day + "','" + lmon_year + ":" + last_day + "')\n"
    line = line + "spawn, 'mv rad_cnts.gif rad_cnts_" + monyear + ".gif'\n"
    line = line + "spawn, 'mv rad_use.gif rad_use_"   + monyear + ".gif'\n"
    line = line + "spawn, 'mv eph_diff.gif eph_diff_" + monyear + ".gif'\n"
    line = line + "spawn, 'mv mon_diff.gif mon_diff_" + monyear + ".gif'\n"
    line = line + "spawn, 'mv per_diff.gif per_diff_" + monyear + ".gif'\n"
    line = line + "spawn, 'mv xper_diff.gif mon_per_diff_" + monyear + ".gif'\n"
    line = line + "print, 'Done " + monyear + "'\n"
    
    line = line + "x = mta_rad('" + last_year + ":001','" + lmon_year + ":001')\n"
    line = line + "spawn, 'mv rad_cnts.gif rad_cnts_" + syear2 + ".gif'\n"
    line = line + "spawn, 'mv rad_use.gif rad_use_"   + syear2 + ".gif'\n"
    line = line + "spawn, 'mv eph_diff.gif eph_diff_" + syear2 + ".gif'\n"
    line = line + "spawn, 'mv mon_diff.gif mon_diff_" + syear2 + ".gif'\n"
    line = line + "spawn, 'mv per_diff.gif per_diff_" + syear2 + ".gif'\n"
    line = line + "spawn, 'mv xper_diff.gif mon_per_diff_" + syear2 + ".gif'\n"
    line = line + "print, 'Done " + syear2 + "'\n"
    
    line = line + "x = mta_rad('" + last_year + ":" + today + "','" + lmon_year + ":" + today + "')\n"
    line = line + "spawn, 'mv rad_cnts.gif rad_cnts_last_one_year.gif'\n"
    line = line + "spawn, 'mv rad_use.gif rad_use_last_one_year.gif'\n"
    line = line + "spawn, 'mv eph_diff.gif eph_diff_last_one_year.gif'\n"
    line = line + "spawn, 'mv mon_diff.gif mon_diff_last_one_year.gif'\n"
    line = line + "spawn, 'mv per_diff.gif per_diff_last_one_year.gif'\n"
    line = line + "spawn, 'mv xper_diff.gif mon_per_diff_last_one_year.gif'\n"
    line = line + "print, 'Done Last One Year\n"
    
    line = line + "x = mta_rad()\n"
    line = line + "spawn, 'mv rad_cnts.gif rad_cnts_all.gif'\n"
    line = line + "spawn, 'mv rad_use.gif rad_use_all.gif'\n"
    line = line + "spawn, 'mv eph_diff.gif eph_diff_all.gif'\n"
    line = line + "spawn, 'mv mon_diff.gif mon_diff_all.gif'\n"
    line = line + "spawn, 'mv per_diff.gif per_diff_all.gif'\n"
    line = line + "print, 'Done all'\n"
    line = line + "retall\n"
    line = line + "exit\n"

#    fo   = open('test.pro', 'w')
    fo   = open('run_rad.pro', 'w')
    fo.write(line)
    fo.close()

#----------------------------------------------------------------------------------------------------------
#-- run_idl_script: run idl script "run_rad.pro"                                                        ---
#----------------------------------------------------------------------------------------------------------

def run_idl_script():
    """
    run idl script "run_rad.pro"
    """
    os.environ['IDL_PATH'] = idl_path
    subprocess.call("idl run_rad.pro",  shell=True)
#    cmd = "idl run_rad.pro"
#    os.system(cmd)


#----------------------------------------------------------------------------------------------------------

if __name__ == "__main__":

    print_run_rad()

#    run_idl_script()
#    os.system("rm -rf ./ccd*.dat")
