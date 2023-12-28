import calendar
import logging
import sys
from datetime import datetime

from netCDF4 import num2date

__author__ = 'Trond Kristiansen'
__email__ = 'me@trondkristiansen.com'
__created__ = datetime(2015, 8, 11)
__modified__ = datetime(2018, 4, 25)
__version__ = "1.5"
__status__ = "Development, modified on 19.04.2018 by A.Barthel, " \
             "refactored 25.04.2018 by Trond Kristiansen"


# Methods for returning list of months and days for the given time-step.

def create_list_of_months(confM2R, currentyear):

    if currentyear == confM2R.startdate.year:
        IDS = [confM2R.startdate.month + m for m in range(13 - confM2R.startdate.month)]
        # months from start month to end of first year

    elif currentyear == confM2R.enddate.year:
        IDS = [1 + m for m in range(int(confM2R.enddate.month))]
        # months from first month to last month of last year

    elif confM2R.startdate.year < currentyear < confM2R.enddate.year:
        # months from first month to last month of last year
        IDS = [1 + m for m in range(12)]

    if confM2R.startdate.year == confM2R.enddate.year:
        # months from first month to last month of last year
        IDS = [confM2R.startdate.month + m for m in range(confM2R.enddate.month - confM2R.startdate.month + 1)]

    if confM2R.isclimatology:
        IDS = [i + 1 for i in range(12)]

    if not IDS:
        logging.info("Unable to generate IDS for time looping: main.py -> func createIDS")
        sys.exit()

    return IDS


def create_list_of_days(confM2R, year, month, first_run):
    days = []
    if confM2R.time_frequency_inputdata == 'day':
        daystep = 7
        if daystep > 1:
            logging.info("WARNING!")
            logging.info("----------------------------------------------------------------------")
            logging.info("You are only using every {} days of input data! (model2roms.py)".format(daystep))
            logging.info("----------------------------------------------------------------------")

        # Regularly we want all days in each month
        ndays = int(calendar.monthrange(year, month)[1])
        days = [d + 1 for d in range(0, ndays, daystep)]
        if first_run:
            for day in days:
                if confM2R.start_day>day:
                    days.remove(day)
                
        # Exceptions:
        # We start in the first month after day one
        if month == confM2R.startdate.month and year == confM2R.startdate.year and confM2R.startdate.day > 1:
            days = [i for i in range(confM2R.startdate.day, ndays, daystep)]

        # We finish in the last month before last day of month
        if month == confM2R.enddate.month and year == confM2R.enddate.year:
            days = [i + 1 for i in range(0, confM2R.enddate.day, daystep)]

        # We start and end on different days in the same month         
        if confM2R.startdate.month == confM2R.enddate.month and confM2R.startdate.year == confM2R.enddate.year:
            days = [i for i in range(confM2R.startdate.day, confM2R.enddate.day, daystep)]
    
    if confM2R.time_frequency_inputdata == '5days':
        
        for dd in confM2R.time_object[:]:
            dd_date = num2date(dd, units=confM2R.time_object.units, calendar=confM2R.time_object.calendar)

            if dd_date.year==year and dd_date.month==month:
                days.append(dd_date.day)
     
        if first_run:
            days_tmp=days.copy()
            for day in days:
                if day < confM2R.start_day:
                    days_tmp.remove(day)
                    print("Removing day {} from days {}".format(day, days_tmp))
            days=days_tmp 
        return days

    if confM2R.time_frequency_inputdata == 'month':
        days = [15]

    return days
