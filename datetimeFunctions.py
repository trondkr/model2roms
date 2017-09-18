from netCDF4 import Dataset, datetime, date2num,num2date
import os
import calendar

__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@imr.no'
__created__ = datetime(2015, 8, 11)
__modified__ = datetime(2015, 8, 11)
__version__ = "1.5"
__status__ = "Development, modified on 11.08.2015"

# Methods for returning list of months and days for the given time-step.

def createListOfMonths(currentyear,startdate,enddate,isClimatology):

    print currentyear, startdate.year, enddate.year
    if currentyear == startdate.year:
        IDS=[startdate.month+m for m in xrange(13-startdate.month)]
        # months from start month to end of first year
       
    elif currentyear == enddate.year:
        IDS=[1+m for m in xrange(int(enddate.month))]
        # months from first month to last month of last year
        
    elif startdate.year < currentyear < enddate.year:
        # months from first month to last month of last year
        IDS = [1+m for m in xrange(12)]
        
    if startdate.year == enddate.year:
        # months from first month to last month of last year
        IDS = [startdate.month + m for m in xrange(enddate.month - startdate.month)]
        
    print "Months for year %s : %s"%(currentyear,IDS)

    if isClimatology==True:
        IDS = [i+1 for i in xrange(12)]
      
    if not IDS:
        print "Unable to generate IDS for time looping: main.py -> func createIDS"
        sys.exit()
  
    return IDS

def createListOfDays(year,month,startdate,enddate,isClimatology,timeFrequencyOfInputData):
    days=[]
    if timeFrequencyOfInputData == 'day':
        dayStep=7
        if dayStep > 1:
            print "WARNING!"
            print "----------------------------------------------------------------------"
            print "You are only using every %s days of input data! (model2roms.py)"%(dayStep)
            print "----------------------------------------------------------------------"
        
        # Regulary we want all days in each month                
        ndays = int(calendar.monthrange(year, month)[1])
        days = [d+1 for d in xrange(0,ndays,dayStep)]
        
        # Exceptions:
        # We start in the first month after day one
        if ( month == startdate.month and  year == startdate.year and startdate.day > 1): 
            days = [i for i in xrange(startdate.day,ndays,dayStep)]
       
        # We finish in the last month before last day of month
        if ( month == enddate.month and  year == enddate.year ): 
            days = [i+1 for i in xrange(0,enddate.day,dayStep)]
        
        # We start and end on different days in the same month         
        if ( startdate.month == enddate.month and  startdate.year == enddate.year ): 
            days = [i for i in xrange(startdate.day,enddate.day,dayStep)]
        print "days", days

    if timeFrequencyOfInputData == 'month':
        days = [15]

    return days