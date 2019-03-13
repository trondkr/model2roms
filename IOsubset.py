import numpy as np
from datetime import datetime

__author__ = 'Trond Kristiansen'
__email__ = 'me@trondkristiansen.com'
__created__ = datetime(2010, 1, 7)
__modified__ = datetime(2013, 7, 1)
__version__ = "0.1"
__status__ = "Development, modified on 07.01.2010, 01.07.2013, 11.03.2014"


def findSubsetIndices(grdMODEL, min_lat, max_lat, min_lon, max_lon):
    """
    Get the indices that covers the new grid, and enables us to only store a subset of
    the large input grid. This routine checks wether you cross the 0 longitude line, so
    that start is in the West (negative longitudes) and end point is in the East (positive longitudes).
    if that is the case, the extraction of data is split at the 0 line so that we first get all
    data in the Western region before all the data in the Eastern region is taken. Then we combine the two arrays
    so that West comes before East. This is not the case in an array where the longitude values goes from 0 to 360. In such an
    array (as SODA is), the Westernmost data comes after the Easternmost. Therefore we need to concentate
    the Westernomst data into an array with the Easternomst data in the correct order.
    
    Input: Input longitudes has to range from -180 to 180
    
    Returns: list of indices, boolean value that is True of splitting occurs (crosses 0 long)
    Indices are returned so that the Westernomst comes first
    
    Trond Kristiansen, 22.07.2009, 08.01.2010, 01.07.2013, 11.03.2014
    """

    if min_lon < 0 and max_lon > 0:
        splitExtract = True;
        Turns = 2
        grdMODEL.splitExtract = splitExtract
    else:
        splitExtract = False;
        Turns = 1
        grdMODEL.splitExtract = splitExtract
        grdMODEL.lon = np.where(grdMODEL.lon > 180, grdMODEL.lon - 360, grdMODEL.lon)

    # Array to store the results returned from the function
    res = np.zeros((Turns, 4), dtype=np.float)

    lats = grdMODEL.lat[:, 0]
    lons = grdMODEL.lon[0, :]

    for k in range(Turns):

        if k == 0 and splitExtract == True:
            minLon = min_lon;
            maxLon = 0
            minLon = minLon + 360
            maxLon = maxLon + 360
        elif k == 1 and splitExtract == True:
            minLon = 0;
            maxLon = max_lon
        else:
            minLon = min_lon;
            maxLon = max_lon

        distances1 = []
        distances2 = []
        indices = []
        index = 1
        for point in lats:
            s1 = max_lat - point  # (vector subtract)
            s2 = min_lat - point  # (vector subtract)
            distances1.append((np.dot(s1, s1), point, index))
            distances2.append((np.dot(s2, s2), point, index - 1))
            index = index + 1

        distances1.sort()
        distances2.sort()
        indices.append(distances1[0])
        indices.append(distances2[0])

        distances1 = []
        distances2 = []
        index = 1

        for point in lons:
            s1 = maxLon - point  # (vector subtract)
            s2 = minLon - point  # (vector subtract)
            distances1.append((np.dot(s1, s1), point, index))
            distances2.append((np.dot(s2, s2), point, index - 1))
            index = index + 1

        distances1.sort()
        distances2.sort()
        indices.append(distances1[0])
        indices.append(distances2[0])

        # Save final product: max_lat_indices,min_lat_indices,max_lon_indices,min_lon_indices
        minJ = indices[1][2]
        maxJ = indices[0][2]
        minI = indices[3][2]
        maxI = indices[2][2]

        res[k, 0] = minI;
        res[k, 1] = maxI;
        res[k, 2] = minJ;
        res[k, 3] = maxJ;

    # Save final product: max_lat_indices,min_lat_indices,max_lon_indices,min_lon_indices
    grdMODEL.indices = res


def checkDomain(grdMODEL, grdROMS):
    lonCHECK = False
    latCHECK = False
    if (grdMODEL.lon.min() <= grdROMS.lon_rho.min() and grdMODEL.lon.max() >= grdROMS.lon_rho.max()):
        lonCHECK = True
    if (grdMODEL.lat.min() <= grdROMS.lat_rho.min() and grdMODEL.lat.max() >= grdROMS.lat_rho.max()):
        latCHECK = True

    print("\n--------------------------")
    print(("---> Area output files  : (longitude=%3.2f,latitude=%3.2f) to (longitude=%3.2f,latitude=%3.2f)" % (
    grdROMS.lon_rho.min(), grdROMS.lat_rho.min(), grdROMS.lon_rho.max(), grdROMS.lat_rho.max())))
    print(("---> Area forcing files : (longitude=%3.2f,latitude=%3.2f) to (longitude=%3.2f,latitude=%3.2f)" % (
    grdMODEL.lon.min(), grdMODEL.lat.min(), grdMODEL.lon.max(), grdMODEL.lat.max())))

    if latCHECK is True and lonCHECK is True:
        print("Domain check passed: Input domain data covers output domain")
    else:
        print("WARNING: Your input domain is smaller or not overlaying your output domain")
        print("IOsubset.py: EXIT")
        # sys.exit()
    print("\n--------------------------")


def organizeSplit(grdMODEL, grdROMS):
    if grdMODEL.splitExtract is True:
        print("\nThe subset crosses the 0 degrees longitude line.")
        print("Need to split the extraction of data into Western and")
        print("Eastern region. The two regions will be concatenated again.\n")

        """Save longitude and latitude values so that only variables are needed to be extracted after first loop"""
        lon1 = grdMODEL.lon[int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
               int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]
        lat1 = grdMODEL.lat[int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
               int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]

        lon2 = grdMODEL.lon[int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
               int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])]
        lat2 = grdMODEL.lat[int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
               int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])]

        # Convert the positive numbers above 180 to negative values
        lon1 = np.where(lon1 > 180, lon1 - 360, lon1)

        grdMODEL.lon = np.concatenate((lon1, lon2), axis=1)
        grdMODEL.lat = np.concatenate((lat1, lat2), axis=1)

        print(("Subset extracted for domain from West (%s) to East (%s)\n" % (grdMODEL.lon.min(), grdMODEL.lon.max())))
    else:
        grdMODEL.lon = grdMODEL.lon[int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                       int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]
        grdMODEL.lat = grdMODEL.lat[int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                       int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]

    checkDomain(grdMODEL, grdROMS)
