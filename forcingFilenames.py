from datetime import datetime

from netCDF4 import MFDataset, date2index, num2date


# Main functions called from model2roms
def get_filename(confM2R, year, month, day, defaultvar):
    """Function returning the filename. 
    All ocean_indata_types should appear here.
    The defaultvar is used in get_time() and in convert_MODEL2ROMS() for grid initialization."""
    if confM2R.ocean_indata_type == "SODA3":
        if defaultvar is None:
            defaultvar = "salinity"
        filenamein = getSODA3filename(confM2R, year, month, day, defaultvar)
    if confM2R.ocean_indata_type == "SODA3_5DAY":
        if defaultvar is None:
            defaultvar = "salinity"
        filenamein = getSODA3_5DAYfilename(confM2R, year, month, day, defaultvar)
    if confM2R.ocean_indata_type == "NORESM":
        if defaultvar is None:
            defaultvar = "grid"
        filenamein = getNORESMfilename(confM2R, year, month, defaultvar)
    if confM2R.ocean_indata_type == "GLORYS":
        if defaultvar is None:
            defaultvar = "So"
        filenamein = get_GLORYS_filename(confM2R, year, month, defaultvar)
    if confM2R.ocean_indata_type == "IPSL":
        if defaultvar is None:
            defaultvar = "so"
        (filenamein, junk) = get_IPSL_filename(confM2R, year, defaultvar)
    if confM2R.ocean_indata_type == "ACCESS":
        if defaultvar is None:
            defaultvar = "so"
        (filenamein, junk) = get_ACCESS_filename(confM2R, year, defaultvar)
    return filenamein

def get_filename_and_yr(confM2R, year, month, day, defaultvar):
    """Function returning filename and start year of file, which is used for time indexing in get_3d_data().
    Only ocean_indata_types with files grouped per couple of years should appear here.
    They also need to appear in get_filename()."""
    if confM2R.ocean_indata_type == "IPSL":
        if defaultvar is None:
            defaultvar = "so"
        (filenamein, file_start_year) = get_IPSL_filename(confM2R, year, defaultvar)
    elif confM2R.ocean_indata_type == "ACCESS":
        if defaultvar is None:
            defaultvar = "so"
        (filenamein, file_start_year) = get_ACCESS_filename(confM2R, year, defaultvar)
    else:
        raise Exception(f"[get_filename_and_yr]: ocean_indata_type {confM2R.ocean_indata_type} not implemented yet.")
        # add new ocean_indata_types in this function (if needed) AND in get_filename() (always).
    return (filenamein, file_start_year)


# private functions called from within module
def get_GLORYS_filename(confM2R, year, month, varname):
    if confM2R.use_zarr:
        return f"{confM2R.ocean_forcing_path}{varname.lower()}"
    return f"{confM2R.ocean_forcing_path}glorys_total.nc"
#    return "{}{}/CMEMS_{}_monthly_MERCATOR_{}-{}.nc".format(
#        confM2R.ocean_forcing_path,
#        varname.lower(),
#        varname.capitalize(),
#        year,
#        str(month).zfill(2),
#    )

def get_IPSL_filename(confM2R, year, myvar):
    if myvar in ['clt', 'pr', 'rsds', 'rsus']: # 'uas', 'vas'
        freq_str = '3hr'
        grid_yr_str = 'gr_201501010130-210012312230'
        file_start_year = 2015
    elif myvar in ['hurs', 'tas']:
        freq_str = 'day'
        grid_yr_str = 'gr_20150101-21001231'
        file_start_year = 2015
    elif myvar in ['psl']:
        freq_str = '6hrPlev'
        grid_yr_str = 'gr_201501010300-210012312100'
        file_start_year = 2015
    elif myvar in ['so', 'thetao', 'uo', 'vo', 'zos']:
        freq_str = 'Omon'
        grid_yr_str = 'gn_201501-210012'
    else:
        raise Exception(f"myvar {str(myvar)} not implemented in get_IPSL_filename()")

    start_yr = 2015  # start year of file, same for all cases above

    filename = (confM2R.ocean_forcing_path
                + myvar
                + "_" + freq_str
                + "_IPSL-CM6A-LR"
                + "_ssp245_r14i1p1f1"
                + "_" + grid_yr_str
                + ".nc"
                )
    return (filename, start_yr)

def get_ACCESS_filename(confM2R, year, myvar):
    if myvar in ['clt', 'pr', 'rsds', 'rsus']:
        freq_str = '3hr'
        [yr_str, start_yr] = fn_10_year_str_atm(year, time_str_start='0130', time_str_end='2230')
    elif myvar in ['tas', 'uas', 'vas']:
        freq_str = '3hr'
        [yr_str, start_yr] = fn_10_year_str_atm(year, time_str_start='0300', time_str_end='0000')
    elif myvar in ['hurs']:
        freq_str = 'day'
        [yr_str, start_yr] = fn_50_year_str(year) # also correct for psl day option
    elif myvar in ['psl']:
        freq_str = '6hrPlevPt'        # or day (files also downloaded as alternative)
        [yr_str, start_yr] = fn_10_year_str_atm(year, time_str_start='0600', time_str_end='0000')
    elif myvar in ['so', 'thetao', 'uo', 'vo']:
        freq_str = 'Omon'
        [yr_str, start_yr] = fn_10_year_str_ocean(year)
    elif myvar in ['zos']:
        freq_str = 'Omon'
        yr_str = '201501-210012'
        start_yr = 2015
    else:
        raise Exception(f"myvar {str(myvar)} not implemented in get_ACCESS_filename()")

    filename = (confM2R.ocean_forcing_path
                + myvar
                + "_" + freq_str
                + "_ACCESS-CM2"
                + "_ssp245_r4i1p1f1"
                + "_gn"
                + "_" + yr_str
                + ".nc"
                )
    return (filename, start_yr)

def fn_50_year_str(year):
    """Useful to get part of filename (fn),
    if files are split in 1 file per 50 year. Based on ACCESS model.
    Also returns start year of the file, which is needed for some vars. """
    if year > 2014 and year < 2065:
        yr_string = '20150101-20641231'
        file_start_year = 2015
    elif year < 2101:
        yr_string = '20650101-21001231'
        file_start_year = 2065
    return (yr_string, file_start_year)

def fn_10_year_str_ocean(year):
    """Useful to get part of filename (fn),
    if files are split in 1 file per 5 years. Based on ACCESS model.
    Version for ocean variables i.e. string contains year & month,
    e.g. '201501-202412'.
    Also returns start year of the file, which is needed for some vars. """
    # round down to nnn5 for start year:
    if year%10 >=5:  # year ends in 5,6,7,8 or 9
        file_start_year = year - year%5            # down to nnn5
    if year%10 < 5:  # year ends in 0,1,2,3 or 4
        file_start_year = year - year%5 - 5        # down to nnn0 and go back 5 more
    file_end_year = min(file_start_year + 9, 2100) # up to nn(n+1)4 or 2100
    yr_string = str(file_start_year) + '01-' + str(file_end_year) + '12'
    return (yr_string, file_start_year)

def fn_10_year_str_atm(year, time_str_start, time_str_end):
    """Useful to get part of filename (fn),
    if files are split in 1 file per 5 years. Based on ACCESS model.
    Version for atmospheric variables i.e. string contains year & month & day & time,
    e.g. '201501010300-202401010000' or '201501010130-202412312230'.
    Also returns start year of the file, which is needed for some vars."""
    # round down to nnn5 for start year:
    if year%10 >=5:  # year ends in 5,6,7,8 or 9
        file_start_year = year - year%5            # down to nnn5
    if year%10 < 5:  # year ends in 0,1,2,3 or 4
        file_start_year = year - year%5 - 5        # down to nnn0 and go back 5 more
    file_end_year = min(file_start_year + 9, 2100) # up to nn(n+1)4 or 2100
    if time_str_end == "0000":
        # file includes New Year's Eve 00:00 so last day of file shifts from e.g. 2034-12-31 to 2035-01-01
        yr_string = str(file_start_year)+'0101'+time_str_start + '-' + str(file_end_year + 1)+'0101'+time_str_end
    else:
        # normal situation; end day of last year is 31 Dec
        yr_string = str(file_start_year)+'0101'+time_str_start + '-' + str(file_end_year)+'1231'+time_str_end
    return (yr_string, file_start_year)

# working but not used/needed:
# def fn_5_year_str(year):
#     """Useful to get part of filename (fn),
#     if files are split in 1 file per 5 years.
#     Also returns start year of the file, which is needed for some vars. """
#     if year%10 < 5:    # year ends in 0,1,2,3 or 4
#         file_start_year = year - year%10    # round down to nnn0
#         file_end_year = file_start_year + 4
#     elif year%10 >=5:  # year ends in 5,6,7,8 or 9
#         file_start_year = year - year%5     # round down to nnn5
#         file_end_year = file_start_year + 4
#     yr_string = str(file_start_year) + '01-' + str(file_end_year) + '12'
#     return (yr_string, file_start_year)


def getNORESMfilename(confM2R, year, month, myvar):
    if myvar == "grid":
        filename = confM2R.ocean_forcing_path + "grid_gx1v6.nc"
    else:
        if myvar in ["iage", "uvel", "vvel", "aice", "hi", "hs"]:
            filename = (
                confM2R.ocean_forcing_path
                + "ICE/NRCP45AERCN_f19_g16_CLE_01.cice.h."
                + str(year)
                + ".nc"
            )

        elif myvar in ["dissic", "talk", "po4", "no3", "si", "o2"]:
            filename = (
                confM2R.ocean_forcing_path
                + "BCG_NRCP85BPRPEX_01.micom.hbgcmlvl.2006-2050.nc"
            )

        elif myvar in ["templvl", "salnlvl", "sealv", "uvellvl", "vvellvl"]:
            if myvar in ["salnlvl", "templvl"]:
                filename = (
                    confM2R.ocean_forcing_path + "TS_NRCP85BPRPEX_01.micom.2006-2100.nc"
                )
            else:
                filename = (
                    confM2R.ocean_forcing_path
                    + "VEL_NRCP85BPRPEX_01.micom.hmlvl.2006-2100.nc"
                )

        elif myvar in ["TAUX", "TAUY", "U10"]:    # atmospheric variables
            filename = (
                confM2R.atmospheric_forcing_path
                + f"NRCP45AERCN_f19_g16_CLE_01.cam2.h5.{year}-{month:02d}-{day:02d}-00000.nc"
            )
    return filename


def getSODA3filename(confM2R, year, month, day, myvar):
    if myvar in ["cn", "hi", "hs"]:
        return "{}soda{}_mn_ice_reg_{}.nc".format(
            confM2R.ocean_forcing_path, confM2R.soda_version, str(year)
        )
    else:
        return "{}soda{}_mn_ocean_reg_{}.nc".format(
            confM2R.ocean_forcing_path, confM2R.soda_version, str(year)
        )


def getSODA3_5DAYfilename(confM2R, year, month, day, myvar):
    if len(confM2R.time_object) == 0:
        mcdf = MFDataset(confM2R.ocean_forcing_path + "*.nc")
        confM2R.time_object = mcdf.variables["time"]

    index = date2index(
        datetime(year, month, day, 0, 0),
        confM2R.time_object,
        calendar=confM2R.time_object.calendar,
        select="nearest",
    )
    seldate = num2date(
        confM2R.time_object[index],
        units=confM2R.time_object.units,
        calendar=confM2R.time_object.calendar,
    )

    if myvar in ["cn", "hi", "hs"]:
        return "{}soda{}_5dy_ocean_ice_{:04}_{:02}_{:02}.nc".format(
            confM2R.ocean_forcing_path,
            confM2R.soda_version,
            seldate.year,
            seldate.month,
            seldate.day,
        )
    else:
        return "{}soda{}_5dy_ocean_reg_{:04}_{:02}_{:02}.nc".format(
            confM2R.ocean_forcing_path,
            confM2R.soda_version,
            seldate.year,
            seldate.month,
            seldate.day,
        )
