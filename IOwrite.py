
from netCDF4 import Dataset


def open_output(grdROMS,ntimes):
        
     f1 = Dataset(outfile, mode='w', format='NETCDF3_CLASSIC')
     
     # Define dimensions
     f1.createDimension('lon', len(grdROMS.lon))
     f1.createDimension('lat', len(grdROMS.lat))
     f1.createDimension('time', 1)
     f1.createDimension('s_rho', grdROMS.s_rho)
     f1.createDimension('s_w', grdROMS.s_w)
     
     # Define coordinate variables
     v = f1.createVariable('lon', 'd', ('lon',))
     v.long_name = 'Longitude'
     v.units = 'degrees east'
     v[:] = grdROMS.lon
     
     
     v = f1.createVariable('lat', 'd', ('lat',))
     v.long_name = 'Latitude'
     v.units = 'degrees north'
     v[:] = grdROMS.lat
     
     v = f1.createVariable('ocean_time', 'd', ('ocean_time',))
     v.long_name = 'ocean_time'
     #v.units = v0.units
     v[:] = grdROMS.ocean_time
     
     v = f1.createVariable('s_rho', 'd', ('s_rho',))
     v0 = grdROMs.s_rho
     v[:] = v0[:]
     
     varTZYX=['temp']
     
  
     
     for var in varTZYX:
         print var
         v0 = grdROMS.t
         v1 = f1.createVariable(var, 'f', ('time', 's_rho', 'lat', 'lon'))
         try:
             v1.long_name = 'Temperature test'
             v1.units = 'degrees Celsius'
         except:
             pass
         v1._FillValue = -9.99E-33
         
         for l in xrange(ntimes):
             for k in xrange(grdROMS.s_rho):
                 print "k, l = ", k, l
                 F0 = np.array(v0[l,k,:,:])
                
               
                 v1[l,k,:,:] = F0.astype('float')
     
     f1.close()


