SnowDataFinal
=============

import gdal
import numpy as np
import numpy.ma as ma
import glob

year = '2009'
satellite = 'MOD10A1'

files = np.sort(glob.glob('files/data/MODIS_Snow_Data/%s.A%d*.*hdf'%(satellite,year)))

filename = 'files/data/MODIS_Snow_Data/MOD10A1.A2009001.h09v05.005.2009009120443.hdf'
# Specify one file at this point to show subdatasets
g = gdal.Open(filename)
# g is now a GDAL dataset of the filename
subdatasets = g.GetSubDatasets()

for fname, name in subdatasets:
    print name
    print "\t", fname
    
# From viewing of the subdatasets, it is clear that Fractional_Snow_Cover and Snow_Spatial_QA are wanted

selected_layers = [ "Fractional_Snow_Cover", "Snow_Spatial_QA" ]
# lists are created here with the layer names that have been selected
data = {}
# empty dictionary to store the data form these layers

file_template = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_Snow_500m:%s'
# Template for GDAL filenames where the filename and data layer are substituted where %s occurs.

for i, layer in enumerate ( selected_layers ):
    this_file = file_template % ( filename, layer )
    g = gdal.Open ( this_file )
    data[layer] = g.ReadAsArray() 
# loop counter initialised, files are the index and the layers are wanted for each file

snow = data['Fractional_Snow_Cover']
qc = data['Snow_Spatial_QA']
# The data layers are in seperate numpy arrays
# SnowQA - Good Quality observations are assinged values of 0 and other quality is 1
# snow - Pixel values between 0-100 are classed as snow, this represents a percentage snow cover
# Cloud has a value of 250 in snow
# Tolerance is not considered. The best possible quality of data is wanted.

qc = qc & 1
# where qc bit 1 = 0 is the good data so we find bit 0
snowm = np.ma.array ( snow, mask=qc )
# masked array showing False for all 'good' data values

WRAPPING INTO A FUNCTION:

def read_MODIS_snow(filename):
    g = gdal.Open(filename)
    subdatasets = g.GetSubDatasets()
    selected_layers = [ "Fractional_Snow_Cover", "Snow_Spatial_QA" ]
    data = {}
    file_template = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_Snow_500m:%s'
    for i, layer in enumerate ( selected_layers ):
        this_file = file_template % ( filename, layer )
        g = gdal.Open ( this_file )
        data[layer] = g.ReadAsArray() 
    snow = data['Fractional_Snow_Cover']
    qc = data['Snow_Spatial_QA']
    qc = qc & 1
    snowm = np.ma.array ( snow, mask=qc )
    
    return snow
    
LOOP FOR ALL IMAGES    
# Define year and satellite 
year = 2009
satellite = 'MOD10A1'
# sort through all files
files = np.sort(glob.glob('files/data/MODIS_Snow_Data/%s.A%d*.*hdf'%(satellite,year)))

snow = []
for f in files:
    snow.append(read_MODIS_snow(f))
snow = np.ma.array(snow)
# masked array created of snow cover and qc information for all images












    
    
