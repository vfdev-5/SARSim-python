#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      vfomin
#
# Created:     12/06/2015
# Copyright:   (c) vfomin 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# Python
import csv

# Numpy
import numpy as np

def read(filename, delimiter):
    """
    Text data reader and returns numpy array
    """
    data = []
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter)
        for row in reader:
            rowValues = np.array(row).astype(np.uint8)
            data.append(rowValues)
    reader = None
    return np.array(data)


if __name__ == '__main__':

    filename_1 = "C:\\VFomin_folder\\PISE_project\\SAR-GMTI\\Data\\CSK_RAW_0340_B001_Re_0_0_1024_1024.dat"
    out = read(filename_1,'\t')
    print out.shape



    # For binary files
##    import struct
##    filename = "C:\\Users\\vfomin\\Downloads\\CSK_RAW_Hawaii\\raw\\CSK20140105.raw"
##    with open(filename, 'rb') as f:
##        for i in range(10):
##            byte = f.read(1)
##            print struct.unpack('B', byte)