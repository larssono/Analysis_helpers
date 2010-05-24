# Time-stamp: <20090616 15:09:31 CB\al653@al653rhel-102d.cb.bscb.cornell.edu>
#
# Python code to read fMRI data from AFNI's BRIK format.
# Andre Lehovich, June 2009.
#
# The file README.attributes distributed with AFNI has documentation
# on the BRIK header format.
#
# Example invocation:
# from afni import brik
# x = brik("e004833_r1_po_psc_at3+tlrc")
# x.header['DATASET_DIMENSIONS']
# pylab.figure(); pylab.imshow(x.data[:,:,20,1])
# pylab.figure(); pylab.plot(x.data[27,32,25,:])

import numpy

# Dictionary translating header field BRICK_TYPES into numpy dtypes.
# Will throw a KeyError if an unimplemented type is used.
bricktypes = {
    0 : "i1",
    1 : "i2",
    3 : "f4"
    }

# Dictionary translating header field BYTEORDER_STRING into numpy dtypes.
byteorders = {
    "MSB_FIRST" : ">", # Big endian (e.g. PPC Macs)
    "LSB_FIRST" : "<", # Little endian (e.g. x86)
    }

class brik:
    """Class to represent an AFNI data BRIK"""
    def __init__(self, basename):
        self.basename = basename
        self.header = read_header(basename + ".HEAD")
        self.data = read_data(basename + ".BRIK", self.header)
        
def read_header(filename):
    """Read the header from an AFNI data BRIK"""
    result = dict()
    f = open(filename, "r")
    while 1:
        line = f.readline()
        if not line: break
        words = line.split()
        if not ((len(words) == 3) and (words[0] == "type")): continue
        type = words[2]
        name = f.readline().split()[2]
        count = int(f.readline().split()[2])
        if (type == 'string-attribute'):
            # AFNI header strings begin with ' and end with ~
            val = f.readline().strip()[1:-1]
        elif (type == 'float-attribute'):
            val = numpy.fromfile(f, dtype=numpy.float, count=count, sep=' ')
        elif (type == 'integer-attribute'):
            val = numpy.fromfile(f, dtype=numpy.int, count=count, sep=' ')
        else:
            print >>sys.stderr, "ERROR: Unknown type:", type, "in", filename
            raise "UnknownType"
        result[name] = val
    return result

def read_data(filename, header):
    """Read the raw data from an AFNI data BRIK"""
    if [a for a in header['BRICK_TYPES'] if a != header['BRICK_TYPES'][0]]:
        print >>sys.stderr, "ERROR: different types of sub-bricks in", filename
        raise "DifferentTypes"
    dtype= byteorders[header['BYTEORDER_STRING']]
    dtype += bricktypes[header['BRICK_TYPES'][0]]
    data = numpy.fromfile(file=filename, dtype=dtype)
    x,y,z = header['DATASET_DIMENSIONS'][0:3]
    nval = header['DATASET_RANK'][1]
    return data.reshape(x, y, z, nval, order='F')
