'''
Derived from freemol.mpeg_encode

Copyright (c) Schrodinger, Inc.
'''


def validate():
    '''
    Return the path to the "mpeg_encode" executable or None if not found
    '''
    from shutil import which
    return which('mpeg_encode') or which('mpeg_encode.exe')


INPUT_TEMPLATE = """
# parameter file template with lots of comments to assist you
#
# you can use this as a template, copying it to a separate file then modifying
# the copy
#
#
# any line beginning with '#' is a comment
#
# no line should be longer than 255 characters
#
#
# general format of each line is:
#	<option> <spaces and/or tabs> <value>
#
# lines can generally be in any order
#
# an exception is the option 'INPUT' which must be followed by input
# files in the order in which they must appear, followed by 'END_INPUT'
#
# Also, if you use the `command` method of generating input file names,
# the command will only be executed in the INPUT_DIR if INPUT_DIR preceeds
# the INPUT parameter.
#
# <option> MUST be in UPPER CASE
#

#PATTERN		IPBBBPBBBB
#PATTERN		IPBBPBBPBB
#PATTERN		IPBPBPBPBP
PATTERN		    IPPPPPPPPP

FORCE_ENCODE_LAST_FRAME

OUTPUT		%s

# mpeg_encode really only accepts 3 different file formats, but using a
# conversion statement it can effectively handle ANY file format
#
# You must specify the type of the input files.  The choices are:
#    YUV, PPM, JMOVIE, Y, JPEG, PNM
#	(must be upper case)
#
BASE_FILE_FORMAT	PPM

#
# if YUV format (or using parallel version), must provide width and height
# YUV_SIZE	widthxheight
# this option is ignored if BASE_FILE_FORMAT is not YUV and you're running
# on just one machine
#
# YUV_SIZE	352x240

# If you are using YUV, there are different supported file formats.
# EYUV or UCB are the same as previous versions of this encoder.
# (All the Y's, then U's then V's, in 4:2:0 subsampling.)
# Other formats, such as Abekas, Phillips, or a general format are
# permissible, the general format is a string of Y's, U's, and V's
# to specify the file order.

#INPUT_FORMAT UCB

# the conversion statement
#
# Each occurrence of '*' will be replaced by the input file
#
# e.g., if you have a bunch of GIF files, then this might be:
#	INPUT_CONVERT	giftoppm *
#
# e.g., if you have a bunch of files like a.Y a.U a.V, etc., then:
#	INPUT_CONVERT	cat *.Y *.U *.V
#
# e.g., if you are grabbing from laser disc you might have something like
#	INPUT_CONVERT	goto frame *; grabppm
# 'INPUT_CONVERT *' means the files are already in the base file format
#
INPUT_CONVERT	*

# number of frames in a GOP.
#
# since each GOP must have at least one I-frame, the encoder will find the
# the first I-frame after GOP_SIZE frames to start the next GOP
#
# later, will add more flexible GOP signalling
#
GOP_SIZE	10

# number of slices in a frame
#
# 1 is a good number.  another possibility is the number of macroblock rows
# (which is the height divided by 16)
#
SLICES_PER_FRAME	1

# directory to get all input files from (makes this file easier to read)
INPUT_DIR 	%s

# There are a bunch of ways to specify the input files.
# from a simple one-per-line listing, to the following 
# way of numbering them.  See the manual for more information.
INPUT
# '*' is replaced by the numbers 01, 02, 03, 04
# if I instead do [01-11], it would be 01, 02, ..., 09, 10, 11
# if I instead do [1-11], it would be 1, 2, 3, ..., 9, 10, 11
# if I instead do [1-11+3], it would be 1, 4, 7, 10
# the program assumes none of your input files has a name ending in ']'
# if you do, too bad!!!
#
#
%s*.ppm 	[%04d-%04d]

# can have more files here if you want...there is no limit on the number
# of files
END_INPUT

# Many of the remaining options have to do with the motion search and qscale

# FULL or HALF -- must be upper case
PIXEL		HALF

# means +/- this many pixels for both P and B frame searches
# specify two numbers if you wish to serc different ranges in the two.
RANGE		16

# this must be one of {EXHAUSTIVE, SUBSAMPLE, LOGARITHMIC}
PSEARCH_ALG	LOGARITHMIC

# this must be one of {SIMPLE, CROSS2, EXHAUSTIVE}
#
# note that EXHAUSTIVE is really, really, really slow
#
BSEARCH_ALG	CROSS2

#
# these specify the q-scale for I, P, and B frames
# (values must be between 1 and 31)
# These are the Qscale values for the entire frame in variable bit-rate
# mode, and starting points (but not important) for constant bit rate
#
IQSCALE		%d
PQSCALE		%d
BQSCALE		%d

# this must be ORIGINAL or DECODED
REFERENCE_FRAME	DECODED

# for parallel parameters see parallel.param in the exmaples subdirectory

# if you want constant bit-rate mode, specify it as follows (number is bits/sec):
#BIT_RATE  8000000

# To specify the buffer size (327680 is default, measused in bits, for 16bit words)
#BUFFER_SIZE 327680

# The frame rate is the number of frames/second (legal values:
# 23.976, 24, 25, 29.97, 30, 50 ,59.94, 60
FRAME_RATE 30

# There are many more options, see the users manual for examples....
# ASPECT_RATIO, USER_DATA, GAMMA, IQTABLE, etc.


"""


def input(filename, path, prefix, first, last, quality=12):
    '''
    Return a parameter file (file contents).
    '''
    return INPUT_TEMPLATE % (
        filename,
        path,
        prefix,
        first,
        last,
        min(quality, 31),
        min(quality + 1, 31),
        min(quality + 2, 31),
    )


def run(paraminput):
    '''
    Run "mpeg_encode"

    :param paraminput: Parameter file contents
    :return: (stdout, stderr) bytes tuple
    '''
    from subprocess import Popen, PIPE
    mpeg_encode_exe = validate()

    if not mpeg_encode_exe:
        raise FileNotFoundError('Cannot find "mpeg_encode"')

    if not isinstance(paraminput, bytes):
        paraminput = paraminput.encode()

    p = Popen([mpeg_encode_exe], stdin=PIPE, stdout=PIPE, stderr=PIPE)

    return p.communicate(paraminput)
