#!/usr/bin/env python

"""xvg2hdf.py: Takes a path to a *.xvg and *.sle file and outputs a *.h5 HDF5
database file in the same directory."""

# Import built-in modules.
import sys
import os.path
import os
from time import time
# Import third-party modules.
from tables import openFile, IsDescription, UInt8Col, FloatCol
from numpy import array

def xvg_gen_maker(xvg_file):
    """Makes a generator function that generates individual parsings of one
    complete snapshot of atoms. This method is used to process the very large
    *.xvg files without having to load the whole file into memory."""
    # Internal index used for naming snapshots.
    i = 1
    coord_line = False
    # Initiate list that holds atomic coordinate lines.
    coord_list = []
    # Read line-by-line to minimize resource usage.
    line = xvg_file.readline()
    while line:
        # If the line consists of a hash and 80 dashes check the next line.
        if line.strip() == '#'+80*'-':
            line = xvg_file.readline()
            # If the next line consists of a hash and 80 dashes too then the
            # line below this one should contain atomic coordinates. Switch
            # the coord_line boolean to "True", read the next line, and go to
            # the beginning of the loop.
            if line.strip() == '#'+80*'-':
                coord_line = True
                line = xvg_file.readline()
                continue
            # If the next line does not consist of a hash and 80 dashes then
            # the next line is not an atomic coordinate. Switch the coord_line
            # boolean to "False". If the previous line was an atomic coordinate
            # then the boolean is changed. If the previous line was not an
            # atomic coordinate, then the boolean is not changed (it was
            # "False" and is set to "False" again).  Read the next line and go
            # to the beginning of the loop.
            else:
                coord_line = False
                line = xvg_file.readline()
                continue
        # If the coord_line boolean is "True" append the line to the
        # coord_list.
        if coord_line:
            coord_list.append(line)
            # Read ahead one line to see if it is the end of the atomic
            # coordinates list. If so, return snapshot index and coord_list,
            # reset the coord_list to empty, and increase the index by one,
            # then go to the beginning of the loop. NOTICE: the delimiter for
            # the end of the atomic coordinates list is one line that consists
            # of a hash and 80 dashes. Reaching the end of the line, we do not
            # read another line from the *.xvg file. We go directly to the
            # beginning of the loop. Since the next line is not a hash and 80
            # dashes the above code sets coord_line to "False". Also the *.xvg
            # files abruptly end without a nice hash and 80 dashes at the end
            # of the file so we also test for the EOF. In python EOF is
            # represented by the return of a zero length string from the
            # readline() function.
            line = xvg_file.readline()
            if (line.strip() == '#'+80*'-') or (len(line) == 0):
                yield (i, coord_list)
                coord_list = []
                i += 1
            continue
        # If the line is neither a hash and 80 dashes, nor is coord_line "True"
        # then the line must be some kind of other output that we're not
        # interested in (for now). Read another line.
        line = xvg_file.readline()

def sle_parser(sle_file):
    sle_list = sle_file.readlines()
    box_size = sle_list[1].split()
    print box_size
    return (float(box_size[3]), float(box_size[4]), float(box_size[5]))

class SimParams(IsDescription):
    """Defines a HDF5 table that encapsulates the parameters used for the whole
    MD trajectory simulation."""
    # Double-precision MD simulation box size.
    box_x = FloatCol()
    box_y = FloatCol()
    box_z = FloatCol()

class Snapshot(IsDescription):
    """Defines a HDF5 table that encapsulates a single snapshot of an MD
    simulation."""
    # Unsigned byte.
    atomic_number = UInt8Col(pos=0)
    # Double-precision xyz atom position.
    pos = FloatCol(shape=3, pos=1)
    # Double-precision px, py ,pz atom momentum.
    momentum = FloatCol(shape=3, pos=2)

if __name__ == "__main__":
    xvg_file = sys.argv[1]
    # Properly parse input file path
    pathname, full_filename = os.path.split(xvg_file)
    filename, file_ext = os.path.splitext(full_filename)

    os.chdir(pathname)
    # Make a generator out of the xvg file.
    xvg = xvg_gen_maker(open(full_filename, 'r'))
    # Start timing total conversion time.
    start_time = time()

    # Open a HDF5 file for writing
    h5_file = openFile(filename + '.h5', mode='w', title=filename)
    # Create a new group under '/' (root)
    calc_meta_group = h5_file.createGroup('/', 'calculation_input',
        "MD Calculation Input")
    md_parameters = h5_file.createTable('/calculation_input', 'md_params',
            SimParams)
    parameter_values = md_parameters.row
    box_params = sle_parser(open(filename + '.sle', 'r'))
    parameter_values['box_x'] = box_params[0]
    parameter_values['box_y'] = box_params[1]
    parameter_values['box_z'] = box_params[2]
    parameter_values.append()
    md_parameters.flush()

    simulation_group = h5_file.createGroup('/', 'simulation',
        "Full MD simulation (all snapshots)")

    for id_num, snapshot_list in xvg:
        atom_num = []
        pos_coords = []
        momentum = []

        for line in snapshot_list:
            columns = line.split()
            atom_num.append(columns[0])
            pos_coords.append(columns[1])
            momentum.append(columns[2])

        atom_num = atom_num[::3]
        pos_array = array(pos_coords).reshape(-1,3)
        momentum_array = array(momentum).reshape(-1,3)

        snapshot_name = "SS_%010i" % (id_num)

        snapshot = h5_file.createTable('/simulation', snapshot_name, Snapshot)
        atom = snapshot.row

        for i in range(len(atom_num)):
            atom['atomic_number'] = atom_num[i]
            atom['pos'] = pos_array[i]
            atom['momentum'] = momentum_array[i]
            atom.append()

        snapshot.flush()

    # Close (and flush) the file
    h5_file.close()

    # Stop timing total computation time and output stats.
    end_time = time()
    tot_time = end_time-start_time
    print 80*'#'
    print "Run time (m): ", tot_time/60.
    print 80*'#'
