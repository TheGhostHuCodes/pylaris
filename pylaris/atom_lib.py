#!/usr/bin/env python

"""atom_lib.py: A collection of functions and classes that assist in reading
files, holding atom positions, and calculating atom distances."""

# import third-party modules.
from numpy import sqrt, arccos
from periodictable import elements

################################################################################
# Dictionaries
################################################################################
at_sym = [element.symbol for element in elements]
at_num = [element.number for element in elements]
at_wt = [element.mass for element in elements]

# Takes the atomic symbol and looks up the atomic number
atomic_sym2num = {key: value for (key, value) in zip(at_sym, at_num)}
# Takes the atomic symbol and looks up the atomic weight
atomic_sym2wt = {key: value for (key, value) in zip(at_sym, at_wt)}
# Takes the atomic number and looks up the atomic symbol
atomic_num2sym = {key: value for (key, value) in zip(at_num, at_sym)}
# Takes the atomic number and looks up the atomic weight
atomic_num2wt= {key: value for (key, value) in zip(at_num, at_wt)}

################################################################################
# Classes
################################################################################
class Atom():
    """Class that holds information about atoms."""
    def __init__(self, element, x_pos, y_pos, z_pos):
        self.element = element
        self.x_pos = x_pos
        self.y_pos = y_pos
        self.z_pos = z_pos

    def __sub__(self, other):
        """Overload the subtraction operator to return the Euclidean distance
        between them."""
        i = self.element
        j = other.element
        r_ij = sqrt((self.x_pos-other.x_pos)**2 +\
                (self.y_pos-other.y_pos)**2 + (self.z_pos-other.z_pos)**2)
        return r_ij

class Pair():
    """Class that holds information about atom pairs. given two "Atom"-type
    instances (atom_i, atom_j) this class calculates the pair distance between
    them and saves the atom name as well."""
    def __init__(self, atom_i, atom_j):
        # Point local variable to input "Atom"-type instances.
        self.atom_i = atom_i
        self.atom_j = atom_j
        # Calculate Euclidean distance between the two atoms using the operator
        # overload function of "Atom"-type instances.
        self.pair_dist = self.atom_i-self.atom_j

class Triples():
    """Class that holds information about atom triples. Given three "Atom"-type
    instances (atom_a, atom_b, atom_c) this class calculates pair distances
    between all three atoms (using the Pair class), and angles between all
    three atoms (using the law_of_cosines function, located inside this
    function)."""
    def __init__(self, atom_a, atom_b, atom_c):
        # Point local variable to input "Atom"-type instances.
        self.atom_a = atom_a
        self.atom_b = atom_b
        self.atom_c = atom_c
        # Enumerate all pairs from the three atoms as "Pair"-type instances.
        self.pair_ab = Pair(atom_a, atom_b)
        self.pair_ac = Pair(atom_a, atom_c)
        self.pair_bc = Pair(atom_b, atom_c)

        def law_of_cosines(adj_1, adj_2, opp):
            """Given the length of three sides of a triangle (adj_1, adj_2,
            opp), calculate the angle between sides "adj_1" and "adj_2". This
            is the same as the angle that is opposite side "opp"."""
            return arccos((adj_1**2 + adj_2**2 - opp**2)/(2*adj_1*adj_2))

        # Calculate angles between atoms. For atom_a the angle with atom_a as
        # the central atom is ang_bc, etc. The angle should <= pi radians by
        # nature of the way it is calculated.
        self.ang_ab = law_of_cosines(self.pair_ac.pair_dist,
                self.pair_bc.pair_dist, self.pair_ab.pair_dist)
        self.ang_ac = law_of_cosines(self.pair_ab.pair_dist,
                self.pair_bc.pair_dist, self.pair_ac.pair_dist)
        self.ang_bc = law_of_cosines(self.pair_ab.pair_dist,
                self.pair_ac.pair_dist, self.pair_bc.pair_dist)

################################################################################
# Functions
################################################################################
def read_xyz(filename):
    """Reads a *.xyz file and returns a list of Atom class instances from
    the file."""
    input_file = open(filename, 'r').readlines()
    atom_list = []
    # Read in the number of atoms from the *xyz file for error checking.
    num_of_atoms = int(input_file[0])
    # Read in the header if there is any and just print it.
    header_line = input_file[1].strip()
    if header_line != "":
        print "HEADER: ", header_line
    for line in input_file[2:]:
        # For an *.xyz file column[0] = atom symbol, column[1] = x_pos,
        # column[2] = y_pos, column[3] = z_pos
        columns = line.split()
        atom_instance = Atom(columns[0], float(columns[1]),
                float(columns[2]), float(columns[3]))
        atom_list.append(atom_instance)
    if len(atom_list) == num_of_atoms:
        print "*.xyz story checks out."
    return atom_list
