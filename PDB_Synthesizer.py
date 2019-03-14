#!/usr/bin/env python
import numpy as np
import MDAnalysis
import warnings
import subprocess, os, StringIO, re, datetime, time, sys
import optparse

def concat(u1, u2):
    """Concatenate two MDAnalysis universes

    Parameters
    ----------
    u1 : MDAnalysis universe
        First universe to be put together


    u2 : MDAnalysis universe
        Second universe to be put together

    Returns
    -------
    u_comb : MDAnalysis universe
        Combined universes

    """
    u_comb = MDAnalysis.Merge(u1, u2)
    p1 = u_comb.select_atoms("bynum 1:%s" %str(len(u1)))
    p2 = u_comb.select_atoms("bynum %s:%s" %(str(len(u1) + 1), str(len(u1) + len(u2))))
    p1.segments.segids = "A"
    p2.segments.segids = "B"
    #p2.residues.set_resid(p2.residues.resids + p1.residues.resids[-1])
    p2.residues.resids = (p2.residues.resids + p1.residues.resids[-1])
    return u_comb

#def construct(entries, write_to, offset_vectors, rotation_axes, angles):
def construct(entries, write_to, link_settings):
    """Construct a linear polymer of a given composition based on the monomer sequence provided

    Parameters
    ----------
    entries : list of strings
        The ordered sequence of monomers to be joined together

    write_to: string
        The path to which output pdb will be written to

    link_settings: list of LinkSetting objects
        Defines how to place the next residue with two translations
        Translation 1:
            vector from link_setting.connect_from to link_setting.connect_to
        Translation 2:
            vector from link_setting.bond_ref_from to link_setting.bond_ref_to


    Returns
    -------
        None

    """
    assert len(link_settings) == len(entries) - 1
    #assert len(rotation_axes) == len(entries) -1
    #assert len(angles) == len(entries) -1
    #assert the length of link_settings is one unit less than number of monomers entried
    for i in range(len(entries)):
        monomer = MDAnalysis.Universe(os.path.abspath('{}.pdb'.format(entries[i])))

        #First entry -> create polymer
        if i == 0:
            polymer = monomer
            continue

        position_N = polymer.atoms.positions[np.where(polymer.atoms.names == "N" ), ][0][-1]
        # compute offset vector (offset) np.array([x,y,z])
        #
        # print N
        # mer.atoms.positions = np.sum([mer.atoms.positions, N], axis = 0) + 1 <- BAD
        monomer.atoms.positions += link_settings[i - 1].compute_offset(monomer, polymer)
        polymer = concat(polymer.atoms, monomer.atoms)
    polymer.segments.segids = "A"
    polymer.atoms.write(write_to)
    return

class LinkSetting():
    """
    This class defines how two residues (A and B) are connected
    We need two offsets to decide where to place residue B

    connect_offset:
        This is the vector between the terminal atoms

    bond_offset:
        This is the vector that defines the new bond between the residues

    connect_from:
        This is the name of the atom that terminates A

    connect_to:
        This is the name of the atom that begins B

    bond_ref_from:
        This is the name of the atom that begins the reference bond to compute bond_offset

    bond_ref_to:
        This is the name of the atom that ends the reference bond to compute bond_offset


    This class has 4 ways to be constructed:
        1)
            Define:
                connect_from
                connect_to
                bond_ref_from
                bond_ref_to
        2)
            Define:
                connect_from
                connect_to
                bond_offset
        3)
            Define:
                bond_ref_from
                bond_ref_to
                connect_offset
        4)
            Define:
                connect_offset
                bond_offset

    """

    def __init__(
        self,
        connect_from=None,
        connect_to=None,
        bond_ref_from=None,
        bond_ref_to=None,
        connect_offset=None,
        bond_offset=None
    ):
        self.connect_from = connect_from
        self.connect_to = connect_to
        self.bond_ref_from = bond_ref_from
        self.bond_ref_to = bond_ref_to
        self.connect_offset = connect_offset
        self.bond_offset = bond_offset

    def compute_offset(self, residueA, residueB):
        if (self.connect_offset is None):
            position_0 = self.get_position(residueA, self.connect_from)
            position_1 = self.get_position(residueB, self.connect_to)
            connect_offset = position_1 - position_0
        else:
            connect_offset = self.connect_offset

        if (self.bond_offset is None):
            position_2 = self.get_position(residueA, self.bond_ref_from)
            position_3 = self.get_position(residueA, self.bond_ref_to)
            bond_offset = position_2 - position_3
        else:
            bond_offset = self.bond_offset

        return connect_offset + bond_offset

    def get_position(self, residue, name):
        return residue.atoms.positions[np.where(residue.atoms.names == name ), ][0][-1]



"""
Example usage:

construct(["MB1"]*5  , "./MB1x5.pdb")
"""
