"""AutoClickChem is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

AutoClickChem is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Copyright 2011 Jacob D. Durrant. If you have any questions, comments, or
suggestions, please don't hesitate to contact me at jdurrant [at] ucsd [dot] edu.

The latest version of AutoClickChem can be downloaded from 
http://sourceforge.net/projects/autoclickchem/

If you use AutoClickChem in your work, please cite [REFERENCE HERE]"""

import os
import sys
import glob
import pymolecule
import time
import random
import math

class StructureLocateGroups:
    """This class contains functions that aid the identification of reactive chemical groups."""
    
    def index_of_azide(self, pdb):
        """Identifies all the azide groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified azide group.
        
        """
        
        # This will return the indices of the atoms, proximal to distal. R-X-N=N=N
        # So we're looking for three nitrogens without any hydrogens, arranged in a straight line.
            
        AzideRoots = []
        
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
            if atom.element=="N" and len(atom.indecies_of_atoms_connecting)==1: # This is the terminal N
                    neighbor1_index = atom.indecies_of_atoms_connecting[0]
                    neighbor1_atom = pdb.all_atoms[neighbor1_index]
                    if neighbor1_atom.element=="N" and len(neighbor1_atom.indecies_of_atoms_connecting)==2: # middle N
                        indicies = []
                        for index in neighbor1_atom.indecies_of_atoms_connecting:
                            indicies.append(index)
                        indicies.remove(atom_index)
                        neighbor2_index = indicies[0]
                        neighbor2_atom = pdb.all_atoms[neighbor2_index]
                        if neighbor2_atom.element == "N" and len(neighbor1_atom.indecies_of_atoms_connecting)==2: # proximal N
                            # they have to be in a straight line
                            angle = atom.coordinates.angle_between_three_points(neighbor1_atom.coordinates, neighbor2_atom.coordinates) * 180 / math.pi
                            if abs(angle-180) < 10: # so it is in a straight line
                                indicies = []
                                for index in neighbor2_atom.indecies_of_atoms_connecting:
                                    indicies.append(index)
                                indicies.remove(neighbor1_index)
                                neighbor3_index = indicies[0] # X described above
                                neighbor3_atom = pdb.all_atoms[neighbor3_index]
                                
                                if neighbor3_atom.element == "C": # This is so sulfonyl azides are not detected here.
                                    
                                    #Now, to be able to position sp2 hybridized amines correctly (azide => sp2 amine), we need one more atom index.
                                    indicies = neighbor3_atom.indecies_of_atoms_connecting[:]
                                    for index in indicies: # get rid of hydrogens, so that only heteroatoms remain in the list
                                        if pdb.all_atoms[index].element == 'H': indicies.remove(index)
                                        if index == neighbor2_index: indicies.remove(index)
                                    neighbor4_index = indicies[0]
    
                                    list = [neighbor3_index, neighbor2_index, neighbor1_index, atom_index]
                                    if len(indicies)>0: list.append(neighbor4_index) # neighbor4 index put at end for backwards compatibility (logically, should be at beginning)
                                    
                                    AzideRoots.append(list)
                                    
        return AzideRoots
    
    def index_of_sulfonyl_azide(self, pdb):
        """Identifies all the sulfonyl azide groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified sulfonyl azide group.
        
        """
        
        # This will return the indices of the atoms, proximal to distal. R-S-N=N=N
        # So we're looking for three nitrogens without any hydrogens, arranged in a straight line.
            
        SulfonylAzideRoots = []
        
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
            if atom.element=="N" and len(atom.indecies_of_atoms_connecting)==1: # this is the terminal nitrogen of the azide group
                neighbor1_index = atom.indecies_of_atoms_connecting[0]
                neighbor1_atom = pdb.all_atoms[neighbor1_index] # this is the middle nitrogen of the azide group
                if neighbor1_atom.element=="N" and len(neighbor1_atom.indecies_of_atoms_connecting)==2:
                    indicies = neighbor1_atom.indecies_of_atoms_connecting[:]
                    indicies.remove(atom_index)
                    neighbor2_index = indicies[0]
                    neighbor2_atom = pdb.all_atoms[neighbor2_index] # this is the proximal nitrogen of the azide group
                    if neighbor2_atom.element == "N" and len(neighbor1_atom.indecies_of_atoms_connecting)==2:
                        # they have to be in a straight line
                        angle = atom.coordinates.angle_between_three_points(neighbor1_atom.coordinates, neighbor2_atom.coordinates) * 180 / math.pi
                        if abs(angle-180) < 10: # so it is in a straight line
                            indicies = neighbor2_atom.indecies_of_atoms_connecting[:]
                            indicies.remove(neighbor1_index)
                            neighbor3_index = indicies[0]
                            neighbor3_atom = pdb.all_atoms[neighbor3_index]
                            
                            if neighbor3_atom.element == "S" and neighbor3_atom.number_of_neighbors() == 4: # this is the sulfonyl part of the sulfonyl azide
                                # now make sure at least two of the atoms attached to this S are O's
                                oxygens = pdb.number_of_neighors_of_element(neighbor3_index, "O")
                                if oxygens >= 2:
                                        SulfonylAzideRoots.append([neighbor3_index, neighbor2_index, neighbor1_index, atom_index])
                            
        return SulfonylAzideRoots
    
    def index_of_thio_acid(self, pdb):
        """Identifies all the thio acid groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified thio acid group.
        
        """
    
        # This will return the indices of the atoms, proximal to distal. R-C(=O)SH
            
        ThioAcidRoots = []
        
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
            if atom.element=="C" and len(atom.indecies_of_atoms_connecting)==3 and pdb.number_of_neighors_of_element(atom_index, "O") >= 1 and pdb.number_of_neighors_of_element(atom_index, "S") >= 1:
                listconnected = atom.indecies_of_atoms_connecting[:] # make a copy
                for index in range(len(listconnected)):
                    neighbor_index = listconnected[index]
                    if pdb.all_atoms[neighbor_index].element == "O" and pdb.all_atoms[neighbor_index].number_of_neighbors() == 1: # so this is the carbonyl oxygen. delete it.
                        carbonyl_oxygen_index = neighbor_index
                        listconnected[index] = "DELETE"
                    if pdb.all_atoms[neighbor_index].element == "S" and pdb.all_atoms[neighbor_index].number_of_neighbors() == 2 and pdb.number_of_neighors_of_element(neighbor_index, "H") == 1:
                        sulfur_index = neighbor_index
                        listconnected[index] = "DELETE"
                while "DELETE" in listconnected: listconnected.remove("DELETE")
                first_atom_of_R_group_index = listconnected[0]

                if len(listconnected) == 1: # sanity check. Makes sure the connected oxygen is a carbonyl oxygen
                    sulfur_atom = pdb.all_atoms[sulfur_index]
                    listconnected = sulfur_atom.indecies_of_atoms_connecting[:] # make a copy
                    listconnected.remove(atom_index)
                    terminal_hydrogen_index = listconnected[0]
            
                    if len(listconnected) == 1:
                        ThioAcidRoots.append([carbonyl_oxygen_index, atom_index, sulfur_index, terminal_hydrogen_index])
                            
        return ThioAcidRoots
    
    def index_of_alkene(self, pdb): # identify alkenes (olefins)
        """Identifies all the alkene groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified alkene group.
        
        """
        
        #  2      5
        #   \    /
        #   1 = 4
        #  /     \
        # 3       6
        
        # unfortunately, this is going to mess up on things like -C=C-C=C-
        # the middle two will be epoxidated incorrectly.
    
        AlkeneRoots=[]
        
        for first_C_index in pdb.all_atoms:
            first_C = pdb.all_atoms[first_C_index]
            if first_C.element=="C" and first_C.number_of_neighbors() == 3: # you've found some sp2 hybridized carbon
                for second_C_index in first_C.indecies_of_atoms_connecting:
                    second_C = pdb.all_atoms[second_C_index]
                    if second_C.element=="C" and second_C.number_of_neighbors() == 3: # you've found a neighboring sp2 hybridized carbon
                        # now make sure they are not in the same ring. It would be nice if we could include olefins in rings, ut it would be too hard to code. :(
                        if pdb.in_same_ring(first_C_index, second_C_index) == "FALSE":
                            first_connectors = first_C.indecies_of_atoms_connecting[:]
                            second_connectors = second_C.indecies_of_atoms_connecting[:]
    
                            first_connectors.remove(second_C_index)
                            second_connectors.remove(first_C_index)
                            
                            items=[]
                            
                            items.append(first_C_index)
                            for index in first_connectors:
                                items.append(index)
                            items.append(second_C_index)
                            for index in second_connectors:
                                items.append(index)
                                
                            # now we need to make sure it's not in any kind of ring
                            if pdb.in_same_ring(items[0], items[3]) == "FALSE":
                                if pdb.in_same_ring(items[0], items[1]) == "FALSE":
                                    if pdb.in_same_ring(items[3], items[4]) == "FALSE":
    
                                        # now let's make sure not to include vinyl carboxylates, ketones, etc. These show up suprisingly frequently
                                        neighbors_all = first_C.indecies_of_atoms_connecting[:]
                                        neighbors_all.extend(second_C.indecies_of_atoms_connecting)
                                        neighbors_all.remove(first_C_index)
                                        neighbors_all.remove(second_C_index)
                                        
                                        continue_ok = True
                                        for index in neighbors_all:
                                            neighbor_all_atom = pdb.all_atoms[index]
                                            if neighbor_all_atom.element == "O" and neighbor_all_atom.number_of_neighbors() == 1: # so it's a carbonyl oxygen atom
                                                continue_ok = False
                                                break
                                        
                                        if continue_ok == True: AlkeneRoots.append(items)
    
        # every alkene has been identified in duplicate. remove the duplicates
        toremove = []
        for item1_index in range(0,len(AlkeneRoots)-1):
            item1 = AlkeneRoots[item1_index]
            item1_sorted=item1[:]
            item1_sorted.sort()
            for item2_index in range(item1_index+1,len(AlkeneRoots)):
                item2 = AlkeneRoots[item2_index]
                item2_sorted=item2[:]
                item2_sorted.sort()
                if item1_sorted == item2_sorted:
                    if item2 not in toremove: toremove.append(item2)
                    
        for remove_item in toremove:
            AlkeneRoots.remove(remove_item)
        
        return AlkeneRoots
                
    
    def index_of_alkyne(self, pdb): # Identify an alkyne, terminal or internal
        """Identifies all the alkyne groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified alkyne group.
        
        """
        
        # Think first_neighbor - first_C =- second_C - second_neighbor
        
        AlkyneRoots=[]
    
        for first_C_index in pdb.all_atoms:
            first_C = pdb.all_atoms[first_C_index]
    
            #print str(first_C_index) + "\t" + first_C.element + "\t" + str(len(first_C.indecies_of_atoms_connecting))
            #print first_C.indecies_of_atoms_connecting
            
            if first_C.element=="C" and len(first_C.indecies_of_atoms_connecting)==2: # so it finds the one of the sp1 hybridized carbons
    
                neighbor1_index = first_C.indecies_of_atoms_connecting[0]
                neighbor2_index = first_C.indecies_of_atoms_connecting[1]
                neighbor1 = pdb.all_atoms[neighbor1_index]
                neighbor2 = pdb.all_atoms[neighbor2_index]
    
                # one of these neighbors is the other carbon of the alkyne bond. Identify it.
                if neighbor1.element=="C" and len(neighbor1.indecies_of_atoms_connecting)==2: # so neighbor1 is the other carbon of the alkyne bond
                    second_C = neighbor1
                    second_C_index = neighbor1_index
                    do_continue="TRUE"
                elif neighbor2.element=="C" and len(neighbor2.indecies_of_atoms_connecting)==2:
                    second_C = neighbor2 # so neighbor2 is the other carbon of the alkyne bond
                    second_C_index = neighbor2_index
                    do_continue="TRUE"
                else: do_continue="FALSE" # the identified carbon is not connected to another sp1-hybridized carbon. Bond not identified here.
                
                if do_continue == "TRUE": # so we do have an alkyne bond
                    # make a copy of the list of the two atoms bound to first_C
                    first_C_bonded = first_C.indecies_of_atoms_connecting[:]
                    
                    # remove the index of second_C
                    first_C_bonded.remove(second_C_index)
                    
                    # make a copy of the list of the two atoms bound to second_C
                    second_C_bonded = second_C.indecies_of_atoms_connecting[:]
                    
                    # remove the index of first_C
                    second_C_bonded.remove(first_C_index)
                    
                    first_neighbor_index = first_C_bonded[0]
                    first_neighbor = pdb.all_atoms[first_neighbor_index]
                    second_neighbor_index = second_C_bonded[0]
                    second_neighbor = pdb.all_atoms[second_neighbor_index]
                    
                    angle1 = first_neighbor.coordinates.angle_between_three_points(first_C.coordinates, second_C.coordinates) * 180 / math.pi
                    angle2 = first_C.coordinates.angle_between_three_points(second_C.coordinates, second_neighbor.coordinates) * 180 / math.pi
                    
                    if abs(angle1 - 180) < 10 and abs(angle2 - 180) < 10:
                        AlkyneRoots.append([first_neighbor_index, first_C_index, second_C_index, second_neighbor_index])
        
        # now, because there is no directionality to this bond, everything is in duplicate. Remove the duplicates.
        ToRemove = []
        for index1 in range(0,len(AlkyneRoots)-1):
            for index2 in range(index1 + 1, len(AlkyneRoots)):
                list1 = AlkyneRoots[index1]
                list2 = AlkyneRoots[index2]
                if list1[0] == list2[3] and list1[1] == list2[2] and list1[2] == list2[1] and list1[3] == list2[0]:
                    ToRemove.append(list2)
        
        for list in ToRemove:
            AlkyneRoots.remove(list)
        
        return AlkyneRoots
    
    def index_of_alcohol(self, pdb): # Identify an alcohol
        """Identifies all the alcohol groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified alcohol group.
        
        """
    
        # Think R - O - H (1 - 2 - 3)
        
        AlcoholRoots=[]
        
        for oxygen_index in pdb.all_atoms:
            oxygen = pdb.all_atoms[oxygen_index]
            if oxygen.element == "O" and oxygen.number_of_neighbors() == 2:
                #Now, one of those neighbors needs to be a hydrogen
                for hydrogen_index in oxygen.indecies_of_atoms_connecting:
                    hydrogen = pdb.all_atoms[hydrogen_index]
                    if hydrogen.element == "H":
                        neighbors = oxygen.indecies_of_atoms_connecting[:]
                        neighbors.remove(hydrogen_index)
                        other_index = neighbors[0]
                        other_atom = pdb.all_atoms[other_index]
                        
                        if other_atom.element == "C": # We're only interest in alcohols attached to carbons
                            # We need to make sure that none of the atoms bound to the C are double-bonded oxygens. Carboxylic acids are not acceptable.
                            # get a list of all the neighbors
                            neighbors = other_atom.indecies_of_atoms_connecting[:]
                            IsAlcohol = "TRUE"
                            for index in neighbors:
                                if pdb.all_atoms[index].element == "O" and pdb.all_atoms[index].number_of_neighbors() == 1:
                                    IsAlcohol = "FALSE"
                                    
                            if IsAlcohol== "TRUE" : AlcoholRoots.append([other_index, oxygen_index, hydrogen_index])
    
        return AlcoholRoots
    
    def index_of_thiol(self, pdb): # Identify an alcohol
        """Identifies all the thiol groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified thiol group.
        
        """
        
        # Think R - S - H (1 - 2 - 3)
        
        ThiolRoots=[]
        
        for sulfur_index in pdb.all_atoms:
            sulfur = pdb.all_atoms[sulfur_index]
            if sulfur.element == "S" and sulfur.number_of_neighbors() == 2:
                #Now, one of those neighbors needs to be a hydrogen
                for hydrogen_index in sulfur.indecies_of_atoms_connecting:
                    hydrogen = pdb.all_atoms[hydrogen_index]
                    if hydrogen.element == "H":
                        neighbors = sulfur.indecies_of_atoms_connecting[:]
                        neighbors.remove(hydrogen_index)
                        
                        other_index = neighbors[0] # This is the carbon the thiol is attached to (R above)
                        other_atom = pdb.all_atoms[other_index]
                        
                        if other_atom.element == "C": # We're only interest in thiols attached to carbons
                            # We need to make sure that none of the atoms bound to the C are double-bonded oxygens. Thiol acids are not acceptable.
                            # get a list of all the neighbors
                            neighbors = other_atom.indecies_of_atoms_connecting[:]
                            IsThiol = "TRUE"
                            for index in neighbors:
                                if pdb.all_atoms[index].element == "O" and pdb.all_atoms[index].number_of_neighbors() == 1:
                                    IsThiol = "FALSE"
                                    
                            if IsThiol == "TRUE" : ThiolRoots.append([other_index, sulfur_index, hydrogen_index])
    
        return ThiolRoots
    
    def index_of_epoxide(self, pdb): # Identify an epoxide # still needs to be tested # the epoxide cannot be part of a larger ring system
        """Identifies all the epoxide groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified epoxide group.
        
        """
        
        #           O(0)
        #          /   \
        # X(2) - C(1) - C(4) - X(5)
        #        /       \
        #     X(3)        X(6)
    
        EpoxideRoots=[]
        
        for oxygen_index in pdb.all_atoms:
            oxygen = pdb.all_atoms[oxygen_index]
            if oxygen.element == "O" and pdb.number_of_neighors_of_element(oxygen_index, "C") == 2 and oxygen.number_of_neighbors() == 2:
                # now we need to see if the two carbons connected to that O are connected to each other
                carbon1_index = oxygen.indecies_of_atoms_connecting[0]
                carbon2_index = oxygen.indecies_of_atoms_connecting[1]
                carbon1 = pdb.all_atoms[carbon1_index]
                carbon2 = pdb.all_atoms[carbon2_index]
                if carbon2_index in carbon1.indecies_of_atoms_connecting and carbon1_index in carbon2.indecies_of_atoms_connecting: # so they are connected
                    # now need to make sure that carbon1_index and carbon2_index are not in a ring together (except for the epoxide ring itself)
                    temp_pdb = pdb.copy_of()
                    temp_pdb.delete_atom(oxygen_index)
                    if temp_pdb.in_same_ring(carbon1_index, carbon2_index) == "FALSE":
        
                        carbon1_neighbors = carbon1.indecies_of_atoms_connecting[:]
                        carbon2_neighbors = carbon2.indecies_of_atoms_connecting[:]
                        
                        # remove oxygen from list
                        carbon1_neighbors.remove(oxygen_index)
                        carbon2_neighbors.remove(oxygen_index)
                        
                        # remove other carbon you're bound too
                        carbon1_neighbors.remove(carbon2_index)
                        carbon2_neighbors.remove(carbon1_index)
                        
                        carbon1_neighbor1_index = carbon1_neighbors[0]
                        carbon1_neighbor2_index = carbon1_neighbors[1]
                        carbon2_neighbor1_index = carbon2_neighbors[0]
                        carbon2_neighbor2_index = carbon2_neighbors[1]
                                
                        if pdb.in_same_ring(carbon1_index, carbon1_neighbor1_index) == "FALSE":
                            if pdb.in_same_ring(carbon2_index, carbon2_neighbor1_index) == "FALSE":
                                EpoxideRoots.append([oxygen_index, carbon1_index, carbon1_neighbor1_index, carbon1_neighbor2_index, carbon2_index, carbon2_neighbor1_index, carbon2_neighbor2_index])
        
        return EpoxideRoots
    
    def index_of_primary_amine(self, pdb): # Here's the name of the function we'll call to identify primary amines. The function accepts as input a pdb object.
        """Identifies all the primary amine groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified primary amine group.
        
        """
    
        Indices_of_primary_amines = []  # Here's the list that will eventually contain the locations of all the primary amines in the pdb object.
                                        # This is actually a list of lists. Each primary amine looks like this: R - X - N - H(2)
                                        # To identify the amine, it's not enough just to identify the index of the N atom. We also need to identify
                                        # the index of it's neighbor, X, and it's hydrogens as well. So, we'll be adding quadrulets to the Indices_of_primary_amines list, 
                                        # (index_N, index_X, index_H1, index_H2, ...).
    
                                        # For example, say the pdb is H(2) - N - CH(2) - CH(2) - N - H(2), and the atoms have the following indicies:
                                        #                              1,2   3   4,5,6   7,8,9   10  11,12, or, just focusing on the heavy atoms...
                                        #                                    3   4       7       10
                                        # Then we want Indices_of_primary_amines to look like this in the end: [[4,3,1,2], [7,10,11,12]].
    
        for atom_index in pdb.all_atoms: # Now, we're going to look at each of the atoms in the pdb object. For each iteration of the loop, the variable "atom_index" is the 
                                        # index of the atom we're examining.
            atom = pdb.all_atoms[atom_index]         # "atom_index" was only the index of the atom. Now, the "atom" variable points to the atom object itself, rather than just the index.
            if atom.element=="N":                   # Look and see if the atom is a nitrogen.
                nitrogen_index = atom_index     # If it is a nitrogen, then store the index of that nitrogen in the variable "nitrogen_index"
                if pdb.number_of_neighors_of_element(nitrogen_index,"H") >= atom.number_of_neighbors()-1:      # If it's a primary amine, this should be true
                                                                                                         # So ammonium is not a primary amine.
                    # So we are dealing with a primary amine!
    
                    if atom.number_of_neighbors() == 3 or atom.number_of_neighbors() == 4:       # A primary amine should have only three atoms connected to it (assuming it's not charged)
                                                            # This conditional is actually unnecessary, but I wanted to show you how to detect the
                                                            # number of atoms connected to a given atom of interest
    
                        nitrogen_neighbors = atom.indecies_of_atoms_connecting[:]  # Here's a list the indices of all the atoms connected to the nitrogen of the primary amine.
    
                        for neighbor_index in nitrogen_neighbors:               # We're going to look at each neighbor. The "neighbor_index" variable contains the index of the 
                                                                        # neighbor we're looking at.
                            neighbor_atom = pdb.all_atoms[neighbor_index]    # "neighbor_index" contained the index of the neighbor, but the variable "neighbor_atom" points
                                                                    # to the neighbor atom object itself.
                            if neighbor_atom.element != "H":                # So we're looking at the neighbor that is not a hydrogen, atom "X" above.
                                index_x = neighbor_index                # We're going to save the index of this non-H atom in the variable "index_x"
                                nitrogen_neighbors.remove(index_x)      # We're now removing the index of this non-H atom from the list, so that the only indices remaining
                                                                        # in the list are the indices of the hydrogen neighbors of the amine N.
                        if pdb.all_atoms[index_x].element == "C": # Require that X be a carbon.
                            # Recall that we need to create a list identifying the atoms of this amine in the order (index_N, index_X, index_H1, index_H2)
                            amine_id = [nitrogen_index, index_x]    # Because we removed the index of the non-H nitrogen neighbor,
                            for item in nitrogen_neighbors:         # we know that nitrogen_neighbors[0] and nitrogen_neighbors[1]
                                amine_id.append(item)               # correspond to the indices of the two neighboring hydrogens.
                                                                    # Let's call this list a "primary-amine idenfitier."
                            # Having formed this identifying list, let's add this list to our growing list of "primary-amine indentifiers," Indices_of_primary_amines.
                            if len(nitrogen_neighbors) >=2 :
                                Indices_of_primary_amines.append(amine_id)
    
        return Indices_of_primary_amines # Having identified all the primary amines in the pdb, let's return our list of "primary-amine identifiers."
    
    
    def index_of_primary_amine_attached_to_S(self, pdb): # Because sulfonyl azides can be reduced to sulfonyl amines.
        """Identifies all the primary amine groups attached to sulfur atoms in a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified group.
        
        """
        
        Indices_of_primary_amines = []  # Here's the list that will eventually contain the locations of all the primary amines in the pdb object.
        for atom_index in pdb.all_atoms: # Now, we're going to look at each of the atoms in the pdb object. For each iteration of the loop, the variable "atom_index" is the 
            atom = pdb.all_atoms[atom_index]         # "atom_index" was only the index of the atom. Now, the "atom" variable points to the atom object itself, rather than just the index.
            if atom.element=="N":                   # Look and see if the atom is a nitrogen.
                nitrogen_index = atom_index     # If it is a nitrogen, then store the index of that nitrogen in the variable "nitrogen_index"
                if pdb.number_of_neighors_of_element(nitrogen_index,"H") >= atom.number_of_neighbors()-1:      # If it's a primary amine, this should be true
                    if atom.number_of_neighbors() == 3 or atom.number_of_neighbors() == 4:       # A primary amine should have only three atoms connected to it (assuming it's not charged)
                        nitrogen_neighbors = atom.indecies_of_atoms_connecting[:]  # Here's a list the indices of all the atoms connected to the nitrogen of the primary amine.
                        for neighbor_index in nitrogen_neighbors:               # We're going to look at each neighbor. The "neighbor_index" variable contains the index of the 
                            neighbor_atom = pdb.all_atoms[neighbor_index]    # "neighbor_index" contained the index of the neighbor, but the variable "neighbor_atom" points
                            if neighbor_atom.element != "H":                # So we're looking at the neighbor that is not a hydrogen, atom "X" above.
                                index_x = neighbor_index                # We're going to save the index of this non-H atom in the variable "index_x"
                                nitrogen_neighbors.remove(index_x)      # We're now removing the index of this non-H atom from the list, so that the only indices remaining
                        if pdb.all_atoms[index_x].element == "S": # Require that X be a sulfur.
                            # Recall that we need to create a list identifying the atoms of this amine in the order (index_N, index_X, index_H1, index_H2)
                            amine_id = [nitrogen_index, index_x]    # Because we removed the index of the non-H nitrogen neighbor,
                            for item in nitrogen_neighbors:         # we know that nitrogen_neighbors[0] and nitrogen_neighbors[1]
                                amine_id.append(item)               # correspond to the indices of the two neighboring hydrogens.
                            if len(nitrogen_neighbors) >=2 :
                                Indices_of_primary_amines.append(amine_id)
        return Indices_of_primary_amines # Having identified all the primary amines in the pdb, let's return our list of "primary-amine identifiers."
    
    
    def index_of_secondary_amine(self, pdb): 
        """Identifies all the secondary amine groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified secondary amine group.
        
        """
        
        Indices_of_secondary_amines = []  
        for atom_index in pdb.all_atoms: 
            atom = pdb.all_atoms[atom_index]         
            if atom.element=="N": # So looking for the nitrogen       
                nitrogen_index = atom_index     
                if pdb.number_of_neighors_of_element(nitrogen_index,"C") == 2 and pdb.number_of_neighors_of_element(nitrogen_index,"H") + pdb.number_of_neighors_of_element(nitrogen_index,"C") == atom.number_of_neighbors() and pdb.number_of_neighors_of_element(nitrogen_index,"H") >0: #so there are two C atoms attached to the nitrogen and at least 1 H
                    amine_id = [nitrogen_index]    
    
                    nitrogen_neighbors = atom.indecies_of_atoms_connecting[:]  # here's alist of all the N neighbors
    
                    # Now add in all the non-H neighbors
                    for neighbor_index in nitrogen_neighbors:               
                        neighbor_atom = pdb.all_atoms[neighbor_index]    
                        if neighbor_atom.element != "H": 
                            amine_id.append(neighbor_index) # You only need the index of one of the non-H neighbors
    
                    # Now add the hydrogen atoms
                    for neighbor_index in nitrogen_neighbors:               
                        neighbor_atom = pdb.all_atoms[neighbor_index]    
                        if neighbor_atom.element == "H": 
                            amine_id.append(neighbor_index) # You only need the index of one of the non-H neighbors
                            
                    # Add the list
                    Indices_of_secondary_amines.append(amine_id)
    
        return Indices_of_secondary_amines 
    
    def index_of_carboxylate(self, pdb): #The function accepts a pdb object as input
        """Identifies all the carboxylate groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified carboxylate group.
        
        """
    
        Indices_of_carboxylates = []
    
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
    
            if atom.element=="C":		   # Look and see if the atom is a carbon
                carbon_index = atom_index
    
                if pdb.number_of_neighors_of_element(carbon_index,"O") == 2:	# If it's a carboxylate, it will have two oxygens attached to it
    
                    if atom.number_of_neighbors() == 3:			   # A carboxylate should have only three atoms connected to it
                        carbon_neighbors = atom.indecies_of_atoms_connecting[:]
    
                        # remove anything that isn't an oxygen
                        for neighbor_index in carbon_neighbors:	 #Look at each neighbor and address right object to each index
                            neighbor_atom = pdb.all_atoms[neighbor_index]
                            if neighbor_atom.element != "O": carbon_neighbors.remove(neighbor_index)
                            
                        # So now you have a list of the indicies of the two oxygens
                        if pdb.number_of_neighors_of_element(carbon_neighbors[0],"H") + pdb.number_of_neighors_of_element(carbon_neighbors[1],"H") <=1: # there should only be one hydrogen between the two oxygens
                            oxy1 = pdb.all_atoms[carbon_neighbors[0]]
                            oxy2 = pdb.all_atoms[carbon_neighbors[1]]
                            
                            if oxy1.number_of_neighbors() + oxy2.number_of_neighbors() <=3: # There should be no more than 3 bonds between the two oxygens
                                hydro_index_1 = pdb.index_of_neighbor_of_element(carbon_neighbors[0],"H")
                                hydro_index_2 = pdb.index_of_neighbor_of_element(carbon_neighbors[1],"H")
                                if hydro_index_1 != -1 and hydro_index_2 == -1: # so carbon_neighbors[0] is the carboxylate alcohol
                                    Indices_of_carboxylates.append([carbon_index,carbon_neighbors[1], carbon_neighbors[0], hydro_index_1])
                                elif hydro_index_1 == -1 and hydro_index_2 != -1: # so carbon_neighbors[1] is the carboxylate alcohol
                                    Indices_of_carboxylates.append([carbon_index,carbon_neighbors[0], carbon_neighbors[1], hydro_index_2])
                                elif oxy1.number_of_neighbors() + oxy2.number_of_neighbors() == 2: # so if it's a ester, no go, but otherwise, order doesn't matter (carboxylate vs. carboxylic acid)
                                    Indices_of_carboxylates.append([carbon_index,carbon_neighbors[0], carbon_neighbors[1], -1])
    
        return Indices_of_carboxylates
    
    def index_of_acylhalide(self, pdb): #The function accepts a pdb object as input
        """Identifies all the acyl halide groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified acyl halide group.
        
        """
    
        Indices_of_acylhalides = []
    
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
    
            if atom.element=="C":
                carbon_index = atom_index
    
                if pdb.number_of_neighors_of_element(carbon_index,"O") == 1:        # If it's an acylhalide, it will have one oxygen attached to it...
    
                    if pdb.number_of_neighors_of_element(carbon_index,"CL") == 1 or pdb.number_of_neighors_of_element(carbon_index,"BR") == 1 or pdb.number_of_neighors_of_element(carbon_index,"I") == 1:       #...and one halide 
    
                        if atom.number_of_neighbors() == 3:
                            carbon_neighbors = atom.indecies_of_atoms_connecting[:]
    
                            for neighbor_index in carbon_neighbors:         #Look at each neighbor and address right object to each index
                                neighbor_atom = pdb.all_atoms[neighbor_index]
    
                                if neighbor_atom.element == "O":
                                    oxygen_index = neighbor_index
                                elif neighbor_atom.element == "CL" or neighbor_atom.element == "BR" or neighbor_atom.element == "I":
                                    halide_index = neighbor_index
    
    
                            acylhalide_id = [carbon_index, oxygen_index, halide_index]
                            Indices_of_acylhalides.append(acylhalide_id)
    
        return Indices_of_acylhalides
    
    def index_of_ester(self, pdb): #The function accepts a pdb object as input
        """Identifies all the ester groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified ester group.
        
        """
    
        Indices_of_esters = []
    
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
    
            if atom.element=="C":                   # Look and see if the atom is a carbon
                carbon_index = atom_index
    
                if pdb.number_of_neighors_of_element(carbon_index,"O") == 2 and pdb.number_of_neighors_of_element(carbon_index,"C") == 1:        # If it's a ester, it will have two oxygens attached to it and one carbon
    
                    if atom.number_of_neighbors() == 3:                       # A ester carbon should have only three atoms connected to it
                        carbon_neighbors = atom.indecies_of_atoms_connecting[:]
                        for neighbor_index in carbon_neighbors:         #Remove non-O neighbors
                            neighbor_atom = pdb.all_atoms[neighbor_index]
                            if neighbor_atom.element != "O":        #non-O neighbors are not included in the list   
                                carbon_neighbors.remove(neighbor_index)
    
                        # So now you have a list of the indicies of the two oxygens
                        if pdb.number_of_neighors_of_element(carbon_neighbors[0],"C") + pdb.number_of_neighors_of_element(carbon_neighbors[1],"C") ==3: # there should only be one carbon between the two oxygens
                            oxy1 = pdb.all_atoms[carbon_neighbors[0]]
                            oxy2 = pdb.all_atoms[carbon_neighbors[1]]
                            
                            oxy1_neighbors = oxy1.indecies_of_atoms_connecting[:]
                            oxy2_neighbors = oxy2.indecies_of_atoms_connecting[:]
                            
                            oxy1_neighbors.remove(carbon_index)
                            oxy2_neighbors.remove(carbon_index)
                            
                            if oxy1.number_of_neighbors() + oxy2.number_of_neighbors() == 3: # There should be 3 bonds between the two oxygens
                                if oxy1.number_of_neighbors() == 1 and oxy2.number_of_neighbors() == 2:
                                    carbonyl_oxygen_index = carbon_neighbors[0]
                                    non_carbonyl_oxygen_index = carbon_neighbors[1]
                                    distal_carbon_index = oxy2_neighbors[0]
                                else:
                                    carbonyl_oxygen_index = carbon_neighbors[1]
                                    non_carbonyl_oxygen_index = carbon_neighbors[0]
                                    distal_carbon_index = oxy1_neighbors[0]
                                
                                # now check to see if it's a cyclical ester. If so, no go.
                                if pdb.in_same_ring(carbon_index, non_carbonyl_oxygen_index) == "FALSE": # so cyclical esters not allowed
                                    Indices_of_esters.append([carbon_index, carbonyl_oxygen_index, non_carbonyl_oxygen_index, distal_carbon_index])
    
        return Indices_of_esters
    
    def index_of_carbonochloridate(self, pdb): #The function accepts a pdb object as input, the name could perhaps be misconceiving,the molecule is Cl-CO-OR
        """Identifies all the carbonochloridate groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified carbonochloridate group.
        
        """
    
        Indices_of_carbonochloridate = []
    
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
    
            if atom.element=="C":
                carbon_index = atom_index
    
                if pdb.number_of_neighors_of_element(carbon_index,"O") == 2:        # If it's a carbonochloridate, it will have one oxygen attached to it...
    
                    if pdb.number_of_neighors_of_element(carbon_index,"CL") == 1:       #...and one chloride 
    
                        if atom.number_of_neighbors() == 3:
                            carbon_neighbors = atom.indecies_of_atoms_connecting[:]
                            
                            # Identify the Cl
                            for neighbor_index in carbon_neighbors:
                                neighbor_atom = pdb.all_atoms[neighbor_index]
                                if neighbor_atom.element == "CL":
                                    chloride_index = neighbor_index
                                    carbon_neighbors.remove(chloride_index)
                                            
                            # Identify carbonyl oxygen (list now contains only oxygens)
                            for neighbor_index in carbon_neighbors:
                                neighbor_atom = pdb.all_atoms[neighbor_index]
                                if neighbor_atom.number_of_neighbors() == 1:
                                    oxygen_index = neighbor_index
                                    carbon_neighbors.remove(oxygen_index)
                                            
                            # The only thing remaining is the ester oxygen
                            ester_oxygen = carbon_neighbors[0]
    
                            carboxylate_Cl_2_id = [carbon_index, oxygen_index, chloride_index]
                            Indices_of_carbonochloridate.append([ester_oxygen, carbon_index, oxygen_index, chloride_index])
    
        return Indices_of_carbonochloridate
    
    def index_of_acid_anhydride(self, pdb): #The function accepts a pdb object as input
        """Identifies all the acid anhydride groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified acid anhydride group.
        
        """
    
        Indices_of_acid_anhydrides = []
    
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
    
            if atom.element == "O" and atom.number_of_neighbors() == 2 and pdb.number_of_neighors_of_element(atom_index,"C") == 2: # the ether oxygen
                ether_oxygen_index = atom_index
                ether_oxygen_neighbors = atom.indecies_of_atoms_connecting[:]
                carbonyl_carbon_1_index = ether_oxygen_neighbors[0]
                carbonyl_carbon_1 = pdb.all_atoms[carbonyl_carbon_1_index]
                
                carbonyl_carbon_2_index = ether_oxygen_neighbors[1]
                carbonyl_carbon_2 = pdb.all_atoms[carbonyl_carbon_2_index]
                
                if carbonyl_carbon_1.number_of_neighbors() == 3 and carbonyl_carbon_2.number_of_neighbors() == 3:
                    if pdb.number_of_neighors_of_element(carbonyl_carbon_1_index,"C") == 1 and pdb.number_of_neighors_of_element(carbonyl_carbon_2_index,"C") == 1 and pdb.number_of_neighors_of_element(carbonyl_carbon_1_index,"O") == 2 and pdb.number_of_neighors_of_element(carbonyl_carbon_2_index,"O") == 2: # make sure the carbonyls are really carbonyls connected to ether oxygen
                        carbonyl_carbon_1_neighbors = carbonyl_carbon_1.indecies_of_atoms_connecting[:]
                        carbonyl_carbon_2_neighbors = carbonyl_carbon_2.indecies_of_atoms_connecting[:]
                        
                        carbonyl_carbon_1_neighbors.remove(ether_oxygen_index)
                        carbonyl_carbon_2_neighbors.remove(ether_oxygen_index)
                        
                        for neighbor_atom_index in carbonyl_carbon_1_neighbors:
                            neighbor_atom = pdb.all_atoms[neighbor_atom_index]
                            if neighbor_atom.element == "C":
                                carbonyl_carbon_1_neighbors.remove(neighbor_atom_index)
                            if neighbor_atom_index == ether_oxygen_index:
                                carbonyl_carbon_1_neighbors.remove(neighbor_atom_index)
                
                        for neighbor_atom_index in carbonyl_carbon_2_neighbors:
                            neighbor_atom = pdb.all_atoms[neighbor_atom_index]
                            if neighbor_atom.element == "C":
                                carbonyl_carbon_2_neighbors.remove(neighbor_atom_index)
                            if neighbor_atom_index == ether_oxygen_index:
                                carbonyl_carbon_2_neighbors.remove(neighbor_atom_index)
    
                        carbonyl_oxygen_1_index = carbonyl_carbon_1_neighbors[0]
                        carbonyl_oxygen_1 = pdb.all_atoms[carbonyl_oxygen_1_index]
                        
                        carbonyl_oxygen_2_index = carbonyl_carbon_2_neighbors[0]
                        carbonyl_oxygen_2 = pdb.all_atoms[carbonyl_oxygen_2_index]
                        
                        # now make sure the carbonyl oxygens are really carbonyl oxygens
                        if carbonyl_oxygen_1.element == "O" and carbonyl_oxygen_2.element == "O" and carbonyl_oxygen_1.number_of_neighbors() == 1 and carbonyl_oxygen_2.number_of_neighbors() == 1:
                            
                            # now make sure it's not a cyclical acid anhydride
                            if pdb.in_same_ring(ether_oxygen_index, carbonyl_carbon_1_index) == "FALSE":
    
                                Indices_of_acid_anhydrides.append([ether_oxygen_index, carbonyl_carbon_1_index, carbonyl_oxygen_1_index, carbonyl_carbon_2_index, carbonyl_oxygen_2_index])
    
        return Indices_of_acid_anhydrides
    
    
    
    def index_of_isocyanate(self, pdb):
        """Identifies all the isocyanate groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified isocyanate group.
        
        """
    
        # This will return the indices of the atoms, proximal to distal. R-X-N=C=O
            
        isocyanateRoots = []
        
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
            if atom.element=="O" and len(atom.indecies_of_atoms_connecting)==1: # This is the terminal O
                    neighbor1_index = atom.indecies_of_atoms_connecting[0]
                    neighbor1_atom = pdb.all_atoms[neighbor1_index]
                    if neighbor1_atom.element=="C" and len(neighbor1_atom.indecies_of_atoms_connecting)==2: # middle C
                        indicies = []
                        for index in neighbor1_atom.indecies_of_atoms_connecting:
                            indicies.append(index)
                        indicies.remove(atom_index)
                        neighbor2_index = indicies[0]
                        neighbor2_atom = pdb.all_atoms[neighbor2_index]
                        if neighbor2_atom.element == "N" and len(neighbor1_atom.indecies_of_atoms_connecting)==2: # proximal N
                            # they have to be in a straight line
                            angle = atom.coordinates.angle_between_three_points(neighbor1_atom.coordinates, neighbor2_atom.coordinates) * 180 / math.pi
                            if abs(angle-180) < 10: # so it is in a straight line
                                indicies = []
                                for index in neighbor2_atom.indecies_of_atoms_connecting:
                                    indicies.append(index)
                                indicies.remove(neighbor1_index)
                                neighbor3_index = indicies[0] # X described above
                                neighbor3_atom = pdb.all_atoms[neighbor3_index]
                                
                                if neighbor3_atom.element == "C": # So X must be a carbon
    
                                    # you need to get the atom attached to the isocyanate as well
                                    isocyanateRoots.append([neighbor3_index, neighbor2_index, neighbor1_index, atom_index])
        #print isocyanateRoots                    
        return isocyanateRoots
    
    
    def index_of_isothiocyanate(self, pdb):
        """Identifies all the isothiocyanate groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified isothiocyanate group.
        
        """
    
        # This will return the indices of the atoms, proximal to distal. R-X-N=C=S
            
        isothiocyanateRoots = []
        
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
            if atom.element=="S" and len(atom.indecies_of_atoms_connecting)==1: # This is the terminal S
                    neighbor1_index = atom.indecies_of_atoms_connecting[0]
                    neighbor1_atom = pdb.all_atoms[neighbor1_index]
                    if neighbor1_atom.element=="C" and len(neighbor1_atom.indecies_of_atoms_connecting)==2: # middle C
                        indicies = []
                        for index in neighbor1_atom.indecies_of_atoms_connecting:
                            indicies.append(index)
                        indicies.remove(atom_index)
                        neighbor2_index = indicies[0]
                        neighbor2_atom = pdb.all_atoms[neighbor2_index]
                        if neighbor2_atom.element == "N" and len(neighbor1_atom.indecies_of_atoms_connecting)==2: # proximal N
                            # they have to be in a straight line
                            angle = atom.coordinates.angle_between_three_points(neighbor1_atom.coordinates, neighbor2_atom.coordinates) * 180 / math.pi
                            if abs(angle-180) < 10: # so it is in a straight line
                                indicies = []
                                for index in neighbor2_atom.indecies_of_atoms_connecting:
                                    indicies.append(index)
                                indicies.remove(neighbor1_index)
                                neighbor3_index = indicies[0] # X described above
                                neighbor3_atom = pdb.all_atoms[neighbor3_index]
                                
                                if neighbor3_atom.element == "C": # So X must be a carbon
    
                                    # you need to get the atom attached to the isothiocyanate as well
                                    isothiocyanateRoots.append([neighbor3_index, neighbor2_index, neighbor1_index, atom_index])
        return isothiocyanateRoots
    
    def index_of_hydrogen(self, pdb): # Identify an hydrogen
        """Identifies all the hydrogen atoms of a specified molecular model, together with the heavy atoms to which they are bound.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified heavy-atom-hydrogen group.
        
        """
        
        # Think R - X - H (1 - 2 - 3)
        
        hydrogenRoots=[]
        
        for hydrogen_index in pdb.all_atoms:
            hydrogen = pdb.all_atoms[hydrogen_index]
            if hydrogen.element == "H" and hydrogen.number_of_neighbors() == 1:
                linker_index = hydrogen.indecies_of_atoms_connecting[0]
                hydrogenRoots.append([linker_index, hydrogen_index])
    
        return hydrogenRoots
    
    def index_of_tertiary_halide(self, pdb):
        """Identifies all the tertiary halide groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified tertiary halide group.
        
        """
    
        tertiary_halideRoots = []
        
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
            if atom.element=="CL" or atom.element=="BR" or atom.element=="I": # This is the halide
                    neighbor1_index = atom.indecies_of_atoms_connecting[0]
                    neighbor1_atom = pdb.all_atoms[neighbor1_index]
                    if neighbor1_atom.element=="C" and len(neighbor1_atom.indecies_of_atoms_connecting)==4: # The halide is bonded to an sp3-hybridized carbon
                        if pdb.number_of_neighors_of_element(neighbor1_index,"H") == 0: # So there are no hydrogens connected to the carbon (tertiary)
                            carbon_neighbors = neighbor1_atom.indecies_of_atoms_connecting[:]
                            carbon_neighbors.remove(atom_index) # So only non-halide neighbors are present
                            tertiary_halideRoots.append([carbon_neighbors[0], carbon_neighbors[1], carbon_neighbors[2], neighbor1_index, atom_index])
        return tertiary_halideRoots
    
    def index_of_secondary_halide(self, pdb):
        """Identifies all the secondary halide groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified secondary halide group.
        
        """
    
        secondary_halideRoots = []
        
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
            if atom.element=="CL" or atom.element=="BR" or atom.element=="I": # This is the halide
                    neighbor1_index = atom.indecies_of_atoms_connecting[0]
                    neighbor1_atom = pdb.all_atoms[neighbor1_index]
                    if neighbor1_atom.element=="C" and len(neighbor1_atom.indecies_of_atoms_connecting)==4: # The halide is bonded to an sp3-hybridized carbon
                        if pdb.number_of_neighors_of_element(neighbor1_index,"H") == 1: # So there is only one hydrogens connected to the carbon (secondary)
                            carbon_neighbors = neighbor1_atom.indecies_of_atoms_connecting[:]
                            carbon_neighbors.remove(atom_index) # So only non-halide neighbors are present
                            secondary_halideRoots.append([carbon_neighbors[0], carbon_neighbors[1], carbon_neighbors[2], neighbor1_index, atom_index])
        return secondary_halideRoots
    
    def index_of_halide_bound_to_sp2_carbon(self, pdb):
        """Identifies all the halide atoms bound to sp2-hybridized carbon atoms in a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified halide-carbon group.
        
        """
    
        halideRoots = []
        
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
            if atom.element=="CL" or atom.element=="BR" or atom.element=="I": # This is the halide
                    neighbor1_index = atom.indecies_of_atoms_connecting[0]
                    neighbor1_atom = pdb.all_atoms[neighbor1_index]
                    if neighbor1_atom.element=="C" and pdb.hybridization(neighbor1_index) == 2 and neighbor1_atom.number_of_neighbors() == 3: # The halide is bonded to an sp2-hybridized carbon
                        carbon_neighbors = neighbor1_atom.indecies_of_atoms_connecting[:]
                        carbon_neighbors.remove(atom_index) # So only non-halide neighbors are present
                        halideRoots.append([carbon_neighbors[0], carbon_neighbors[1], neighbor1_index, atom_index])
        return halideRoots
    
    
    def index_of_primary_halide(self, pdb):
        """Identifies all the primary halide groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified primary halide group.
        
        """
    
        primary_halideRoots = []
        
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
            if atom.element=="CL" or atom.element=="BR" or atom.element=="I": # This is the halide
                    neighbor1_index = atom.indecies_of_atoms_connecting[0]
                    neighbor1_atom = pdb.all_atoms[neighbor1_index]
                    if neighbor1_atom.element=="C" and len(neighbor1_atom.indecies_of_atoms_connecting)==4: # The halide is bonded to an sp3-hybridized carbon
                        if pdb.number_of_neighors_of_element(neighbor1_index,"H") == 2: # So there is only one hydrogens connected to the carbon (primary)
                            carbon_neighbors = neighbor1_atom.indecies_of_atoms_connecting[:]
                            carbon_neighbors.remove(atom_index) # So only non-halide neighbors are present
                            primary_halideRoots.append([carbon_neighbors[0], carbon_neighbors[1], carbon_neighbors[2], neighbor1_index, atom_index])
        return primary_halideRoots
    
    def index_of_tertiary_alcohol(self, pdb):
        """Identifies all the tertiary alcohol groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified tertiary alcohol group.
        
        """
    
        tertiary_alcoholRoots = []
        
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
            if atom.element=="O" and len(atom.indecies_of_atoms_connecting) == 2 and pdb.number_of_neighors_of_element(atom_index,"H") == 1: # This is the alcohol
                    oxygen_neighbors = atom.indecies_of_atoms_connecting[:]
                    hydroxy_hydrogen_index = pdb.index_of_neighbor_of_element(atom_index,"H")
                    oxygen_neighbors.remove(hydroxy_hydrogen_index)        
                    neighbor1_index = oxygen_neighbors[0] # this is the carbon neighbor
                    neighbor1_atom = pdb.all_atoms[neighbor1_index]
                    if neighbor1_atom.element=="C" and len(neighbor1_atom.indecies_of_atoms_connecting)==4: # The alcohol is bonded to an sp3-hybridized carbon
                        if pdb.number_of_neighors_of_element(neighbor1_index,"H") == 0: # So there are no hydrogens connected to the carbon (tertiary)
                            carbon_neighbors = neighbor1_atom.indecies_of_atoms_connecting[:]
                            carbon_neighbors.remove(atom_index) # So only non-alcohol neighbors are present
                            tertiary_alcoholRoots.append([carbon_neighbors[0], carbon_neighbors[1], carbon_neighbors[2], neighbor1_index, atom_index, hydroxy_hydrogen_index])
        return tertiary_alcoholRoots
    
    def index_of_secondary_alcohol(self, pdb):
        """Identifies all the secondary alcohol groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified secondary alcohol group.
        
        """
            
        secondary_alcoholRoots = []
        
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
            if atom.element=="O" and len(atom.indecies_of_atoms_connecting) == 2 and pdb.number_of_neighors_of_element(atom_index,"H") == 1: # This is the alcohol
                    oxygen_neighbors = atom.indecies_of_atoms_connecting[:]
                    hydroxy_hydrogen_index = pdb.index_of_neighbor_of_element(atom_index,"H")
                    oxygen_neighbors.remove(hydroxy_hydrogen_index)        
                    neighbor1_index = oxygen_neighbors[0] # this is the carbon neighbor
                    neighbor1_atom = pdb.all_atoms[neighbor1_index]
                    if neighbor1_atom.element=="C" and len(neighbor1_atom.indecies_of_atoms_connecting)==4: # The alcohol is bonded to an sp3-hybridized carbon
                        if pdb.number_of_neighors_of_element(neighbor1_index,"H") == 1: # So there are no hydrogens connected to the carbon (secondary)
                            carbon_neighbors = neighbor1_atom.indecies_of_atoms_connecting[:]
                            carbon_neighbors.remove(atom_index) # So only non-alcohol neighbors are present
                            secondary_alcoholRoots.append([carbon_neighbors[0], carbon_neighbors[1], carbon_neighbors[2], neighbor1_index, atom_index, hydroxy_hydrogen_index])
        return secondary_alcoholRoots
    
    def index_of_alcohol_bound_to_sp2_carbon(self, pdb):
        """Identifies all the alcohol groups bound to sp2-hybridized carbons in a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified alcohol group.
        
        """
    
        alcoholRoots = []
        
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
            if atom.element=="O" and len(atom.indecies_of_atoms_connecting) == 2 and pdb.number_of_neighors_of_element(atom_index,"H") == 1: # This is the alcohol
                    oxygen_neighbors = atom.indecies_of_atoms_connecting[:]
                    hydroxy_hydrogen_index = pdb.index_of_neighbor_of_element(atom_index,"H")
                    oxygen_neighbors.remove(hydroxy_hydrogen_index)        
                    neighbor1_index = oxygen_neighbors[0] # this is the carbon neighbor
                    neighbor1_atom = pdb.all_atoms[neighbor1_index]
                    if neighbor1_atom.element=="C" and pdb.hybridization(neighbor1_index) == 2 and neighbor1_atom.number_of_neighbors() == 3: # The alcohol is bonded to an sp2-hybridized carbon
                        carbon_neighbors = neighbor1_atom.indecies_of_atoms_connecting[:]
                        carbon_neighbors.remove(atom_index) # So only non-alcohol neighbors are present
                        alcoholRoots.append([carbon_neighbors[0], carbon_neighbors[1], neighbor1_index, atom_index, hydroxy_hydrogen_index])
        return alcoholRoots
    
    def index_of_primary_alcohol(self, pdb):
        """Identifies all the primary alcohol groups of a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified primary alcohol group.
        
        """
    
        primary_alcoholRoots = []
        
        for atom_index in pdb.all_atoms:
            atom = pdb.all_atoms[atom_index]
            if atom.element=="O" and len(atom.indecies_of_atoms_connecting) == 2 and pdb.number_of_neighors_of_element(atom_index,"H") == 1: # This is the alcohol
                    oxygen_neighbors = atom.indecies_of_atoms_connecting[:]
                    hydroxy_hydrogen_index = pdb.index_of_neighbor_of_element(atom_index,"H")
                    oxygen_neighbors.remove(hydroxy_hydrogen_index)        
                    neighbor1_index = oxygen_neighbors[0] # this is the carbon neighbor
                    neighbor1_atom = pdb.all_atoms[neighbor1_index]
                    if neighbor1_atom.element=="C" and len(neighbor1_atom.indecies_of_atoms_connecting)==4: # The alcohol is bonded to an sp3-hybridized carbon
                        if pdb.number_of_neighors_of_element(neighbor1_index,"H") == 2: # So there are no hydrogens connected to the carbon (primary)
                            carbon_neighbors = neighbor1_atom.indecies_of_atoms_connecting[:]
                            carbon_neighbors.remove(atom_index) # So only non-alcohol neighbors are present
                            primary_alcoholRoots.append([carbon_neighbors[0], carbon_neighbors[1], carbon_neighbors[2], neighbor1_index, atom_index, hydroxy_hydrogen_index])
        return primary_alcoholRoots
    
    def index_of_secondary_amine_sp2_in_ring(self, pdb): # Mostly, AutoClickChem can't deal with ringed structures. This is one that is specific enough that it can be delt with in carbonyl chemistry
        """Identifies all the secondary amine groups in (aromatic) rings in a specified molecular model.
        
        Arguments:
        pdb -- A molecular model (pymolecule.Molecule).
        
        Returns:
        A list of lists, where each list contains the indices (int) of the atoms present in an identified secondary amine group.
        
        """
    
        Indices_of_secondary_amines = []  
        for atom_index in pdb.all_atoms: 
            atom = pdb.all_atoms[atom_index]         
            if atom.element=="N": # So looking for the nitrogen
                nitrogen_index = atom_index     
                if pdb.hybridization(nitrogen_index) == 2: #So it's sp2 hybridized
                    if pdb.number_of_neighors_of_element(nitrogen_index,"H") == 1 and atom.number_of_neighbors() == 3:
                        
                        amine_id = [nitrogen_index]    
        
                        nitrogen_neighbors = atom.indecies_of_atoms_connecting[:]  # here's alist of all the N neighbors
        
                        # Now add in all the non-H neighbors
                        for neighbor_index in nitrogen_neighbors:               
                            neighbor_atom = pdb.all_atoms[neighbor_index]    
                            if neighbor_atom.element != "H": 
                                amine_id.append(neighbor_index) # You only need the index of one of the non-H neighbors
        
                        # Now add the hydrogen atoms
                        for neighbor_index in nitrogen_neighbors:               
                            neighbor_atom = pdb.all_atoms[neighbor_index]    
                            if neighbor_atom.element == "H": 
                                amine_id.append(neighbor_index) # You only need the index of one of the non-H neighbors
                        
                        # Now, make sure it's in a ring
                        if pdb.in_same_ring(nitrogen_index, amine_id[1]) == 'TRUE':
                        
                            # Add the list
                            Indices_of_secondary_amines.append(amine_id)
    
        return Indices_of_secondary_amines 

class StericProblemSolver:
    """This module contains functions that are useful for simultaneously manipulating multiple pymolecule.Molecule objects."""
    
    def reduce_steric_hindrance(self, pdbs):
        """Rotates aligned molecular models so as to reduce steric hindrance.
        
        Arguments:
        pdbs -- A list of triplets, where each triplet contains a pymolecule.Molecule object and two points (pymolecule.Point) forming a line segment about which that pymolecule.Molecule object is rotated.
        
        """
        
        # go through each of the triplets
        for t in range(0,2):
            for triplet in pdbs:
                if len(triplet) > 1: # if it's just 1, you won't be rotating
                    # now rotate this pdb and determine the best angle
                    pdb = triplet[0]
                    point1 = triplet[1]
                    point2 = triplet[2]
                    
                    pdb.set_undo_point()
                    best_angle_rad = 0.0
                    best_score = 1000000000.0
                    for angle in range(0, 360, 30):
                        #index = index + 1
                        angle_rad = angle * math.pi / 180
                        
                        # rotate the Molecule
                        pdb.undo()
                        pdb.rotate_molecule_around_a_line(point1, point2, angle_rad)
                        
                        score = self.__reduce_streic_hinrance_score(pdbs)
                        if score < best_score:
                            best_score = score
                            best_angle_rad= angle_rad
                    pdb.undo()
                    pdb.rotate_molecule_around_a_line(point1, point2, best_angle_rad) # so move the branch to the best position.
                    pdb.set_undo_point()
    
    def __reduce_streic_hinrance_score(self, pdbs): # lower scores are better
        """Determines whether or not a group of aligned pdb objects have steric hindrance.
        
        Arguments:
        pdbs -- A list of triplets, where each triplet contains a pymolecule.Molecule object and two points (pymolecule.Point) forming a line segment about which the pymolecule.Molecule object is rotated.
        
        Returns:
        score -- A float, a measure of the degree of steric hindrance. Lower scores correspond to less steric hindrance.
        
        """
        
        score = 0.0
        
        #pdbs is a listof triplets: triplets (Molecule, point1, point2), where point1 and point2 form a line around which Molecule will be rotated.
        for triplet_index1 in range(0,len(pdbs)-1):
            triplet1 = pdbs[triplet_index1]
            triplet1_pdb = triplet1[0]
            for triplet_index2 in range(triplet_index1 + 1, len(pdbs)): # so these two loops compare all pdbs against all pdbs
                triplet2 = pdbs[triplet_index2]
                triplet2_pdb = triplet2[0]
                for atom1_index in triplet1_pdb.all_atoms: # now these two loops do atom-by-atom comparison between the two pdbs selected
                    atom1 = triplet1_pdb.all_atoms[atom1_index]
                    for atom2_index in triplet2_pdb.all_atoms:
                        atom2 = triplet2_pdb.all_atoms[atom2_index]
                        the_dist = atom1.coordinates.distance_to_another_point(atom2.coordinates)
                        if the_dist < 2.5: # so only worried about close contacts
                            score = score + 2.5 - the_dist # so you want lower scores
        return score

class OperatorsReact:
    """This class contains functions that "react" two molecular models according to the rules of click chemistry."""
    
    structure_locate_groups = StructureLocateGroups()
    structure_pdb_functions = StericProblemSolver()
    
    # Here's where the click chemistry is simulated in silico
    
    def react_molecules(self, pdb1, pdb2, allowed_reaction_types=["azide_and_alkyne_to_azole", "epoxide_alcohol_opening", "epoxide_thiol_opening", "chloroformate_and_amine_to_carbamate", "sulfonyl_azide_and_thio_acid", "carboxylate_and_alcohol_to_ester", "carboxylate_and_thiol_to_thioester", "acyl_halide_and_alcohol_to_ester", "acyl_halide_and_thiol_to_thioester", "ester_and_alcohol_to_ester", "ester_and_thiol_to_thioester", "acid_anhydride_and_alcohol_to_ester", "acid_anhydride_and_thiol_to_thioester", "carboxylate_and_amine_to_amide", "acyl_halide_and_amine_to_amide", "ester_and_amine_to_amide", "acid_anhydride_and_amine_to_amide", "isocyanate_and_amine_to_urea", "isothiocyanate_and_amine_to_thiourea", "isocyanate_and_alcohol_to_carbamate", "isothiocyanate_and_alcohol_to_carbamothioate", "isocyanate_and_thiol_to_carbamothioate", "isothiocyanate_and_thiol_to_carbamodithioate", "alkene_to_epoxide", "halide_to_cyanide", "alcohol_to_cyanide", "carboxylate_to_cyanide", "acyl_halide_to_cyanide", "acid_anhydride_to_cyanide", "halide_to_azide", "alcohol_to_azide", "carboxylate_to_azide", "acyl_halide_to_azide", "acid_anhydride_to_azide", "amine_to_azide", "amine_to_isocyanate", "amine_to_isothiocyanate", "azide_to_amine"]): # pdb1 is the one that will remain stable, pdb2 will be moved.
        """Combines two molecular models into products according to the rules of click chemistry.
        
        Arguments:
        pdb1 -- A molecular model (pymolecule.Molecule).
        pdb2 -- A molecular model (pymolecule.Molecule).
        allowed_reaction_types -- A list of strings identifying which click-chemistry reactions are permitted.
        
        Returns:
        A list of pymolecule.Molecule objects corresponding to the possible products.
        
        """
    
        if len(pdb2.all_atoms)==0:
            pdb2_is_empty = True
        else:
            pdb2_is_empty = False
    
        # first we need to identify what functional groups are on each of these pdbs
        
        azides_pdb1 = self.structure_locate_groups.index_of_azide(pdb1)
        azides_pdb2 = self.structure_locate_groups.index_of_azide(pdb2)
        
        alkyne_pdb1 = self.structure_locate_groups.index_of_alkyne(pdb1)
        alkyne_pdb2 = self.structure_locate_groups.index_of_alkyne(pdb2)
    
        sulfonyl_azide_pdb1 = self.structure_locate_groups.index_of_sulfonyl_azide(pdb1)
        sulfonyl_azide_pdb2 = self.structure_locate_groups.index_of_sulfonyl_azide(pdb2)

        thio_acid_pdb1 = self.structure_locate_groups.index_of_thio_acid(pdb1)
        thio_acid_pdb2 = self.structure_locate_groups.index_of_thio_acid(pdb2)

        alkene_pdb1 = self.structure_locate_groups.index_of_alkene(pdb1)
        alkene_pdb2 = self.structure_locate_groups.index_of_alkene(pdb2)
        
        epoxide_pdb1 = self.structure_locate_groups.index_of_epoxide(pdb1)
        epoxide_pdb2 = self.structure_locate_groups.index_of_epoxide(pdb2)
        
        alcohol_pdb1 = self.structure_locate_groups.index_of_alcohol(pdb1)
        alcohol_pdb2 = self.structure_locate_groups.index_of_alcohol(pdb2)
        
        thiol_pdb1 = self.structure_locate_groups.index_of_thiol(pdb1)
        thiol_pdb2 = self.structure_locate_groups.index_of_thiol(pdb2)
        
        primary_amine_pdb1 = self.structure_locate_groups.index_of_primary_amine(pdb1)
        primary_amine_pdb2 = self.structure_locate_groups.index_of_primary_amine(pdb2)
        
        secondary_amine_pdb1 = self.structure_locate_groups.index_of_secondary_amine(pdb1)
        secondary_amine_pdb2 = self.structure_locate_groups.index_of_secondary_amine(pdb2)
    
        carbonochloridate_pdb1 = self.structure_locate_groups.index_of_carbonochloridate(pdb1)
        carbonochloridate_pdb2 = self.structure_locate_groups.index_of_carbonochloridate(pdb2)
    
        carboxylate_pdb1 = self.structure_locate_groups.index_of_carboxylate(pdb1)
        carboxylate_pdb2 = self.structure_locate_groups.index_of_carboxylate(pdb2)
    
        acylhalide_pdb1 = self.structure_locate_groups.index_of_acylhalide(pdb1)
        acylhalide_pdb2 = self.structure_locate_groups.index_of_acylhalide(pdb2)
    
        ester_pdb1 = self.structure_locate_groups.index_of_ester(pdb1)
        ester_pdb2 = self.structure_locate_groups.index_of_ester(pdb2)
    
        acid_anhydride_pdb1 = self.structure_locate_groups.index_of_acid_anhydride(pdb1)
        acid_anhydride_pdb2 = self.structure_locate_groups.index_of_acid_anhydride(pdb2)
    
        primary_halide_pdb1 = self.structure_locate_groups.index_of_primary_halide(pdb1)
        primary_halide_pdb2 = self.structure_locate_groups.index_of_primary_halide(pdb2)
    
        secondary_halide_pdb1 = self.structure_locate_groups.index_of_secondary_halide(pdb1)
        secondary_halide_pdb2 = self.structure_locate_groups.index_of_secondary_halide(pdb2)
        
        halide_bound_to_sp2_carbon_pdb1 = self.structure_locate_groups.index_of_halide_bound_to_sp2_carbon(pdb1)
        halide_bound_to_sp2_carbon_pdb2 = self.structure_locate_groups.index_of_halide_bound_to_sp2_carbon(pdb2)
        
        tertiary_halide_pdb1 = self.structure_locate_groups.index_of_tertiary_halide(pdb1)
        tertiary_halide_pdb2 = self.structure_locate_groups.index_of_tertiary_halide(pdb2)
    
        isocyanate_pdb1 = self.structure_locate_groups.index_of_isocyanate(pdb1)
        isocyanate_pdb2 = self.structure_locate_groups.index_of_isocyanate(pdb2)
    
        isothiocyanate_pdb1 = self.structure_locate_groups.index_of_isothiocyanate(pdb1)
        isothiocyanate_pdb2 = self.structure_locate_groups.index_of_isothiocyanate(pdb2)
    
        primary_alcohol_pdb1 = self.structure_locate_groups.index_of_primary_alcohol(pdb1)
        primary_alcohol_pdb2 = self.structure_locate_groups.index_of_primary_alcohol(pdb2)
    
        secondary_alcohol_pdb1 = self.structure_locate_groups.index_of_secondary_alcohol(pdb1)
        secondary_alcohol_pdb2 = self.structure_locate_groups.index_of_secondary_alcohol(pdb2)
    
        alcohol_bound_to_sp2_carbon_pdb1 = self.structure_locate_groups.index_of_alcohol_bound_to_sp2_carbon(pdb1)
        alcohol_bound_to_sp2_carbon_pdb2 = self.structure_locate_groups.index_of_alcohol_bound_to_sp2_carbon(pdb2)
    
        tertiary_alcohol_pdb1 = self.structure_locate_groups.index_of_tertiary_alcohol(pdb1)
        tertiary_alcohol_pdb2 = self.structure_locate_groups.index_of_tertiary_alcohol(pdb2)
    
        secondary_amine_sp2_in_ring_pdb1  = self.structure_locate_groups.index_of_secondary_amine_sp2_in_ring(pdb1)
        secondary_amine_sp2_in_ring_pdb2  = self.structure_locate_groups.index_of_secondary_amine_sp2_in_ring(pdb2)
    
        primary_amine_attached_to_S_pdb1  = self.structure_locate_groups.index_of_primary_amine_attached_to_S(pdb1)
        primary_amine_attached_to_S_pdb2  = self.structure_locate_groups.index_of_primary_amine_attached_to_S(pdb2)
    
    
        # Now compile a list of all the interactions.
        possible_reactions = []
           
        if "azide_and_alkyne_to_azole" in allowed_reaction_types: 
            
            if len(azides_pdb1) != 0 and len(alkyne_pdb2) !=0:
                for atoms1 in azides_pdb1:
                    for atoms2 in alkyne_pdb2:
                        possible_reactions.append(["AZIDE", atoms1, "ALKYNE", atoms2])
                        
            if len(alkyne_pdb1) !=0 and len(azides_pdb2) != 0:
                for atoms1 in alkyne_pdb1:
                    for atoms2 in azides_pdb2:
                        possible_reactions.append(["ALKYNE", atoms1, "AZIDE", atoms2])
        
        if "sulfonyl_azide_and_thio_acid" in allowed_reaction_types: 
    
            if len(thio_acid_pdb1) != 0 and len(sulfonyl_azide_pdb2) !=0:
                for atoms1 in thio_acid_pdb1:
                    for atoms2 in sulfonyl_azide_pdb2:
                        possible_reactions.append(["THIO_ACID", atoms1, "SULFONYL_AZIDE", atoms2])
                        
            if len(sulfonyl_azide_pdb1) !=0 and len(thio_acid_pdb2) != 0:
                for atoms1 in sulfonyl_azide_pdb1:
                    for atoms2 in thio_acid_pdb2:
                        possible_reactions.append(["SULFONYL_AZIDE", atoms1, "THIO_ACID", atoms2])
    
        if "chloroformate_and_amine_to_carbamate" in allowed_reaction_types: 
    
            if len(primary_amine_pdb1) !=0 and len(carbonochloridate_pdb2) != 0:
                for atoms1 in primary_amine_pdb1:
                    for atoms2 in carbonochloridate_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "CARBONOCHLORIDATE", atoms2])
            
            if len(carbonochloridate_pdb1) != 0 and len(primary_amine_pdb2) !=0:
                for atoms1 in carbonochloridate_pdb1:
                    for atoms2 in primary_amine_pdb2:
                        possible_reactions.append(["CARBONOCHLORIDATE", atoms1, "AMINE", atoms2])
        
            if len(secondary_amine_pdb1) !=0 and len(carbonochloridate_pdb2) != 0:
                for atoms1 in secondary_amine_pdb1:
                    for atoms2 in carbonochloridate_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "CARBONOCHLORIDATE", atoms2])
            
            if len(carbonochloridate_pdb1) != 0 and len(secondary_amine_pdb2) !=0:
                for atoms1 in carbonochloridate_pdb1:
                    for atoms2 in secondary_amine_pdb2:
                        possible_reactions.append(["CARBONOCHLORIDATE", atoms1, "AMINE", atoms2])
    
            if len(secondary_amine_sp2_in_ring_pdb1) !=0 and len(carbonochloridate_pdb2) != 0:
                for atoms1 in secondary_amine_sp2_in_ring_pdb1:
                    for atoms2 in carbonochloridate_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "CARBONOCHLORIDATE", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
            
            if len(carbonochloridate_pdb1) != 0 and len(secondary_amine_sp2_in_ring_pdb2) !=0:
                for atoms1 in carbonochloridate_pdb1:
                    for atoms2 in secondary_amine_sp2_in_ring_pdb2:
                        possible_reactions.append(["CARBONOCHLORIDATE", atoms1, "AMINE", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
    
        if pdb2_is_empty == True: # so only do these reactions with one reactant if the second pdb is empty
    
            if len(primary_amine_pdb1) != 0 : 
                for atoms1 in primary_amine_pdb1:
                    possible_reactions.append(["PRIMARY_AMINE", atoms1]) # this can be azide, isocyanate, or isothiocyanate, so resolve allowed_reaction_types below
        
            if len(primary_amine_attached_to_S_pdb1) != 0 : 
                for atoms1 in primary_amine_attached_to_S_pdb1:
                    possible_reactions.append(["PRIMARY_AMINE_ATTACHED_TO_S", atoms1]) # this can be azide, isocyanate, or isothiocyanate, so resolve allowed_reaction_types below
        
            if "azide_to_amine" in allowed_reaction_types:
        
                if len(azides_pdb1) != 0 : # Azide to amine
                    for atoms1 in azides_pdb1:
                        possible_reactions.append(["AZIDE", atoms1])
            
                if len(sulfonyl_azide_pdb1) != 0 : # Azide to amine
                    for atoms1 in sulfonyl_azide_pdb1:
                        possible_reactions.append(["SULFONYL_AZIDE", atoms1])
        
            if "alkene_to_epoxide" in allowed_reaction_types:
        
                if len(alkene_pdb1) != 0 :
                    for atoms1 in alkene_pdb1:
                        possible_reactions.append(["ALKENE", atoms1])
    
        if "epoxide_alcohol_opening" in allowed_reaction_types:
    
            if len(epoxide_pdb1) !=0 and len(alcohol_pdb2) != 0:
                for atoms1 in epoxide_pdb1:
                    for atoms2 in alcohol_pdb2:
                        possible_reactions.append(["EPOXIDE", atoms1, "ALCOHOL", atoms2, "O"])
        
            if len(alcohol_pdb1) !=0 and len(epoxide_pdb2) != 0:
                for atoms1 in alcohol_pdb1:
                    for atoms2 in epoxide_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "EPOXIDE", atoms2, "O"])
    
        if "epoxide_thiol_opening" in allowed_reaction_types:
    
            if len(epoxide_pdb1) !=0 and len(thiol_pdb2) != 0:
                for atoms1 in epoxide_pdb1:
                    for atoms2 in thiol_pdb2:
                        possible_reactions.append(["EPOXIDE", atoms1, "ALCOHOL", atoms2, "S"]) # not actually an alcohol, but a thiol. However, processed the same way.
        
            if len(thiol_pdb1) !=0 and len(epoxide_pdb2) != 0:
                for atoms1 in thiol_pdb1:
                    for atoms2 in epoxide_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "EPOXIDE", atoms2, "S"]) # not actually an alcohol, but a thiol. However, processed the same way.
    
        if "carboxylate_and_amine_to_amide" in allowed_reaction_types:
    
            if len(primary_amine_pdb1) !=0 and len(carboxylate_pdb2) != 0:
                for atoms1 in primary_amine_pdb1:
                    for atoms2 in carboxylate_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "CARBOXYLATE", atoms2])
            
            if len(carboxylate_pdb1) != 0 and len(primary_amine_pdb2) !=0:
                for atoms1 in carboxylate_pdb1:
                    for atoms2 in primary_amine_pdb2:
                        possible_reactions.append(["CARBOXYLATE", atoms1, "AMINE", atoms2])
        
            if len(secondary_amine_pdb1) !=0 and len(carboxylate_pdb2) != 0:
                for atoms1 in secondary_amine_pdb1:
                    for atoms2 in carboxylate_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "CARBOXYLATE", atoms2])
            
            if len(carboxylate_pdb1) != 0 and len(secondary_amine_pdb2) !=0:
                for atoms1 in carboxylate_pdb1:
                    for atoms2 in secondary_amine_pdb2:
                        possible_reactions.append(["CARBOXYLATE", atoms1, "AMINE", atoms2])
    
            if len(secondary_amine_sp2_in_ring_pdb1) !=0 and len(carboxylate_pdb2) != 0:
                for atoms1 in secondary_amine_sp2_in_ring_pdb1:
                    for atoms2 in carboxylate_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "CARBOXYLATE", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
            
            if len(carboxylate_pdb1) != 0 and len(secondary_amine_sp2_in_ring_pdb2) !=0:
                for atoms1 in carboxylate_pdb1:
                    for atoms2 in secondary_amine_sp2_in_ring_pdb2:
                        possible_reactions.append(["CARBOXYLATE", atoms1, "AMINE", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
    
        if "acyl_halide_and_amine_to_amide" in allowed_reaction_types:
            
            if len(primary_amine_pdb1) !=0 and len(acylhalide_pdb2) != 0:
                for atoms1 in primary_amine_pdb1:
                    for atoms2 in acylhalide_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ACYLHALIDE", atoms2])
            
            if len(acylhalide_pdb1) != 0 and len(primary_amine_pdb2) !=0:
                for atoms1 in acylhalide_pdb1:
                    for atoms2 in primary_amine_pdb2:
                        possible_reactions.append(["ACYLHALIDE", atoms1, "AMINE", atoms2])
            
            if len(secondary_amine_pdb1) !=0 and len(acylhalide_pdb2) != 0:
                for atoms1 in secondary_amine_pdb1:
                    for atoms2 in acylhalide_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ACYLHALIDE", atoms2])
            
            if len(acylhalide_pdb1) != 0 and len(secondary_amine_pdb2) !=0:
                for atoms1 in acylhalide_pdb1:
                    for atoms2 in secondary_amine_pdb2:
                        possible_reactions.append(["ACYLHALIDE", atoms1, "AMINE", atoms2])
    
            if len(secondary_amine_sp2_in_ring_pdb1) !=0 and len(acylhalide_pdb2) != 0:
                for atoms1 in secondary_amine_sp2_in_ring_pdb1:
                    for atoms2 in acylhalide_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ACYLHALIDE", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
            
            if len(acylhalide_pdb1) != 0 and len(secondary_amine_sp2_in_ring_pdb2) !=0:
                for atoms1 in acylhalide_pdb1:
                    for atoms2 in secondary_amine_sp2_in_ring_pdb2:
                        possible_reactions.append(["ACYLHALIDE", atoms1, "AMINE", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
    
        if "ester_and_amine_to_amide" in allowed_reaction_types:
        
            if len(primary_amine_pdb1) !=0 and len(ester_pdb2) != 0:
                for atoms1 in primary_amine_pdb1:
                    for atoms2 in ester_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ESTER", atoms2])
            
            if len(ester_pdb1) != 0 and len(primary_amine_pdb2) !=0:
                for atoms1 in ester_pdb1:
                    for atoms2 in primary_amine_pdb2:
                        possible_reactions.append(["ESTER", atoms1, "AMINE", atoms2])
            
            if len(secondary_amine_pdb1) !=0 and len(ester_pdb2) != 0:
                for atoms1 in secondary_amine_pdb1:
                    for atoms2 in ester_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ESTER", atoms2])
            
            if len(ester_pdb1) != 0 and len(secondary_amine_pdb2) !=0:
                for atoms1 in ester_pdb1:
                    for atoms2 in secondary_amine_pdb2:
                        possible_reactions.append(["ESTER", atoms1, "AMINE", atoms2])
    
            if len(secondary_amine_sp2_in_ring_pdb1) !=0 and len(ester_pdb2) != 0:
                for atoms1 in secondary_amine_sp2_in_ring_pdb1:
                    for atoms2 in ester_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ESTER", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
            
            if len(ester_pdb1) != 0 and len(secondary_amine_sp2_in_ring_pdb2) !=0:
                for atoms1 in ester_pdb1:
                    for atoms2 in secondary_amine_sp2_in_ring_pdb2:
                        possible_reactions.append(["ESTER", atoms1, "AMINE", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
    
        if "acid_anhydride_and_amine_to_amide" in allowed_reaction_types:
    
            if len(primary_amine_pdb1) !=0 and len(acid_anhydride_pdb2) != 0:
                for atoms1 in primary_amine_pdb1:
                    for atoms2 in acid_anhydride_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ACID_ANHYDRIDE", atoms2])
            
            if len(acid_anhydride_pdb1) != 0 and len(primary_amine_pdb2) !=0:
                for atoms1 in acid_anhydride_pdb1:
                    for atoms2 in primary_amine_pdb2:
                        possible_reactions.append(["ACID_ANHYDRIDE", atoms1, "AMINE", atoms2])
            
            if len(secondary_amine_pdb1) !=0 and len(acid_anhydride_pdb2) != 0:
                for atoms1 in secondary_amine_pdb1:
                    for atoms2 in acid_anhydride_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ACID_ANHYDRIDE", atoms2])
            
            if len(acid_anhydride_pdb1) != 0 and len(secondary_amine_pdb2) !=0:
                for atoms1 in acid_anhydride_pdb1:
                    for atoms2 in secondary_amine_pdb2:
                        possible_reactions.append(["ACID_ANHYDRIDE", atoms1, "AMINE", atoms2])
    
            if len(secondary_amine_sp2_in_ring_pdb1) !=0 and len(acid_anhydride_pdb2) != 0:
                for atoms1 in secondary_amine_sp2_in_ring_pdb1:
                    for atoms2 in acid_anhydride_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ACID_ANHYDRIDE", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
            
            if len(acid_anhydride_pdb1) != 0 and len(secondary_amine_sp2_in_ring_pdb2) !=0:
                for atoms1 in acid_anhydride_pdb1:
                    for atoms2 in secondary_amine_sp2_in_ring_pdb2:
                        possible_reactions.append(["ACID_ANHYDRIDE", atoms1, "AMINE", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
    
        if "carboxylate_and_alcohol_to_ester" in allowed_reaction_types:
    
            if len(alcohol_pdb1) !=0 and len(carboxylate_pdb2) != 0:
                for atoms1 in alcohol_pdb1:
                    for atoms2 in carboxylate_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "CARBOXYLATE", atoms2, "O"])
            
            if len(carboxylate_pdb1) != 0 and len(alcohol_pdb2) !=0:
                for atoms1 in carboxylate_pdb1:
                    for atoms2 in alcohol_pdb2:
                        possible_reactions.append(["CARBOXYLATE", atoms1, "ALCOHOL", atoms2, "O"])
                        
        if "carboxylate_and_thiol_to_thioester" in allowed_reaction_types:
                        
            if len(thiol_pdb1) !=0 and len(carboxylate_pdb2) != 0:
                for atoms1 in thiol_pdb1:
                    for atoms2 in carboxylate_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "CARBOXYLATE", atoms2, "S"])
            
            if len(carboxylate_pdb1) != 0 and len(thiol_pdb2) !=0:
                for atoms1 in carboxylate_pdb1:
                    for atoms2 in thiol_pdb2:
                        possible_reactions.append(["CARBOXYLATE", atoms1, "ALCOHOL", atoms2, "S"])
    
        if "acyl_halide_and_alcohol_to_ester" in allowed_reaction_types:
    
            if len(alcohol_pdb1) !=0 and len(acylhalide_pdb2) != 0:
                for atoms1 in alcohol_pdb1:
                    for atoms2 in acylhalide_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "ACYLHALIDE", atoms2, "O"])
            
            if len(acylhalide_pdb1) != 0 and len(alcohol_pdb2) !=0:
                for atoms1 in acylhalide_pdb1:
                    for atoms2 in alcohol_pdb2:
                        possible_reactions.append(["ACYLHALIDE", atoms1, "ALCOHOL", atoms2, "O"])
    
        if "acyl_halide_and_thiol_to_thioester" in allowed_reaction_types:
    
            if len(thiol_pdb1) !=0 and len(acylhalide_pdb2) != 0:
                for atoms1 in thiol_pdb1:
                    for atoms2 in acylhalide_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "ACYLHALIDE", atoms2, "S"])
            
            if len(acylhalide_pdb1) != 0 and len(thiol_pdb2) !=0:
                for atoms1 in acylhalide_pdb1:
                    for atoms2 in thiol_pdb2:
                        possible_reactions.append(["ACYLHALIDE", atoms1, "ALCOHOL", atoms2, "S"])
    
        if "ester_and_alcohol_to_ester" in allowed_reaction_types:
    
            if len(alcohol_pdb1) !=0 and len(ester_pdb2) != 0:
                for atoms1 in alcohol_pdb1:
                    for atoms2 in ester_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "ESTER", atoms2, "O"])
            
            if len(ester_pdb1) != 0 and len(alcohol_pdb2) !=0:
                for atoms1 in ester_pdb1:
                    for atoms2 in alcohol_pdb2:
                        possible_reactions.append(["ESTER", atoms1, "ALCOHOL", atoms2, "O"])
                    
        if "ester_and_thiol_to_thioester" in allowed_reaction_types:
                    
            if len(thiol_pdb1) !=0 and len(ester_pdb2) != 0:
                for atoms1 in thiol_pdb1:
                    for atoms2 in ester_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "ESTER", atoms2, "S"])
            
            if len(ester_pdb1) != 0 and len(thiol_pdb2) !=0:
                for atoms1 in ester_pdb1:
                    for atoms2 in thiol_pdb2:
                        possible_reactions.append(["ESTER", atoms1, "ALCOHOL", atoms2, "S"])
    
        if "acid_anhydride_and_alcohol_to_ester" in allowed_reaction_types:
    
            if len(alcohol_pdb1) !=0 and len(acid_anhydride_pdb2) != 0:
                for atoms1 in alcohol_pdb1:
                    for atoms2 in acid_anhydride_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "ACID_ANHYDRIDE", atoms2, "O"])
            
            if len(acid_anhydride_pdb1) != 0 and len(alcohol_pdb2) !=0:
                for atoms1 in acid_anhydride_pdb1:
                    for atoms2 in alcohol_pdb2:
                        possible_reactions.append(["ACID_ANHYDRIDE", atoms1, "ALCOHOL", atoms2, "O"])
        
        if "acid_anhydride_and_thiol_to_thioester" in allowed_reaction_types:
            
            if len(thiol_pdb1) !=0 and len(acid_anhydride_pdb2) != 0:
                for atoms1 in thiol_pdb1:
                    for atoms2 in acid_anhydride_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "ACID_ANHYDRIDE", atoms2, "S"])
            
            if len(acid_anhydride_pdb1) != 0 and len(thiol_pdb2) !=0:
                for atoms1 in acid_anhydride_pdb1:
                    for atoms2 in thiol_pdb2:
                        possible_reactions.append(["ACID_ANHYDRIDE", atoms1, "ALCOHOL", atoms2, "S"])
    
        if "isocyanate_and_amine_to_urea" in allowed_reaction_types:
    
            if len(isocyanate_pdb1) !=0 and len(primary_amine_pdb2) != 0:
                for atoms1 in isocyanate_pdb1:
                    for atoms2 in primary_amine_pdb2:
                        possible_reactions.append(["ISOCYANATE", atoms1, "AMINE", atoms2])
            
            if len(primary_amine_pdb1) != 0 and len(isocyanate_pdb2) !=0:
                for atoms1 in primary_amine_pdb1:
                    for atoms2 in isocyanate_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ISOCYANATE", atoms2])
        
            if len(isocyanate_pdb1) !=0 and len(secondary_amine_pdb2) != 0:
                for atoms1 in isocyanate_pdb1:
                    for atoms2 in secondary_amine_pdb2:
                        possible_reactions.append(["ISOCYANATE", atoms1, "AMINE", atoms2])
            
            if len(secondary_amine_pdb1) != 0 and len(isocyanate_pdb2) !=0:
                for atoms1 in secondary_amine_pdb1:
                    for atoms2 in isocyanate_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ISOCYANATE", atoms2])
    
            if len(isocyanate_pdb1) !=0 and len(secondary_amine_sp2_in_ring_pdb2) != 0:
                for atoms1 in isocyanate_pdb1:
                    for atoms2 in secondary_amine_sp2_in_ring_pdb2:
                        possible_reactions.append(["ISOCYANATE", atoms1, "AMINE", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
            
            if len(secondary_amine_sp2_in_ring_pdb1) != 0 and len(isocyanate_pdb2) !=0:
                for atoms1 in secondary_amine_sp2_in_ring_pdb1:
                    for atoms2 in isocyanate_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ISOCYANATE", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
    
        if "isothiocyanate_and_amine_to_thiourea" in allowed_reaction_types:
                    
            if len(isothiocyanate_pdb1) !=0 and len(primary_amine_pdb2) != 0:
                for atoms1 in isothiocyanate_pdb1:
                    for atoms2 in primary_amine_pdb2:
                        possible_reactions.append(["ISOCYANATE", atoms1, "AMINE", atoms2])
            
            if len(primary_amine_pdb1) != 0 and len(isothiocyanate_pdb2) !=0:
                for atoms1 in primary_amine_pdb1:
                    for atoms2 in isothiocyanate_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ISOCYANATE", atoms2])
                        
            if len(isothiocyanate_pdb1) !=0 and len(secondary_amine_pdb2) != 0:
                for atoms1 in isothiocyanate_pdb1:
                    for atoms2 in secondary_amine_pdb2:
                        possible_reactions.append(["ISOCYANATE", atoms1, "AMINE", atoms2])
            
            if len(secondary_amine_pdb1) != 0 and len(isothiocyanate_pdb2) !=0:
                for atoms1 in secondary_amine_pdb1:
                    for atoms2 in isothiocyanate_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ISOCYANATE", atoms2])
    
            if len(isothiocyanate_pdb1) !=0 and len(secondary_amine_sp2_in_ring_pdb2) != 0:
                for atoms1 in isothiocyanate_pdb1:
                    for atoms2 in secondary_amine_sp2_in_ring_pdb2:
                        possible_reactions.append(["ISOCYANATE", atoms1, "AMINE", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
            
            if len(secondary_amine_sp2_in_ring_pdb1) != 0 and len(isothiocyanate_pdb2) !=0:
                for atoms1 in secondary_amine_sp2_in_ring_pdb1:
                    for atoms2 in isothiocyanate_pdb2:
                        possible_reactions.append(["AMINE", atoms1, "ISOCYANATE", atoms2, "SECONDARY_AMINE_SP2_IN_RING"])
    
        if "isocyanate_and_alcohol_to_carbamate" in allowed_reaction_types:
                    
            if len(isocyanate_pdb1) !=0 and len(alcohol_pdb2) != 0:
                for atoms1 in isocyanate_pdb1:
                    for atoms2 in alcohol_pdb2:
                        possible_reactions.append(["ISOCYANATE", atoms1, "ALCOHOL", atoms2])
            
            if len(alcohol_pdb1) != 0 and len(isocyanate_pdb2) !=0:
                for atoms1 in alcohol_pdb1:
                    for atoms2 in isocyanate_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "ISOCYANATE", atoms2])
    
        if "isocyanate_and_thiol_to_carbamothioate" in allowed_reaction_types:
    
            if len(isocyanate_pdb1) !=0 and len(thiol_pdb2) != 0:
                for atoms1 in isocyanate_pdb1:
                    for atoms2 in thiol_pdb2:
                        possible_reactions.append(["ISOCYANATE", atoms1, "ALCOHOL", atoms2])
            
            if len(thiol_pdb1) != 0 and len(isocyanate_pdb2) !=0:
                for atoms1 in thiol_pdb1:
                    for atoms2 in isocyanate_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "ISOCYANATE", atoms2])
                        
        if "isothiocyanate_and_alcohol_to_carbamothioate" in allowed_reaction_types:
                        
            if len(isothiocyanate_pdb1) !=0 and len(alcohol_pdb2) != 0:
                for atoms1 in isothiocyanate_pdb1:
                    for atoms2 in alcohol_pdb2:
                        possible_reactions.append(["ISOCYANATE", atoms1, "ALCOHOL", atoms2])
            
            if len(alcohol_pdb1) != 0 and len(isothiocyanate_pdb2) !=0:
                for atoms1 in alcohol_pdb1:
                    for atoms2 in isothiocyanate_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "ISOCYANATE", atoms2])
                    
        if "isothiocyanate_and_thiol_to_carbamodithioate" in allowed_reaction_types:
    
            if len(isothiocyanate_pdb1) !=0 and len(thiol_pdb2) != 0:
                for atoms1 in isothiocyanate_pdb1:
                    for atoms2 in thiol_pdb2:
                        possible_reactions.append(["ISOCYANATE", atoms1, "ALCOHOL", atoms2])
            
            if len(thiol_pdb1) != 0 and len(isothiocyanate_pdb2) !=0:
                for atoms1 in thiol_pdb1:
                    for atoms2 in isothiocyanate_pdb2:
                        possible_reactions.append(["ALCOHOL", atoms1, "ISOCYANATE", atoms2])
    
        if pdb2_is_empty == True:
        
            if len(primary_halide_pdb1) != 0 : 
                for atoms1 in primary_halide_pdb1:
                    possible_reactions.append(["PRIMARY_HALIDE", atoms1]) # this can go to azide or cyanide, so allowed_reaction_types delt with below
        
            if len(secondary_halide_pdb1) != 0 : 
                for atoms1 in secondary_halide_pdb1:
                    possible_reactions.append(["SECONDARY_HALIDE", atoms1]) # this can go to azide or cyanide, so allowed_reaction_types delt with below
        
            if len(halide_bound_to_sp2_carbon_pdb1) != 0 : 
                for atoms1 in halide_bound_to_sp2_carbon_pdb1:
                    possible_reactions.append(["HALIDE_BOUND_TO_SP2_CARBON", atoms1]) # this can go to azide or cyanide, so allowed_reaction_types delt with below
        
            if len(tertiary_halide_pdb1) != 0 : 
                for atoms1 in tertiary_halide_pdb1:
                    possible_reactions.append(["TERTIARY_HALIDE", atoms1]) # this can go to azide or cyanide, so allowed_reaction_types delt with below
                    
            if len(primary_alcohol_pdb1) != 0 : 
                for atoms1 in primary_alcohol_pdb1:
                    possible_reactions.append(["PRIMARY_ALCOHOL", atoms1]) # this can go to azide or cyanide, so allowed_reaction_types delt with below
        
            if len(secondary_alcohol_pdb1) != 0 : 
                for atoms1 in secondary_alcohol_pdb1:
                    possible_reactions.append(["SECONDARY_ALCOHOL", atoms1]) # this can go to azide or cyanide, so allowed_reaction_types delt with below
                    possible_reactions.append(["INVERTABLE_SECONDARY_ALCOHOL", atoms1]) # this can go to azide or cyanide, so allowed_reaction_types delt with below
        
            if len(alcohol_bound_to_sp2_carbon_pdb1) != 0 : 
                for atoms1 in alcohol_bound_to_sp2_carbon_pdb1:
                    possible_reactions.append(["ALCOHOL_BOUND_TO_SP2_CARBON", atoms1]) # this can go to azide or cyanide, so allowed_reaction_types delt with below
         
            if len(tertiary_alcohol_pdb1) != 0 : 
                for atoms1 in tertiary_alcohol_pdb1:
                    possible_reactions.append(["TERTIARY_ALCOHOL", atoms1]) # this can go to azide or cyanide, so allowed_reaction_types delt with below
        
            if len(carboxylate_pdb1) != 0 : 
                for atoms1 in carboxylate_pdb1:
                    if "carboxylate_to_cyanide" in allowed_reaction_types: possible_reactions.append(["carboxylate_cyanide", atoms1])
                    if "carboxylate_to_azide" in allowed_reaction_types: possible_reactions.append(["carboxylate_azide", atoms1])
        
            if len(acylhalide_pdb1) != 0 : 
                for atoms1 in acylhalide_pdb1:
                    if "acyl_halide_to_cyanide" in allowed_reaction_types: possible_reactions.append(["acylhalide_cyanide", atoms1])
                    if "acyl_halide_to_azide" in allowed_reaction_types: possible_reactions.append(["acylhalide_azide", atoms1])
        
            if len(acid_anhydride_pdb1) != 0 : 
                for atoms1 in acid_anhydride_pdb1:
                    if "acid_anhydride_to_azide" in allowed_reaction_types: possible_reactions.append(["acid_anhydride_azide", atoms1])
                    if "acid_anhydride_to_cyanide" in allowed_reaction_types: possible_reactions.append(["acid_anhydride_cyanide", atoms1])
    
        # Now that you've enumerated all possible reactions, pick one of them
        products = []
        for reaction in possible_reactions:
            if len(reaction) < 4:
                if reaction[0] == "PRIMARY_AMINE": 
                    if "amine_to_azide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__primary_amine_to_azide_or_NCO_or_NCS(pdb1_copy, reaction, 0)
                        product.remarks.append("amine to azide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                    
                    if "amine_to_isocyanate" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__primary_amine_to_azide_or_NCO_or_NCS(pdb1_copy, reaction, 1)
                        product.remarks.append("amine to isocyanate")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                    
                    if "amine_to_isothiocyanate" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__primary_amine_to_azide_or_NCO_or_NCS(pdb1_copy, reaction, 2)
                        product.remarks.append("amine to isothiocyanate")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                elif reaction[0] == "PRIMARY_AMINE_ATTACHED_TO_S": 
                    if "amine_to_azide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__primary_amine_to_azide_or_NCO_or_NCS(pdb1_copy, reaction, 0)
                        product.remarks.append("amine to azide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                    if "amine_to_isocyanate" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__primary_amine_to_azide_or_NCO_or_NCS(pdb1_copy, reaction, 1)
                        product.remarks.append("amine to isocyanate")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                    if "amine_to_isothiocyanate" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__primary_amine_to_azide_or_NCO_or_NCS(pdb1_copy, reaction, 2)
                        product.remarks.append("amine to isothiocyanate")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                # primary azide goes to amine
                elif reaction[0] == "AZIDE": 
                    pdb1_copy = pdb1.copy_of()
                    product = self.__azide_to_primary_amine(pdb1_copy, reaction, 0)
                    product.remarks.append("azide to amine")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
                
                    pdb1_copy = pdb1.copy_of()
                    product = self.__azide_to_primary_amine(pdb1_copy, reaction, 1)
                    product.remarks.append("azide to amine")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
    
                elif reaction[0] == "SULFONYL_AZIDE": 
                    pdb1_copy = pdb1.copy_of()
                    product = self.__azide_to_primary_amine(pdb1_copy, reaction, 0)
                    product.remarks.append("azide to amine")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
                
                    pdb1_copy = pdb1.copy_of()
                    product = self.__azide_to_primary_amine(pdb1_copy, reaction, 1)
                    product.remarks.append("azide to amine")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
                
                elif reaction[0] == "PRIMARY_HALIDE": 
    
                    if "halide_to_azide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__tertiary_or_primary_halide_to_azide(pdb1_copy, reaction, 0)
                        product.remarks.append("halide to azide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                    if "halide_to_cyanide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__tertiary_or_primary_halide_to_azide(pdb1_copy, reaction, 1)
                        product.remarks.append("halide to cyanide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                elif reaction[0] == "TERTIARY_HALIDE": 
    
                    if "halide_to_azide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__tertiary_or_primary_halide_to_azide(pdb1_copy, reaction, 0)
                        product.remarks.append("halide to azide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                    if "halide_to_cyanide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__tertiary_or_primary_halide_to_azide(pdb1_copy, reaction, 1)
                        product.remarks.append("halide to cyanide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                elif reaction[0] == "SECONDARY_HALIDE": 
    
                    if "halide_to_azide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__secondary_halide_to_azide(pdb1_copy, reaction, 0)
                        product.remarks.append("halide to azide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                    if "halide_to_cyanide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__secondary_halide_to_azide(pdb1_copy, reaction, 1)
                        product.remarks.append("halide to cyanide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                elif reaction[0] == "HALIDE_BOUND_TO_SP2_CARBON": 
    
                    if "halide_to_azide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__halide_bound_to_sp2_carbon_to_azide(pdb1_copy, reaction, 0)
                        product.remarks.append("halide to azide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                    if "halide_to_cyanide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__halide_bound_to_sp2_carbon_to_azide(pdb1_copy, reaction, 1)
                        product.remarks.append("halide to cyanide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                    
                elif reaction[0] == "PRIMARY_ALCOHOL": 
    
                    if "alcohol_to_azide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__alcohol_to_azide_same_stereochem(pdb1_copy, reaction, 0)
                        product.remarks.append("alcohol to azide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                    if "alcohol_to_cyanide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__alcohol_to_azide_same_stereochem(pdb1_copy, reaction, 1)
                        product.remarks.append("alcohol to cyanide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                elif reaction[0] == "TERTIARY_ALCOHOL": 
    
                    if "alcohol_to_azide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__alcohol_to_azide_same_stereochem(pdb1_copy, reaction, 0)
                        product.remarks.append("alcohol to azide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                    if "alcohol_to_cyanide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__alcohol_to_azide_same_stereochem(pdb1_copy, reaction, 1)
                        product.remarks.append("alcohol to cyanide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                elif reaction[0] == "SECONDARY_ALCOHOL": 
    
                    if "alcohol_to_azide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__alcohol_to_azide_same_stereochem(pdb1_copy, reaction, 0)
                        product.remarks.append("alcohol to azide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                    if "alcohol_to_cyanide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__alcohol_to_azide_same_stereochem(pdb1_copy, reaction, 1)
                        product.remarks.append("alcohol to cyanide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                elif reaction[0] == "ALCOHOL_BOUND_TO_SP2_CARBON": 
    
                    if "alcohol_to_azide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__alcohol_bound_to_sp2_carbon_to_azide(pdb1_copy, reaction, 0)
                        product.remarks.append("alcohol to azide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                    if "alcohol_to_cyanide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__alcohol_bound_to_sp2_carbon_to_azide(pdb1_copy, reaction, 1)
                        product.remarks.append("alcohol to cyanide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                    
                elif reaction[0] == "INVERTABLE_SECONDARY_ALCOHOL": 
    
                    if "alcohol_to_azide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__invertable_secondary_alcohol_to_azide(pdb1_copy, reaction, 0)
                        product.remarks.append("alcohol to azide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                    if "alcohol_to_cyanide" in allowed_reaction_types:
                        pdb1_copy = pdb1.copy_of()
                        product = self.__invertable_secondary_alcohol_to_azide(pdb1_copy, reaction, 1)
                        product.remarks.append("alcohol to cyanide")
                        product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                        products.append(product)
                
                elif reaction[0] == "carboxylate_cyanide": # SECONDARY_AMINE_SP2_IN_RING
                    pdb1_copy = pdb1.copy_of()
                    product = self.__carboxylate_cyanide(pdb1_copy, reaction)
                    product.remarks.append("carboxylate to cyanide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
                    
                elif reaction[0] == "carboxylate_azide": 
                    pdb1_copy = pdb1.copy_of()
                    product = self.__carboxylate_azide(pdb1_copy, reaction)
                    product.remarks.append("carboxylate to azide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
                    
                elif reaction[0] == "acylhalide_azide": 
                    pdb1_copy = pdb1.copy_of()
                    product = self.__acylhalide_azide(pdb1_copy, reaction)
                    product.remarks.append("carboxylate to cyanide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
                
                elif reaction[0] == "acylhalide_cyanide": 
                    pdb1_copy = pdb1.copy_of()
                    product = self.__acylhalide_cyanide(pdb1_copy, reaction)
                    product.remarks.append("acyl halide to cyanide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
                
                elif reaction[0] == "acid_anhydride_azide": 
                    pdb1_copy = pdb1.copy_of()
                    product = self.__acid_anhydride_azide(pdb1_copy, reaction, 0)
                    product.remarks.append("acid anhydride to azide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
                
                    pdb1_copy = pdb1.copy_of()
                    product = self.__acid_anhydride_azide(pdb1_copy, reaction, 1)
                    product.remarks.append("acid anhydride to azide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
                
                elif reaction[0] == "acid_anhydride_cyanide": 
                    pdb1_copy = pdb1.copy_of()
                    product = self.__acid_anhydride_cyanide(pdb1_copy, reaction, 0)
                    product.remarks.append("acid anhydride to cyanide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
                
                    pdb1_copy = pdb1.copy_of()
                    product = self.__acid_anhydride_cyanide(pdb1_copy, reaction, 1)
                    product.remarks.append("acid anhydride to cyanide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
    
                # alkene goes to an epoxide
                elif reaction[0] == "ALKENE":
                    pdb1_copy = pdb1.copy_of()
                    product = self.__alkene_to_epoxide(pdb1_copy, reaction, 0)
                    product.remarks.append("alkene to epoxide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    product = self.__alkene_to_epoxide(pdb1_copy, reaction, 1)
                    product.remarks.append("alkene to epoxide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    product = self.__alkene_to_epoxide(pdb1_copy, reaction, 2)
                    product.remarks.append("alkene to epoxide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    product = self.__alkene_to_epoxide(pdb1_copy, reaction, 3)
                    product.remarks.append("alkene to epoxide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename)
                    products.append(product)
        
            elif len(reaction) >= 4:
                if (reaction[0] == "EPOXIDE" and reaction[2] == "ALCOHOL") or (reaction[0] == "ALCOHOL" and reaction[2] == "EPOXIDE"): # this actually accounts for alcohol and thiol # This has been updated
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__epoxide_alcohol(pdb1_copy, pdb2_copy, reaction, 0, 0)
                    if "O" in reaction:
                        product.remarks.append("epoxide opening by alcohol")
                    else:
                        product.remarks.append("epoxide opening by thiol")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
                    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__epoxide_alcohol(pdb1_copy, pdb2_copy, reaction, 1, 0)
                    if "O" in reaction:
                        product.remarks.append("epoxide opening by alcohol")
                    else:
                        product.remarks.append("epoxide opening by thiol")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
                    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__epoxide_alcohol(pdb1_copy, pdb2_copy, reaction, 2, 0)
                    if "O" in reaction:
                        product.remarks.append("epoxide opening by alcohol")
                    else:
                        product.remarks.append("epoxide opening by thiol")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
                    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__epoxide_alcohol(pdb1_copy, pdb2_copy, reaction, 3, 0)
                    if "O" in reaction:
                        product.remarks.append("epoxide opening by alcohol")
                    else:
                        product.remarks.append("epoxide opening by thiol")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
                    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__epoxide_alcohol(pdb1_copy, pdb2_copy, reaction, 0, 1)
                    if "O" in reaction:
                        product.remarks.append("epoxide opening by alcohol")
                    else:
                        product.remarks.append("epoxide opening by thiol")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
                    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__epoxide_alcohol(pdb1_copy, pdb2_copy, reaction, 1, 1)
                    if "O" in reaction:
                        product.remarks.append("epoxide opening by alcohol")
                    else:
                        product.remarks.append("epoxide opening by thiol")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
                    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__epoxide_alcohol(pdb1_copy, pdb2_copy, reaction, 2, 1)
                    if "O" in reaction:
                        product.remarks.append("epoxide opening by alcohol")
                    else:
                        product.remarks.append("epoxide opening by thiol")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
                    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__epoxide_alcohol(pdb1_copy, pdb2_copy, reaction, 3, 1)
                    if "O" in reaction:
                        product.remarks.append("epoxide opening by alcohol")
                    else:
                        product.remarks.append("epoxide opening by thiol")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                # sulfonyl_azide - thio acid reactions
                elif (reaction[0] == "SULFONYL_AZIDE" and reaction[2] == "THIO_ACID") or (reaction[0] == "THIO_ACID" and reaction[2] == "SULFONYL_AZIDE"): # this one is updated
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__sulfonyl_azide_thio_acid(pdb1_copy, pdb2_copy, reaction, 0)
                    product.remarks.append("sulfonyl azide + thio acid => ???")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__sulfonyl_azide_thio_acid(pdb1_copy, pdb2_copy, reaction, 1)
                    product.remarks.append("sulfonyl azide + thio acid => ???")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                # Azide alkyne Huisgen cycloaddition reactions
                elif (reaction[0] == "ALKYNE" and reaction[2] == "AZIDE") or (reaction[0] == "AZIDE" and reaction[2] == "ALKYNE"): # This one is updated.
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__alkyne_azide(pdb1_copy, pdb2_copy, reaction, 0)
                    product.remarks.append("alkyne + azide => 1,2,3-triazole")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__alkyne_azide(pdb1_copy, pdb2_copy, reaction, 1)
                    product.remarks.append("alkyne + azide => 1,2,3-triazole")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                # amine and carbonochloridate
                elif (reaction[0] == "AMINE" and reaction[2] == "CARBONOCHLORIDATE") or (reaction[0] == "CARBONOCHLORIDATE" and reaction[2] == "AMINE"): # This one is updated.
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__carbonochloridate_amine(pdb1_copy, pdb2_copy, reaction, 0)
                    product.remarks.append("carbonochloridate + amine => carbamate")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__carbonochloridate_amine(pdb1_copy, pdb2_copy, reaction, 1)
                    product.remarks.append("carbonochloridate + amine => carbamate")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                # amine and carboxylate
                elif (reaction[0] == "AMINE" and reaction[2] == "CARBOXYLATE") or (reaction[0] == "CARBOXYLATE" and reaction[2] == "AMINE"): # This one is updated.
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__carboxylate_amine(pdb1_copy, pdb2_copy, reaction, 0)
                    product.remarks.append("carboxylate + amine => amide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__carboxylate_amine(pdb1_copy, pdb2_copy, reaction, 1)
                    product.remarks.append("carboxylate + amine => amide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                # amine and acylhalide
                elif (reaction[0] == "AMINE" and reaction[2] == "ACYLHALIDE") or (reaction[0] == "ACYLHALIDE" and reaction[2] == "AMINE"): # This one is updated.
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__acylhalide_amine(pdb1_copy, pdb2_copy, reaction, 0)
                    product.remarks.append("acyl halide + amine => amide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__acylhalide_amine(pdb1_copy, pdb2_copy, reaction, 1)
                    product.remarks.append("acyl halide + amine => amide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                # amine and ester
                elif (reaction[0] == "AMINE" and reaction[2] == "ESTER") or (reaction[0] == "ESTER" and reaction[2] == "AMINE"): # This one is updated.
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__ester_amine(pdb1_copy, pdb2_copy, reaction, 0)
                    product.remarks.append("ester + amine => amide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__ester_amine(pdb1_copy, pdb2_copy, reaction, 1)
                    product.remarks.append("ester + amine => amide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                # amine and acid_anhydride
                elif (reaction[0] == "AMINE" and reaction[2] == "ACID_ANHYDRIDE") or (reaction[0] == "ACID_ANHYDRIDE" and reaction[2] == "AMINE"): # This one is updated.
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__acid_anhydride_amine(pdb1_copy, pdb2_copy, reaction, 0, 0)
                    product.remarks.append("acid anhydride + amine => amide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__acid_anhydride_amine(pdb1_copy, pdb2_copy, reaction, 1, 0)
                    product.remarks.append("acid anhydride + amine => amide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__acid_anhydride_amine(pdb1_copy, pdb2_copy, reaction, 0, 1)
                    product.remarks.append("acid anhydride + amine => amide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__acid_anhydride_amine(pdb1_copy, pdb2_copy, reaction, 1, 1)
                    product.remarks.append("acid anhydride + amine => amide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                elif (reaction[0] == "ALCOHOL" and reaction[2] == "CARBOXYLATE") or (reaction[0] == "CARBOXYLATE" and reaction[2] == "ALCOHOL"): # This is for both alcohols and thiols
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__carboxylate_alcohol(pdb1_copy, pdb2_copy, reaction)
                    if "O" in reaction:
                        product.remarks.append("carboxylate + alcohol => ester")
                    else:
                        product.remarks.append("carboxylate + thiol => thioester")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                elif (reaction[0] == "ALCOHOL" and reaction[2] == "ACYLHALIDE") or (reaction[0] == "ACYLHALIDE" and reaction[2] == "ALCOHOL"): # This is for both alcohols and thiols
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__acylhalide_alcohol(pdb1_copy, pdb2_copy, reaction)
                    if "O" in reaction:
                        product.remarks.append("acyl halide + alcohol => ester")
                    else:
                        product.remarks.append("acyl halide + thiol => thioester")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                elif (reaction[0] == "ALCOHOL" and reaction[2] == "ESTER") or (reaction[0] == "ESTER" and reaction[2] == "ALCOHOL"): # This is for both alcohols and thiols
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__ester_alcohol(pdb1_copy, pdb2_copy, reaction)
                    if "O" in reaction:
                        product.remarks.append("ester + alcohol => ester")
                    else:
                        product.remarks.append("ester + thiol => thioester")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                elif (reaction[0] == "ALCOHOL" and reaction[2] == "ACID_ANHYDRIDE") or (reaction[0] == "ACID_ANHYDRIDE" and reaction[2] == "ALCOHOL"): # This is for both alcohols and thiols
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__acid_anhydride_alcohol(pdb1_copy, pdb2_copy, reaction, 0)
                    if "O" in reaction:
                        product.remarks.append("acid anhydride + alcohol => ester")
                    else:
                        product.remarks.append("acid anhydride + thiol => thioester")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__acid_anhydride_alcohol(pdb1_copy, pdb2_copy, reaction, 1)
                    if "O" in reaction:
                        product.remarks.append("acid anhydride + alcohol => ester")
                    else:
                        product.remarks.append("acid anhydride + thiol => thioester")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                elif (reaction[0] == "ISOCYANATE" and reaction[2] == "AMINE") or (reaction[0] == "AMINE" and reaction[2] == "ISOCYANATE"): # This works for isocyanates, isothiocyanates, primary and secondary amines
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.isocyanate_amine(pdb1_copy, pdb2_copy, reaction, 0, 0)
                    product.remarks.append("isocyanate + amine => carbamide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.isocyanate_amine(pdb1_copy, pdb2_copy, reaction, 1, 0)
                    product.remarks.append("isocyanate + amine => carbamide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.isocyanate_amine(pdb1_copy, pdb2_copy, reaction, 0, 1)
                    product.remarks.append("isocyanate + amine => carbamide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.isocyanate_amine(pdb1_copy, pdb2_copy, reaction, 1, 1)
                    product.remarks.append("isocyanate + amine => carbamide")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                elif (reaction[0] == "ISOCYANATE" and reaction[2] == "ALCOHOL") or (reaction[0] == "ALCOHOL" and reaction[2] == "ISOCYANATE"): # This works for isocyanates, isothiocyanates, alcohols, thiols.
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__isocyanate_alcohol(pdb1_copy, pdb2_copy, reaction, 0)
                    product.remarks.append("isocyanate + amine => carbamate")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
    
                    pdb1_copy = pdb1.copy_of()
                    pdb2_copy = pdb2.copy_of()
                    product = self.__isocyanate_alcohol(pdb1_copy, pdb2_copy, reaction, 1)
                    product.remarks.append("isocyanate + amine => carbamate")
                    product.remarks.append("SOURCE FILES: " + pdb1_copy.filename + "; " + pdb2_copy.filename)
                    products.append(product)
        
        return products
    
    def __epoxide_alcohol(self, pdb1, pdb2, reaction, chiral_choice, hydroxyl_side_choice): # this actually accounts for alcohol and thiol # This has been updated
        """Simulates the reaction between an epoxide and an alcohol/thiol.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), either an epoxide or alcohol/thiol.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        chiral_choice -- As multiple chiral centers are possible, this int dictates which product should be made.
        hydroxyl_side_choice -- As the hydroxyl oxygen atom can open to either side of the epoxide, this int specifies which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        # In some other reactions, I tried to preserve the position of a number of atoms so I could develop a crossover operator like that used in AutoGrow
        # I realize now that such a crossover could lead to compounds that could not be formulated with click chemistry, so I'm abandoning these efforts.
    
        #            O
        #          /   \
        #    R1 - C  -  C - R3
        #        /       \
        #      R2         R4
    
        # do choices about chirality and such need to be made
        if chiral_choice == -1: chiral_choice = random.randrange(0,4)
        if hydroxyl_side_choice == -1: hydroxyl_side_choice = random.randrange(0,2)
    
        # first, identify which is the epoxide and which is the alcohol
        if reaction[0] == "EPOXIDE":
            epoxide = pdb1
            epoxide_indices = reaction[1]
            alcohol = pdb2
            alcohol_indices =  reaction[3]
        else:
            epoxide = pdb2
            epoxide_indices = reaction[3]
            alcohol = pdb1
            alcohol_indices =  reaction[1]
    
        # let's make the names easy        
        carbon_one_index = epoxide_indices[1]
        carbon_one_neighbor_one_index = epoxide_indices[2]
        carbon_one_neighbor_two_index = epoxide_indices[3]
        
        carbon_two_index = epoxide_indices[4]
        carbon_two_neighbor_one_index = epoxide_indices[5]
        carbon_two_neighbor_two_index = epoxide_indices[6]
        
        # now we first need to define the various fragments
        carbon_one_frag_one = epoxide.get_branch(carbon_one_index, carbon_one_neighbor_one_index)
        carbon_one_frag_two = epoxide.get_branch(carbon_one_index, carbon_one_neighbor_two_index)
        carbon_two_frag_one = epoxide.get_branch(carbon_two_index, carbon_two_neighbor_one_index)
        carbon_two_frag_two = epoxide.get_branch(carbon_two_index, carbon_two_neighbor_two_index)
    
        # Now, randomly pick which fragments will be R1 and R2, and which will be R3 and R4. This effectively randomizes the chirality at both centers.
        if chiral_choice == 0:
            R1 = carbon_one_frag_one
            R2 = carbon_one_frag_two
            R3 = carbon_two_frag_one
            R4 = carbon_two_frag_two
            
            R1_aim = epoxide_indices[2]
            R2_aim = epoxide_indices[3]
            R3_aim = epoxide_indices[5]
            R4_aim = epoxide_indices[6]
            
        elif chiral_choice == 1:
            R2 = carbon_one_frag_one
            R1 = carbon_one_frag_two
            R3 = carbon_two_frag_one
            R4 = carbon_two_frag_two
    
            R2_aim = epoxide_indices[2]
            R1_aim = epoxide_indices[3]
            R3_aim = epoxide_indices[5]
            R4_aim = epoxide_indices[6]
    
        elif chiral_choice == 2:
            R1 = carbon_one_frag_one
            R2 = carbon_one_frag_two
            R4 = carbon_two_frag_one
            R3 = carbon_two_frag_two
    
            R1_aim = epoxide_indices[2]
            R2_aim = epoxide_indices[3]
            R4_aim = epoxide_indices[5]
            R3_aim = epoxide_indices[6]
            
        elif chiral_choice == 3:
            R2 = carbon_one_frag_one
            R1 = carbon_one_frag_two
            R4 = carbon_two_frag_one
            R3 = carbon_two_frag_two
            
            R2_aim = epoxide_indices[2]
            R1_aim = epoxide_indices[3]
            R4_aim = epoxide_indices[5]
            R3_aim = epoxide_indices[6]
    
        # okay, now let's load in the intermediate
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb('./intermediates/intermediate3.pdb')
        
        # okay, now let's load in the hydroxyl group
        hydroxyl = pymolecule.Molecule()
        hydroxyl.load_pdb('./intermediates/intermediate5.pdb')
        
        # okay, according to the Sharpless paper "A Few Good Reactions" this can go one of four ways (hydroxyl to right or left)
        intermediate_R1_pivot = 1           # always the same
        intermediate_R1_aim = 7             # always the same
        intermediate_R2_pivot = 1           # always the same
        intermediate_R3_pivot = 2           # always the same
        intermediate_R3_aim = 3             # always the same
        intermediate_R4_pivot = 2           # always the same
    
        if hydroxyl_side_choice == 0:
            intermediate_hydroxyl_pivot = 1
            intermediate_hydroxyl_aim = 8
            intermediate_R2_aim = 6
            intermediate_nucleophile_pivot = 2
            intermediate_nucleophile_aim = 4
            intermediate_R4_aim = 5
        elif hydroxyl_side_choice == 1:
            intermediate_hydroxyl_pivot = 2
            intermediate_hydroxyl_aim = 5
            intermediate_R2_aim = 8
            intermediate_nucleophile_pivot = 1
            intermediate_nucleophile_aim = 6
            intermediate_R4_aim = 4
    
        # Okay, start aligning the fragments to the intermediate.
        # Position R1
        tethers = [[intermediate_R1_pivot, carbon_one_index], [intermediate_R1_aim, R1_aim]]
        R1 = intermediate.align_another_molecule_to_this_one(R1, tethers)
        
        # Position R2
        tethers = [[intermediate_R2_pivot, carbon_one_index], [intermediate_R2_aim, R2_aim]] 
        R2 = intermediate.align_another_molecule_to_this_one(R2, tethers)
        
        # Position R3
        tethers = [[intermediate_R3_pivot, carbon_two_index], [intermediate_R3_aim, R3_aim]] 
        R3 = intermediate.align_another_molecule_to_this_one(R3, tethers)
        
        # Position R4
        tethers = [[intermediate_R4_pivot, carbon_two_index], [intermediate_R4_aim, R4_aim]] 
        R4 = intermediate.align_another_molecule_to_this_one(R4, tethers)
        
        # position the nucleophile (the alcohol)
        
        # first load in an intermediate
        nucleophile_heteroatom = reaction[4]
        ether_intermediate = pymolecule.Molecule()
        if nucleophile_heteroatom == "O":
            ether_intermediate.load_pdb('./intermediates/intermediate4.pdb')
        elif nucleophile_heteroatom == "S":
            ether_intermediate.load_pdb('./intermediates/intermediate7.pdb')
        
        # position the intermediate
        tethers = [[intermediate_nucleophile_pivot, 1], [intermediate_nucleophile_aim, 2]] 
        ether_intermediate = intermediate.align_another_molecule_to_this_one(ether_intermediate, tethers)
    
        # position the nucleophile (the alcohol)
        tethers = [[2, alcohol_indices[1]], [3, alcohol_indices[0]], [1, alcohol_indices[2]]] 
        alcohol = ether_intermediate.align_another_molecule_to_this_one(alcohol, tethers)
        
        # position the hydroxyl group
        tethers = [[intermediate_hydroxyl_pivot, 1], [intermediate_hydroxyl_aim, 2]] 
        hydroxyl = intermediate.align_another_molecule_to_this_one(hydroxyl, tethers)
    
        # Now minimize steric problems
        item1 = [R1, intermediate.all_atoms[intermediate_R1_pivot].coordinates, intermediate.all_atoms[intermediate_R1_aim].coordinates]
        item2 = [R2, intermediate.all_atoms[intermediate_R2_pivot].coordinates, intermediate.all_atoms[intermediate_R2_aim].coordinates]
        item3 = [R3, intermediate.all_atoms[intermediate_R3_pivot].coordinates, intermediate.all_atoms[intermediate_R3_aim].coordinates]
        item4 = [R4, intermediate.all_atoms[intermediate_R4_pivot].coordinates, intermediate.all_atoms[intermediate_R4_aim].coordinates]
        item5 = [alcohol, alcohol.all_atoms[alcohol_indices[1]].coordinates, alcohol.all_atoms[alcohol_indices[2]].coordinates] # right?
        item6 = [hydroxyl, hydroxyl.all_atoms[1].coordinates, hydroxyl.all_atoms[2].coordinates] 
        item7 = [intermediate]
        
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        thelist.append(item3)
        thelist.append(item4)
        thelist.append(item5)
        thelist.append(item6)
        thelist.append(item7)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        # now delete some overlapping atoms
        intermediate.delete_atom(7)
        intermediate.delete_atom(3)
        intermediate.delete_atom(4)
        intermediate.delete_atom(5)
        intermediate.delete_atom(6)
        intermediate.delete_atom(8)
        R1.delete_atom(carbon_one_index)
        R2.delete_atom(carbon_one_index)
        R3.delete_atom(carbon_two_index)
        R4.delete_atom(carbon_two_index)
        hydroxyl.delete_atom(1)
        alcohol.delete_atom(alcohol_indices[2])
        
        intermediate.change_residue('FR1')
        R1.change_residue('FR2')
        R2.change_residue('FR3')
        R3.change_residue('FR4')
        R4.change_residue('FR5')
        alcohol.change_residue('FR6')
        hydroxyl.change_residue('FR7')
        
        build = intermediate.merge_with_another_molecule(R1)
        build = build.merge_with_another_molecule(R2)
        build = build.merge_with_another_molecule(R3)
        build = build.merge_with_another_molecule(R4)
        build = build.merge_with_another_molecule(alcohol)
        build = build.merge_with_another_molecule(hydroxyl)
        
        return build
    
    def __alkene_to_epoxide(self, pdb1, reaction, chiral_choice):
        """Simulates the oxidation of an alkene to an epoxide.
        
        Arguments:
        pdb1 -- A molecular model (pymolecule.Molecule) containing an alkene.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        chiral_choice -- As multiple chiral centers are possible, this int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        # In some other reactions, I tried to preserve the position of a number of atoms so I could develop a crossover operator like that used in AutoGrow
        # I realize now that such a crossover could lead to compounds that could not be formulated with click chemistry, so I'm abandoning these efforts.
    
        #  R1        R3
        #   \       /
        #    C1 = C2
        #   /       \
        # R2         R4
    
        # choose the chirality
        if chiral_choice == -1: chiral_choice = random.randrange(0,4)
    
        # first, identify the alkene
        alkene = pdb1
        alkene_indices = reaction[1]
            
        # let's make the names easy        
        carbon_one_index = alkene_indices[0]
        carbon_one_neighbor_one_index = alkene_indices[1]
        carbon_one_neighbor_two_index = alkene_indices[2]
        
        carbon_two_index = alkene_indices[3]
        carbon_two_neighbor_one_index = alkene_indices[4]
        carbon_two_neighbor_two_index = alkene_indices[5]
        
        # now we first need to define the various fragments
        carbon_one_frag_one = alkene.get_branch(carbon_one_index, carbon_one_neighbor_one_index)
        carbon_one_frag_two = alkene.get_branch(carbon_one_index, carbon_one_neighbor_two_index)
        carbon_two_frag_one = alkene.get_branch(carbon_two_index, carbon_two_neighbor_one_index)
        carbon_two_frag_two = alkene.get_branch(carbon_two_index, carbon_two_neighbor_two_index)
    
        # Now, randomly pick which fragments will be R1 and R2, and which will be R3 and R4
        if chiral_choice == 0:
            R1 = carbon_one_frag_one
            R2 = carbon_one_frag_two
            R3 = carbon_two_frag_one
            R4 = carbon_two_frag_two
            
            R1_aim = alkene_indices[1]
            R2_aim = alkene_indices[2]
            R3_aim = alkene_indices[4]
            R4_aim = alkene_indices[5]
            
        elif chiral_choice == 1:
            R2 = carbon_one_frag_one
            R1 = carbon_one_frag_two
            R3 = carbon_two_frag_one
            R4 = carbon_two_frag_two
    
            R2_aim = alkene_indices[1]
            R1_aim = alkene_indices[2]
            R3_aim = alkene_indices[4]
            R4_aim = alkene_indices[5]
    
        elif chiral_choice == 2:
            R1 = carbon_one_frag_one
            R2 = carbon_one_frag_two
            R4 = carbon_two_frag_one
            R3 = carbon_two_frag_two
    
            R1_aim = alkene_indices[1]
            R2_aim = alkene_indices[2]
            R4_aim = alkene_indices[4]
            R3_aim = alkene_indices[5]
            
        elif chiral_choice == 3:
            R2 = carbon_one_frag_one
            R1 = carbon_one_frag_two
            R4 = carbon_two_frag_one
            R3 = carbon_two_frag_two
            
            R2_aim = alkene_indices[1]
            R1_aim = alkene_indices[2]
            R4_aim = alkene_indices[4]
            R3_aim = alkene_indices[5]
    
        # Now load in the intermediate (a generic epoxide)
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb("./intermediates/generic_epoxide.pdb")
        
        # Now align the fragments to the intermediate model
        tethers = [[1,carbon_one_index], [5,R1_aim]]
        R1 = intermediate.align_another_molecule_to_this_one(R1, tethers)
        
        tethers = [[1,carbon_one_index], [4,R2_aim]]
        R2 = intermediate.align_another_molecule_to_this_one(R2, tethers)
    
        tethers = [[3,carbon_two_index], [6,R3_aim]]
        R3 = intermediate.align_another_molecule_to_this_one(R3, tethers)
        
        tethers = [[3,carbon_two_index], [7,R4_aim]]
        R4 = intermediate.align_another_molecule_to_this_one(R4, tethers)
        
        # Now minimize steric problems
        item1 = [R1, R1.all_atoms[carbon_one_index].coordinates, R1.all_atoms[R1_aim].coordinates]
        item2 = [R2, R2.all_atoms[carbon_one_index].coordinates, R2.all_atoms[R2_aim].coordinates]
        item3 = [R3, R3.all_atoms[carbon_two_index].coordinates, R3.all_atoms[R3_aim].coordinates]
        item4 = [R4, R4.all_atoms[carbon_two_index].coordinates, R4.all_atoms[R4_aim].coordinates]
        item5 = [intermediate]
        
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        thelist.append(item3)
        thelist.append(item4)
        thelist.append(item5)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        # Now some atoms will be deleted
        intermediate.delete_atom(5)
        intermediate.delete_atom(4)
        intermediate.delete_atom(6)
        intermediate.delete_atom(7)
        intermediate.delete_atom(1)
        intermediate.delete_atom(3)
        R1.delete_atom(carbon_one_index)
        R3.delete_atom(carbon_two_index)
        
        intermediate.change_residue('FR1')
        R1.change_residue('FR2')
        R2.change_residue('FR3')
        R3.change_residue('FR4')
        R4.change_residue('FR5')
        
        # Now merge the pdbs
        build = intermediate.merge_with_another_molecule(R1)
        build = build.merge_with_another_molecule(R2)
        build = build.merge_with_another_molecule(R3)
        build = build.merge_with_another_molecule(R4)
    
        return build
            
    def __sulfonyl_azide_thio_acid(self, pdb1, pdb2, reaction, cis_or_trans): # sulfonyl_azide - thio acid reactions
        """Simulates the reaction between a sulfonyl azide and a thio acid.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), either a sulfonyl azide or a thio acid.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        cis_or_trans -- As either cis or trans products are possible, this int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        # two choices are available! cis and trans around amide bond!
        if cis_or_trans == -1: cis_or_trans = random.randrange(0,2)
        
        # First, identify which is which
        if reaction[0] == "SULFONYL_AZIDE":
            sulfyl_azide = pdb1
            sulfyl_azide_S_index = reaction[1][0]
            sulfyl_azide_N_proximal_index = reaction[1][1]
            sulfyl_azide_N_middle_index = reaction[1][2]
            sulfyl_azide_N_distal_index = reaction[1][3]
            
            thio_acid = pdb2
            thio_acid_carbonyl_oxygen_index = reaction[3][0]
            thio_acid_carbonyl_carbon_index = reaction[3][1]
            thio_acid_carbonyl_S_index = reaction[3][2]
            thio_acid_carbonyl_H_index = reaction[3][3]
        
        else:
            
            sulfyl_azide = pdb2
            sulfyl_azide_S_index = reaction[3][0]
            sulfyl_azide_N_proximal_index = reaction[3][1]
            sulfyl_azide_N_middle_index = reaction[3][2]
            sulfyl_azide_N_distal_index = reaction[3][3]
            
            thio_acid = pdb1
            thio_acid_carbonyl_oxygen_index = reaction[1][0]
            thio_acid_carbonyl_carbon_index = reaction[1][1]
            thio_acid_carbonyl_S_index = reaction[1][2]
            thio_acid_carbonyl_H_index = reaction[1][3]
            
            
        # load intermediate
        intermediate = pymolecule.Molecule()
        if cis_or_trans == 0:
            intermediate.load_pdb("./intermediates/intermediate2.pdb")
        else:
            intermediate.load_pdb("./intermediates/intermediate2_isomer.pdb")
        
        # now move the sulfyl_azide
        tethers = [[3, sulfyl_azide_N_proximal_index], [4, sulfyl_azide_S_index]] 
        sulfyl_azide = intermediate.align_another_molecule_to_this_one(sulfyl_azide, tethers)
    
        # now move the thio acid in place
        tethers = [[1,thio_acid_carbonyl_carbon_index], [3,thio_acid_carbonyl_S_index], [2,thio_acid_carbonyl_oxygen_index]]
        thio_acid = intermediate.align_another_molecule_to_this_one(thio_acid, tethers)
    
        # Now delete some of the atoms
        sulfyl_azide.delete_atom(sulfyl_azide_N_middle_index)
        sulfyl_azide.delete_atom(sulfyl_azide_N_distal_index)
        thio_acid.delete_atom(thio_acid_carbonyl_oxygen_index)
        thio_acid.delete_atom(thio_acid_carbonyl_carbon_index)
        thio_acid.delete_atom(thio_acid_carbonyl_S_index)
        thio_acid.delete_atom(thio_acid_carbonyl_H_index)
        intermediate.delete_atom(4)
    
        # Now rotate the SO(2)-R part of the molecule to minimize steric hindrance. The R group coming of the amide cannot be rotated.
        item1 = [sulfyl_azide, sulfyl_azide.all_atoms[sulfyl_azide_N_proximal_index].coordinates, sulfyl_azide.all_atoms[sulfyl_azide_S_index].coordinates]
        item2 = [thio_acid]
        item3 = [intermediate]
        
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        thelist.append(item3)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        # Now delete some more of the atoms
        sulfyl_azide.delete_atom(sulfyl_azide_N_proximal_index)
    
        intermediate.change_residue('FR1')
        sulfyl_azide.change_residue('FR2')
        thio_acid.change_residue('FR3')
    
        # Now combine all the pdbs into one
        build = intermediate.merge_with_another_molecule(sulfyl_azide)
        build = build.merge_with_another_molecule(thio_acid)
        
        return build
            
    def __alkyne_azide(self, pdb1, pdb2, reaction, orientation_choice): # Azide alkyne Huisgen cycloaddition reactions
        """Simulates the reaction between an alkyne and an azide.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), either an alkyne or an azide.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        orientation_choice -- As two products are possible, this int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        # Pick between two configurations
        if orientation_choice == -1: orientation_choice = random.randrange(0,2)
        
        # Identify which is the alkyne and which is the azide
        if reaction[0] == "ALKYNE":
            alkyne_pdb = pdb1
            alkyne_first_neighbor_index = reaction[1][0]
            alkyne_first_C_index = reaction[1][1]
            alkyne_second_C_index = reaction[1][2]
            alkyne_second_neighbor_index = reaction[1][3]
            
            azide_pdb = pdb2
            azide_neighbor_index = reaction[3][0]
            azide_proximal_N_index = reaction[3][1]
            azide_middle_N_index = reaction[3][2]
            azide_distal_N_index = reaction[3][3]
        else:
            alkyne_pdb = pdb2
            alkyne_first_neighbor_index = reaction[3][0]
            alkyne_first_C_index = reaction[3][1]
            alkyne_second_C_index = reaction[3][2]
            alkyne_second_neighbor_index = reaction[3][3]
    
            azide_pdb = pdb1
            azide_neighbor_index = reaction[1][0]
            azide_proximal_N_index = reaction[1][1]
            azide_middle_N_index = reaction[1][2]
            azide_distal_N_index = reaction[1][3]
            
        # Now, get the alkyne fragments
        alkyne_frag1 = alkyne_pdb.get_branch(alkyne_first_C_index, alkyne_first_neighbor_index)
        alkyne_frag2 = alkyne_pdb.get_branch(alkyne_second_C_index, alkyne_second_neighbor_index)
            
        # Now, load in the intermediate
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb('./intermediates/intermediate1.pdb')
    
        # Now add the alkyne fragments to the intermediate. There are two possible orientations.
        if orientation_choice == 0:
            tethers = [[1, alkyne_first_C_index], [8, alkyne_first_neighbor_index]]
            alkyne_frag1 = intermediate.align_another_molecule_to_this_one(alkyne_frag1, tethers)
    
            tethers = [[2, alkyne_second_C_index], [7, alkyne_second_neighbor_index]]
            alkyne_frag2 = intermediate.align_another_molecule_to_this_one(alkyne_frag2, tethers)
    
        else: 
    
            tethers = [[2, alkyne_first_C_index], [7, alkyne_first_neighbor_index]]
            alkyne_frag1 = intermediate.align_another_molecule_to_this_one(alkyne_frag1, tethers)
    
            tethers = [[1, alkyne_second_C_index], [8, alkyne_second_neighbor_index]]
            alkyne_frag2 = intermediate.align_another_molecule_to_this_one(alkyne_frag2, tethers)
    
        # Now add the azide
        tethers = [[3, azide_proximal_N_index], [6, azide_neighbor_index]]
        azide_pdb = intermediate.align_another_molecule_to_this_one(azide_pdb, tethers)
        
        # Now delete a few extra atoms
        intermediate.delete_atom(6)
        intermediate.delete_atom(7)
        intermediate.delete_atom(8)
    
        azide_pdb.delete_atom(azide_middle_N_index)
        azide_pdb.delete_atom(azide_distal_N_index)
    
        # Now try to reduce steric hindrance
        item1 = [azide_pdb, azide_pdb.all_atoms[azide_proximal_N_index].coordinates, azide_pdb.all_atoms[azide_neighbor_index].coordinates]
        item2 = [alkyne_frag1, alkyne_frag1.all_atoms[alkyne_first_C_index].coordinates, alkyne_frag1.all_atoms[alkyne_first_neighbor_index].coordinates]
        item3 = [alkyne_frag2, alkyne_frag2.all_atoms[alkyne_second_C_index].coordinates, alkyne_frag2.all_atoms[alkyne_second_neighbor_index].coordinates]
        item4 = [intermediate]
        
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        thelist.append(item3)
        thelist.append(item4)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        # Delete a few more atoms
    
        azide_pdb.delete_atom(azide_proximal_N_index)
        alkyne_frag1.delete_atom(alkyne_first_C_index)
        alkyne_frag2.delete_atom(alkyne_second_C_index)
    
        intermediate.change_residue('FR1')
        alkyne_frag1.change_residue('FR2')
        alkyne_frag2.change_residue('FR3')
        azide_pdb.change_residue('FR4')
    
        # Now merge all pdbs into one
        
        build = intermediate.merge_with_another_molecule(alkyne_frag1)
        build = build.merge_with_another_molecule(alkyne_frag2)
        build = build.merge_with_another_molecule(azide_pdb)
        
        return build    
    
    def __primary_amine_to_azide_or_NCO_or_NCS(self, pdb1, reaction, choice):
        """Simulates the conversion of a primary amine to an azide, NCO, or NCS.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), containing a primary amine.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        choice -- This int dictates which product should be made (azide, NCO, or NCS).
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        if choice == -1 : choice = random.randrange(0,3)
        
        # rename the atoms of the amine to make it easier
        amine_X_index = reaction[1][1]
        amine_N_index = reaction[1][0]
    
        azide_X_index = 1
        azide_N_proximal_index = 2
        
        # Load in the azide
        azide = pymolecule.Molecule()
        if choice == 0:
            azide.load_pdb('./intermediates/intermediate9.pdb')
        elif choice == 1:
            azide.load_pdb('./intermediates/N-C-O.pdb')
        elif choice == 2:
            azide.load_pdb('./intermediates/N-C-S.pdb')
        
        # Position the azide
        tethers = [[amine_X_index, azide_X_index], [amine_N_index, azide_N_proximal_index]]
        azide = pdb1.align_another_molecule_to_this_one(azide, tethers)
    
        # Now delete hydrogen atoms on amine
        for t in range(2,len(reaction[1])):
            pdb1.delete_atom(reaction[1][t])
        
        # delete nitrogen on amine
        pdb1.delete_atom(amine_N_index)
        
        # Now rotate the azide
        item1 = [azide, azide.all_atoms[azide_X_index].coordinates, azide.all_atoms[azide_N_proximal_index].coordinates]
        item2 = [pdb1]
        
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        # delete X on azide
        azide.delete_atom(azide_X_index)
    
        pdb1.change_residue('FR1')
        azide.change_residue('FR2')
        
        # merge the pdbs
        build = pdb1.merge_with_another_molecule(azide)
        
        return build
    
    def __azide_to_primary_amine(self, pdb1, reaction, charge_choice):
        """Simulates the conversion of an azide to a primary amine.
        
        Arguments:
        pdb1 -- The molecular model (pymolecule.Molecule), containing an azide.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        charge_choice -- The primary amine can be charged or neutral. This int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        if charge_choice == -1 : charge_choice = random.randrange(0,2)
        
        # get easy names of azide
        azide_X_index = reaction[1][0]
        azide_N_proximal_index = reaction[1][1]
        azide_N_middle_index = reaction[1][2]
        azide_N_distal_index = reaction[1][3]
        if len(reaction[1]) == 5:
            azide_X_neighbor_index = reaction[1][4]
        else:
            azide_X_neighbor_index=-999
        
        # load in a primary amine model. Both NH(2) and NH(3) (neither of which is chiral, even in computational terms)
        amine = pymolecule.Molecule()
        sp2 = 'NO'
        if charge_choice == 0:
            # Now you need to determine whether this amine will be bound to an sp3 or sp2 hybridized carbon
            if pdb1.hybridization(azide_X_index) == 3:
                amine.load_pdb('./intermediates/intermediate10.pdb')
            else: # so it must be bound to an sp2 hybridized atom
                amine.load_pdb('./intermediates/sp2_primary_amine.pdb')
                sp2 = 'YES'
        else:
            amine.load_pdb('./intermediates/charged_amine.pdb')
            # I wonder if it is possible to have a charged amine next to a sp2 hybridized carbon. You might want to look into that...
            
        # get easy names of amine
        amine_X_index = 1
        amine_N_index = 2
        
        # Position the amino group
        tethers = [[azide_X_index, amine_X_index], [azide_N_proximal_index, amine_N_index]]
        
        # If it's an sp2 hybridized amine, add an extra tether
        if sp2 == 'YES' and azide_X_neighbor_index != -999: tethers.append([azide_X_neighbor_index,5])
        
        amine = pdb1.align_another_molecule_to_this_one(amine, tethers)
    
        # delete azide
        pdb1.delete_atom(azide_N_proximal_index)
        pdb1.delete_atom(azide_N_middle_index)
        pdb1.delete_atom(azide_N_distal_index)
    
        #delete extra C from amine
        if sp2 == 'YES': amine.delete_atom(5)
        
        if sp2 == 'NO':
            # rotate amine to minimize steric hindrance
            item1 = [amine, amine.all_atoms[amine_X_index].coordinates, amine.all_atoms[amine_N_index].coordinates]
            item2 = [pdb1]
        
            thelist = []
            thelist.append(item1)
            thelist.append(item2)
            self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        amine.delete_atom(amine_X_index)
    
        pdb1.change_residue('FR1')
        amine.change_residue('FR2')
        
        # Now merge the Molecules
        build = pdb1.merge_with_another_molecule(amine)
    
        return build
        
    
    def __carbonochloridate_amine(self, pdb1, pdb2, reaction, cis_or_trans_choice): # Note that pdb2 can be a primary or secondary amine
        """Simulates the reaction between a carbonochloridate and an amine.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), containing either a carbonochloridate or an amine.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        cis_or_trans_choice -- Either cis or trans products can be generated. This int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        # pdb1 contains the first pdb data structure
        # pdb2 contains the second pdb data structure
        # reaction is a list of lists. It looks like this:
        #   ["NAME OF REACTIVE GROUP ON PDB1", INDICES OF THE ATOMS OF THAT PDB1 GROUP, "NAME OF REACTIVE GROUP ON PDB2", INDICES OF THE ATOMS OF THAT PDB2 GROUP, ]
        #   Here's a sample:
        #       ['CARBONOCHLORIDATE', [2, 1, 5, 6], 'AMINE', [3, 4, 11, 12]]
        # The order of the indices of the atoms of each group are given by the function defined in self.structure_locate_groups.py
    
        if cis_or_trans_choice == -1 : cis_or_trans_choice = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
        
        # You need to figure out which pdb is the carbonochloridate and which is the primary amine
        if reaction[0] == "CARBONOCHLORIDATE": # pdb1 is the carbonochloridate
            carbonochloridate_pdb = pdb1
            carbonochloridate_indecies = reaction[1]
        
            amine_pdb = pdb2
            amine_indecies = reaction[3]
            
        else: # pdb2 is the carbonochloridate
            carbonochloridate_pdb = pdb2
            carbonochloridate_indecies = reaction[3]
    
            amine_pdb = pdb1
            amine_indecies = reaction[1]
             
        # get easy names of the carbonochloridate, so we don't have to refer to the reaction[] list
        carbonochloridate_ester_oxygen = carbonochloridate_indecies[0]
        carbonochloridate_ester_carybonyl_carbon = carbonochloridate_indecies[1]
        carbonochloridate_ester_carybonyl_oxygen = carbonochloridate_indecies[2]
        carbonochloridate_ester_chloride = carbonochloridate_indecies[3]
    
        # get easy names for the amine atoms
        amine_nitrogen = amine_indecies[0]
        amine_Rs = []
        amine_hydrogens = []
        for t in range(1,len(amine_indecies)): # starting at index 1, so skipping amine_nitrogen
            index = amine_indecies[t]
            atom = amine_pdb.all_atoms[index]
            if atom.element == "H":
                amine_hydrogens.append(index)
            else :
                amine_Rs.append(index)
        
        # If it's a primary amine, they'll be only one element in amine_Rs list. Move one of the hydrogens to that list to complete it
        while len(amine_Rs) == 1:
            amine_Rs.append(amine_hydrogens[0])
            amine_hydrogens.pop(0)
            
        # I guess this would technically accept an amonium ion. Perhaps that should be added to the self.structure_locate_groups.py file...?
            
        # Now, load in the intermediate. This is as a model carbamate, the product of this reaction.
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb('./intermediates/carbamate.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_esther_oxygen_index = 3
        intermediate_carbonyl_carbon_index = 1
        intermediate_carbonyl_oxygen_index = 2
        intermediate_nitrogen_index = 4
        intermediate_amine_R1_index = 5
        intermediate_amine_R2_index = 6
        
        # Now, posiition the carbonochloridate onto the model carbamate
        tethers = [[intermediate_carbonyl_carbon_index, carbonochloridate_ester_carybonyl_carbon], [intermediate_carbonyl_oxygen_index, carbonochloridate_ester_carybonyl_oxygen], [intermediate_nitrogen_index, carbonochloridate_ester_chloride]]
        carbonochloridate_pdb = intermediate.align_another_molecule_to_this_one(carbonochloridate_pdb, tethers)
        # So what does thae above do? First, it looks at the first pair in the tethers list: [intermediate_carbonyl_carbon_index, carbonochloridate_ester_carybonyl_carbon]
        # It moves pdb2 onto pdb1 such that the atoms corresponding to these two indices are right on top of each other. This point is a "pivot point."
        # Then it rotates pdb2 about that pivot point so as to minimize the distance between the atoms of the subsequent pairs.
    
        # Now the nitrogen on the amine is going to be going from SP3 to SP2. We need to disect this amine by removing it's branches
        # cis_or_trans_choice was either specified by the user or randomly set to either 0 or 1.
        # Use this varialbe in deciding which branch is which.
        if cis_or_trans_choice == 0 :
            amine_branch_R1_index = amine_Rs[0]
            amine_branch_R2_index = amine_Rs[1]
            
        else :
            amine_branch_R1_index = amine_Rs[1]
            amine_branch_R2_index = amine_Rs[0]
        
        # Is the nitrogen an SP2 hybridized nitrogen in a ring? If so, it can be processed, even though ring nitrogens generally can't
        sp2_nitrogen_in_ring = 'FALSE'
        if len(reaction) == 5:
            if reaction[4] == 'SECONDARY_AMINE_SP2_IN_RING':
                sp2_nitrogen_in_ring = 'TRUE'
    
        
        # First, consider if it is not in a ring
        if sp2_nitrogen_in_ring == 'FALSE':
            # Make sure the amide is not in a loop. These must be discarded, unless they are SP2 hybridized
            if amine_pdb.in_same_ring(amine_nitrogen, amine_branch_R2_index) == "TRUE":
                return pymolecule.Molecule() # return a blank Molecule
                
            amine_branch_1 = amine_pdb.get_branch(amine_nitrogen, amine_branch_R1_index)
            amine_branch_2 = amine_pdb.get_branch(amine_nitrogen, amine_branch_R2_index)
            
            # Now we need to position the amine branches onto the model carbamate
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R1_index, amine_branch_R1_index]]
            amine_branch_1 = intermediate.align_another_molecule_to_this_one(amine_branch_1, tethers)
        
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R2_index, amine_branch_R2_index]]
            amine_branch_2 = intermediate.align_another_molecule_to_this_one(amine_branch_2, tethers)
        
            # Now delete a few extra atoms. We won't even need to save the intermediate
            carbonochloridate_pdb.delete_atom(carbonochloridate_ester_chloride)
        
            # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
            # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
            # Note that the carbonochloridate is fixed and cannot rotate. The two branches of the amine can, though.
            item1 = [amine_branch_1, amine_branch_1.all_atoms[amine_nitrogen].coordinates, amine_branch_1.all_atoms[amine_branch_R1_index].coordinates]
            item2 = [amine_branch_2, amine_branch_2.all_atoms[amine_nitrogen].coordinates, amine_branch_2.all_atoms[amine_branch_R2_index].coordinates]
            item3 = [carbonochloridate_pdb] # These won't rotate, so the coordinates are omitted
            item4 = [intermediate]
            
            thelist = []
            thelist.append(item1)
            thelist.append(item2)
            thelist.append(item3)
            thelist.append(item4)
            self.structure_pdb_functions.reduce_steric_hindrance(thelist)
            
            # Now delete a few more atoms.
            amine_branch_1.delete_atom(amine_nitrogen)
    
            carbonochloridate_pdb.change_residue('FR1')
            amine_branch_1.change_residue('FR2')
            amine_branch_2.change_residue('FR3')
    
            # Now merge all pdbs into one
            build = carbonochloridate_pdb.merge_with_another_molecule(amine_branch_1)
            build = build.merge_with_another_molecule(amine_branch_2)
    
        elif sp2_nitrogen_in_ring == 'TRUE': # Now deal with the special case where the amine is sp2 in a ring
            
            # Now we need to position the amine onto the model carbamate
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R1_index, amine_branch_R1_index], [intermediate_amine_R2_index, amine_branch_R2_index]]
            amine_pdb = intermediate.align_another_molecule_to_this_one(amine_pdb, tethers)
        
            # Now delete a few extra atoms. We won't even need to save the intermediate
            carbonochloridate_pdb.delete_atom(carbonochloridate_ester_chloride)
            amine_pdb.delete_atom(amine_hydrogens[0])
    
            carbonochloridate_pdb.change_residue('FR1')
            amine_pdb.change_residue('FR2')
    
            # Now merge all pdbs into one
            build = carbonochloridate_pdb.merge_with_another_molecule(amine_pdb)
        
        return build    
        
    def __carboxylate_amine(self, pdb1, pdb2, reaction, cis_or_trans_choice): # Note that pdb2 can be a primary or secondary amine
        """Simulates the reaction between a carboxylate and an amine.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), either an carboxylate or an amine.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        cis_or_trans_choice -- Both cis and trans products can be formed. This int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        # pdb1 contains the first pdb data structure
        # pdb2 contains the second pdb data structure
        # reaction is a list of lists. It looks like this:
        #   ["NAME OF REACTIVE GROUP ON PDB1", INDICES OF THE ATOMS OF THAT PDB1 GROUP, "NAME OF REACTIVE GROUP ON PDB2", INDICES OF THE ATOMS OF THAT PDB2 GROUP, ]
        #   Here's a sample:
        #       ['CARBOXYLATE', [2, 1, 5, 6], 'AMINE', [3, 4, 11, 12]]
        # The order of the indices of the atoms of each group are given by the function defined in self.structure_locate_groups.py
    
        if cis_or_trans_choice == -1 : cis_or_trans_choice = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
        
        # You need to figure out which pdb is the carbonochloridate and which is the primary amine
        if reaction[0] == "CARBOXYLATE": # pdb1 is the carbonochloridate
            carboxylate_pdb = pdb1
            carboxylate_indecies = reaction[1]
        
            amine_pdb = pdb2
            amine_indecies = reaction[3]
            
        else: # pdb2 is the carbonochloridate
            carboxylate_pdb = pdb2
            carboxylate_indecies = reaction[3]
    
            amine_pdb = pdb1
            amine_indecies = reaction[1]
             
        # get easy names of the carboxylate, so we don't have to refer to the reaction[] list
        carboxylate_carybonyl_carbon = carboxylate_indecies[0]
        carboxylate_carybonyl_oxygen = carboxylate_indecies[1]
        carboxylate_hydroxyl_oxygen = carboxylate_indecies[2]
        carboxylate_hydroxyl_hydrogen = carboxylate_indecies[3]
    
        # get easy names for the amine atoms
        amine_nitrogen = amine_indecies[0]
        amine_Rs = []
        amine_hydrogens = []
        for t in range(1,len(amine_indecies)): # starting at index 1, so skipping amine_nitrogen
            index = amine_indecies[t]
            atom = amine_pdb.all_atoms[index]
            if atom.element == "H":
                amine_hydrogens.append(index)
            else :
                amine_Rs.append(index)
        
        # If it's a primary amine, they'll be only one element in amine_Rs list. Move one of the hydrogens to that list to complete it
        while len(amine_Rs) == 1:
            amine_Rs.append(amine_hydrogens[0])
            amine_hydrogens.pop(0)
            
        # I guess this would technically accept an amonium ion. Perhaps that should be added to the self.structure_locate_groups.py file...?
            
        # Now, load in the intermediate. This is as a model amide, the product of this reaction.
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb('./intermediates/amide.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_oppose_carbon_index = 3
        intermediate_carbonyl_carbon_index = 1
        intermediate_carbonyl_oxygen_index = 2
        intermediate_nitrogen_index = 4
        intermediate_amine_R1_index = 5
        intermediate_amine_R2_index = 6
        
        # Now, position the carboxylate onto the model amide
        tethers = [[intermediate_carbonyl_carbon_index, carboxylate_carybonyl_carbon], [intermediate_nitrogen_index, carboxylate_hydroxyl_oxygen], [intermediate_carbonyl_oxygen_index, carboxylate_carybonyl_oxygen]]
        carboxylate_pdb = intermediate.align_another_molecule_to_this_one(carboxylate_pdb, tethers)
    
        # So what does thae above do? First, it looks at the first pair in the tethers list: [intermediate_carbonyl_carbon_index, carboxylate_ester_carybonyl_carbon]
        # It moves pdb2 onto pdb1 such that the atoms corresponding to these two indices are right on top of each other. This point is a "pivot point."
        # Then it rotates pdb2 about that pivot point so as to minimize the distance between the atoms of the subsequent pairs.
    
        # Now the nitrogen on the amine is going to be going from SP3 to SP2. We need to disect this amine by removing it's branches
        # cis_or_trans_choice was either specified by the user or randomly set to either 0 or 1.
        # Use this varialbe in deciding which branch is which.
        if cis_or_trans_choice == 0 :
            amine_branch_R1_index = amine_Rs[0]
            amine_branch_R2_index = amine_Rs[1]
            
        else :
            amine_branch_R1_index = amine_Rs[1]
            amine_branch_R2_index = amine_Rs[0]
    
        # Is the nitrogen an SP2 hybridized nitrogen in a ring? If so, it can be processed, even though ring nitrogens generally can't
        sp2_nitrogen_in_ring = 'FALSE'
        if len(reaction) == 5:
            if reaction[4] == 'SECONDARY_AMINE_SP2_IN_RING':
                sp2_nitrogen_in_ring = 'TRUE'
                
        if sp2_nitrogen_in_ring == 'FALSE':
            
            # Make sure the amide is not in a loop. These must be discarded
            if amine_pdb.in_same_ring(amine_nitrogen, amine_branch_R2_index) == "TRUE":
                return pymolecule.Molecule() # return a blank Molecule
        
            amine_branch_1 = amine_pdb.get_branch(amine_nitrogen, amine_branch_R1_index)
            amine_branch_2 = amine_pdb.get_branch(amine_nitrogen, amine_branch_R2_index)
            
            # Now we need to position the amine branches onto the model amide
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R1_index, amine_branch_R1_index]]
            amine_branch_1 = intermediate.align_another_molecule_to_this_one(amine_branch_1, tethers)
        
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R2_index, amine_branch_R2_index]]
            amine_branch_2 = intermediate.align_another_molecule_to_this_one(amine_branch_2, tethers)
        
            # Now delete a few more atoms.
            carboxylate_pdb.delete_atom(carboxylate_hydroxyl_oxygen)
            if carboxylate_hydroxyl_hydrogen!=-1: carboxylate_pdb.delete_atom(carboxylate_hydroxyl_hydrogen)
        
            # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
            # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
            # Note that the carboxylate is fixed and cannot rotate. The two branches of the amine can, though.
            item1 = [amine_branch_1, amine_branch_1.all_atoms[amine_nitrogen].coordinates, amine_branch_1.all_atoms[amine_branch_R1_index].coordinates]
            item2 = [amine_branch_2, amine_branch_2.all_atoms[amine_nitrogen].coordinates, amine_branch_2.all_atoms[amine_branch_R2_index].coordinates]
            item3 = [carboxylate_pdb] # These won't rotate, so the coordinates are omitted
            item4 = [intermediate]
            
            thelist = []
            thelist.append(item1)
            thelist.append(item2)
            thelist.append(item3)
            thelist.append(item4)
            self.structure_pdb_functions.reduce_steric_hindrance(thelist)
        
            # Now delete a few more atoms.
            amine_branch_1.delete_atom(amine_nitrogen)
            
            carboxylate_pdb.change_residue('FR1')
            amine_branch_1.change_residue('FR2')
            amine_branch_2.change_residue('FR3')
            
            # Now merge all pdbs into one
            build = carboxylate_pdb.merge_with_another_molecule(amine_branch_1)
            build = build.merge_with_another_molecule(amine_branch_2)
        elif sp2_nitrogen_in_ring == 'TRUE':
        
            # Now we need to position the amine branches onto the model amide
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R1_index, amine_branch_R1_index], [intermediate_amine_R2_index, amine_branch_R2_index]]
            amine_pdb = intermediate.align_another_molecule_to_this_one(amine_pdb, tethers)
        
            # Now delete a few more atoms.
            carboxylate_pdb.delete_atom(carboxylate_hydroxyl_oxygen)
            if carboxylate_hydroxyl_hydrogen!=-1: carboxylate_pdb.delete_atom(carboxylate_hydroxyl_hydrogen)
        
            # Now delete a few more atoms.
            amine_pdb.delete_atom(amine_hydrogens[0])
            
            carboxylate_pdb.change_residue('FR1')
            amine_pdb.change_residue('FR2')
            
            # Now merge all pdbs into one
            build = carboxylate_pdb.merge_with_another_molecule(amine_pdb)
            
        return build    
        
    def __acylhalide_amine(self, pdb1, pdb2, reaction, cis_or_trans_choice): # Note that pdb2 can be a primary or secondary amine
        """Simulates the reaction between an acyl halide and an amine.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), containing either an acyl halide or an amine.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        cis_or_trans_choice -- Both cis and trans products can be formed. This int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        # pdb1 contains the first pdb data structure
        # pdb2 contains the second pdb data structure
        # reaction is a list of lists. It looks like this:
        #   ["NAME OF REACTIVE GROUP ON PDB1", INDICES OF THE ATOMS OF THAT PDB1 GROUP, "NAME OF REACTIVE GROUP ON PDB2", INDICES OF THE ATOMS OF THAT PDB2 GROUP, ]
        #   Here's a sample:
        #       ['acylhalide', [2, 1, 5, 6], 'AMINE', [3, 4, 11, 12]]
        # The order of the indices of the atoms of each group are given by the function defined in self.structure_locate_groups.py
    
        if cis_or_trans_choice == -1 : cis_or_trans_choice = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
        
        # You need to figure out which pdb is the carbonochloridate and which is the primary amine
        if reaction[0] == "ACYLHALIDE": # pdb1 is the carbonochloridate
            acylhalide_pdb = pdb1
            acylhalide_indecies = reaction[1]
        
            amine_pdb = pdb2
            amine_indecies = reaction[3]
            
        else: # pdb2 is the carbonochloridate
            acylhalide_pdb = pdb2
            acylhalide_indecies = reaction[3]
    
            amine_pdb = pdb1
            amine_indecies = reaction[1]
             
        # get easy names of the acylhalide, so we don't have to refer to the reaction[] list
        acylhalide_carybonyl_carbon = acylhalide_indecies[0]
        acylhalide_carybonyl_oxygen = acylhalide_indecies[1]
        acylhalide_chloride = acylhalide_indecies[2]
    
        # get easy names for the amine atoms
        amine_nitrogen = amine_indecies[0]
        amine_Rs = []
        amine_hydrogens = []
        for t in range(1,len(amine_indecies)): # starting at index 1, so skipping amine_nitrogen
            index = amine_indecies[t]
            atom = amine_pdb.all_atoms[index]
            if atom.element == "H":
                amine_hydrogens.append(index)
            else :
                amine_Rs.append(index)
        
        # If it's a primary amine, they'll be only one element in amine_Rs list. Move one of the hydrogens to that list to complete it
        while len(amine_Rs) == 1:
            amine_Rs.append(amine_hydrogens[0])
            amine_hydrogens.pop(0)
            
        # I guess this would technically accept an amonium ion. Perhaps that should be added to the self.structure_locate_groups.py file...?
            
        # Now, load in the intermediate. This is as a model amide, the product of this reaction.
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb('./intermediates/amide.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_oppose_carbon_index = 3
        intermediate_carbonyl_carbon_index = 1
        intermediate_carbonyl_oxygen_index = 2
        intermediate_nitrogen_index = 4
        intermediate_amine_R1_index = 5
        intermediate_amine_R2_index = 6
        
        # Now, posiition the acylhalide onto the model amide
        tethers = [[intermediate_carbonyl_carbon_index, acylhalide_carybonyl_carbon], [intermediate_carbonyl_oxygen_index, acylhalide_carybonyl_oxygen], [intermediate_nitrogen_index, acylhalide_chloride]]
        acylhalide_pdb = intermediate.align_another_molecule_to_this_one(acylhalide_pdb, tethers)
        # So what does thae above do? First, it looks at the first pair in the tethers list: [intermediate_carbonyl_carbon_index, acylhalide_ester_carybonyl_carbon]
        # It moves pdb2 onto pdb1 such that the atoms corresponding to these two indices are right on top of each other. This point is a "pivot point."
        # Then it rotates pdb2 about that pivot point so as to minimize the distance between the atoms of the subsequent pairs.
    
        # Now the nitrogen on the amine is going to be going from SP3 to SP2. We need to disect this amine by removing it's branches
        # cis_or_trans_choice was either specified by the user or randomly set to either 0 or 1.
        # Use this varialbe in deciding which branch is which.
        if cis_or_trans_choice == 0 :
            amine_branch_R1_index = amine_Rs[0]
            amine_branch_R2_index = amine_Rs[1]
            
        else :
            amine_branch_R1_index = amine_Rs[1]
            amine_branch_R2_index = amine_Rs[0]
    
        # Is the nitrogen an SP2 hybridized nitrogen in a ring? If so, it can be processed, even though ring nitrogens generally can't
        sp2_nitrogen_in_ring = 'FALSE'
        if len(reaction) == 5:
            if reaction[4] == 'SECONDARY_AMINE_SP2_IN_RING':
                sp2_nitrogen_in_ring = 'TRUE'
                
        if sp2_nitrogen_in_ring == 'FALSE':
            # Make sure the amide is not in a loop. These must be discarded
            if amine_pdb.in_same_ring(amine_nitrogen, amine_branch_R2_index) == "TRUE":
                return pymolecule.Molecule() # return a blank Molecule
        
            amine_branch_1 = amine_pdb.get_branch(amine_nitrogen, amine_branch_R1_index)
            amine_branch_2 = amine_pdb.get_branch(amine_nitrogen, amine_branch_R2_index)
            
            # Now we need to position the amine branches onto the model amide
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R1_index, amine_branch_R1_index]]
            amine_branch_1 = intermediate.align_another_molecule_to_this_one(amine_branch_1, tethers)
        
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R2_index, amine_branch_R2_index]]
            amine_branch_2 = intermediate.align_another_molecule_to_this_one(amine_branch_2, tethers)
           
            # Now delete a few more atoms.
            acylhalide_pdb.delete_atom(acylhalide_chloride)
        
            # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
            # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
            # Note that the acylhalide is fixed and cannot rotate. The two branches of the amine can, though.
            item1 = [amine_branch_1, amine_branch_1.all_atoms[amine_nitrogen].coordinates, amine_branch_1.all_atoms[amine_branch_R1_index].coordinates]
            item2 = [amine_branch_2, amine_branch_2.all_atoms[amine_nitrogen].coordinates, amine_branch_2.all_atoms[amine_branch_R2_index].coordinates]
            item3 = [acylhalide_pdb] # These won't rotate, so the coordinates are omitted
            item4 = [intermediate]
            
            thelist = []
            thelist.append(item1)
            thelist.append(item2)
            thelist.append(item3)
            thelist.append(item4)
            self.structure_pdb_functions.reduce_steric_hindrance(thelist)
        
            # Now delete a few more atoms.
            amine_branch_1.delete_atom(amine_nitrogen)
            
            acylhalide_pdb.change_residue('FR1')
            amine_branch_1.change_residue('FR2')
            amine_branch_2.change_residue('FR3')
        
            # Now merge all pdbs into one
            build = acylhalide_pdb.merge_with_another_molecule(amine_branch_1)
            build = build.merge_with_another_molecule(amine_branch_2)
    
        elif sp2_nitrogen_in_ring == 'TRUE':
    
            # Now we need to position the amine branches onto the model amide
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R1_index, amine_branch_R1_index], [intermediate_amine_R2_index, amine_branch_R2_index]]
            amine_pdb = intermediate.align_another_molecule_to_this_one(amine_pdb, tethers)
        
            # Now delete a few more atoms.
            acylhalide_pdb.delete_atom(acylhalide_chloride)
            
            # Now delete a few more atoms.
            amine_pdb.delete_atom(amine_hydrogens[0])
        
            acylhalide_pdb.change_residue('FR1')
            amine_pdb.change_residue('FR2')
        
            # Now merge all pdbs into one
            build = acylhalide_pdb.merge_with_another_molecule(amine_pdb)
        
        return build    
        
    def __ester_amine(self, pdb1, pdb2, reaction, cis_or_trans_choice): # Note that pdb2 can be a primary or secondary amine
        """Simulates the reaction between an ester and an amine.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), containing either an ester or an amine.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        cis_or_trans_choice -- Both cis and trans products can be formed. This int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        # pdb1 contains the first pdb data structure
        # pdb2 contains the second pdb data structure
        # reaction is a list of lists. It looks like this:
        #   ["NAME OF REACTIVE GROUP ON PDB1", INDICES OF THE ATOMS OF THAT PDB1 GROUP, "NAME OF REACTIVE GROUP ON PDB2", INDICES OF THE ATOMS OF THAT PDB2 GROUP, ]
        #   Here's a sample:
        #       ['ESTER', [2, 1, 5, 6], 'AMINE', [3, 4, 11, 12]]
        # The order of the indices of the atoms of each group are given by the function defined in self.structure_locate_groups.py
    
        if cis_or_trans_choice == -1 : cis_or_trans_choice = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
        
        # You need to figure out which pdb is the carbonochloridate and which is the primary amine
        if reaction[0] == "ESTER": # pdb1 is the carbonochloridate
            ester_pdb = pdb1
            ester_indecies = reaction[1]
        
            amine_pdb = pdb2
            amine_indecies = reaction[3]
            
        else: # pdb2 is the carbonochloridate
            ester_pdb = pdb2
            ester_indecies = reaction[3]
    
            amine_pdb = pdb1
            amine_indecies = reaction[1]
             
        # get easy names of the ester, so we don't have to refer to the reaction[] list
        ester_carybonyl_carbon = ester_indecies[0]
        ester_carybonyl_oxygen = ester_indecies[1]
        ester_oxygen = ester_indecies[2]
        ester_carbon_beyond_oxygen = ester_indecies[3]
    
        # get easy names for the amine atoms
        amine_nitrogen = amine_indecies[0]
        amine_Rs = []
        amine_hydrogens = []
        for t in range(1,len(amine_indecies)): # starting at index 1, so skipping amine_nitrogen
            index = amine_indecies[t]
            atom = amine_pdb.all_atoms[index]
            if atom.element == "H":
                amine_hydrogens.append(index)
            else :
                amine_Rs.append(index)
        
        # If it's a primary amine, they'll be only one element in amine_Rs list. Move one of the hydrogens to that list to complete it
        while len(amine_Rs) == 1:
            amine_Rs.append(amine_hydrogens[0])
            amine_hydrogens.pop(0)
            
        # I guess this would technically accept an amonium ion. Perhaps that should be added to the self.structure_locate_groups.py file...?
            
        # Now, load in the intermediate. This is as a model amide, the product of this reaction.
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb('./intermediates/amide.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_oppose_carbon_index = 3
        intermediate_carbonyl_carbon_index = 1
        intermediate_carbonyl_oxygen_index = 2
        intermediate_nitrogen_index = 4
        intermediate_amine_R1_index = 5
        intermediate_amine_R2_index = 6
        
        # Now, we need to discard the alcohol portion of the ester
        ester_pdb = ester_pdb.get_branch(ester_oxygen, ester_carybonyl_carbon)
        
        # Now, posiition the ester onto the model amide
        tethers = [[intermediate_carbonyl_carbon_index, ester_carybonyl_carbon], [intermediate_carbonyl_oxygen_index, ester_carybonyl_oxygen], [intermediate_nitrogen_index, ester_oxygen]]
        ester_pdb = intermediate.align_another_molecule_to_this_one(ester_pdb, tethers)
        # So what does thae above do? First, it looks at the first pair in the tethers list: [intermediate_carbonyl_carbon_index, ester_ester_carybonyl_carbon]
        # It moves pdb2 onto pdb1 such that the atoms corresponding to these two indices are right on top of each other. This point is a "pivot point."
        # Then it rotates pdb2 about that pivot point so as to minimize the distance between the atoms of the subsequent pairs.
    
        # Now the nitrogen on the amine is going to be going from SP3 to SP2. We need to disect this amine by removing it's branches
        # cis_or_trans_choice was either specified by the user or randomly set to either 0 or 1.
        # Use this varialbe in deciding which branch is which.
        if cis_or_trans_choice == 0 :
            amine_branch_R1_index = amine_Rs[0]
            amine_branch_R2_index = amine_Rs[1]
            
        else :
            amine_branch_R1_index = amine_Rs[1]
            amine_branch_R2_index = amine_Rs[0]
    
        # Is the nitrogen an SP2 hybridized nitrogen in a ring? If so, it can be processed, even though ring nitrogens generally can't
        sp2_nitrogen_in_ring = 'FALSE'
        if len(reaction) == 5:
            if reaction[4] == 'SECONDARY_AMINE_SP2_IN_RING':
                sp2_nitrogen_in_ring = 'TRUE'
                
        if sp2_nitrogen_in_ring == 'FALSE':
            # Make sure the amide is not in a loop. These must be discarded
            if amine_pdb.in_same_ring(amine_nitrogen, amine_branch_R2_index) == "TRUE":
                return pymolecule.Molecule() # return a blank Molecule
                
            amine_branch_1 = amine_pdb.get_branch(amine_nitrogen, amine_branch_R1_index)
            amine_branch_2 = amine_pdb.get_branch(amine_nitrogen, amine_branch_R2_index)
            
            # Now we need to position the amine branches onto the model amide
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R1_index, amine_branch_R1_index]]
            amine_branch_1 = intermediate.align_another_molecule_to_this_one(amine_branch_1, tethers)
        
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R2_index, amine_branch_R2_index]]
            amine_branch_2 = intermediate.align_another_molecule_to_this_one(amine_branch_2, tethers)
        
            # Now delete a few more atoms.
            ester_pdb.delete_atom(ester_oxygen)
        
            another_ester_rotate_indexes = ester_pdb.all_atoms[ester_carybonyl_carbon].indecies_of_atoms_connecting
            another_ester_rotate_indexes.remove(ester_carybonyl_oxygen)
            another_ester_rotate_index = another_ester_rotate_indexes[0]
        
            ester_pdb.delete_atom(ester_carybonyl_oxygen)
            
            # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
            # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
            # Note that the ester is fixed and cannot rotate. The two branches of the amine can, though.
            item1 = [amine_branch_1, amine_branch_1.all_atoms[amine_nitrogen].coordinates, amine_branch_1.all_atoms[amine_branch_R1_index].coordinates]
            item2 = [amine_branch_2, amine_branch_2.all_atoms[amine_nitrogen].coordinates, amine_branch_2.all_atoms[amine_branch_R2_index].coordinates]
            item3 = [ester_pdb, ester_pdb.all_atoms[ester_carybonyl_carbon].coordinates, ester_pdb.all_atoms[another_ester_rotate_index].coordinates]
            item4 = [intermediate] # These won't rotate, so the coordinates are omitted
            
            thelist = []
            thelist.append(item1)
            thelist.append(item2)
            thelist.append(item3)
            thelist.append(item4)
            self.structure_pdb_functions.reduce_steric_hindrance(thelist)
            
            # Now delete a few more atoms.
            amine_branch_1.delete_atom(amine_nitrogen)
            ester_pdb.delete_atom(ester_carybonyl_carbon)
            intermediate.delete_atom(intermediate_oppose_carbon_index)
            intermediate.delete_atom(intermediate_amine_R1_index)
            intermediate.delete_atom(intermediate_amine_R2_index)
            intermediate.delete_atom(intermediate_nitrogen_index)
    
            intermediate.change_residue('FR1')
            ester_pdb.change_residue('FR2')
            amine_branch_1.change_residue('FR3')
            amine_branch_2.change_residue('FR4')
        
            # Now merge all pdbs into one
            build = intermediate.merge_with_another_molecule(ester_pdb)
            build = build.merge_with_another_molecule(amine_branch_1)
            build = build.merge_with_another_molecule(amine_branch_2)
    
        elif sp2_nitrogen_in_ring == 'TRUE':
            
            # Now we need to position the amine branches onto the model amide
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R1_index, amine_branch_R1_index], [intermediate_amine_R2_index, amine_branch_R2_index]]
            amine_pdb = intermediate.align_another_molecule_to_this_one(amine_pdb, tethers)
        
            # Now delete a few more atoms.
            ester_pdb.delete_atom(ester_oxygen)
        
            # Now delete a few more atoms.
            amine_pdb.delete_atom(amine_hydrogens[0])
        
            ester_pdb.change_residue('FR1')
            amine_pdb.change_residue('FR2')
            
            # Now merge all pdbs into one
            build = ester_pdb.merge_with_another_molecule(amine_pdb)
            
        return build    
        
    def __acid_anhydride_amine(self, pdb1, pdb2, reaction, cis_or_trans_choice, side_to_keep_choice): # Note that pdb2 can be a primary or secondary amine
        """Simulates the reaction between an acid anhydride and an amine.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), containing either an acid anhydride or an amine.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        cis_or_trans_choice -- Both cis and trans products can be formed. This int dictates which product should be made.
        side_to_keep_choice -- The acid anhydride will be split. One half of the molecule will be discarded, and the other half will be incorporated into the product. This int determines which half to keep.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
    
        # pdb1 contains the first pdb data structure
        # pdb2 contains the second pdb data structure
        # reaction is a list of lists. It looks like this:
        #   ["NAME OF REACTIVE GROUP ON PDB1", INDICES OF THE ATOMS OF THAT PDB1 GROUP, "NAME OF REACTIVE GROUP ON PDB2", INDICES OF THE ATOMS OF THAT PDB2 GROUP, ]
        #   Here's a sample:
        #       ['ACID_ANHYDRIDE', [2, 1, 5, 6], 'AMINE', [3, 4, 11, 12]]
        # The order of the indices of the atoms of each group are given by the function defined in self.structure_locate_groups.py
    
        if cis_or_trans_choice == -1 : cis_or_trans_choice = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
        if side_to_keep_choice == -1 : side_to_keep_choice = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
        
        # You need to figure out which pdb is the carbonochloridate and which is the primary amine
        if reaction[0] == "ACID_ANHYDRIDE": # pdb1 is the carbonochloridate
            acid_anhydride_pdb = pdb1
            acid_anhydride_indecies = reaction[1]
        
            amine_pdb = pdb2
            amine_indecies = reaction[3]
            
        else: # pdb2 is the carbonochloridate
            acid_anhydride_pdb = pdb2
            acid_anhydride_indecies = reaction[3]
    
            amine_pdb = pdb1
            amine_indecies = reaction[1]
             
        # get easy names of the acid_anhydride, so we don't have to refer to the reaction[] list
        acid_anhydride_ether_oxgen = acid_anhydride_indecies[0]
        acid_anhydride_carybonyl_carbon_1 = acid_anhydride_indecies[1]
        acid_anhydride_carybonyl_oxygen_1 = acid_anhydride_indecies[2]
        acid_anhydride_carybonyl_carbon_2 = acid_anhydride_indecies[3]
        acid_anhydride_carybonyl_oxygen_2 = acid_anhydride_indecies[4]
    
        # get easy names for the amine atoms
        amine_nitrogen = amine_indecies[0]
        amine_Rs = []
        amine_hydrogens = []
        for t in range(1,len(amine_indecies)): # starting at index 1, so skipping amine_nitrogen
            index = amine_indecies[t]
            atom = amine_pdb.all_atoms[index]
            if atom.element == "H":
                amine_hydrogens.append(index)
            else :
                amine_Rs.append(index)
        
        # If it's a primary amine, they'll be only one element in amine_Rs list. Move one of the hydrogens to that list to complete it
        while len(amine_Rs) == 1:
            amine_Rs.append(amine_hydrogens[0])
            amine_hydrogens.pop(0)
            
        # I guess this would technically accept an amonium ion. Perhaps that should be added to the self.structure_locate_groups.py file...?
            
        # Now, load in the intermediate. This is as a model amide, the product of this reaction.
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb('./intermediates/amide.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_oppose_carbon_index = 3
        intermediate_carbonyl_carbon_index = 1
        intermediate_carbonyl_oxygen_index = 2
        intermediate_nitrogen_index = 4
        intermediate_amine_R1_index = 5
        intermediate_amine_R2_index = 6
        
        # Now, we need to pick which side of the acid anhydride you'll keep
        if side_to_keep_choice == 0 :
            acid_anhydride_pdb = acid_anhydride_pdb.get_branch(acid_anhydride_ether_oxgen, acid_anhydride_carybonyl_carbon_1)
            acid_anhydride_carybonyl_carbon = acid_anhydride_carybonyl_carbon_1
            acid_anhydride_carybonyl_oxygen = acid_anhydride_carybonyl_oxygen_1
    
        else :
            acid_anhydride_pdb = acid_anhydride_pdb.get_branch(acid_anhydride_ether_oxgen, acid_anhydride_carybonyl_carbon_2)
            acid_anhydride_carybonyl_carbon = acid_anhydride_carybonyl_carbon_2
            acid_anhydride_carybonyl_oxygen = acid_anhydride_carybonyl_oxygen_2
            
        
        # Now, posiition the acid_anhydride onto the model amide
        tethers = [[intermediate_carbonyl_carbon_index, acid_anhydride_carybonyl_carbon], [intermediate_carbonyl_oxygen_index, acid_anhydride_carybonyl_oxygen], [intermediate_nitrogen_index, acid_anhydride_ether_oxgen]]
        acid_anhydride_pdb = intermediate.align_another_molecule_to_this_one(acid_anhydride_pdb, tethers)
        # So what does thae above do? First, it looks at the first pair in the tethers list: [intermediate_carbonyl_carbon_index, acid_anhydride_acid_anhydride_carybonyl_carbon]
        # It moves pdb2 onto pdb1 such that the atoms corresponding to these two indices are right on top of each other. This point is a "pivot point."
        # Then it rotates pdb2 about that pivot point so as to minimize the distance between the atoms of the subsequent pairs.
    
        # Now the nitrogen on the amine is going to be going from SP3 to SP2. We need to disect this amine by removing it's branches
        # cis_or_trans_choice was either specified by the user or randomly set to either 0 or 1.
        # Use this varialbe in deciding which branch is which.
        if cis_or_trans_choice == 0 :
            amine_branch_R1_index = amine_Rs[0]
            amine_branch_R2_index = amine_Rs[1]
            
        else :
            amine_branch_R1_index = amine_Rs[1]
            amine_branch_R2_index = amine_Rs[0]
    
        # Is the nitrogen an SP2 hybridized nitrogen in a ring? If so, it can be processed, even though ring nitrogens generally can't
        sp2_nitrogen_in_ring = 'FALSE'
        if len(reaction) == 5:
            if reaction[4] == 'SECONDARY_AMINE_SP2_IN_RING':
                sp2_nitrogen_in_ring = 'TRUE'
    
        if sp2_nitrogen_in_ring == 'FALSE':
            # Make sure the amide is not in a loop. These must be discarded
            if amine_pdb.in_same_ring(amine_nitrogen, amine_branch_R2_index) == "TRUE":
                return pymolecule.Molecule() # return a blank Molecule
        
            amine_branch_1 = amine_pdb.get_branch(amine_nitrogen, amine_branch_R1_index)
            amine_branch_2 = amine_pdb.get_branch(amine_nitrogen, amine_branch_R2_index)
            
            # Now we need to position the amine branches onto the model amide
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R1_index, amine_branch_R1_index]]
            amine_branch_1 = intermediate.align_another_molecule_to_this_one(amine_branch_1, tethers)
        
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R2_index, amine_branch_R2_index]]
            amine_branch_2 = intermediate.align_another_molecule_to_this_one(amine_branch_2, tethers)
        
            # Now delete a few more atoms.
            acid_anhydride_pdb.delete_atom(acid_anhydride_ether_oxgen)
        
            # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
            # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
            # Note that the acid_anhydride is fixed and cannot rotate. The two branches of the amine can, though.
            item1 = [amine_branch_1, amine_branch_1.all_atoms[amine_nitrogen].coordinates, amine_branch_1.all_atoms[amine_branch_R1_index].coordinates]
            item2 = [amine_branch_2, amine_branch_2.all_atoms[amine_nitrogen].coordinates, amine_branch_2.all_atoms[amine_branch_R2_index].coordinates]
            item3 = [acid_anhydride_pdb] # These won't rotate, so the coordinates are omitted
            item4 = [intermediate]
            
            thelist = []
            thelist.append(item1)
            thelist.append(item2)
            thelist.append(item3)
            thelist.append(item4)
            self.structure_pdb_functions.reduce_steric_hindrance(thelist)
        
            # Now delete a few more atoms.
            amine_branch_1.delete_atom(amine_nitrogen)
        
            acid_anhydride_pdb.change_residue('FR1')
            amine_branch_1.change_residue('FR2')
            amine_branch_2.change_residue('FR3')
            
            # Now merge all pdbs into one
            build = acid_anhydride_pdb.merge_with_another_molecule(amine_branch_1)
            build = build.merge_with_another_molecule(amine_branch_2)
        elif sp2_nitrogen_in_ring == 'TRUE':
            
            # Now we need to position the amine branches onto the model amide
            tethers = [[intermediate_nitrogen_index, amine_nitrogen], [intermediate_amine_R1_index, amine_branch_R1_index], [intermediate_amine_R2_index, amine_branch_R2_index]]
            amine_pdb = intermediate.align_another_molecule_to_this_one(amine_pdb, tethers)
        
            # Now delete a few more atoms.
            acid_anhydride_pdb.delete_atom(acid_anhydride_ether_oxgen)
        
            # Now delete a few more atoms.
            amine_pdb.delete_atom(amine_hydrogens[0])
        
            acid_anhydride_pdb.change_residue('FR1')
            amine_pdb.change_residue('FR2')
        
            # Now merge all pdbs into one
            build = acid_anhydride_pdb.merge_with_another_molecule(amine_pdb)
            
        return build    
    
    def __carboxylate_alcohol(self, pdb1, pdb2, reaction): # This covers thiols and alcohols
        """Simulates the reaction between a carboxylate and an alcohol.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), containing either a carboxylate or an alcohol.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
    
        # pdb1 contains the first pdb data structure
        # pdb2 contains the second pdb data structure
        # reaction is a list of lists. It looks like this:
        #   ["NAME OF REACTIVE GROUP ON PDB1", INDICES OF THE ATOMS OF THAT PDB1 GROUP, "NAME OF REACTIVE GROUP ON PDB2", INDICES OF THE ATOMS OF THAT PDB2 GROUP, ]
        #   Here's a sample:
        #       ['CARBOXYLATE', [2, 1, 5, 6], 'ALCOHOL', [3, 4, 11, 12], 'O']
        # The order of the indices of the atoms of each group are given by the function defined in self.structure_locate_groups.py
    
        # You need to figure out which pdb is the carbonochloridate and which is the primary alcohol
        if reaction[0] == "CARBOXYLATE": # pdb1 is the carbonochloridate
            carboxylate_pdb = pdb1
            carboxylate_indecies = reaction[1]
        
            alcohol_pdb = pdb2
            alcohol_indecies = reaction[3]
            
        else: # pdb2 is the carbonochloridate
            carboxylate_pdb = pdb2
            carboxylate_indecies = reaction[3]
    
            alcohol_pdb = pdb1
            alcohol_indecies = reaction[1]
             
        # get easy names of the carboxylate, so we don't have to refer to the reaction[] list
        carboxylate_carybonyl_carbon = carboxylate_indecies[0]
        carboxylate_carybonyl_oxygen = carboxylate_indecies[1]
        carboxylate_hydroxyl_oxygen = carboxylate_indecies[2]
        carboxylate_hydroxyl_hydrogen = carboxylate_indecies[3]
    
        # get easy names for the alcohol atoms
        alcohol_neighbor = alcohol_indecies[0]
        alcohol_hydroxyl_oxygen = alcohol_indecies[1]
        alcohol_hydroxyl_hydrogen = alcohol_indecies[2]
            
        # Now, load in the intermediate. This is a model carboxylic acid or thiol acid, depending on whether we're dealing with an alcohol or a thiol.
        intermediate = pymolecule.Molecule()
        if (reaction[4] == "O"):
            intermediate.load_pdb('./intermediates/carboxylic_acid.pdb')
        else:
            intermediate.load_pdb('./intermediates/thio_acid.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_oppose_carbon_index = 4
        intermediate_carbonyl_carbon_index = 1
        intermediate_carbonyl_oxygen_index = 2
        intermediate_hydroxyl_oxygen_index = 3
        
        # Now, posiition the carboxylate onto the model acid
        tethers = [[intermediate_carbonyl_carbon_index, carboxylate_carybonyl_carbon], [intermediate_carbonyl_oxygen_index, carboxylate_carybonyl_oxygen], [intermediate_hydroxyl_oxygen_index, carboxylate_hydroxyl_oxygen]]
        carboxylate_pdb = intermediate.align_another_molecule_to_this_one(carboxylate_pdb, tethers)
    
        # Now we need to position the alcohol/thiol onto the model acid
        tethers = [[intermediate_hydroxyl_oxygen_index, alcohol_hydroxyl_oxygen], [intermediate_carbonyl_carbon_index, alcohol_hydroxyl_hydrogen]]
        alcohol_pdb = intermediate.align_another_molecule_to_this_one(alcohol_pdb, tethers)
    
        # Now delete a few atoms.
        carboxylate_pdb.delete_atom(carboxylate_hydroxyl_oxygen)
        if carboxylate_hydroxyl_hydrogen!=-1: carboxylate_pdb.delete_atom(carboxylate_hydroxyl_hydrogen)
    
        # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
        # Note that the carboxylate is fixed and cannot rotate. The two branches of the alcohol can, though.
        item1 = [alcohol_pdb, alcohol_pdb.all_atoms[alcohol_hydroxyl_oxygen].coordinates, alcohol_pdb.all_atoms[alcohol_hydroxyl_hydrogen].coordinates]
        item2 = [carboxylate_pdb] # These won't rotate, so the coordinates are omitted
        item3 = [intermediate]
        
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        thelist.append(item3)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        # Now delete a few more atoms.
        alcohol_pdb.delete_atom(alcohol_hydroxyl_hydrogen)
    
        carboxylate_pdb.change_residue('FR1')
        alcohol_pdb.change_residue('FR2')
        
        # Now merge all pdbs into one
        build = carboxylate_pdb.merge_with_another_molecule(alcohol_pdb)
        
        return build    
    
    def __acylhalide_alcohol(self, pdb1, pdb2, reaction): # This covers thiols and alcohols
        """Simulates the reaction between an acyl halide and an alcohol.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), containing either an acyl halide or an alcohopl.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
    
        # pdb1 contains the first pdb data structure
        # pdb2 contains the second pdb data structure
        # reaction is a list of lists. It looks like this:
        #   ["NAME OF REACTIVE GROUP ON PDB1", INDICES OF THE ATOMS OF THAT PDB1 GROUP, "NAME OF REACTIVE GROUP ON PDB2", INDICES OF THE ATOMS OF THAT PDB2 GROUP, ]
        #   Here's a sample:
        #       ['acylhalide', [2, 1, 5, 6], 'ALCOHOL', [3, 4, 11, 12], 'O']
        # The order of the indices of the atoms of each group are given by the function defined in self.structure_locate_groups.py
    
        # You need to figure out which pdb is the carbonochloridate and which is the primary alcohol
        if reaction[0] == "ACYLHALIDE": # pdb1 is the carbonochloridate
            acylhalide_pdb = pdb1
            acylhalide_indecies = reaction[1]
        
            alcohol_pdb = pdb2
            alcohol_indecies = reaction[3]
            
        else: # pdb2 is the carbonochloridate
            acylhalide_pdb = pdb2
            acylhalide_indecies = reaction[3]
    
            alcohol_pdb = pdb1
            alcohol_indecies = reaction[1]
             
        # get easy names of the acylhalide, so we don't have to refer to the reaction[] list
        acylhalide_carybonyl_carbon = acylhalide_indecies[0]
        acylhalide_carybonyl_oxygen = acylhalide_indecies[1]
        acylhalide_chloride = acylhalide_indecies[2]
    
        # get easy names for the alcohol atoms
        alcohol_neighbor = alcohol_indecies[0]
        alcohol_hydroxyl_oxygen = alcohol_indecies[1]
        alcohol_hydroxyl_hydrogen = alcohol_indecies[2]
            
        # Now, load in the intermediate. This is a model carboxylic acid or thiol acid, depending on whether we're dealing with an alcohol or a thiol.
        intermediate = pymolecule.Molecule()
        if (reaction[4] == "O"):
            intermediate.load_pdb('./intermediates/carboxylic_acid.pdb')
        else:
            intermediate.load_pdb('./intermediates/thio_acid.pdb')
            
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_oppose_carbon_index = 4
        intermediate_carbonyl_carbon_index = 1
        intermediate_carbonyl_oxygen_index = 2
        intermediate_hydroxyl_oxygen_index = 3
        
        # Now, posiition the acylhalide onto the model acid
        tethers = [[intermediate_carbonyl_carbon_index, acylhalide_carybonyl_carbon], [intermediate_carbonyl_oxygen_index, acylhalide_carybonyl_oxygen], [intermediate_hydroxyl_oxygen_index, acylhalide_chloride]]
        acylhalide_pdb = intermediate.align_another_molecule_to_this_one(acylhalide_pdb, tethers)
    
        # Now we need to position the alcohol/thiol onto the model acid
        tethers = [[intermediate_hydroxyl_oxygen_index, alcohol_hydroxyl_oxygen], [intermediate_carbonyl_carbon_index, alcohol_hydroxyl_hydrogen]]
        alcohol_pdb = intermediate.align_another_molecule_to_this_one(alcohol_pdb, tethers)
    
        # Now delete a few atoms.
        acylhalide_pdb.delete_atom(acylhalide_chloride)
    
        # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
        # Note that the acylhalide is fixed and cannot rotate. The two branches of the alcohol can, though.
        item1 = [alcohol_pdb, alcohol_pdb.all_atoms[alcohol_hydroxyl_oxygen].coordinates, alcohol_pdb.all_atoms[alcohol_hydroxyl_hydrogen].coordinates]
        item2 = [acylhalide_pdb] # These won't rotate, so the coordinates are omitted
        item3 = [intermediate]
        
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        thelist.append(item3)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        # Now delete a few more atoms.
        alcohol_pdb.delete_atom(alcohol_hydroxyl_hydrogen)
    
        acylhalide_pdb.change_residue('FR1')
        alcohol_pdb.change_residue('FR2')
        
        # Now merge all pdbs into one
        build = acylhalide_pdb.merge_with_another_molecule(alcohol_pdb)
        
        return build
    
    def __ester_alcohol(self, pdb1, pdb2, reaction): # Note that pdb2 can be a primary or secondary or tertiary alcohol
        """Simulates the reaction between an ester and an alcohol.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), containing either an ester or an alcohol.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
    
        # pdb1 contains the first pdb data structure
        # pdb2 contains the second pdb data structure
        # reaction is a list of lists. It looks like this:
        #   ["NAME OF REACTIVE GROUP ON PDB1", INDICES OF THE ATOMS OF THAT PDB1 GROUP, "NAME OF REACTIVE GROUP ON PDB2", INDICES OF THE ATOMS OF THAT PDB2 GROUP, ]
        #   Here's a sample:
        #       ['ESTER', [2, 1, 5, 6], 'ALCOHOL', [3, 4, 11, 12]]
        # The order of the indices of the atoms of each group are given by the function defined in self.structure_locate_groups.py
        
        # You need to figure out which pdb is the carbonochloridate and which is the primary alcohol
        if reaction[0] == "ESTER": # pdb1 is the carbonochloridate
            ester_pdb = pdb1
            ester_indecies = reaction[1]
        
            alcohol_pdb = pdb2
            alcohol_indecies = reaction[3]
            
        else: # pdb2 is the carbonochloridate
            ester_pdb = pdb2
            ester_indecies = reaction[3]
    
            alcohol_pdb = pdb1
            alcohol_indecies = reaction[1]
             
        # get easy names of the ester, so we don't have to refer to the reaction[] list
        ester_carybonyl_carbon = ester_indecies[0]
        ester_carybonyl_oxygen = ester_indecies[1]
        ester_oxygen = ester_indecies[2]
        ester_carbon_beyond_oxygen = ester_indecies[3]
    
        # get easy names for the alcohol atoms
        alcohol_neighbor = alcohol_indecies[0]
        alcohol_hydroxyl_oxygen = alcohol_indecies[1]
        alcohol_hydroxyl_hydrogen = alcohol_indecies[2]
            
        # Now, load in the intermediate. This is a model carboxylic acid or thiol acid, depending on whether we're dealing with an alcohol or a thiol.
        intermediate = pymolecule.Molecule()
        if (reaction[4] == "O"):
            intermediate.load_pdb('./intermediates/carboxylic_acid.pdb')
        else:
            intermediate.load_pdb('./intermediates/thio_acid.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_oppose_carbon_index = 4
        intermediate_carbonyl_carbon_index = 1
        intermediate_carbonyl_oxygen_index = 2
        intermediate_hydroxyl_oxygen_index = 3
        
        # Now, we need to discard the alcohol portion of the ester
        ester_pdb = ester_pdb.get_branch(ester_oxygen, ester_carybonyl_carbon)
        
        # Now, posiition the ester onto the model acid
        tethers = [[intermediate_carbonyl_carbon_index, ester_carybonyl_carbon], [intermediate_carbonyl_oxygen_index, ester_carybonyl_oxygen], [intermediate_hydroxyl_oxygen_index, ester_oxygen]]
        ester_pdb = intermediate.align_another_molecule_to_this_one(ester_pdb, tethers)
    
        # Now we need to position the alcohol onto the model amide
        tethers = [[intermediate_hydroxyl_oxygen_index, alcohol_hydroxyl_oxygen], [intermediate_carbonyl_carbon_index, alcohol_hydroxyl_hydrogen]]
        alcohol_pdb = intermediate.align_another_molecule_to_this_one(alcohol_pdb, tethers)
    
        # Now delete a few atoms.
        ester_pdb.delete_atom(ester_oxygen)
    
        # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
        # Note that the ester is fixed and cannot rotate. The two branches of the alcohol can, though.
        item1 = [alcohol_pdb, alcohol_pdb.all_atoms[alcohol_hydroxyl_oxygen].coordinates, alcohol_pdb.all_atoms[alcohol_hydroxyl_hydrogen].coordinates]
        item2 = [ester_pdb] # These won't rotate, so the coordinates are omitted
        item3 = [intermediate]
        
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        thelist.append(item3)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        # Now delete a few more atoms.
        alcohol_pdb.delete_atom(alcohol_hydroxyl_hydrogen)
    
        ester_pdb.change_residue('FR1')
        alcohol_pdb.change_residue('FR2')
        
        # Now merge all pdbs into one
        build = ester_pdb.merge_with_another_molecule(alcohol_pdb)
        
        return build    
    
    def __acid_anhydride_alcohol(self, pdb1, pdb2, reaction, side_to_keep_choice): # Note that pdb2 can be a primary or secondary alcohol
        """Simulates the reaction between an acid anhydride and an alcohol.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), containing either an acid anhydride or an alcohol.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        side_to_keep_choice -- The acid anhydride will be split. One half of the molecule will be discarded, and the other half will be incorporated into the product. This int determines which half to keep.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
    
        # pdb1 contains the first pdb data structure
        # pdb2 contains the second pdb data structure
        # reaction is a list of lists. It looks like this:
        #   ["NAME OF REACTIVE GROUP ON PDB1", INDICES OF THE ATOMS OF THAT PDB1 GROUP, "NAME OF REACTIVE GROUP ON PDB2", INDICES OF THE ATOMS OF THAT PDB2 GROUP, ]
        #   Here's a sample:
        #       ['ACID_ANHYDRIDE', [2, 1, 5, 6], 'ALCOHOL', [3, 4, 11, 12]]
        # The order of the indices of the atoms of each group are given by the function defined in self.structure_locate_groups.py
    
        if side_to_keep_choice == -1 : side_to_keep_choice = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
        
        # You need to figure out which pdb is the carbonochloridate and which is the primary alcohol
        if reaction[0] == "ACID_ANHYDRIDE": # pdb1 is the carbonochloridate
            acid_anhydride_pdb = pdb1
            acid_anhydride_indecies = reaction[1]
        
            alcohol_pdb = pdb2
            alcohol_indecies = reaction[3]
            
        else: # pdb2 is the carbonochloridate
            acid_anhydride_pdb = pdb2
            acid_anhydride_indecies = reaction[3]
    
            alcohol_pdb = pdb1
            alcohol_indecies = reaction[1]
             
        # get easy names of the acid_anhydride, so we don't have to refer to the reaction[] list
        acid_anhydride_ether_oxgen = acid_anhydride_indecies[0]
        acid_anhydride_carybonyl_carbon_1 = acid_anhydride_indecies[1]
        acid_anhydride_carybonyl_oxygen_1 = acid_anhydride_indecies[2]
        acid_anhydride_carybonyl_carbon_2 = acid_anhydride_indecies[3]
        acid_anhydride_carybonyl_oxygen_2 = acid_anhydride_indecies[4]
        
        # get easy names for the alcohol atoms
        alcohol_neighbor = alcohol_indecies[0]
        alcohol_hydroxyl_oxygen = alcohol_indecies[1]
        alcohol_hydroxyl_hydrogen = alcohol_indecies[2]
            
        # Now, load in the intermediate. This is a model carboxylic acid or thiol acid, depending on whether we're dealing with an alcohol or a thiol.
        intermediate = pymolecule.Molecule()
        if (reaction[4] == "O"):
            intermediate.load_pdb('./intermediates/carboxylic_acid.pdb')
        else:
            intermediate.load_pdb('./intermediates/thio_acid.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_oppose_carbon_index = 4
        intermediate_carbonyl_carbon_index = 1
        intermediate_carbonyl_oxygen_index = 2
        intermediate_hydroxyl_oxygen_index = 3
    
        # Now, we need to pick which side of the acid anhydride you'll keep
        
        if side_to_keep_choice == 0 :
            acid_anhydride_pdb = acid_anhydride_pdb.get_branch(acid_anhydride_ether_oxgen, acid_anhydride_carybonyl_carbon_1)
            acid_anhydride_carybonyl_carbon = acid_anhydride_carybonyl_carbon_1
            acid_anhydride_carybonyl_oxygen = acid_anhydride_carybonyl_oxygen_1
    
        else :
            acid_anhydride_pdb = acid_anhydride_pdb.get_branch(acid_anhydride_ether_oxgen, acid_anhydride_carybonyl_carbon_2)
            acid_anhydride_carybonyl_carbon = acid_anhydride_carybonyl_carbon_2
            acid_anhydride_carybonyl_oxygen = acid_anhydride_carybonyl_oxygen_2
            
        
        # Now, posiition the acid_anhydride onto the model acid
        tethers = [[intermediate_carbonyl_carbon_index, acid_anhydride_carybonyl_carbon], [intermediate_carbonyl_oxygen_index, acid_anhydride_carybonyl_oxygen], [intermediate_hydroxyl_oxygen_index, acid_anhydride_ether_oxgen]]
        acid_anhydride_pdb = intermediate.align_another_molecule_to_this_one(acid_anhydride_pdb, tethers)
        # So what does thae above do? First, it looks at the first pair in the tethers list: [intermediate_carbonyl_carbon_index, acid_anhydride_acid_anhydride_carybonyl_carbon]
        # It moves pdb2 onto pdb1 such that the atoms corresponding to these two indices are right on top of each other. This point is a "pivot point."
        # Then it rotates pdb2 about that pivot point so as to minimize the distance between the atoms of the subsequent pairs.
    
        # Now we need to position the alcohol onto the model amide
        tethers = [[intermediate_hydroxyl_oxygen_index, alcohol_hydroxyl_oxygen], [intermediate_carbonyl_carbon_index, alcohol_hydroxyl_hydrogen]]
        alcohol_pdb = intermediate.align_another_molecule_to_this_one(alcohol_pdb, tethers)
    
        # Now delete a few more atoms.
        acid_anhydride_pdb.delete_atom(acid_anhydride_ether_oxgen)
    
        # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
        # Note that the acid_anhydride is fixed and cannot rotate. The two branches of the alcohol can, though.
        item1 = [alcohol_pdb, alcohol_pdb.all_atoms[alcohol_hydroxyl_oxygen].coordinates, alcohol_pdb.all_atoms[alcohol_hydroxyl_hydrogen].coordinates]
        item2 = [acid_anhydride_pdb] # These won't rotate, so the coordinates are omitted
        item3 = [intermediate]
        
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        thelist.append(item3)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        # Now delete a few more atoms.
        alcohol_pdb.delete_atom(alcohol_hydroxyl_hydrogen)
    
        acid_anhydride_pdb.change_residue('FR1')
        alcohol_pdb.change_residue('FR2')
        
        # Now merge all pdbs into one
        build = acid_anhydride_pdb.merge_with_another_molecule(alcohol_pdb)
        
        return build    
    
    def isocyanate_amine(self, pdb1, pdb2, reaction, amide_cis_or_trans_choice, nucleophile_amine_cis_or_trans_choice): # Note that pdb2 can be a primary or secondary amine. Note that this also works for isothiocyanates
        """Simulates the reaction between an isocyanate and an amine.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), containing either an isocyanate or an amine.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        amide_cis_or_trans_choice -- Both cis and trans amide products can be formed. This int dictates which product should be made.
        nucleophile_amine_cis_or_trans_choice -- Both cis and trans nucleophilic-amine products can be formed. This int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        # pdb1 contains the first pdb data structure
        # pdb2 contains the second pdb data structure
        # reaction is a list of lists. It looks like this:
        #   ["NAME OF REACTIVE GROUP ON PDB1", INDICES OF THE ATOMS OF THAT PDB1 GROUP, "NAME OF REACTIVE GROUP ON PDB2", INDICES OF THE ATOMS OF THAT PDB2 GROUP, ]
        #   Here's a sample:
        #       ['ISOCYANATE', [2, 1, 5, 6], 'AMINE', [3, 4, 11, 12]]
        # The order of the indices of the atoms of each group are given by the function defined in self.structure_locate_groups.py
    
        if nucleophile_amine_cis_or_trans_choice == -1 : nucleophile_amine_cis_or_trans_choice = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
        if amide_cis_or_trans_choice == -1 : amide_cis_or_trans_choice = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
            
        # You need to figure out which pdb is the carbonochloridate and which is the primary amine
        if reaction[0] == "ISOCYANATE" or reaction[0] == "ISOTHIOCYANATE": # pdb1 is the isocyanate or isothiocyanate
            isocyanate_pdb = pdb1
            isocyanate_indecies = reaction[1]
        
            amine_pdb = pdb2
            amine_indecies = reaction[3]
            
        else: # pdb2 is the carbonochloridate
            isocyanate_pdb = pdb2
            isocyanate_indecies = reaction[3]
    
            amine_pdb = pdb1
            amine_indecies = reaction[1]
             
        # get easy names of the isocyanate, so we don't have to refer to the reaction[] list
        isocyanate_neighbor_carbon = isocyanate_indecies[0]
        isocyanate_proximal_nitrogen = isocyanate_indecies[1]
        isocyanate_medial_carbon = isocyanate_indecies[2]
        isocyanate_distal_oxygen = isocyanate_indecies[3] # or sulfur
    
        # get easy names for the amine atoms
        amine_nitrogen = amine_indecies[0]
        amine_Rs = []
        amine_hydrogens = []
        for t in range(1,len(amine_indecies)): # starting at index 1, so skipping amine_nitrogen
            index = amine_indecies[t]
            atom = amine_pdb.all_atoms[index]
            if atom.element == "H":
                amine_hydrogens.append(index)
            else :
                amine_Rs.append(index)
        
        # If it's a primary amine, they'll be only one element in amine_Rs list. Move one of the hydrogens to that list to complete it
        while len(amine_Rs) == 1:
            amine_Rs.append(amine_hydrogens[0])
            amine_hydrogens.pop(0)
            
        # I guess this would technically accept an amonium ion. Perhaps that should be added to the self.structure_locate_groups.py file...?
            
        # Now, load in the intermediate. This is as a model urea (or thio urea), the product of this reaction.
        # But the intermediate can be a urea or a thiourea. Figure out which one.
        element = isocyanate_pdb.all_atoms[isocyanate_distal_oxygen].element
        
        # Now, the amide could be cis or trans
        intermediate = pymolecule.Molecule()
        if amide_cis_or_trans_choice == 0: # so cis
            if element == "O": # so urea
                intermediate.load_pdb('./intermediates/urea_cis.pdb')
            else: # so thiourea
                intermediate.load_pdb('./intermediates/thiourea_cis.pdb')
        else: # so trans
            if element == "O": # so urea
                intermediate.load_pdb('./intermediates/urea_trans.pdb')
            else: # so thiourea
                intermediate.load_pdb('./intermediates/thiourea_trans.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
    
        intermediate_carbonyl_carbon_index = 3
        intermediate_carbonyl_oxygen_index = 5 # or sulfur
        intermediate_nitrogen_nucleophile_index = 4
        intermediate_amine_R1_index = 7
        intermediate_amine_R2_index = 8
        intermediate_oppose_nitrogen_index = 2
        intermediate_oppose_nitrogen_neighbor_index = 1
        
        # Now, posiition the isocyanate onto the model amide
        tethers = [[intermediate_oppose_nitrogen_index, isocyanate_proximal_nitrogen], [intermediate_oppose_nitrogen_neighbor_index, isocyanate_neighbor_carbon]]
        isocyanate_pdb = intermediate.align_another_molecule_to_this_one(isocyanate_pdb, tethers)
        # So what does thae above do? First, it looks at the first pair in the tethers list: [intermediate_carbonyl_carbon_index, isocyanate_ester_carybonyl_carbon]
        # It moves pdb2 onto pdb1 such that the atoms corresponding to these two indices are right on top of each other. This point is a "pivot point."
        # Then it rotates pdb2 about that pivot point so as to minimize the distance between the atoms of the subsequent pairs.
    
        # Now the nitrogen on the amine is going to be going from SP3 to SP2. We need to disect this amine by removing it's branches
        # cis_or_trans_choice was either specified by the user or randomly set to either 0 or 1.
        # Use this varialbe in deciding which branch is which.
        if nucleophile_amine_cis_or_trans_choice  == 0 :
            amine_branch_R1_index = amine_Rs[0]
            amine_branch_R2_index = amine_Rs[1]
            
        else :
            amine_branch_R1_index = amine_Rs[1]
            amine_branch_R2_index = amine_Rs[0]
    
        # Is the nitrogen an SP2 hybridized nitrogen in a ring? If so, it can be processed, even though ring nitrogens generally can't
        sp2_nitrogen_in_ring = 'FALSE'
        if len(reaction) == 5:
            if reaction[4] == 'SECONDARY_AMINE_SP2_IN_RING':
                sp2_nitrogen_in_ring = 'TRUE'
    
        if sp2_nitrogen_in_ring == 'FALSE':
            # Make sure the amide is not in a loop. These must be discarded
            if amine_pdb.in_same_ring(amine_nitrogen, amine_branch_R2_index) == "TRUE":
                return pymolecule.Molecule() # return a blank Molecule
        
            amine_branch_1 = amine_pdb.get_branch(amine_nitrogen, amine_branch_R1_index)
            amine_branch_2 = amine_pdb.get_branch(amine_nitrogen, amine_branch_R2_index)
            
            # Now we need to position the amine branches onto the model amide
            tethers = [[intermediate_nitrogen_nucleophile_index, amine_nitrogen], [intermediate_amine_R1_index, amine_branch_R1_index]]
            amine_branch_1 = intermediate.align_another_molecule_to_this_one(amine_branch_1, tethers)
        
            tethers = [[intermediate_nitrogen_nucleophile_index, amine_nitrogen], [intermediate_amine_R2_index, amine_branch_R2_index]]
            amine_branch_2 = intermediate.align_another_molecule_to_this_one(amine_branch_2, tethers)
        
            # Now delete a few more atoms.
            isocyanate_pdb.delete_atom(isocyanate_medial_carbon)
            isocyanate_pdb.delete_atom(isocyanate_distal_oxygen)
        
            # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
            # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
            # Note that the isocyanate is fixed and cannot rotate. The two branches of the amine can, though.
            item1 = [amine_branch_1, amine_branch_1.all_atoms[amine_nitrogen].coordinates, amine_branch_1.all_atoms[amine_branch_R1_index].coordinates]
            item2 = [amine_branch_2, amine_branch_2.all_atoms[amine_nitrogen].coordinates, amine_branch_2.all_atoms[amine_branch_R2_index].coordinates]
            item3 = [isocyanate_pdb, isocyanate_pdb.all_atoms[isocyanate_neighbor_carbon].coordinates, isocyanate_pdb.all_atoms[isocyanate_proximal_nitrogen].coordinates] 
            item4 = [intermediate]
                
            thelist = []
            thelist.append(item1)
            thelist.append(item2)
            thelist.append(item3)
            thelist.append(item4)
            self.structure_pdb_functions.reduce_steric_hindrance(thelist)
        
            # Now delete a few more atoms.
            amine_branch_1.delete_atom(amine_nitrogen)
            isocyanate_pdb.delete_atom(isocyanate_proximal_nitrogen)
            intermediate.delete_atom(intermediate_nitrogen_nucleophile_index)
            intermediate.delete_atom(intermediate_amine_R1_index)
            intermediate.delete_atom(intermediate_amine_R2_index)
        
            intermediate.delete_atom(intermediate_oppose_nitrogen_neighbor_index)
            
            intermediate.change_residue('FR1')
            isocyanate_pdb.change_residue('FR2')
            amine_branch_1.change_residue('FR3')
            amine_branch_2.change_residue('FR4')
            
            # Now merge all pdbs into one
            build = intermediate.merge_with_another_molecule(isocyanate_pdb)
            build = build.merge_with_another_molecule(amine_branch_1)
            build = build.merge_with_another_molecule(amine_branch_2)
        elif sp2_nitrogen_in_ring == 'TRUE':
        
            # Now we need to position the amine branches onto the model amide
            tethers = [[intermediate_nitrogen_nucleophile_index, amine_nitrogen], [intermediate_amine_R1_index, amine_branch_R1_index], [intermediate_amine_R2_index, amine_branch_R2_index]]
            amine_pdb = intermediate.align_another_molecule_to_this_one(amine_pdb, tethers)
        
            # Now delete a few more atoms.
            isocyanate_pdb.delete_atom(isocyanate_medial_carbon)
            isocyanate_pdb.delete_atom(isocyanate_distal_oxygen)        
            
            # Now delete a few more atoms.
            amine_pdb.delete_atom(amine_hydrogens[0])
            isocyanate_pdb.delete_atom(isocyanate_proximal_nitrogen)
            intermediate.delete_atom(intermediate_nitrogen_nucleophile_index)
            intermediate.delete_atom(intermediate_amine_R1_index)
            intermediate.delete_atom(intermediate_amine_R2_index)
            intermediate.delete_atom(intermediate_oppose_nitrogen_neighbor_index)
            
            intermediate.change_residue('FR1')
            isocyanate_pdb.change_residue('FR2')
            amine_pdb.change_residue('FR3')
            
            # Now merge all pdbs into one
            build = intermediate.merge_with_another_molecule(isocyanate_pdb)
            build = build.merge_with_another_molecule(amine_pdb)
            
        return build
    
    def __isocyanate_alcohol(self, pdb1, pdb2, reaction, amide_cis_or_trans_choice): # This works for isocyanates, isothiocyanates, alcohols, thiols.
        """Simulates the reaction between an isocyanate and an alcohol.
        
        Arguments:
        pdb1 -- The first molecular model (pymolecule.Molecule), containing either an isocyanate or an alcohol.
        pdb2 -- The other molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        amide_cis_or_trans_choice -- Both cis and trans amide products can be formed. This int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
            # pdb1 contains the first pdb data structure
        # pdb2 contains the second pdb data structure
        # reaction is a list of lists. It looks like this:
        #   ["NAME OF REACTIVE GROUP ON PDB1", INDICES OF THE ATOMS OF THAT PDB1 GROUP, "NAME OF REACTIVE GROUP ON PDB2", INDICES OF THE ATOMS OF THAT PDB2 GROUP, ]
        #   Here's a sample:
        #       ['ISOCYANATE', [2, 1, 5, 6], 'ALCOHOL', [3, 4, 11, 12]]
        # The order of the indices of the atoms of each group are given by the function defined in self.structure_locate_groups.py
    
        if amide_cis_or_trans_choice == -1 : amide_cis_or_trans_choice = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
            
        # You need to figure out which pdb is the isocyanate and which is the primary alcohol
        if reaction[0] == "ISOCYANATE" or reaction[0] == "ISOTHIOCYANATE": # pdb1 is the isocyanate or isothiocyanate
            isocyanate_pdb = pdb1
            isocyanate_indecies = reaction[1]
        
            alcohol_pdb = pdb2
            alcohol_indecies = reaction[3]
            
        else: # pdb2 is the isocyanate
            isocyanate_pdb = pdb2
            isocyanate_indecies = reaction[3]
    
            alcohol_pdb = pdb1
            alcohol_indecies = reaction[1]
             
        # get easy names of the isocyanate, so we don't have to refer to the reaction[] list
        isocyanate_neighbor_carbon = isocyanate_indecies[0]
        isocyanate_proximal_nitrogen = isocyanate_indecies[1]
        isocyanate_medial_carbon = isocyanate_indecies[2]
        isocyanate_distal_oxygen = isocyanate_indecies[3] # or sulfur
    
        # get easy names for the alcohol atoms
        alcohol_hydroxyl_oxygen = alcohol_indecies[1]
        alcohol_hydroxyl_hydrogen = alcohol_indecies[2]
    
            
        # Now, load in the intermediate. This is as a model urea (or thio urea), the product of this reaction.
        # But the intermediate can be a carbamate, carbamothioate, or carbamodithioate. Figure out which one.
        isocyanate_element = isocyanate_pdb.all_atoms[isocyanate_distal_oxygen].element
        alcohol_element = alcohol_pdb.all_atoms[alcohol_hydroxyl_oxygen].element
        
        # Now, the amide could be cis or trans
        intermediate = pymolecule.Molecule()
    
        if amide_cis_or_trans_choice == 0: # so cis
            if isocyanate_element == "O": # so carbonyl
                if alcohol_element == "O": # so alcohol
                    # cis, carbonyl, alcohol
                    intermediate.load_pdb('./intermediates/carbamate_2_cis.pdb')
                else: # so thiol
                    # cis, carbonyl, thiol
                    intermediate.load_pdb('./intermediates/carbamothioate_2_cis.pdb')
            else: # so thiocarbonyl
                if alcohol_element == "O": # so alcohol
                    # cis, thiocarbonyl, alcohol
                    intermediate.load_pdb('./intermediates/carbamothioate_1_cis.pdb')
                else: # so thiol
                    # cis, thiocarbonyl, thiol
                    intermediate.load_pdb('./intermediates/carbamodithioate_cis.pdb')
        else: # so trans
            if isocyanate_element == "O": # so carbonyl
                if alcohol_element == "O": # so alcohol
                    # trans, carbonyl, alcohol
                    intermediate.load_pdb('./intermediates/carbamate_2_trans.pdb')
                else: # so thiol
                    # trans, carbonyl, thiol
                    intermediate.load_pdb('./intermediates/carbamothioate_2_trans.pdb')
            else: # so thiocarbonyl
                if alcohol_element == "O": # so alcohol
                    # trans, thiocarbonyl, alcohol
                    intermediate.load_pdb('./intermediates/carbamothioate_1_trans.pdb')
                else: # so thiol
                    # trans, thiocarbonyl, thiol
                    intermediate.load_pdb('./intermediates/carbamodithioate_trans.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
    
        intermediate_carbonyl_carbon_index = 3
        intermediate_carbonyl_oxygen_index = 5 # or sulfur
        intermediate_oxygen_nucleophile_index = 4 # or sulfur
        intermediate_oppose_nitrogen_index = 2
        intermediate_oppose_nitrogen_neighbor_index = 1
        
        # Now, posiition the isocyanate onto the model amide
        tethers = [[intermediate_oppose_nitrogen_index, isocyanate_proximal_nitrogen], [intermediate_oppose_nitrogen_neighbor_index, isocyanate_neighbor_carbon]]
        isocyanate_pdb = intermediate.align_another_molecule_to_this_one(isocyanate_pdb, tethers)
        # So what does thae above do? First, it looks at the first pair in the tethers list: [intermediate_carbonyl_carbon_index, isocyanate_ester_carybonyl_carbon]
        # It moves pdb2 onto pdb1 such that the atoms corresponding to these two indices are right on top of each other. This point is a "pivot point."
        # Then it rotates pdb2 about that pivot point so as to minimize the distance between the atoms of the subsequent pairs.
    
        # Now we need to position the alcohol branches onto the model amide
        tethers = [[intermediate_oxygen_nucleophile_index, alcohol_hydroxyl_oxygen], [intermediate_carbonyl_carbon_index, alcohol_hydroxyl_hydrogen]]
        alcohol_pdb = intermediate.align_another_molecule_to_this_one(alcohol_pdb, tethers)
    
        # Now delete a few more atoms.
        isocyanate_pdb.delete_atom(isocyanate_medial_carbon)
        isocyanate_pdb.delete_atom(isocyanate_distal_oxygen)
    
        # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
        # Note that the isocyanate is fixed and cannot rotate. The two branches of the alcohol can, though.
        item1 = [alcohol_pdb, alcohol_pdb.all_atoms[alcohol_hydroxyl_oxygen].coordinates, alcohol_pdb.all_atoms[alcohol_hydroxyl_hydrogen].coordinates]
        item2 = [isocyanate_pdb, isocyanate_pdb.all_atoms[isocyanate_neighbor_carbon].coordinates, isocyanate_pdb.all_atoms[isocyanate_proximal_nitrogen].coordinates] # These won't rotate, so the coordinates are omitted
        item3 = [intermediate]
    
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        thelist.append(item3)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        # Now delete a few more atoms.
        alcohol_pdb.delete_atom(alcohol_hydroxyl_hydrogen)
        isocyanate_pdb.delete_atom(isocyanate_proximal_nitrogen)
        intermediate.delete_atom(intermediate_oxygen_nucleophile_index)
        intermediate.delete_atom(intermediate_oppose_nitrogen_neighbor_index)
    
        intermediate.change_residue('FR1')
        isocyanate_pdb.change_residue('FR2')
        alcohol_pdb.change_residue('FR3')
        
        # Now merge all pdbs into one
        build = intermediate.merge_with_another_molecule(isocyanate_pdb)
        build = build.merge_with_another_molecule(alcohol_pdb)
        
        return build
    
    
    def __tertiary_or_primary_halide_to_azide(self, pdb1, reaction, cyanide_or_azide): # This works for tertiary_or_primary_halides. These are the ones that have no inversions of their stereocenters.
        """Simulates the conversion of a tertiary or primary halide to an azide or cyanide.
        
        Arguments:
        pdb1 -- The molecular model (pymolecule.Molecule), containing a tertiary or primary halide.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        cyanide_or_azide -- Either azide or cyanide products are possible. This int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        if cyanide_or_azide == -1 : cyanide_or_azide = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
    
        tertiary_or_primary_halide_pdb = pdb1
        tertiary_or_primary_halide_indecies = reaction[1]
    
        # get easy names of the tertiary_or_primary_halide, so we don't have to refer to the reaction[] list
        tertiary_or_primary_halide_index = tertiary_or_primary_halide_indecies[4]
        tertiary_or_primary_halide_proximal_carbon_index = tertiary_or_primary_halide_indecies[3]
            
        # Now, load in the intermediate. This is as a model azide
        intermediate = pymolecule.Molecule()
        if cyanide_or_azide == 0:
            intermediate.load_pdb('./intermediates/intermediate9.pdb')
        elif cyanide_or_azide == 1:
            intermediate.load_pdb('./intermediates/cyanide.pdb')
        
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
    
        intermediate_carbon_index = 1
        intermediate_first_nitrogen_index = 2
    
        # Now, posiition the tertiary_or_primary_halide onto the model amide
        tethers = [[tertiary_or_primary_halide_proximal_carbon_index, intermediate_carbon_index], [tertiary_or_primary_halide_index, intermediate_first_nitrogen_index]]
        intermediate = tertiary_or_primary_halide_pdb.align_another_molecule_to_this_one(intermediate, tethers)
        
        # Now delete a few more atoms.
        intermediate.delete_atom(intermediate_carbon_index)
        tertiary_or_primary_halide_pdb.delete_atom(tertiary_or_primary_halide_index)
    
        # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
        # Note that the tertiary_or_primary_halide is fixed and cannot rotate. The two branches of the alcohol can, though.
        item1 = [intermediate, intermediate.all_atoms[intermediate_first_nitrogen_index].coordinates, tertiary_or_primary_halide_pdb.all_atoms[tertiary_or_primary_halide_proximal_carbon_index].coordinates] # These won't rotate, so the coordinates are omitted
        item2 = [tertiary_or_primary_halide_pdb]
    
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        tertiary_or_primary_halide_pdb.change_residue('FR1')
        intermediate.change_residue('FR2')
        
        # Now merge all pdbs into one
        build = tertiary_or_primary_halide_pdb.merge_with_another_molecule(intermediate)
        
        return build
    
    def __secondary_halide_to_azide(self, pdb1, reaction, cyanide_or_azide): # This works for secondary_halides. These are the ones that do have inversions of their stereocenters.
        """Simulates the conversion of a secondary halide to an azide or cyanide.
        
        Arguments:
        pdb1 -- The molecular model (pymolecule.Molecule), containing a tertiary or primary halide.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        cyanide_or_azide -- Either azide or cyanide products are possible. This int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        if cyanide_or_azide == -1 : cyanide_or_azide = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
    
        secondary_halide_pdb = pdb1
        secondary_halide_indecies = reaction[1]
    
        # get easy names of the secondary_halide, so we don't have to refer to the reaction[] list
        secondary_halide_index = secondary_halide_indecies[4]
        secondary_halide_proximal_carbon_index = secondary_halide_indecies[3]
        
        for t in range(0,3):
            if pdb1.all_atoms[secondary_halide_indecies[t]].element == "H":
                secondary_halide_hydrogen_index = secondary_halide_indecies[t]
    
        # Now, load in the intermediate. This is as a model azide
        intermediate = pymolecule.Molecule()
        if cyanide_or_azide == 0:
            intermediate.load_pdb('./intermediates/intermediate9.pdb')
        elif cyanide_or_azide == 1:
            intermediate.load_pdb('./intermediates/cyanide.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_carbon_index = 1
        intermediate_first_nitrogen_index = 2
    
        # Now, posiition the secondary_halide onto the model amide
        tethers = [[secondary_halide_proximal_carbon_index, intermediate_carbon_index], [secondary_halide_hydrogen_index, intermediate_first_nitrogen_index]]
        intermediate = secondary_halide_pdb.align_another_molecule_to_this_one(intermediate, tethers)
        
        # Now get the hydrogen fragment
        hydrogen_frag = secondary_halide_pdb.get_branch(secondary_halide_proximal_carbon_index, secondary_halide_hydrogen_index)
        
        # Now position the hydrogen frag
        tethers = [[secondary_halide_proximal_carbon_index, secondary_halide_proximal_carbon_index], [secondary_halide_index, secondary_halide_hydrogen_index]]
        hydrogen_frag = secondary_halide_pdb.align_another_molecule_to_this_one(hydrogen_frag, tethers)
    
        # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
        # Note that the secondary_halide is fixed and cannot rotate. The two branches of the alcohol can, though.
        item1 = [intermediate, intermediate.all_atoms[intermediate_first_nitrogen_index].coordinates, intermediate.all_atoms[intermediate_carbon_index].coordinates] # These won't rotate, so the coordinates are omitted
        item2 = [secondary_halide_pdb]
        item3 = [hydrogen_frag]
    
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        thelist.append(item3)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        # Now delete a few more atoms.
        intermediate.delete_atom(intermediate_carbon_index)
        hydrogen_frag.delete_atom(secondary_halide_proximal_carbon_index)
        secondary_halide_pdb.delete_atom(secondary_halide_index)
        secondary_halide_pdb.delete_atom(secondary_halide_hydrogen_index)
        
        secondary_halide_pdb.change_residue('FR1')
        intermediate.change_residue('FR2')
        hydrogen_frag.change_residue('FR3')
        
        # Now merge all pdbs into one
        build = secondary_halide_pdb.merge_with_another_molecule(intermediate)
        build = build.merge_with_another_molecule(hydrogen_frag)
        
        return build
    
    def __halide_bound_to_sp2_carbon_to_azide(self, pdb1, reaction, cyanide_or_azide): # This works for secondary_halides. These are the ones that do have inversions of their stereocenters.
        """Simulates the conversion of a halide bound to an sp2-hybridized carbon, to an azide or cyanide.
        
        Arguments:
        pdb1 -- The molecular model (pymolecule.Molecule), containing a tertiary or primary halide.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        cyanide_or_azide -- Either azide or cyanide products are possible. This int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        if cyanide_or_azide == -1 : cyanide_or_azide = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
    
        halide_bound_to_sp2_carbon_pdb = pdb1
        halide_bound_to_sp2_carbon_indecies = reaction[1]
    
        # get easy names of the halide_bound_to_sp2_carbon, so we don't have to refer to the reaction[] list
        halide_bound_to_sp2_carbon_index = halide_bound_to_sp2_carbon_indecies[3]
        halide_bound_to_sp2_carbon_proximal_carbon_index = halide_bound_to_sp2_carbon_indecies[2]
        
        for t in range(0,2):
            if pdb1.all_atoms[halide_bound_to_sp2_carbon_indecies[t]].element == "H":
                halide_bound_to_sp2_carbon_hydrogen_index = halide_bound_to_sp2_carbon_indecies[t]
    
        # Now, load in the intermediate. This is as a model azide
        intermediate = pymolecule.Molecule()
        if cyanide_or_azide == 0:
            intermediate.load_pdb('./intermediates/intermediate9.pdb')
        elif cyanide_or_azide == 1:
            intermediate.load_pdb('./intermediates/cyanide.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_carbon_index = 1
        intermediate_first_nitrogen_index = 2
    
        # Now, posiition the halide_bound_to_sp2_carbon onto the model amide
        tethers = [[halide_bound_to_sp2_carbon_proximal_carbon_index, intermediate_carbon_index], [halide_bound_to_sp2_carbon_index, intermediate_first_nitrogen_index]]
        intermediate = halide_bound_to_sp2_carbon_pdb.align_another_molecule_to_this_one(intermediate, tethers)
    
        # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
        # Note that the halide_bound_to_sp2_carbon is fixed and cannot rotate. The two branches of the alcohol can, though.
        item1 = [intermediate, intermediate.all_atoms[intermediate_first_nitrogen_index].coordinates, intermediate.all_atoms[intermediate_carbon_index].coordinates] # These won't rotate, so the coordinates are omitted
        item2 = [halide_bound_to_sp2_carbon_pdb]
    
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        # Now delete a few more atoms.
        intermediate.delete_atom(intermediate_carbon_index)
        halide_bound_to_sp2_carbon_pdb.delete_atom(halide_bound_to_sp2_carbon_index)
         
        halide_bound_to_sp2_carbon_pdb.change_residue('FR1')
        intermediate.change_residue('FR2')
         
        # Now merge all pdbs into one
        build = halide_bound_to_sp2_carbon_pdb.merge_with_another_molecule(intermediate)
        
        return build
    
    def __alcohol_to_azide_same_stereochem(self, pdb1, reaction, cyanide_or_azide): # Any alcohol connected to an SP3-hybridized carbon, primary, secondary, or tertiary, can become an azide with no change in stereochemistry
        """Simulates the conversion of an alcohol to an azide or cyanide, maintaining the same stereochemistry.
        
        Arguments:
        pdb1 -- The molecular model (pymolecule.Molecule), containing an alcohol.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        cyanide_or_azide -- Either azide or cyanide products are possible. This int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        if cyanide_or_azide == -1 : cyanide_or_azide = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
    
        alcohol_pdb = pdb1
        alcohol_indecies = reaction[1]
    
        # get easy names of the alcohol, so we don't have to refer to the reaction[] list
        alcohol_hydroxyl_hydrogen_index = alcohol_indecies[5]
        alcohol_hydroxyl_oxygen_index = alcohol_indecies[4]
        alcohol_proximal_carbon_index = alcohol_indecies[3]
    
        # Let's just delete the hydrogen early on.
        alcohol_pdb.delete_atom(alcohol_hydroxyl_hydrogen_index)
    
        # Now, load in the intermediate. This is as a model azide
        intermediate = pymolecule.Molecule()
        if cyanide_or_azide == 0:
            intermediate.load_pdb('./intermediates/intermediate9.pdb')
        elif cyanide_or_azide == 1:
            intermediate.load_pdb('./intermediates/cyanide.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
    
        intermediate_carbon_index = 1
        intermediate_first_nitrogen_index = 2
    
        # Now, posiition the alcohol onto the model amide
        tethers = [[alcohol_proximal_carbon_index, intermediate_carbon_index], [alcohol_hydroxyl_oxygen_index, intermediate_first_nitrogen_index]]
        intermediate = alcohol_pdb.align_another_molecule_to_this_one(intermediate, tethers)
        
        # Now delete a few more atoms.
        intermediate.delete_atom(intermediate_carbon_index)
        alcohol_pdb.delete_atom(alcohol_hydroxyl_oxygen_index)
    
        # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
        # Note that the alcohol is fixed and cannot rotate. The two branches of the alcohol can, though.
        item1 = [intermediate, intermediate.all_atoms[intermediate_first_nitrogen_index].coordinates, alcohol_pdb.all_atoms[alcohol_proximal_carbon_index].coordinates] # These won't rotate, so the coordinates are omitted
        item2 = [alcohol_pdb]
    
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        alcohol_pdb.change_residue('FR1')
        intermediate.change_residue('FR2')
        
        # Now merge all pdbs into one
        build = alcohol_pdb.merge_with_another_molecule(intermediate)
        
        return build
    
    def __alcohol_bound_to_sp2_carbon_to_azide(self, pdb1, reaction, cyanide_or_azide): # Any alcohol connected to an SP3-hybridized carbon, primary, secondary, or tertiary, can become an azide with no change in stereochemistry
        """Simulates the conversion of an alcohol bound to an sp2-hybridized carbon, to an azide or cyanide.
        
        Arguments:
        pdb1 -- The molecular model (pymolecule.Molecule), containing an alcohol bond to an sp2-hybridized carbon.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        cyanide_or_azide -- Either azide or cyanide products are possible. This int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        if cyanide_or_azide == -1 : cyanide_or_azide = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
    
        alcohol_pdb = pdb1
        alcohol_indecies = reaction[1]
    
        # get easy names of the alcohol, so we don't have to refer to the reaction[] list
        alcohol_hydroxyl_hydrogen_index = alcohol_indecies[4]
        alcohol_hydroxyl_oxygen_index = alcohol_indecies[3]
        alcohol_proximal_carbon_index = alcohol_indecies[2]
    
        # Let's just delete the hydrogen early on.
        alcohol_pdb.delete_atom(alcohol_hydroxyl_hydrogen_index)
    
        # Now, load in the intermediate. This is as a model azide
        intermediate = pymolecule.Molecule()
        if cyanide_or_azide == 0:
            intermediate.load_pdb('./intermediates/intermediate9.pdb')
        elif cyanide_or_azide == 1:
            intermediate.load_pdb('./intermediates/cyanide.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
    
        intermediate_carbon_index = 1
        intermediate_first_nitrogen_index = 2
    
        # Now, posiition the alcohol onto the model amide
        tethers = [[alcohol_proximal_carbon_index, intermediate_carbon_index], [alcohol_hydroxyl_oxygen_index, intermediate_first_nitrogen_index]]
        intermediate = alcohol_pdb.align_another_molecule_to_this_one(intermediate, tethers)
        
        # Now delete a few more atoms.
        intermediate.delete_atom(intermediate_carbon_index)
        alcohol_pdb.delete_atom(alcohol_hydroxyl_oxygen_index)
    
        # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
        # Note that the alcohol is fixed and cannot rotate. The two branches of the alcohol can, though.
        item1 = [intermediate, intermediate.all_atoms[intermediate_first_nitrogen_index].coordinates, alcohol_pdb.all_atoms[alcohol_proximal_carbon_index].coordinates] # These won't rotate, so the coordinates are omitted
        item2 = [alcohol_pdb]
    
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        alcohol_pdb.change_residue('FR1')
        intermediate.change_residue('FR2')
    
        # Now merge all pdbs into one
        build = alcohol_pdb.merge_with_another_molecule(intermediate)
        
        return build
    
    def __invertable_secondary_alcohol_to_azide(self, pdb1, reaction, cyanide_or_azide): # Secondary alcohols can become azides with or without inversion of stereochemistry. This function handles the case when there is inversion of stereochemistry.
        """Simulates the conversion of a secondary alcohol to an azide or cyanide, with stereochemistry inverted.
        
        Arguments:
        pdb1 -- The molecular model (pymolecule.Molecule), containing a secondary alcohol.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        cyanide_or_azide -- Either azide or cyanide products are possible. This int dictates which product should be made.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        if cyanide_or_azide == -1 : cyanide_or_azide = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
    
        invertable_secondary_alcohol_pdb = pdb1
        invertable_secondary_alcohol_indecies = reaction[1]
    
        # get easy names of the invertable_secondary_alcohol, so we don't have to refer to the reaction[] list
        invertable_hydroxyl_hydrogen_index = invertable_secondary_alcohol_indecies[5]
        invertable_hydroxyl_oxygen_index = invertable_secondary_alcohol_indecies[4]
        invertable_secondary_alcohol_proximal_carbon_index = invertable_secondary_alcohol_indecies[3]
        
        # let's just delete the hydroxyl hydrogen right off the bat
        invertable_secondary_alcohol_pdb.delete_atom(invertable_hydroxyl_hydrogen_index)
        
        for t in range(0,3):
            if pdb1.all_atoms[invertable_secondary_alcohol_indecies[t]].element == "H":
                invertable_secondary_alcohol_hydrogen_index = invertable_secondary_alcohol_indecies[t]
    
        # Now, load in the intermediate. This is as a model azide
        intermediate = pymolecule.Molecule()
        if cyanide_or_azide == 0:
            intermediate.load_pdb('./intermediates/intermediate9.pdb')
        elif cyanide_or_azide == 1:
            intermediate.load_pdb('./intermediates/cyanide.pdb')
    
        # Now, let's define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_carbon_index = 1
        intermediate_first_nitrogen_index = 2
    
        # Now, posiition the invertable_secondary_alcohol onto the model amide
        tethers = [[invertable_secondary_alcohol_proximal_carbon_index, intermediate_carbon_index], [invertable_secondary_alcohol_hydrogen_index, intermediate_first_nitrogen_index]]
        intermediate = invertable_secondary_alcohol_pdb.align_another_molecule_to_this_one(intermediate, tethers)
        
        # Now get the hydrogen fragment
        hydrogen_frag = invertable_secondary_alcohol_pdb.get_branch(invertable_secondary_alcohol_proximal_carbon_index, invertable_secondary_alcohol_hydrogen_index)
        
        # Now position the hydrogen frag
        tethers = [[invertable_secondary_alcohol_proximal_carbon_index, invertable_secondary_alcohol_proximal_carbon_index], [invertable_hydroxyl_oxygen_index, invertable_secondary_alcohol_hydrogen_index]]
        hydrogen_frag = invertable_secondary_alcohol_pdb.align_another_molecule_to_this_one(hydrogen_frag, tethers)
    
        # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
        # Note that the invertable_secondary_alcohol is fixed and cannot rotate. The two branches of the alcohol can, though.
        item1 = [intermediate, intermediate.all_atoms[intermediate_first_nitrogen_index].coordinates, intermediate.all_atoms[intermediate_carbon_index].coordinates] # These won't rotate, so the coordinates are omitted
        item2 = [invertable_secondary_alcohol_pdb]
        item3 = [hydrogen_frag]
    
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        thelist.append(item3)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        # Now delete a few more atoms.
        intermediate.delete_atom(intermediate_carbon_index)
        hydrogen_frag.delete_atom(invertable_secondary_alcohol_proximal_carbon_index)
        invertable_secondary_alcohol_pdb.delete_atom(invertable_hydroxyl_oxygen_index)
        invertable_secondary_alcohol_pdb.delete_atom(invertable_secondary_alcohol_hydrogen_index)
        
        invertable_secondary_alcohol_pdb.change_residue('FR1')
        intermediate.change_residue('FR2')
        hydrogen_frag.change_residue('FR3')
        
        # Now merge all pdbs into one
        build = invertable_secondary_alcohol_pdb.merge_with_another_molecule(intermediate)
        build = build.merge_with_another_molecule(hydrogen_frag)
        
        return build
    
    def __carboxylate_cyanide(self, pdb1, reaction): # Note only one pdb as input
        """Simulates the conversion of a carboxylate to a cyanide.
        
        Arguments:
        pdb1 -- The molecular model (pymolecule.Molecule), containing a carboxylate.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        # pdb1 is the carboxyalate
        carboxylate_pdb = pdb1
        carboxylate_indecies = reaction[1]
    
        # get easy names of the carboxylate
        carboxylate_carbonyl_carbon = carboxylate_indecies[0]
        carboxylate_carbonyl_oxygen = carboxylate_indecies[1]
        carboxylate_oxygen = carboxylate_indecies[2]
        carboxylate_oxygen_hydrogen = carboxylate_indecies[3]
    
        # load in the intermediate, cyanide is thought to be in the solution
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb('./intermediates/cyanide.pdb')
    
        # define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_carbon_index = 1
        intermediate_carbon_carbon_index = 2
        intermediate_carbon_carbon_nitrogen_index = 3
    
        # position the carboxylate onto the cyanide model
        tethers = [[intermediate_carbon_index, carboxylate_carbonyl_carbon], [intermediate_carbon_carbon_index, carboxylate_oxygen]]
        carboxylate_pdb = intermediate.align_another_molecule_to_this_one(carboxylate_pdb, tethers)
    
        # delete a few extra atoms
        carboxylate_pdb.delete_atom(carboxylate_oxygen)
        if carboxylate_oxygen_hydrogen!=-1: carboxylate_pdb.delete_atom(carboxylate_oxygen_hydrogen)
    
        # rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        item1 = [intermediate, intermediate.all_atoms[intermediate_carbon_index].coordinates, intermediate.all_atoms[intermediate_carbon_carbon_index].coordinates]
        item2 = [carboxylate_pdb] # These won't rotate, so the coordinates are omitted
    
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        carboxylate_pdb.delete_atom(carboxylate_carbonyl_carbon)
    
        carboxylate_pdb.change_residue('FR1')
        intermediate.change_residue('FR2')
    
        # merge all pdbs into one
        build = carboxylate_pdb.merge_with_another_molecule(intermediate)
    
        return build
    
    def __carboxylate_azide(self, pdb1, reaction): # Note only one pdb as input
        """Simulates the conversion of a carboxylate to an azide.
        
        Arguments:
        pdb1 -- The molecular model (pymolecule.Molecule), containing a carboxylate.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        # pdb1 is the carboxyalate
        carboxylate_pdb = pdb1
        carboxylate_indecies = reaction[1]
    
        # get easy names of the carboxylate
        carboxylate_carbonyl_carbon = carboxylate_indecies[0]
        carboxylate_carbonyl_oxygen = carboxylate_indecies[1]
        carboxylate_oxygen = carboxylate_indecies[2]
        carboxylate_oxygen_hydrogen = carboxylate_indecies[3]
    
        # load in the intermediate, azide is thought to be in the solution
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb('./intermediates/intermediate9.pdb')
    
        # define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_carbon_index = 1
        intermediate_carbon_nitrogen_index = 2
        intermediate_carbon_nitrogen_nitrogen_index = 3
        intermediate_carbon_nitrogen_nitrogen_nitrogen_index = 4
    
        # position the carboxylate onto the azide model
        tethers = [[intermediate_carbon_index, carboxylate_carbonyl_carbon], [intermediate_carbon_nitrogen_index, carboxylate_oxygen]]
        carboxylate_pdb = intermediate.align_another_molecule_to_this_one(carboxylate_pdb, tethers)
    
        # delete a few extra atoms
        carboxylate_pdb.delete_atom(carboxylate_oxygen)
        if carboxylate_oxygen_hydrogen!=-1: carboxylate_pdb.delete_atom(carboxylate_oxygen_hydrogen)
    
        # rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        item1 = [intermediate, intermediate.all_atoms[intermediate_carbon_index].coordinates, intermediate.all_atoms[intermediate_carbon_nitrogen_index].coordinates]
        item2 = [carboxylate_pdb] # These won't rotate, so the coordinates are omitted
    
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
        
        carboxylate_pdb.delete_atom(carboxylate_carbonyl_carbon)
    
        carboxylate_pdb.change_residue('FR1')
        intermediate.change_residue('FR2')
    
        # merge all pdbs into one
        build = carboxylate_pdb.merge_with_another_molecule(intermediate)
    
        return build
    
    def __acylhalide_azide(self, pdb1, reaction): # Note only one pdb as input
        """Simulates the conversion of an acyl halide to an azide.
        
        Arguments:
        pdb1 -- The molecular model (pymolecule.Molecule), containing an acyl halide.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        # pdb1 is the carboxyalate_Cl
        acylhalide_pdb = pdb1
        acylhalide_indecies = reaction[1]
                
        # get easy names of the carboxylate
        acylhalide_carbonyl_carbon = acylhalide_indecies[0]
        acylhalide_carbonyl_oxygen = acylhalide_indecies[1]
        acylhalide_chloride = acylhalide_indecies[2]
    
        # load in the intermediate, azide is thought to be in the solution
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb('./intermediates/intermediate9.pdb')
    
        # define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_carbon_index = 1
        intermediate_carbon_nitrogen_index = 2 
        intermediate_carbon_nitrogen_nitrogen_index = 3
        intermediate_carbon_nitrogen_nitrogen_nitrogen_index = 4
    
        # position the carboxylate onto the azide model
        tethers = [[intermediate_carbon_index, acylhalide_carbonyl_carbon], [intermediate_carbon_nitrogen_index, acylhalide_chloride]]
        acylhalide_pdb = intermediate.align_another_molecule_to_this_one(acylhalide_pdb, tethers)
    
        # delete a few extra atoms
        acylhalide_pdb.delete_atom(acylhalide_chloride)
        
        # rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        item1 = [intermediate, intermediate.all_atoms[intermediate_carbon_index].coordinates, intermediate.all_atoms[intermediate_carbon_nitrogen_index].coordinates]
        item2 = [acylhalide_pdb] # These won't rotate, so the coordinates are omitted
        
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
        
        acylhalide_pdb.delete_atom(acylhalide_carbonyl_carbon)
        
        acylhalide_pdb.change_residue('FR1')
        intermediate.change_residue('FR2')
        
        # merge all pdbs into one
        build = acylhalide_pdb.merge_with_another_molecule(intermediate)
    
        return build
    
    def __acylhalide_cyanide(self, pdb1, reaction): # Note only one pdb as input
        """Simulates the conversion of an acyl halide to a cyanide.
        
        Arguments:
        pdb1 -- The molecular model (pymolecule.Molecule), containing an acyl halide.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        # pdb1 is the carboxyalate
        acylhalide_pdb = pdb1
        acylhalide_indecies = reaction[1]
    
        # get easy names of the acylhalide
        acylhalide_carbonyl_carbon = acylhalide_indecies[0]
        acylhalide_carbonyl_oxygen = acylhalide_indecies[1]
        acylhalide_chloride = acylhalide_indecies[2]
    
        # load in the intermediate, cyanide is thought to be in the solution
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb('./intermediates/cyanide.pdb')
    
        # define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_carbon_index = 1
        intermediate_carbon_carbon_index = 2
        intermediate_carbon_carbon_nitrogen_index = 3
    
        # position the acylhalide onto the cyanide model
        tethers = [[intermediate_carbon_index, acylhalide_carbonyl_carbon], [intermediate_carbon_carbon_index, acylhalide_chloride]]
        acylhalide_pdb = intermediate.align_another_molecule_to_this_one(acylhalide_pdb, tethers)
    
        # delete a few extra atoms
        acylhalide_pdb.delete_atom(acylhalide_chloride)
    
        # rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        item1 = [intermediate, intermediate.all_atoms[intermediate_carbon_index].coordinates, intermediate.all_atoms[intermediate_carbon_carbon_index].coordinates]
        item2 = [acylhalide_pdb] # These won't rotate, so the coordinates are omitted
    
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        acylhalide_pdb.delete_atom(acylhalide_carbonyl_carbon)
    
        acylhalide_pdb.change_residue('FR1')
        intermediate.change_residue('FR2')
    
        # merge all pdbs into one
        build = acylhalide_pdb.merge_with_another_molecule(intermediate)
    
        return build
    
    def __acid_anhydride_azide(self, pdb1, reaction, side_to_keep_choice): # Note that pdb2 can be a primary or secondary azide
        """Simulates the conversion of an acid anhydride to an azide.
        
        Arguments:
        pdb1 -- The molecular model (pymolecule.Molecule), containing an acid anhydride.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        side_to_keep_choice -- The acid anhydride will be split. One half of the molecule will be discarded, and the other half will be incorporated into the product. This int determines which half to keep.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
        
        if side_to_keep_choice == -1 : side_to_keep_choice = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
    
        acid_anhydride_pdb = pdb1
        acid_anhydride_indecies = reaction[1]
            
        # get easy names of the acid_anhydride, so we don't have to refer to the reaction[] list
        acid_anhydride_ether_oxgen = acid_anhydride_indecies[0]
        acid_anhydride_carybonyl_carbon_1 = acid_anhydride_indecies[1]
        acid_anhydride_carybonyl_oxygen_1 = acid_anhydride_indecies[2]
        acid_anhydride_carybonyl_carbon_2 = acid_anhydride_indecies[3]
        acid_anhydride_carybonyl_oxygen_2 = acid_anhydride_indecies[4]
    
        # load in the intermediate, azide is thought to be in the solution
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb('./intermediates/intermediate9.pdb')
    
        # define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_carbon_index = 1
        intermediate_carbon_nitrogen_index = 2 
        intermediate_carbon_nitrogen_nitrogen_index = 3
        intermediate_carbon_nitrogen_nitrogen_nitrogen_index = 4
    
        # Now, we need to pick which side of the acid anhydride you'll keep
        if side_to_keep_choice == 0 :
            acid_anhydride_pdb = acid_anhydride_pdb.get_branch(acid_anhydride_ether_oxgen, acid_anhydride_carybonyl_carbon_1)
            acid_anhydride_carybonyl_carbon = acid_anhydride_carybonyl_carbon_1
            acid_anhydride_carybonyl_oxygen = acid_anhydride_carybonyl_oxygen_1
    
        else :
            acid_anhydride_pdb = acid_anhydride_pdb.get_branch(acid_anhydride_ether_oxgen, acid_anhydride_carybonyl_carbon_2)
            acid_anhydride_carybonyl_carbon = acid_anhydride_carybonyl_carbon_2
            acid_anhydride_carybonyl_oxygen = acid_anhydride_carybonyl_oxygen_2
            
        # Now, posiition the acid_anhydride onto the model amide
        tethers = [[intermediate_carbon_index, acid_anhydride_carybonyl_carbon], [intermediate_carbon_nitrogen_index, acid_anhydride_ether_oxgen]]
        acid_anhydride_pdb = intermediate.align_another_molecule_to_this_one(acid_anhydride_pdb, tethers)
    
        # Now delete a few more atoms.
        acid_anhydride_pdb.delete_atom(acid_anhydride_ether_oxgen)
        acid_anhydride_pdb.delete_atom(acid_anhydride_carybonyl_carbon)
    
        # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
        # Note that the acid_anhydride is fixed and cannot rotate. The two branches of the azide can, though.
        item1 = [intermediate, intermediate.all_atoms[intermediate_carbon_index].coordinates, intermediate.all_atoms[intermediate_carbon_nitrogen_index].coordinates]
        item2 = [acid_anhydride_pdb] # These won't rotate, so the coordinates are omitted
        
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        acid_anhydride_pdb.change_residue('FR1')
        intermediate.change_residue('FR2')
    
        # Now merge all pdbs into one
        build = acid_anhydride_pdb.merge_with_another_molecule(intermediate)
        
        return build    
    
    def __acid_anhydride_cyanide(self, pdb1, reaction, side_to_keep_choice):
        """Simulates the conversion of an acid anhydride to a cyanide.
        
        Arguments:
        pdb1 -- The molecular model (pymolecule.Molecule), containing an acid anhydride.
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        side_to_keep_choice -- The acid anhydride will be split. One half of the molecule will be discarded, and the other half will be incorporated into the product. This int determines which half to keep.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
    
        if side_to_keep_choice == -1 : side_to_keep_choice = random.randrange(0,2) # so if the user doesn't specify a value, it randomly picks between 0 and 1
    
        acid_anhydride_pdb = pdb1
        acid_anhydride_indecies = reaction[1]
            
        # get easy names of the acid_anhydride, so we don't have to refer to the reaction[] list
        acid_anhydride_ether_oxgen = acid_anhydride_indecies[0]
        acid_anhydride_carybonyl_carbon_1 = acid_anhydride_indecies[1]
        acid_anhydride_carybonyl_oxygen_1 = acid_anhydride_indecies[2]
        acid_anhydride_carybonyl_carbon_2 = acid_anhydride_indecies[3]
        acid_anhydride_carybonyl_oxygen_2 = acid_anhydride_indecies[4]
    
        # load in the intermediate, cyanide is thought to be in the solution
        intermediate = pymolecule.Molecule()
        intermediate.load_pdb('./intermediates/cyanide.pdb')
    
        # define variable names for the indices of the intermediate (the model), so we don't have to remember the numbers
        intermediate_carbon_index = 1
        intermediate_carbon_carbon_index = 2
        intermediate_carbon_carbon_nitrogen_index = 3
    
        # Now, we need to pick which side of the acid anhydride you'll keep
        if side_to_keep_choice == 0 :
            acid_anhydride_pdb = acid_anhydride_pdb.get_branch(acid_anhydride_ether_oxgen, acid_anhydride_carybonyl_carbon_1)
            acid_anhydride_carybonyl_carbon = acid_anhydride_carybonyl_carbon_1
            acid_anhydride_carybonyl_oxygen = acid_anhydride_carybonyl_oxygen_1
    
        else :
            acid_anhydride_pdb = acid_anhydride_pdb.get_branch(acid_anhydride_ether_oxgen, acid_anhydride_carybonyl_carbon_2)
            acid_anhydride_carybonyl_carbon = acid_anhydride_carybonyl_carbon_2
            acid_anhydride_carybonyl_oxygen = acid_anhydride_carybonyl_oxygen_2
            
        
        # Now, posiition the acid_anhydride onto the model amide
        tethers = [[intermediate_carbon_index, acid_anhydride_carybonyl_carbon], [intermediate_carbon_carbon_index, acid_anhydride_ether_oxgen]]
        acid_anhydride_pdb = intermediate.align_another_molecule_to_this_one(acid_anhydride_pdb, tethers)
    
        # Now delete a few more atoms.
        acid_anhydride_pdb.delete_atom(acid_anhydride_ether_oxgen)
        acid_anhydride_pdb.delete_atom(acid_anhydride_carybonyl_carbon)
    
        # Now, we need to rotate all these fragments coming off the intermediate so as to minimize steric hindrance
        # We define a series of items (triplets), containing a given pdb and two coordinates forming a line about which the pdb can be rotated
        # Note that the acid_anhydride is fixed and cannot rotate. The two branches of the cyanide can, though.
        item1 = [intermediate, intermediate.all_atoms[intermediate_carbon_index].coordinates, intermediate.all_atoms[intermediate_carbon_carbon_index].coordinates]
        item2 = [acid_anhydride_pdb] # These won't rotate, so the coordinates are omitted
        
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        acid_anhydride_pdb.change_residue('FR1')
        intermediate.change_residue('FR2')
    
        # Now merge all pdbs into one
        build = acid_anhydride_pdb.merge_with_another_molecule(intermediate)
        
        return build    
    
    def __autogrow(self, pdb1, pdb2, reaction):
        """Replicates the fragment-addition protocol of AutoGrow 1.0, replacing hydrogen atoms with molecular fragments.
        
        Arguments:
        pdb1 -- A molecular model (pymolecule.Molecule).
        pdb2 -- A molecular model (pymolecule.Molecule).
        reaction -- A list identifying the reactive chemical groups and which atoms will participate in the reaction.
        
        Returns:
        A pymolecule.Molecule model of the product.
        
        """
    
        # This simulates autogrow, linking fragments by hydrogens.
    
        # You need to figure out which pdb is the carbonochloridate and which is the primary alcohol
        linker_index_1 = reaction[1][0]
        linker_1 = pdb1.all_atoms[linker_index_1]
        hydrogen_index_1 = reaction[1][1]
        hydrogen_1 = pdb1.all_atoms[hydrogen_index_1]
    
        linker_index_2 = reaction[3][0]
        linker_2 = pdb2.all_atoms[linker_index_2]
        hydrogen_index_2 = reaction[3][1]
        hydrogen_2 = pdb2.all_atoms[hydrogen_index_2]
    
        # Now figure out where pdb2 needs to be moved
        # Get the vector between the first X-H.
        vec = hydrogen_1.coordinates.subtract(linker_1.coordinates)
        vec = vec.normalized()
        vec = vec.scalar_multiply(pdb1.bond_length(linker_1.element, linker_2.element)) # all Molecule objects have bond_length function, so pdb1 is just as good as any.
        new_coor = linker_1.coordinates.add(vec)
        hydrogen_1.coordinates = new_coor
    
        # Now rotatione pdb2 to align it with the bond to be formed
        tethers = [[hydrogen_index_1, linker_index_2], [linker_index_1, hydrogen_index_2]]
        pdb2 = pdb1.align_another_molecule_to_this_one(pdb2, tethers)
    
        # Now delete the hydrogens
        pdb1.delete_atom(hydrogen_index_1)
        pdb2.delete_atom(hydrogen_index_2)
        
        # Now rotate the fragments to try to minimize steric hindrance
        item1 = [pdb2, linker_1.coordinates, linker_2.coordinates]
        item2 = [pdb1]
    
        thelist = []
        thelist.append(item1)
        thelist.append(item2)
        self.structure_pdb_functions.reduce_steric_hindrance(thelist)
    
        pdb1.change_residue('FR1')
        pdb2.change_residue('FR2')
    
        # Now merge these two pdbs
        build = pdb1.merge_with_another_molecule(pdb2)
    
        return build
    
class AutoClickChem:
    """The main AutoClickChem class."""
    
    version = '1.0.0'
    operators_react = OperatorsReact()

    def __init__(self):
        """Initialize an AutoClickChem object."""
        
        # always check to see if intermediate directory exists. If not, make it and populate it with files.
        if not os.path.exists("./intermediates"):
            os.mkdir("./intermediates")
            self.__make_intermediates()

    def __make_intermediates(self):
        """Generates PDB files containing models used to construct click-chemistry products."""
        
        f = open('./intermediates/amide.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -0.215   3.275  -0.011  1.00  0.00           C" + "\n")
        f.write("HETATM    2  O2  UNK     0      -0.002   4.458  -0.301  1.00  0.00           O" + "\n")
        f.write("HETATM    3  C3  UNK     0       0.940   2.360   0.292  1.00  0.00           C" + "\n")
        f.write("HETATM    4  N4  UNK     0      -1.532   2.795   0.033  1.00  0.00           N" + "\n")
        f.write("HETATM    5  C5  UNK     0      -1.796   1.422   0.368  1.00  0.00           C" + "\n")
        f.write("HETATM    6  C6  UNK     0      -2.623   3.680  -0.259  1.00  0.00           C" + "\n")
        f.close()
        f = open('./intermediates/carbamate_2_cis.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -6.353   3.159   0.167  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2  UNK     0      -5.179   2.357  -0.020  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3  UNK     0      -3.916   2.952   0.004  1.00  0.00           C" + "\n")
        f.write("HETATM    4  O4  UNK     0      -3.801   4.291   0.206  1.00  0.00           O" + "\n")
        f.write("HETATM    5  O5  UNK     0      -2.904   2.261  -0.157  1.00  0.00           O" + "\n")
        f.write("HETATM    6  H4  UNK     0      -5.264   1.377  -0.168  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/carbamate_2_trans.pdb','w')
        f.write("HETATM    1  C1          0      -6.358   3.168  -0.001  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2          0      -5.179   2.351  -0.000  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3          0      -3.917   2.947  -0.002  1.00  0.00           C" + "\n")
        f.write("HETATM    4  O4          0      -2.774   2.154  -0.001  1.00  0.00           O" + "\n")
        f.write("HETATM    5  O5          0      -3.815   4.179  -0.003  1.00  0.00           O" + "\n")
        f.write("HETATM    6  H4          0      -5.261   1.359   0.001  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/carbamate.pdb','w')
        f.write("HETATM    1  C1          0      -0.223   3.281  -0.013  1.00  0.00           C" + "\n")
        f.write("HETATM    2  O2          0      -0.008   4.463  -0.302  1.00  0.00           O" + "\n")
        f.write("HETATM    3  O3          0       0.839   2.428   0.268  1.00  0.00           O" + "\n")
        f.write("HETATM    4  N4          0      -1.534   2.801   0.031  1.00  0.00           N" + "\n")
        f.write("HETATM    5  C5          0      -1.784   1.429   0.367  1.00  0.00           C" + "\n")
        f.write("HETATM    6  C6          0      -2.630   3.680  -0.259  1.00  0.00           C" + "\n")
        f.close()
        f = open('./intermediates/carbamodithioate_cis.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -6.338   3.150   0.315  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2  UNK     0      -5.171   2.395  -0.037  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3  UNK     0      -3.911   2.996   0.005  1.00  0.00           C" + "\n")
        f.write("HETATM    4  S4  UNK     0      -3.759   4.641   0.480  1.00  0.00           S" + "\n")
        f.write("HETATM    5  S5  UNK     0      -2.638   2.171  -0.380  1.00  0.00           S" + "\n")
        f.write("HETATM    6  H4  UNK     0      -5.258   1.443  -0.312  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/carbamodithioate_trans.pdb','w')
        f.write("HETATM    1  C1          0      -6.358   3.187  -0.000  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2          0      -5.170   2.381  -0.000  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3          0      -3.913   2.989  -0.002  1.00  0.00           C" + "\n")
        f.write("HETATM    4  S4          0      -2.449   1.993  -0.002  1.00  0.00           S" + "\n")
        f.write("HETATM    5  S5          0      -3.797   4.550  -0.004  1.00  0.00           S" + "\n")
        f.write("HETATM    6  H4          0      -5.243   1.389   0.001  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/carbamothioate_1_cis.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -6.354   3.179   0.169  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2  UNK     0      -5.172   2.389  -0.020  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3  UNK     0      -3.916   2.998  -0.008  1.00  0.00           C" + "\n")
        f.write("HETATM    4  O4  UNK     0      -3.814   4.339   0.183  1.00  0.00           O" + "\n")
        f.write("HETATM    5  S5  UNK     0      -2.627   2.137  -0.214  1.00  0.00           S" + "\n")
        f.write("HETATM    6  H4  UNK     0      -5.247   1.407  -0.160  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/carbamothioate_1_trans.pdb','w')
        f.write("HETATM    1  C1          0      -6.360   3.188  -0.000  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2          0      -5.172   2.382  -0.000  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3          0      -3.915   2.991  -0.002  1.00  0.00           C" + "\n")
        f.write("HETATM    4  O4          0      -2.765   2.209  -0.002  1.00  0.00           O" + "\n")
        f.write("HETATM    5  S5          0      -3.801   4.552  -0.004  1.00  0.00           S" + "\n")
        f.write("HETATM    6  H4          0      -5.244   1.390   0.001  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/carbamothioate_2_cis.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -6.351   3.201   0.200  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2  UNK     0      -5.163   2.429  -0.024  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3  UNK     0      -3.911   3.045   0.019  1.00  0.00           C" + "\n")
        f.write("HETATM    4  S4  UNK     0      -3.792   4.729   0.341  1.00  0.00           S" + "\n")
        f.write("HETATM    5  O5  UNK     0      -2.888   2.378  -0.175  1.00  0.00           O" + "\n")
        f.write("HETATM    6  H4  UNK     0      -5.232   1.454  -0.211  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/carbamothioate_2_trans.pdb','w')
        f.write("HETATM    1  C1          0      -6.359   3.212  -0.000  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2          0      -5.164   2.418  -0.000  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3          0      -3.913   3.039  -0.002  1.00  0.00           C" + "\n")
        f.write("HETATM    4  S4          0      -2.438   2.058  -0.002  1.00  0.00           S" + "\n")
        f.write("HETATM    5  O5          0      -3.835   4.272  -0.004  1.00  0.00           O" + "\n")
        f.write("HETATM    6  H4          0      -5.226   1.425   0.001  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/carboxylic_acid.pdb','w')
        f.write("HETATM    1  C1          0      -2.797   0.397  -0.101  1.00  0.00           C" + "\n")
        f.write("HETATM    2  O2          0      -2.340   1.197   0.724  1.00  0.00           O" + "\n")
        f.write("HETATM    3  O3          0      -2.066   0.083  -1.244  1.00  0.00           O" + "\n")
        f.write("HETATM    4  C4          0      -4.140  -0.235   0.132  1.00  0.00           C" + "\n")
        f.close()
        f = open('./intermediates/charged_amine.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -5.340   2.851   0.001  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2  UNK     0      -3.870   2.852   0.000  1.00  0.00           N" + "\n")
        f.write("HETATM    3  H1  UNK     0      -3.527   2.367  -0.841  1.00  0.00           H" + "\n")
        f.write("HETATM    4  H2  UNK     0      -3.527   3.823   0.001  1.00  0.00           H" + "\n")
        f.write("HETATM    5  H1  UNK     0      -3.526   2.366   0.841  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/cyanide.pdb','w')
        f.write("HETATM    1  C1          0       1.243   1.263   0.364  1.00  0.00           C" + "\n")
        f.write("HETATM    2  C2          0       0.014   1.575   1.094  1.00  0.00           C" + "\n")
        f.write("HETATM    3  N3          0      -0.953   1.821   1.668  1.00  0.00           N" + "\n")
        f.close()
        f = open('./intermediates/generic_epoxide.pdb','w')
        f.write("HETATM    1  C1          0      -2.998   3.356   0.117  1.00  0.00           C" + "\n")
        f.write("HETATM    2  O2          0      -2.352   3.636  -1.158  1.00  0.00           O" + "\n")
        f.write("HETATM    3  C3          0      -1.505   3.454   0.013  1.00  0.00           C" + "\n")
        f.write("HETATM    4  C1          0      -3.598   1.941   0.244  1.00  0.00           C" + "\n")
        f.write("HETATM    5  C5          0      -3.735   4.555   0.749  1.00  0.00           C" + "\n")
        f.write("HETATM    6  C6          0      -0.850   4.745   0.547  1.00  0.00           C" + "\n")
        f.write("HETATM    7  C7          0      -0.714   2.130   0.042  1.00  0.00           C" + "\n")
        f.close()
        f = open('./intermediates/intermediate10.pdb','w')
        f.write("HETATM    1  C1          0      -5.348   2.852   0.001  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2          0      -3.878   2.852  -0.000  1.00  0.00           N" + "\n")
        f.write("HETATM    3  H1          0      -3.534   3.337   0.841  1.00  0.00           H" + "\n")
        f.write("HETATM    4  H2          0      -3.535   3.337  -0.841  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/intermediate1.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -8.108   6.038  -0.000  1.00  0.00           C" + "\n")
        f.write("HETATM    2  C2  UNK     0      -6.995   5.337   0.001  1.00  0.00           C" + "\n")
        f.write("HETATM    3  N3  UNK     0      -7.382   4.039   0.001  1.00  0.00           N" + "\n")
        f.write("HETATM    4  N4  UNK     0      -8.737   3.987   0.001  1.00  0.00           N" + "\n")
        f.write("HETATM    5  N5  UNK     0      -9.186   5.132   0.000  1.00  0.00           N" + "\n")
        f.write("HETATM    6  C6  UNK     0      -6.504   2.903   0.002  1.00  0.00           C" + "\n")
        f.write("HETATM    7  C7  UNK     0      -5.589   5.867   0.001  1.00  0.00           C" + "\n")
        f.write("HETATM    8  C8  UNK     0      -8.212   7.536  -0.001  1.00  0.00           C" + "\n")
        f.close()
        f = open('./intermediates/intermediate2_isomer.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -6.186   7.659   0.926  1.00  0.00           C" + "\n")
        f.write("HETATM    2  O2  UNK     0      -6.507   6.465   0.931  1.00  0.00           O" + "\n")
        f.write("HETATM    3  N4  UNK     0      -4.839   8.020   0.855  1.00  0.00           N" + "\n")
        f.write("HETATM    4  S4  UNK     0      -4.449   9.473   0.849  1.00  0.00           S" + "\n")
        f.write("HETATM    5  H4  UNK     0      -4.138   7.317   0.808  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/intermediate2.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -6.207   7.596   0.961  1.00  0.00           C" + "\n")
        f.write("HETATM    2  O2  UNK     0      -7.136   8.411   1.025  1.00  0.00           O" + "\n")
        f.write("HETATM    3  N4  UNK     0      -4.845   7.993   0.863  1.00  0.00           N" + "\n")
        f.write("HETATM    4  S4  UNK     0      -4.361   9.643   0.825  1.00  0.00           S" + "\n")
        f.write("HETATM    5  H4  UNK     0      -4.143   7.290   0.816  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/intermediate3.pdb','w')
        f.write("HETATM    1  C1          0      -4.349   2.095  -0.336  1.00  0.00           C" + "\n")
        f.write("HETATM    2  C2          0      -2.902   2.051  -0.914  1.00  0.00           C" + "\n")
        f.write("HETATM    3  C3          0      -2.648   3.300  -1.784  1.00  0.00           C" + "\n")
        f.write("HETATM    4  C4          0      -2.722   0.787  -1.780  1.00  0.00           C" + "\n")
        f.write("HETATM    5  C5          0      -1.877   2.023   0.240  1.00  0.00           C" + "\n")
        f.write("HETATM    6  C6          0      -5.374   2.123  -1.489  1.00  0.00           C" + "\n")
        f.write("HETATM    7  C7          0      -4.603   0.846   0.535  1.00  0.00           C" + "\n")
        f.write("HETATM    8  C8          0      -4.529   3.359   0.531  1.00  0.00           C" + "\n")
        f.close()
        f = open('./intermediates/intermediate4.pdb','w')
        f.write("HETATM    1  C1          0      -6.126   3.494   0.001  1.00  0.00           C" + "\n")
        f.write("HETATM    2  O2          0      -5.006   4.384  -0.001  1.00  0.00           O" + "\n")
        f.write("HETATM    3  C3          0      -3.794   3.625  -0.002  1.00  0.00           C" + "\n")
        f.close()
        f = open('./intermediates/intermediate5.pdb','w')
        f.write("HETATM    1  C1          0      -6.114   3.496   0.001  1.00  0.00           C" + "\n")
        f.write("HETATM    2  O2          0      -4.984   4.372  -0.001  1.00  0.00           O" + "\n")
        f.write("HETATM    3  H3          0      -4.151   3.836  -0.002  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/intermediate6.pdb','w')
        f.write("HETATM    1  C1          0      -6.412   3.387   0.002  1.00  0.00           C" + "\n")
        f.write("HETATM    2  S2          0      -4.999   4.519  -0.001  1.00  0.00           S" + "\n")
        f.write("HETATM    3  H3          0      -3.835   3.797  -0.002  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/intermediate7.pdb','w')
        f.write("HETATM    1  C1          0      -6.418   3.389  -0.001  1.00  0.00           C" + "\n")
        f.write("HETATM    2  S2          0      -5.006   4.522   0.000  1.00  0.00           S" + "\n")
        f.write("HETATM    3  C3          0      -3.467   3.570  -0.002  1.00  0.00           C" + "\n")
        f.close()
        f = open('./intermediates/intermediate8.pdb','w')
        f.write("HETATM    1  C1          0      -6.421   3.387   0.000  1.00  0.00           C" + "\n")
        f.write("HETATM    2  H2          0      -5.563   4.075   0.000  1.00  0.00           H" + "\n")
        f.close()
        f = open('./intermediates/intermediate9.pdb','w')
        f.write("HETATM    1  C1          0      -5.522   4.870  -0.000  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2          0      -4.094   5.008   0.000  1.00  0.00           N" + "\n")
        f.write("HETATM    3  N3          0      -3.409   4.048   0.000  1.00  0.00           N" + "\n")
        f.write("HETATM    4  N4          0      -2.725   3.089   0.001  1.00  0.00           N" + "\n")
        f.close()
        f = open('./intermediates/N-C-O.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -5.498   4.861   0.000  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2  UNK     0      -4.073   4.986   0.001  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3  UNK     0      -3.373   3.969   0.001  1.00  0.00           C" + "\n")
        f.write("HETATM    4  O4  UNK     0      -2.693   2.980   0.002  1.00  0.00           O" + "\n")
        f.close()
        f = open('./intermediates/N-C-S.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -5.472   4.848   0.001  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2  UNK     0      -4.042   4.946   0.001  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3  UNK     0      -3.355   3.919   0.002  1.00  0.00           C" + "\n")
        f.write("HETATM    4  S4  UNK     0      -2.504   2.649   0.002  1.00  0.00           S" + "\n")
        f.close()
        f = open('./intermediates/sp2_primary_amine.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -5.290   2.852   0.001  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2  UNK     0      -3.894   2.852   0.000  1.00  0.00           N" + "\n")
        f.write("HETATM    3  H1  UNK     0      -3.398   2.242  -0.609  1.00  0.00           H" + "\n")
        f.write("HETATM    4  H2  UNK     0      -3.397   3.461   0.609  1.00  0.00           H" + "\n")
        f.write("HETATM    5  C1  UNK     0      -6.040   1.933  -0.918  1.00  0.00           C" + "\n")
        f.close()
        f = open('./intermediates/thio_acid.pdb','w')
        f.write("HETATM    1  C1          0      -2.739   0.349   0.055  1.00  0.00           C" + "\n")
        f.write("HETATM    2  O2          0      -2.356   1.183   0.883  1.00  0.00           O" + "\n")
        f.write("HETATM    3  S3          0      -1.644  -0.195  -1.227  1.00  0.00           S" + "\n")
        f.write("HETATM    4  C4          0      -4.133  -0.204   0.137  1.00  0.00           C" + "\n")
        f.close()
        f = open('./intermediates/thiourea_cis.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -5.614   1.049   0.232  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2  UNK     0      -5.242   2.440   0.051  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3  UNK     0      -3.872   2.978   0.012  1.00  0.00           C" + "\n")
        f.write("HETATM    4  N4  UNK     0      -2.588   2.267   0.144  1.00  0.00           N" + "\n")
        f.write("HETATM    5  S5  UNK     0      -3.768   4.525  -0.198  1.00  0.00           S" + "\n")
        f.write("HETATM    6  H4  UNK     0      -5.988   3.094  -0.059  1.00  0.00           H" + "\n")
        f.write("HETATM    7  C7  UNK     0      -1.391   3.062   0.066  1.00  0.00           C" + "\n")
        f.write("HETATM    8  C8  UNK     0      -2.411   0.841   0.345  1.00  0.00           C" + "\n")
        f.close()
        f = open('./intermediates/thiourea_trans.pdb','w')
        f.write("HETATM    1  C1          0      -6.332   3.224  -0.054  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2          0      -5.139   2.423  -0.021  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3          0      -3.878   3.050   0.055  1.00  0.00           C" + "\n")
        f.write("HETATM    4  N4          0      -2.682   2.299   0.091  1.00  0.00           N" + "\n")
        f.write("HETATM    5  S5          0      -3.803   4.612   0.101  1.00  0.00           S" + "\n")
        f.write("HETATM    6  H4          0      -5.213   1.432  -0.051  1.00  0.00           H" + "\n")
        f.write("HETATM    7  C7          0      -1.419   2.982   0.168  1.00  0.00           C" + "\n")
        f.write("HETATM    8  C8          0      -2.706   0.862   0.051  1.00  0.00           C" + "\n")
        f.close()
        f = open('./intermediates/urea_cis.pdb','w')
        f.write("HETATM    1  C1  UNK     0      -5.576   0.982  -0.196  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2  UNK     0      -5.233   2.385  -0.054  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3  UNK     0      -3.876   2.943   0.070  1.00  0.00           C" + "\n")
        f.write("HETATM    4  N4  UNK     0      -2.579   2.244   0.081  1.00  0.00           N" + "\n")
        f.write("HETATM    5  O5  UNK     0      -3.818   4.173   0.182  1.00  0.00           O" + "\n")
        f.write("HETATM    6  H4  UNK     0      -5.992   3.034  -0.039  1.00  0.00           H" + "\n")
        f.write("HETATM    7  C7  UNK     0      -1.402   3.059   0.217  1.00  0.00           C" + "\n")
        f.write("HETATM    8  C8  UNK     0      -2.372   0.813  -0.034  1.00  0.00           C" + "\n")
        f.close()
        f = open('./intermediates/urea_trans.pdb','w')
        f.write("HETATM    1  C1          0      -6.333   3.197  -0.055  1.00  0.00           C" + "\n")
        f.write("HETATM    2  N2          0      -5.141   2.394  -0.022  1.00  0.00           N" + "\n")
        f.write("HETATM    3  C3          0      -3.880   3.021   0.054  1.00  0.00           C" + "\n")
        f.write("HETATM    4  N4          0      -2.684   2.271   0.090  1.00  0.00           N" + "\n")
        f.write("HETATM    5  O5          0      -3.820   4.255   0.091  1.00  0.00           O" + "\n")
        f.write("HETATM    6  H4          0      -5.215   1.404  -0.052  1.00  0.00           H" + "\n")
        f.write("HETATM    7  C7          0      -1.421   2.955   0.167  1.00  0.00           C" + "\n")
        f.write("HETATM    8  C8          0      -2.707   0.834   0.049  1.00  0.00           C" + "\n")
        f.close()

    def __get_pdb_files(self, loc):
        """Generates a list of PDB files.
        
        Arguments:
        loc -- A string describing the directory where PDB files can be found.
        
        Returns:
        A list of the identified PDB files (string).
        
        """
        
        files = []
        if os.path.isdir(loc): # so it's a directory, go through the directory and find all the pdb files
            if loc[-1:]!=os.sep: loc = loc + os.sep # so add a / to the end of the directory
            files.extend(glob.glob(loc + '*.pdb'))
            files.extend(glob.glob(loc + '*.PDB'))
        else: # so it's a file
            if loc[-3:]=="PDB" or loc[-3:]=="pdb":
                files.append(loc)
            else:
                __log_file_output("The file " + loc + " does not have the PDB file extention and so cannot be used.",log)
        return files
    
    def __log_file_output(self, text, file=""):
        """Logs program output, either to the screen or to a log file.
        
        Arguments:
        text -- A string, the text to log.
        file -- An optional file, to which the text is saved.
        
        """
        
        if self.log_file == "":
            print text
        else: # so write to a file instead
            file.write(text + "\n")

    def __extended_help(self, log, kinds_of_reactions):
        """Provide extended help if requested."""
        
        self.__log_file_output("AutoClickChem generates two lists of PDB files, called reactants1 and reactants2.",log)
        self.__log_file_output("The program tries to combine the compounds of the first list with the compounds",log)
        self.__log_file_output("of the second list using the reactions of click chemistry. To build these two",log)
        self.__log_file_output("lists, use the tags -reactants1 and -reactants2 from the command line. Optionally,",log)
        self.__log_file_output("you can also tell AutoClickChem where to save PDB files of the products using the",log)
        self.__log_file_output("-output_dir tag. For example,\n",log)
        self.__log_file_output("python autoclickchem.py -reactants1 ../azides/myazide.pdb -reactants2 ../alkynes/myalkyne.pdb -output_dir ./output/\n",log)
        self.__log_file_output("The -reactants1 and -reactants2 tags can also specify directories, in which case",log)
        self.__log_file_output("AutoClickChem will add all the PDB files in the specified directory to the respective",log)
        self.__log_file_output("list. For example,\n",log)
        self.__log_file_output("python autoclickchem.py -reactants1 ../azides/ -reactants2 ../alkynes/ -output_dir ./output/\n",log)
        self.__log_file_output("The program is not limited to one -reactants1 or -reactants2 tag. For example, the",log)
        self.__log_file_output("following is valid:\n",log)
        self.__log_file_output("python autoclickchem.py -reactants1 ../azides/myazide.pdb -reactants1 ../azides2/ -reactants2 ../alkynes/myalkyne.pdb -reactants2 ../alkynes2/ -output_dir ./output/\n",log)
        
        self.__log_file_output("Some of the AutoClickChem reactions (e.g., alkene to epoxide) require only one type of reactant,",log)
        self.__log_file_output("not two. To access these reactions, use only the reactants1 tag. For example,\n",log)
        self.__log_file_output("python autoclickchem.py -reactants1 ../alkenes/ -output_dir ./my_epoxides/\n",log)
        
        self.__log_file_output("Note that, unlike the online version of AutoClickChem, the command-line version only accepts",log)
        self.__log_file_output("PDB files as input. We recommend using OpenBabel to convert your small-molecule models",log)
        self.__log_file_output("from other formats into the PDB format.\n",log)
        self.__log_file_output("If you would like to limit the number of pairs of compounds that are reacted together,",log)
        self.__log_file_output("the max_reactions tag can be used:\n",log)
        self.__log_file_output("python autoclickchem.py -reactants1 ../azides/ -reactants2 ../alkynes/ -output_dir ./output/ -max_reactions 10\n",log)
        self.__log_file_output("The following are special tags that refer to each of the AutoClickChem reactions:\n",log)
        
        for reaction in kinds_of_reactions:
            self.__log_file_output(reaction,log)
        
        self.__log_file_output("",log)
        self.__log_file_output("By default, AutoClickChem performs all reactions possible. Specific reactions can be turned",log)
        self.__log_file_output("on or off by by placing a + or -, respectively, before the appropriate tag. For example,",log)
        self.__log_file_output("if you want AutoClickChem to avoid opening epoxides, the following tags can be used:\n",log)
        self.__log_file_output("python autoclickchem.py -reactants1 ../azides/ -reactants2 ../alkynes/ -output_dir ./output/ -epoxide_alcohol_opening -epoxide_thiol_opening\n",log)
        self.__log_file_output("The tag \"all_reactions\" is also valid. If you'd like to perform only the",log)
        self.__log_file_output("azide_and_alkyne_to_azole reaction, for example, the following tags can be used:\n",log)
        self.__log_file_output("python autoclickchem.py -reactants1 ../azides/ -reactants2 ../alkynes/ -output_dir ./output/ -all_reactions +azide_and_alkyne_to_azole\n",log)
        self.__log_file_output("I hope you enjoy using AutoClickChem!\n",log)

    def start(self):
        """The function called to start AutoClickChem."""
        
        args = sys.argv[:]
        
        # before anything else, look for an input file and add those parameters to the commandline parameters
        for index in range(1,len(args)):
            arg = args[index]
            if arg.lower()=="-input_file": # you can also pass arguments through an input file. used mostly for the server application.
                filename = args[index + 1]
        
                args[index] = ""
                args[index + 1] = ""
        
                f = open(filename, 'r')
                cmnds = ""
                for line in f.readlines():
                    cmnds = cmnds + line.strip() + " "
                f.close()
                while "\n" in cmnds: cmnds = cmnds.replace("\n"," ")
                while "\t" in cmnds: cmnds = cmnds.replace("\t"," ")
                while "  " in cmnds: cmnds = cmnds.replace("  "," ")
                cmnds = cmnds.strip()
                cmnds = cmnds.split(" ")
                
                # so notice that parameters in the input file are placed before parameters specified from the commandline
                # This essentially makes it so that parameters specified from the command line will always take precidence
                # over ones specified in the input file.
                
                for index in range(len(cmnds)-1,-1,-1):
                    args.insert(1,cmnds[index])
                    
        self.log_file = ""
        
        # next, look to see if there's a log_file specified or if extended help is requested
        show_help = False
        for index in range(1,len(args)):
            arg = args[index]
            if arg.lower()=="-help":
                show_help = True
            if arg.lower()=="-log_file": # this is a hidden optional command used for the server application.
                self.log_file = args[index + 1]
                args[index] = ""
                args[index + 1] = ""
        
        if self.log_file != "":
            log = open(self.log_file, 'w')
        else:
            log = ""
        
        kinds_of_reactions = ["azide_and_alkyne_to_azole", "epoxide_alcohol_opening", "epoxide_thiol_opening", "chloroformate_and_amine_to_carbamate", "sulfonyl_azide_and_thio_acid", "carboxylate_and_alcohol_to_ester", "carboxylate_and_thiol_to_thioester", "acyl_halide_and_alcohol_to_ester", "acyl_halide_and_thiol_to_thioester", "ester_and_alcohol_to_ester", "ester_and_thiol_to_thioester", "acid_anhydride_and_alcohol_to_ester", "acid_anhydride_and_thiol_to_thioester", "carboxylate_and_amine_to_amide", "acyl_halide_and_amine_to_amide", "ester_and_amine_to_amide", "acid_anhydride_and_amine_to_amide", "isocyanate_and_amine_to_urea", "isothiocyanate_and_amine_to_thiourea", "isocyanate_and_alcohol_to_carbamate", "isothiocyanate_and_alcohol_to_carbamothioate", "isocyanate_and_thiol_to_carbamothioate", "isothiocyanate_and_thiol_to_carbamodithioate", "alkene_to_epoxide", "halide_to_cyanide", "alcohol_to_cyanide", "carboxylate_to_cyanide", "acyl_halide_to_cyanide", "acid_anhydride_to_cyanide", "halide_to_azide", "alcohol_to_azide", "carboxylate_to_azide", "acyl_halide_to_azide", "acid_anhydride_to_azide", "amine_to_azide", "amine_to_isocyanate", "amine_to_isothiocyanate", "azide_to_amine"]
        
        self.__log_file_output("\nAutoClickChem " + self.version,log)
        self.__log_file_output("\nIf you use AutoClickChem in your research, please cite the following reference:",log)
        self.__log_file_output("  {REFERENCE HERE}\n",log)
        
        if show_help == True:
            self.__extended_help(log, kinds_of_reactions)
            sys.exit(0)
        else: self.__log_file_output("Use the -help tag to see an extended description of how to use AutoClickChem.\n",log)
        
        # start by enabling all reactions
        reactions_to_perform = kinds_of_reactions[:]
        
        # Get reactants
        reactants1 = []
        reactants2 = []
        
        output_dir = "." + os.sep
        required_input1 = False
        #required_input2 = False
        max_reactions = -1
        make_zip = False
        
        # First, you need to see if there are multiple instances where the user has specified the output_dir tag. Delete all but the last instance, so as not to create multiple unused directories
        # find all the places -output_dir is specified
        index_of_output_dir = []
        for index in range(len(args)):
            #print "*", args[index]
            if args[index].lower() == "-output_dir": index_of_output_dir.append(index)
        
        if len(index_of_output_dir) > 1: # see if there are more than one
            index_of_output_dir.pop(-1) # if so, remove the last one
            # go through the rest and remove them from the args list, as well as the one after them (the name of the output dir specified)
            for index in index_of_output_dir:
                args[index] = ""
                args[index + 1] = ""
        
        for index in range(1,len(args)):
            arg = args[index]
            if arg.lower()=="-reactants1":
                loc = args[index + 1]
                reactants1.extend(self.__get_pdb_files(loc)) # so it's a directory
                required_input1 = True
                args[index] = ""
                args[index + 1] = ""
            elif arg.lower()=="-reactants2":
                loc = args[index + 1]
                reactants2.extend(self.__get_pdb_files(loc)) # so it's a directory
                #required_input2 = True
                args[index] = ""
                args[index + 1] = ""
            elif arg.lower()=="-output_dir":
                
                output_dir = args[index + 1]
                if output_dir[-1:] != os.sep: output_dir = output_dir + os.sep
                
                # if the output directory doesn't exist, create it
                if not os.path.exists(output_dir):
                    self.__log_file_output("The specified output directory \"" + output_dir + "\" doesn't exist, so I'll create it...\n",log)
                    os.mkdir(output_dir)
                
                args[index] = ""
                args[index + 1] = ""
            elif arg.lower()=="-max_reactions":
                max_reactions = int(args[index + 1])
                args[index] = ""
                args[index + 1] = ""
            elif arg.lower()=="-all_reactions":
                reactions_to_perform = []
                args[index] = ""
            elif arg.lower()=="+all_reactions":
                reactions_to_perform = kinds_of_reactions[:]
                args[index] = ""
            elif arg.lower()=="-make_zip": # this is a hidden optional command used for the server application.
                make_zip = True
                args[index] = ""
            elif arg[:1]=="+": # so it starts with a plus
                tag = arg[1:].lower()
                if tag in kinds_of_reactions: # so it's a valid tag
                    if tag not in reactions_to_perform: # so it's not already in the list of reactions that will be performed
                        reactions_to_perform.append(tag) # add it to the list
                        args[index] = ""
            elif arg[:1]=="-": # so it starts with a minus
                tag = arg[1:].lower()
                if tag in reactions_to_perform:
                    while tag in reactions_to_perform: # just in case for some reason it appears more than once in the list
                        reactions_to_perform.remove(tag)
                        args[index] = ""
        
        if required_input1 == False: # or required_input2 == False:
            self.__log_file_output("Error! Required input not specified. At the very least, the -reactants1 tag must\nbe specified.\n",log)
            self.__log_file_output("Example: python autoclickchem.py -reactants1 ../azides/ -reactants2 ../alkynes/\n",log)
            
            sys.exit()
        
        # Now give the user a warning if there were some tags that weren't used
        args[0]=""
        while "" in args: args.remove("")
        if len(args) > 0:
            toprint = "Warning: The following parameters were not recognized: "
            for arg in args:
                toprint = toprint + arg + " "
            self.__log_file_output(toprint,log)
            self.__log_file_output("",log)
        
        # now let the user know which reactants have been selected
        self.__log_file_output("Reactants in the group \"Reactants1\":",log)
        for reactant in reactants1:
            self.__log_file_output("   " + reactant,log)
        self.__log_file_output("",log)
        
        if len(reactants2) > 0:
            self.__log_file_output("Reactants in the group \"Reactants2\":",log)
            for reactant in reactants2:
                self.__log_file_output("   " + reactant,log)
            self.__log_file_output("",log)
        
        self.__log_file_output("Files will be output to " + output_dir + "\n",log)
        
        self.__log_file_output("The following reactions will be allowed: ",log)
        
        #for reaction in reactions_to_perform:
        self.__log_file_output(", ".join(reactions_to_perform),log)
        
        self.__log_file_output("",log)
        
        reactions = []
        # first, consider all the reactions where there's just one reactant
        for react1 in reactants1:
            reactions.append([react1])
        for react1 in reactants2:
            reactions.append([react1])
        # now get all combinations of reactant1 x reactant2
        for react1 in reactants1:
            for react2 in reactants2:
                reactions.append([react1,react2])
        
        index = 0
        count = 1
        total = len(reactions)
        if max_reactions != -1 and max_reactions < total: total = max_reactions
        self.__log_file_output(str(total) + " pairs of compounds will be reacted..." ,log)
        
        for reaction in reactions:
        
            if count > total: 
                print "Breaking with some reactions not yet performed..."
                break
        
            file1 = reaction[0]
            pdb1 = pymolecule.Molecule()
            pdb1.load_pdb(file1)
            file1 = os.path.basename(file1)
        
            file2 = "empty"
            pdb2 = pymolecule.Molecule()
            
            if len(reaction) > 1:
                file2 = reaction[1]
                pdb2.load_pdb(file2)
                file2 = os.path.basename(file2)
                self.__log_file_output("\nReacting " + file1 + " with " + file2 + " ... (Pair " + str(count) + " of " + str(total) + ")",log)
            else:
                self.__log_file_output("\nReacting " + file1 + " ... (Pair " + str(count) + " of " + str(total) + ")",log)
            
            pdb_list = self.operators_react.react_molecules(pdb1, pdb2, reactions_to_perform)
            if len(pdb_list) == 0:
                self.__log_file_output("No click reactions possible!",log)
                total = total - 1
            else:
                count = count + 1 # so only add one to the count if some reaction was possible
                
            for pdb in pdb_list:
                index = index + 1
                
                printout = pdb.save_pdb(output_dir + file1 + "." + file2 + "." + str(index) + ".pdb")
                self.__log_file_output(printout, log)
        
            if os.path.exists(output_dir + 'stop_job'): # this is only used for the server application
                self.__log_file_output("\nStopping AutoClickChem prematurely, as requested...",log)
                break
        
        
        if make_zip == True: # This is only used for the server application
            t = os.system('cd ' + output_dir + ';zip products.zip *.pdb')
        
        self.__log_file_output("",log)
        self.__log_file_output("Program done!",log)
        self.__log_file_output("",log)
        
        if self.log_file != "": log.close()


if __name__ == '__main__':
    runit = AutoClickChem()
    runit.start()
