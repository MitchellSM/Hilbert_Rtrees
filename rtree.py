"""
Author: Mitchell Sulz-Martin
ID: 001185643
Date: Feburary 1, 2019
Desc: Implimentation of R-tree
"""

import numpy as np
objList = []

class RealObject:
    """ Object with lower left (x1, y1) and upper-right (x2, y2) coordinates. """
    def __init__(self, ll, ur):
        """ Initialize an (non-)empty polygon/point object """
        self.lowleft = ll
        self.upright = ur

    

class Node: 
    """ Non-leaf node object. """
    def __init__(self, children = [], ll = None, ur = None):
        """Initialize an (non-)empty non-leaf node"""
        self.page = children # List of Children 
        self.lowleft = ll # Lower-left coord for mbr (x1, y1)
        self.upright = ur # Upper-right coord for mbr (x2, y2)
        self.height = 0 # Height of node within the tree
    
   
    def setup(self):
        """ Predefined list of objects to insert into the tree """
        # m1 (2, 13)(4, 18)
        objList.append(RealObject((2, 13), (4, 18)))
        # m2 (3, 9)(5, 14)
        objList.append(RealObject((3, 9), (5, 14)))
        # m3 (12, 2)(14, 7)
        objList.append(RealObject((12, 2), (14, 7)))
        # m4 (7, 2)(8, 10)
        objList.append(RealObject((7, 2), (8, 10)))
        # m5 (16, 2)(22, 9)
        objList.append(RealObject((16, 2), (22, 9)))
        # m7 (20, 7)(25, 14)
        objList.append(RealObject((20, 7), (25, 14)))
        # m8 (26, 5)(28, 8)
        objList.append(RealObject((26, 5), (28, 8)))
        # m9 (18, 16)(25, 18)
        objList.append(RealObject((18, 16), (25, 18)))
        # p5 (12, 14)
        objList.append(RealObject((12, 14), (12, 14)))
        # p9 (27, 3)
        objList.append(RealObject((27, 3), (27, 3)))
        return 
    
    def new_mbr(self):
        """sets the minimum bounding rectangle for a node"""
        x1, y1, x2, y2 = np.inf, np.inf, 0, 0
        if(len(self.page) > 0):
            for c in self.page:
                x1 = c.lowleft[0] if c.lowleft[0] < x1 else x1
                y1 = c.lowleft[1] if c.lowleft[1] < y1 else y1
                x2 = c.upright[0] if c.upright[0] < x2 else x2
                y2 = c.upright[1] if c.upright[1] < y2 else y2
            self.lowleft = (x1, y1)
            self.upright(x2, y2)
        else:
            self.lowleft = (np.NINF, np.NINF)
            self.upright = (np.inf, np.inf)
        return
    
    def printObjectData(self):
        """ Prints all predefined objects """
        for i in objList:
            print(i.type, ",", i.lowleft, ",", i.upright)
        return
    
class RTree:
    """ Implementation of RTree Methods and Functions """
    
    def __init__(self):
        """ Initilize an empty R-tree """
        self.root = Node()
        self.root.new_mbr() # (x1, y1)(x2, y2) = (-inf, -inf)(inf, inf)
        self.totalHeight = 0 # Total Height of the current tree
        self.maxCap = 4 # Node bucket maximum cap
        self.minCap = 1 # Node bucket minimum cap
    
    def build(self, node = Node(), idx=0, height=0):
        """ Build tree from a predefined list of objects """
        return
    
    #(0, 1, 2, 3) -> (x1, y1) (x2, y2)
    def intersect(self, node1, node2):
        """ Bool Check for node MBR intersection """
        if((node1.lowleft[0] < node2.lowleft[0] < node1.upright[0] or
            node1.lowleft[0] < node2.upright[0] < node1.upright[1]) and
           (node1.lowleft[1] < node2.lowleft[1] < node1.upright[0] or
            node1.lowleft[1] < node2.upright[1] < node1.upright[1])):
            return True
        else:
            return False
    
    def enlargementValue(self, node, leaf):
        """ Calculate the size of the enlargement of an internal node """
        xL = abs(node.lowleft[0] - leaf.lowleft[0])
        yL = abs(node.lowleft[1] - leaf.lowleft[1])
        xH = abs(node.upright[0] - leaf.upright[0])
        yH = abs(node.upright[1] - leaf.upright[1])
        
        return

    def insert():
        return
    
    def search():
        return
    
    def delete():
        return
    
    def overflow():
        return
    
    def underflow():
        return
