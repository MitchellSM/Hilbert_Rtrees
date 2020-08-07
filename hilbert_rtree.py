"""
@NAME: Hilbert R-Trees
@AUTHOR: Mitchell Sulz-Martin
@DESC: An implementation of Hilbert R-Trees for 2D objects
@LIMIT: Complexity of implementing a function to map
the coordinates of a N-dimensional object to a hilbert value is beyond
the scope of this project.
I limited the dimensions of real objects to 2D.
For higher dimensions than 2 a new function for such mappings will need to be implemented
"""

import time
import sys
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
sys.setrecursionlimit(1000000)
totalTime = []
pointSearch = []
querySearch = []
insertTime = []
buildTime = []

#used to locate the next quadrant for point(x,y) during hilbert value mapping
HMAP = np.array([
    [(0, 3), (1, 0), (3, 1), (2, 0)],
    [(2, 1), (1, 1), (3, 1), (0, 2)],
    [(2, 2), (3, 1), (1, 2), (0, 1)],
    [(0, 0), (3, 2), (1, 3), (2, 3)]])


class RObj:
    """Container for real objects"""
    def __init__(self, data):
        self.mbr = data #coordinates of the object
        self.hilbert_value = None #Hilbert Value of object

    #The order of a hilbert curve defines the size of grid which it will fill
    #(i.e. order=16, grid = 16*16)
    #As the order decreases, it will "generate" a hilbert curve of that order
    #within the quadrant specified by the hMap

    def calc_h_val(self, center, order=16):
        """ Return hilbert value of a coordinate"""
        _x, _y = center
        curr, h_val = 0, 0
        for i in range(order-1, -1, -1):
            h_val <<= 2
            x_quad = 1 if _x & (1 << i) else 0
            y_quad = 1 if _y & (1 << i) else 0
            _q, curr = HMAP[curr, ((x_quad <<1)| y_quad)]
            h_val |= _q
        return h_val

class Node:
    """A node contained within a hilbert R tree"""
    def __init__(self):
        self.hilbert_value = 0 #largest hilbert value of its children, used as a key
        self.page = []
        self.level = 0
        self.mbr = []
        #coords for bounding box using in graphing (xLow, yLow, Xhigh, yHigh)

    def set_up(self, lvl=None):
        """assigns a nodes level, and hilbert value"""
        self.hilbert_value = self.page[len(self.page)-1].hilbert_value
        if lvl is not None:
            self.level = lvl
        return

    def new_mbr(self):
        """sets the minimum bounding rectangle for a node"""
        _xl, _yl, _xh, _yh = np.inf, np.inf, 0, 0
        for i in self.page:
            _xl = i.mbr[0] if i.mbr[0] < _xl else _xl
            _yl = i.mbr[1] if i.mbr[1] < _yl else _yl
            _xh = i.mbr[2] if i.mbr[2] > _xh else _xh
            _yh = i.mbr[3] if i.mbr[3] > _yh else _yh
        self.mbr = [_xl, _yl, _xh, _yh]
        return


class HilbertRtree:
    """The structure and methods for creating a hilbert R tree"""
    def __init__(self):
        self.obj_data = [] #contains the RObj, used for packing the tree
        self.root = None
        self.height = None
        self.cap = 15

    ##Helper Functions

    def new_data(self, size, maximum, data=None):
        """Creates RObj from a list of coordinates"""
        centerpoint = lambda i: ((int((i[0]+i[2])/2)), (int((i[1]+i[3])/2)))
        if data is None:
            data = []
            datalow = zip(np.random.randint(0, maximum*2, size), 
                          np.random.randint(0, maximum*2, size))
            for i in datalow:
                xhighdata = np.random.randint(i[0]+1, maximum*2+i[0])
                yhighdata = np.random.randint(i[1]+1, maximum*2+i[1])
                data.append((i[0], i[1], xhighdata, yhighdata)) 
        for i in data:
            #print(i)
            obj = RObj(i)
            obj.hilbert_value = obj.calc_h_val(centerpoint(i))
            self.obj_data.append(obj)
        self.obj_data.sort(key=lambda x: x.hilbert_value)
        return

    def print_tree_new(self, root):
        """Prints tree left to right"""
        for i in root.page:
            if isinstance(i, Node):
                _str = "page: " + str(i.hilbert_value) + " " + str(i.mbr)
                if i.level < 1:
                    _str += " " + str(i.level)
                    print(end="     ")
                else:
                    if i.level == 1:
                        print(end="---")
                    if i.level == 2:
                        print(end="--")
                    if i.level == 3:
                        print(end="-")
                    _str += " @ level: " + str(i.level)
                    
                print(_str)
                self.print_tree_new(i)
            else:
                print("    -", i.hilbert_value, i.mbr)
        print()
        return

    ##Tree Operations:

    def build_tree(self, size, maximum, test=None):
        """build_trees a hilbert tree with bottom up construction"""
        self.root = Node()
        node = Node()
        self.new_data(size, maximum, test)
        self.pack_tree(self.root, self.obj_data, node)
        self.cap += 1
        return

    def pack_tree(self, root, obj, node, idx=0, height=0):
        """packs tree using a bottom up construction"""
        while idx < len(obj):
            node.page = obj[idx : idx+self.cap]
            node.set_up(height)
            node.new_mbr()
            self.root.page.append(node)
            self.root.set_up(height+1)
            node = Node()
            idx += self.cap

        if len(root.page) > self.cap:
            root.set_up(height)
            root.new_mbr()
            node = Node()
            parents = self.root
            self.root = Node()
            self.pack_tree(self.root, parents.page, node, 0, height+1)
        self.root.new_mbr()
        self.height = self.root.level
        return

    def find_parent(self, root, k, height, level):
        """Finds the parents of a specific node"""
        if height > level+1:
            if k.hilbert_value > root.page[-1].hilbert_value:
                return self.find_parent(root.page[-1], k, height-1, level)
            else:
                for i in root.page:
                    if i.hilbert_value >= k.hilbert_value:
                        return self.find_parent(i, k, height-1, level)
        return root
    

    def insert_node(self, root, k, height=0):
        """inserts a real object into the tree"""
        leaf = self.find_parent(root, k, height, -1)
        if k.hilbert_value >= leaf.page[-1].hilbert_value:
            leaf.page.insert(len(leaf.page), k)
            leaf.set_up()
            leaf.new_mbr()
        else:
            for i in range(0, len(leaf.page)):
                if leaf.page[i].hilbert_value >= k.hilbert_value:
                    leaf.page.insert(i, k)
                    leaf.set_up()
                    leaf.new_mbr()
                    break
        if len(leaf.page) > self.cap:
            self.split_node(leaf)
            
        return


    def split_node(self, leaf):
        """split node into two by taking lowest elements"""
        if leaf != self.root:
            parent = self.find_parent(self.root, leaf, self.height, leaf.level)
            idx = parent.page.index(leaf)
            if len(parent.page)-1 > idx:
                if len(parent.page[idx+1].page) < self.cap:
                    parent.page[idx+1].page.insert(0, parent.page[idx].page[-1])
                    del parent.page[idx].page[-1]
                    
                    parent.page[idx].set_up()
                    parent.page[idx].new_mbr()
                    parent.page[idx+1].new_mbr()
                    parent.set_up()
                    parent.new_mbr()
                else:
                    self.split_leaf(parent, leaf, idx)
                    
            elif idx > 0:
                if len(parent.page[idx-1].page) < self.cap:
                    parent.page[idx-1].page.append(parent.page[idx].page[0])
                    del parent.page[idx].page[0]
                    parent.page[idx].set_up()
                    parent.page[idx].new_mbr()
                    parent.page[idx-1].new_mbr()
                    parent.set_up()
                    parent.new_mbr()
                else:
                    self.split_leaf(parent, leaf, idx)
            else:
                self.split_leaf(parent, leaf, idx)

            if len(parent.page) > self.cap:
                self.split_node(parent)
                return
        else:
            #splitting root node
            child = self.root
            self.root, node = Node(), Node()
            self.pack_tree(self.root, child.page, node, 0, self.height)
        return
    
    def split_leaf(self, parent, leaf, idx):
        new_leaf = Node()
        for i in range(0, int(len(leaf.page)/2)):
            new_leaf.page.append(leaf.page[i])
        del leaf.page[:(int(len(leaf.page)/2))]
        new_leaf.set_up(leaf.level)
        new_leaf.new_mbr()
        leaf.new_mbr()
        parent.page.insert(idx, new_leaf)
        parent.set_up()
        parent.new_mbr()
        return
    
    def intersect(self, w, n):
        """checks if w instersects n"""
        if ((n[0] < w[0] < n[2] or n[0] < w[2] < n[2]) and 
            (n[1] < w[1] < n[3] or n[1] < w[3] < n[3])):
            return True
        else:
            return False
        
    def search_tree_query(self, root, w, height, results):
        """searches tree for all nodes which overlap a bounding box"""
        for i in root.page:
            if self.intersect(w.mbr, i.mbr):
                if height > 0:
                    self.search_tree_query(i, w, height-1, results)
                else:
                    results.append(i)
        return results
    
    def search_tree_point(self, root, k):
        """searchs tree for the closes node to a point"""
        parent = self.find_parent(root, k, self.height, -1)
        if k.hilbert_value > parent.page[-1].hilbert_value:
            return parent.page[-1]
        else:
            for i in parent.page:
                if i.hilbert_value >= k.hilbert_value:
                    return i

    def delete_node(self, root, k, height):
        """deletes a node from  the tree"""
        parent = self.find_parent(root, k, height, -1)        
        for i in range(0, len(parent.page)):            
            if parent.page[i].hilbert_value == k.hilbert_value:                
                del parent.page[i]                
                if len(parent.page) < 2:                    
                    gparent = self.find_parent(root, parent, height, 0)                    
                    for j in gparent.page:
                        if (j.hilbert_value != parent.hilbert_value 
                            and len(j.page) < self.cap):                            
                            j.page.insert(0, parent.page[0])
                            del gparent.page[gparent.page.index(parent)]                            
                            j.new_mbr()
                            return
                    return
                else:
                    parent.set_up()
                    parent.new_mbr()
                    gparent = self.find_parent(root, parent, height, 0)
                    gparent.set_up()
                    gparent.new_mbr()
                    return
        return

    ##Testing/Graphing Functions:

    def small_test(self):
        """Method for testing a tree"""
        results = []
        maximum = 10 
        self.cap = 2 

        #build a tree using points from testpoints array
        print("---------------------------------")
        print("     PACKING ALGORITHM TEST")
        print("---------------------------------")
        self.build_tree(25, maximum)
        self.print_tree_new(self.root)
        
        #delete testing:
        print("---------------------------------")
        print("    DELETION ALGORITHM TEST")
        print("---------------------------------")
        print("Deleting nodes:")
        for i in self.obj_data[0:1]:
            print("  -", i.hilbert_value)
            self.delete_node(self.root, i, self.height)
        print("\n=======================================")
        self.print_tree_new(self.root)
  
        #search query and point search test:
        print("---------------------------------")
        print("QUERY/POINT SEARCH ALGORITHM TEST")
        print("---------------------------------")
        self.obj_data = []
        self.new_data(2, sqrt(maximum))
        for i in self.obj_data:    
            print("Query Searching for:", i.hilbert_value, i.mbr, "\n   Over laping nodes:")
            
            query = self.search_tree_query(self.root, i, self.height, results)
            if query:
                for j in query:
                    print("   -", j.hilbert_value, j.mbr)
            else:
                print("   - no nodes overlap with", i.hilbert_value)
            query, results = [], []
            print("\nPoint Searching for:", i.hilbert_value)
            point = self.search_tree_point(self.root, i)
            print("    -", point.hilbert_value)
            print("\n=======================================")
        self.print_tree_new(self.root)
        
        #insertion testing:
        print("---------------------------------")
        print("    INSERTION ALGORITHM TEST")
        print("---------------------------------")
        self.obj_data = []
        self.new_data(5, maximum+2)
        print("Inserting:")
        for i in self.obj_data:
            print("  -", i.hilbert_value, i.mbr)
            self.insert_node(self.root, i, self.height)
        print("\n=======================================")
        self.print_tree_new(self.root)
    
    def build_test(self, size):
        global buildTime
        t_1 = time.time()
        self.build_tree(size*10, 100)
        t_2 = time.time()
        buildTime.append((t_2 - t_1))
        return
        
    def insertion_test(self, size):
        global insertTime, sizeArr
        self.build_tree(size*10, 100)
        self.obj_data = []
        self.new_data(size, 100)
        t_3 = time.time()
        for i in self.obj_data:
            self.insert_node(self.root, i, self.height)
        t_4 = time.time()
        insertTime.append((t_4-t_3))
        return
    
    def query_test(self, size, cap):
        global querySearch
        self.cap = cap
        results = []
        self.build_tree(size*10, 100)
        self.obj_data = []
        self.new_data(size, 100)
        t_5 = time.time()
        for i in self.obj_data:
            self.search_tree_query(self.root, i, self.height, results)
        t_6 = time.time()
        querySearch.append((t_6-t_5))
        return
        
HRT = HilbertRtree()
HRT.small_test()
ctr = 100
sizeArr = []

while ctr <= 1000:
    cap = int(sqrt(ctr))
    HRT1 = HilbertRtree()
    HRT2 = HilbertRtree()
    HRT3 = HilbertRtree()
    HRT1.build_test(ctr)
    HRT2.insertion_test(ctr)
    HRT3.query_test(ctr, cap) 
    sizeArr.append(ctr)
    ctr += 100    
    
plt.figure(1)
plt.plot(buildTime, sizeArr, 'g-')
plt.plot(buildTime, sizeArr, 'mo')
plt.title("Packing time")
plt.ylabel("Number of nodes")
plt.xlabel("Time (s)")
plt.show()

plt.figure(2)
plt.plot(querySearch, sizeArr, 'y-')
plt.plot(querySearch, sizeArr, 'go')
plt.title("Query Search time")
plt.ylabel("Number of nodes")
plt.xlabel("Time (s)")
plt.show()

plt.figure(1)
plt.plot(insertTime, sizeArr, 'm-')
plt.plot(insertTime, sizeArr, 'ro')
plt.title("Insertion time")
plt.ylabel("Number of nodes")
plt.xlabel("Time (s)")
plt.show()


    

