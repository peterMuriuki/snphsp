"""The goal here is to create an implementation of the gauss elimination algorithm
for finding unkowns in a system of linear equations

The solution will then hence forth be used to solve the solvit problem that happens 
in the leisure page of the standard newspaper
B + E + D + F = 17
A +B + D + G = 26
C + f + E + A = 16
H + J + G + H =18
B + A + C + H = 20
E + B + F + J = 19
D + D + E + G = 20
F + G + A + H = 18
"""
import unittest

class Matrix(object):
    def __init__(self, object, columns=0,filler=0):
        """:parameters: a list of objects, alist of lists of objects
        a list of tuples of objects"""
        self.data = []
        self.rows = len(self.data)
        self.columns = len(self.data[0])
        
        if isinstance(object, int) and isinstance(object, int):
            # we have dimensions for row;
            self.data = [[filler] * columns] * object
        elif isinstance(object, list) or isinstance(object, tuple):
            #first we consider the 2 dimensional matrix, columns should have the same no. of elements
            if type(object[0][0]) == list or type(object[0][0] == tuple)
                self.data = list(object)
                for el in self.data:
                    el = list(el)
            # first we  consider a one dimensional sequence of data
            if type(object[0][0]) == int:
                self.data = [element for element in object if isinstance(element, int) else raise Exception("Expected int data, but got {}".format(type(element)))]
            
        return
    
    def __repr__(self):
        for row in self.data:
            templ = [' '.join(row) for row in self.data]
        return '\n'.join(templ)
    
    def __eq__(self, other):
        if not isinstance(other, Matrix):
            return False
        if self.rows != other.rows or self.columns != other.columns:
            return False
        # we can now compare element by element
    
class MatrixTests(unittest.Testcase):
    def setUp(self):
        """"""
        pass
    
    def tearDown(self):
        """"""
        pass
    
    def test_matrix_creation_using_row_and_column_input_as_int(self):
        sample = Matrix(1,1,0)
        self.assert