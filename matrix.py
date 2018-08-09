"""The goal here is to create an implementation of 
    1. inverse_method
    2. the gauss elimination algorithm: 
        reducing a matrix to row echelon form
        performing back substitution on the returned row echelon matrix
for finding unkowns in a system of linear equations

The solution will then hence forth be used to solve the solvit problem that happens 
in the leisure page of the standard newspaper
:sample:
B + E + D + F = 17
A + B + D + G = 26
C + f + E + A = 16
H + J + G + H =18
B + A + C + H = 20
E + B + F + J = 19
D + D + E + G = 20
F + G + A + H = 18
"""
import unittest, random
import numpy as np
import numpy.matlib as mt
from copy import deepcopy
from fractions import Fraction
import unittest

def translate(eqstring, solutions):
    solutions = [int(i) for i in solutions]
    eqstring = [char.upper() for char in eqstring]
    eqlist = [map(str.strip, char.split("+")) for char in eqstring]
    temp = []
    for charlist in eqlist:
        temp.extend(charlist)
    unknownset = set(temp)

    #issue 1: the system does not translate to a square matrix
    translated_matrix  = []
    for string in eqstring:
        coeff = []
        for unknown in unknownset:
            coeff.append(string.count(unknown))
        translated_matrix.append(coeff)
    return translated_matrix, list(unknownset), solutions


def add_list(list1, list2):
    if len(list1) == len(list2):
        ans = []
        for i in range(len(list1)):
            ans.append(list1[i] + list2[i])
        return ans
    else: raise Exception()
        
def squarify(A, x, b):
    #make A into a square matrix with order len(b) * len(b)
    diff = len(x) - len(b)
    for i in range(diff):
        first_rand = random.randint(0, len(A) - 1)
        second_rand = random.randint(0, len(A) - 1)
        A.append(add_list(A[first_rand],A[second_rand]))
        b.append(b[first_rand] + b[second_rand])
    return A, x, b
    
def inverse_method(A, x, b):
    """linear system format: Ax = b"""
    A, x, b = squarify(A, x, b)
    try:
        inverse = np.linalg.inv(A)
        sols = inverse * b
        sols_dic  = {}
        for i in range(len(sols)):
            sols_dic[x[i]] = sols[i]
        return sols_dic
    except np.linalg.linalg.LinAlgError as err:
        print(err.message) #usually singular matrix not invertible
    return A, x, b

A, x, b = squarify(translated_matrix, unknownset, solutions)


#Gaussian eliminations
def echelon(A, x, b):
    #reducing to row echelon format
    x = list(x)
    try:
        n , m = len(A), len(A[0])
    except IndexError:
        print("not a Matrix, seems like a vector")
    # for a matrix A with order n x m ; n rows and m columns
    rows, columns = 0, 0
    while rows < n and columns < m:
        # get the max for certain columns for all rows
        temp = [abs(A[i][columns]) for i in range(rows, n)]
        big_ind = np.argmax(temp) #found pivot
        big_ind = big_ind + len(A) - len(temp)
        big_value = A[big_ind][columns]
        if A[big_ind][columns] == 0:
            #no pivot, skip column
            columns += 1
        else:
            #swap big_ind with current rows **** DONT FORGET ABOUT THE b ****
            if big_ind != rows:
                A[big_ind], A[rows] = A[rows], A[big_ind]
                b[big_ind], b[rows] = b[rows], b[big_ind]
            #optimize full big_ind row by dividing with big-value use fractions
            for index in range(m):
                A[rows][index] = Fraction(A[rows][index], big_value)
            b[rows] = Fraction(b[rows], big_value)
            #for all rows below pivot: and b
            for index in range(rows+1, n):
                fact = deepcopy(A[index][columns])
                if A[index][columns] < 0 and A[rows][columns]*fact < 0 or A[index][columns] > 0 and A[rows][columns]*fact > 0:
                    sign = -1
                else:
                    sign = 1
                
                if fact != 0:
                    for idx in range(columns, m):
                        A[index][idx] = A[index][idx] + A[rows][idx] * fact *sign
                    b[index] = b[index] + b[rows] * fact * sign
                    #subtract A[rows][column] * 1 - full rows
            columns += 1
            rows += 1
    return A, x, b

echelon(A, x, b)

def back_substitution(A, x, b):
    """Back substitution"""
    # seems like the puzzle leaves 2 variables to be guessed- for this we have no other choice but brute force
    # populate dict with unkowns and initialize with zero
    d = {}
    for _x in list(x):
        d[_x] = 0
    cset_dict = {} # dict to preserve the state for corr last check
    
    def corr(key, iter1, iter2):
        #key is between 0 and 10 excluding both
        _in = key < 10 and key > 0
        if not _in:
            return False
        #key is int and not float
        _int = Fraction(key, 1).denominator == 1
        if not _int:
            return False
        #key is not the same for the current iterations
        if cset_dict.get((iter1, iter2)) is not None:
            #means that we have already started recording for this iter
            _found = key in cset_dict[(iter1, iter2)]
            if not _found:
                cset_dict[(iter1, iter2)].add(key)
        else:
            cset_dict[(iter1, iter2)] = set()
            cset_dict[(iter1, iter2)].add(key)
            _found = False
        return _in and _int and not _found

    check = False
    try:
        for d[x[len(x) -2]] in range(1, 10):
            for d[x[len(x) -1]] in range(1, 10):
                if d[x[len(x) - 2]] != d[x[len(x) - 1]]:
                    for row in range(len(A) - 1, -1, -1):
                        #find the last pivot: will be the first non zero value in current row
                        def pivot(A):
                            """Given a reduced row it will return first non zero"""
                            for index, value in enumerate(A[row]):
                                if value != 0:
                                    return value, index
                            return 0, len(A[row])-1
                        pivvalue, pivindex = pivot(A)
                        if pivvalue == 0:
                            continue #? what if there does not exist a pivot in this row
                        total = 0
                        for i in range(len(x)):
                            if i != pivindex:
                                total += A[row][i] * d[x[i]]
                        temp = (b[pivindex] - total)/pivvalue
                        if corr(temp, d[x[len(x) -2]], d[x[len(x) -1]]):
                            d[x[pivindex]] = temp
                            if pivindex == 0:
                                check = True
                                raise StopIteration("Found our matches")
                        else:
                            break
    except StopIteration as err:
        pass
    if check:
        return d
    else:
        return "Could not converge, no values of unknowns found within range(1,10)"
    
def validate(_input):
    """format :E + B + F + J = 19"""
    #^\S{1}\s*[+]\s*\S{1}\s*[+]\s*\S{1}\s*[+]\s*\S{1}\s*=\s*\d+\s*$
    pattern = r'^\s*\S{1}\s*[+]\s*\S{1}\s*[+]\s*\S{1}\s*[+]\s*\S{1}\s*=\s*\d+\s*$'
    if _input.lower() == "quit" or _input.lower() == "q":
        raise  Exception("Program stopping")
    return re.search(pattern, _input)
    
def separate(_input):
    """sample input = :E + B + F + J = 19"""
    temp = _input.split("=")
    return [char.strip() for char in temp]

def input():
    """defines the input data structure and form"""
    # am thinking using the command line and filling each linear system in a linear line
    eqstring, solutions = [], []
    print("type in the equations below: sample: E + B + F + J = 19")
    eqstring = []
    for i in range(1,9):
        while True:
            uinput = input("#{}. eq:".format(i))
            if validates(uinput):
                res = separate(uinput)[0]
                eqstring.append(res[0])
                solutions.append(res[1])
                break
            else:
                print("seems like something went wrong, please retype that:")
    #here should have A, x, and b
    return eqstring, solutions
            
if __name__ == '__main__':
    while True:
        try:
            eqstring, solutions = input() #put in while loop
            A, x, b = translate(eqstring, solutions)
            A, x, b = echelon(A,x,b)
            d = back_substitution(A,x,b)
            print(d)
            #the end
        except Exception as error:
            if error.args[0] == "Program stopping":
                print(error)
                break
            else:
                raise(error)

class MatrixTests(unittest.TestCase):
    
    def test_echelon_function_simple(self):
        A = [[2, 1, 1], 
             [1, 2, 1],
             [1, 1, 2]]
        b = [1, 1, 1]
        Asol = [
                [Fraction(1, 1), Fraction(1,2), Fraction(1,2)],
                [Fraction(0, 1), Fraction(1, 1), Fraction(1,3)],
                [Fraction(0, 1), Fraction(0,1), Fraction(1,1)]
        ]
        bsol = [Fraction(1, 2), Fraction(1, 3), Fraction(1, 4)]
        self.assertListEqual(Asol, echelon(A, [], b)[0])
        self.assertListEqual(bsol, echelon(A, [], b)[2])
    
    def test_echelon_function_with_sample_problem_matrix(self):
        A =[[1, 1, 1, 0, 0, 0, 0, 0, 1],
            [1, 0, 0, 0, 0, 1, 0, 1, 1], 
            [0, 1, 1, 1, 0, 1, 0, 0, 0], 
            [0, 0, 0, 0, 2, 0, 1, 1, 0], 
            [0, 0, 0, 1, 1, 1, 0, 0, 1], 
            [0, 1, 1, 0, 0, 0, 1, 0, 1],
            [2, 0, 1, 0, 0, 0, 0, 1, 0], 
            [0, 1, 0, 0, 1, 1, 0, 1, 0]]
        b = [17, 26, 16, 18, 20, 19, 20, 18]
        x = ["D","F","E","C","H","A","J","G","B"]
        Asol = [
            [Fraction(1, 1), Fraction(0, 1), Fraction(1, 2), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 2), Fraction(0, 1)], 
            [Fraction(0, 1), Fraction(1, 1), Fraction(1, 1), Fraction(1, 1), Fraction(0, 1), Fraction(1, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)], 
            [Fraction(0, 1), Fraction(0, 1), Fraction(1, 1), Fraction(1, 1), Fraction(-1, 1), Fraction(0, 1), Fraction(0, 1), Fraction(-1, 1), Fraction(0, 1)],
            [Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 1), Fraction(1, 1), Fraction(1, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 1)], 
            [Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 1), Fraction(0, 1), Fraction(1, 2), Fraction(1, 2), Fraction(0, 1)], 
            [Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 1), Fraction(0, 1), Fraction(2, 1), Fraction(-3, 1)], 
            [Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 1), Fraction(-1, 1), Fraction(4, 1)], 
            [Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)]]
        bsol = [Fraction(10, 1), Fraction(16, 1), Fraction(-2, 1), Fraction(20, 1), Fraction(9, 1),
                Fraction(0, 1), Fraction(28, 1), Fraction(0, 1)]
        res = echelon(A, x, b)
        self.assertListEqual(res[0], Asol)
        self.assertListEqual(res[2], bsol)