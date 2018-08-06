"""The goal here is to create an implementation of 
    1. inverse_method
    2. the gauss elimination algorithm
for finding unkowns in a system of linear equations

The solution will then hence forth be used to solve the solvit problem that happens 
in the leisure page of the standard newspaper
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

eqstring = ["B + E + D + F","A +B + D + G ","C + F + E + A",
            "H + J + G + H","B + A + C + H","E + B + F + J ","D + D + E + G", "F + G + A + H"]
solutions = [17, 26, 16, 18, 20, 19, 20, 18]
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
def gauss(A, x, b):
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
        temp = [A[i][columns] for i in range(rows, n)]
        big_ind = np.argmax(temp) #found pivot
        big_value = temp[big_ind]
        big_ind = big_ind + len(A) - len(temp) 
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
            for i in range(len(A)):
                print(str([float(j) for j in A[i]]), end="\t")
#                 print(x[i], end="\t")
                print(str(b[i]))
                print(x)
            print("*" * 100)

gauss(A, [], b)
def prettify_to_file(A,b):
    with open('sample.txt', 'w') as file:
        for i in range(len(A)):
            file.write(str([float(j) for j in A[i]]))
            file.write(str(b[i]))
    return
prettify_to_file(A, b)

def back_substitution(A, x, b):
    """Back substitution"""
    # seems like the puzzle leaves 2 variables to be guessed- for this we have no other choice but brute force
    #populate dict with unkowns and initialize with zero
    d = {}
    for _x in list(x):
        d[_x] = 0
    while d[x[len(x) -2]] < 10:
        while d[x[len(x) -1]] < 10:
            if d[x[len(x) - 2]] != d[x[len(x) - 1]]:
                for row in range(len(A) - 1, 0):
                    #find the last pivot: will be the first non zero value in current row
                    def pivot(A):
                        for index in enumerate(A[len(A) - 1]):
                            if val != 0:
                                return value, index
                    pivvalue, pivindex = pivot(A)
                    total = 0
                    for i in range(len(x)):
                        if i != pivindex:
                            total = A[row][i] * d[x[i]]
                    temp = (b[pivindex] - total)/pivvalue
                    if temp > 10:
                        d[x[len[x] - 1]] += 1
                        continue
                    else:
                        d[x[pivindex]] = temp
        d[x[len[x] - 2]] += 1
    return d

class MatrixTests(unittest.TestCase):
    def setup(self):
        pass
    def tearDown(self):
        pass
    
"""
A = [   
    [1, 0, 1, 0, 0, 1, 1, 0, 0], 
    [1, 0, 0, 0, 1, 1, 0, 0, 1], 
    [0, 0, 1, 0, 0, 0, 1, 1, 1], 
    [0, 2, 0, 1, 1, 0, 0, 0, 0], 
    [1, 1, 0, 0, 0, 0, 0, 1, 1], 
    [1, 0, 1, 1, 0, 0, 1, 0, 0], 
    [0, 0, 0, 0, 1, 2, 1, 0, 0], 
    [0, 1, 1, 0, 1, 0, 0, 0, 1],
    [1, 2, 1, 0, 1, 0, 0, 1, 2]
]
b = [17, 26, 16, 18, 20, 19, 20, 18, 38]

#Gaussian eliminations
def gauss(A, x, b):
    #reducing to row echelon format
    try:
        n , m = len(A), len(A[0])
    except IndexError:
        print("not a Matrix, seems like a vector")
    # for a matrix A with order n x m ; n rows and m columns
    rows, columns = 0, 0
    while rows < n and columns < m:
        # get the max for certain columns for all rows
        temp = [A[i][columns] for i in range(rows, n)]
        big_ind, big_value = max(enumerate(temp), key=operator.itemgetter(1))
        # big_ind = np.argmax(temp) #found pivot
        # big_value = temp[big_ind]
        big_ind = big_ind + len(A) - len(temp) 
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
            for i in range(len(A)):
                print(str([float(j) for j in A[i]]), end="\t")
                print(str(b[i]))
            print("*" * 100)
gauss(A, [], b)

"""