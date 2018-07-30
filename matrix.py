"""The goal here is to create an implementation of 
    1. inverse_method
    2. the gauss elimination algorithm
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
import unittest, random
import numpy as np
import numpy.matlib as mt

eqstring = ["B + E + D + F","A +B + D + G ","C + f + E + A",
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

# print(translated_matrix)

def add_list(list1, list2):
    if len(list1) == len(list2):
        ans = []
        for i in range(len(list1)):
            ans.append(list1[i] + list2[i])
        return ans
    else: raise Exception()
    
def inverse_method(A, x, b):
    """linear system format: Ax = b"""
    #make A into a square matrix with order len(b) * len(b)
    diff = len(x) - len(A)
    for i in range(diff):
        first_rand = random.randint(0, len(A) - 1)
        second_rand = random.randint(0, len(A) - 1)
        A.append(add_list(A[first_rand],A[second_rand]))
        b.append(b[first_rand] + b[second_rand])
    try:
        inverse = np.linalg.inv(A)
        sols = inverse * b
        sols_dic  = {}
        for i in range(len(sols)):
            sols_dic[x[i]] = sols[i]
        return sols_dic
    except np.linalg.linalg.LinAlgError as err:
        print("err.message") #usually singular matrix not invertible
    return A, x, b

A, x, b = inverse_method(translated_matrix, unknownset, solutions)
x=list(x) # for debug purposes
for i in range(len(x)):
    print(A[i], x[i], b[i])