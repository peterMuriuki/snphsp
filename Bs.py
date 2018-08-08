from fractions import Fraction
a, b, c, d, e, f, g, h, j= 0,0,0,0,0,0,0,0,0
cset_dict = {}
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
    for h in range(1,11):
      for j in range(1,11):
        if j != h:
            temp = Fraction(79,3) -2*j-Fraction(7,3)*h
            if corr(temp, h, j):
              g = temp
            else: continue
            temp=-Fraction(5,2)-Fraction(1,2)*g+Fraction(1,2)*h+j
            if corr(temp, h, j):
              f=temp
            else: continue
            temp=25-h-2*j
            if corr(temp, h, j):
              e=temp
            else: continue
            temp=2+f
            if corr(temp, h, j):
              d=temp
            else: continue
            temp=16-j-h-g
            if corr(temp, h, j):
              c = temp
            else: continue
            temp=9-Fraction(1,2)*d-Fraction(1,2)*e
            if corr(temp, h, j):
              b =temp
            else: continue
            temp=17-c-f-g
            if corr(temp, h, j):
              a = temp
              check = True
              raise Exception()
except Exception:
    pass
    
if check:
  print(a,b,c,d,e,f,g,h,j)
else:
  print("not converge")

print(a,b,c,d,e,f,g,h,j)
