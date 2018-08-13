# SNPSolvit Hard Puzzle
[![Build Status](https://travis-ci.org/p-netm/snphsp.svg?branch=master)](https://travis-ci.org/p-netm/snphsp)

## installation
copy the git url and clone it, or download the files
create virtualenv 
	virtualenv venv

## Execution
	python matrix.py 
	
### Sample
##### (A sample problem, sample input, sample output)
### Sample Problem

-|-|-|-|
D|G|E|T|22
B|A|H|D|23
E|C|G|C|14
H|J|F|A|22
16|26|11|28

### Sample Input

>type in the equations below: sample: E + B + F + J = 19(q/quit)
>eq. #1 : D + G + E + J = 22
>eq. #2 : B + A + H + D = 23
>eq. #3 : E+C+G+C=14
>eq. #4 : h+j+f+a=22
>eq. #5 : d+B+E+h = 16
>eq. #6 : G+A+C+J=26
>eq. #7 : E+H+G+f=11
>eq. #8 : J+D+C+A=28

### Sample output

```python
{'B': 1, 'C': 6, 'D': 8, 'H': 5, 'A': 9, 'G': 4, 'E': 2, 'J': 8, 'F': 0}
```