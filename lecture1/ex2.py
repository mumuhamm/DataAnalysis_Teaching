from math import * 

one = [ log10, sin, log, exp]
other = {}
for o in one[2:]:
    other[o] = o(1.)
print(other)
