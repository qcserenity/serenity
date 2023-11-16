#!/usr/bin/python3
import sys
"""
@author Jan P. Unsleber
@date   Dec 01, 2015

A script to convert some of the raw C6ab data from Grimme et al to C++ code.

The filename is taken as the first argument.

The data should be in a plain file, ordered like this:

               3.0267, 1, 1, 0.91180, 0.91180, 
               2.0835, 2, 1, 0.00000, 0.91180, 
               1.5583, 2, 2, 0.00000, 0.00000, 
               38.9448, 3, 1, 0.00000, 0.91180, 
               22.1508, 3, 2, 0.00000, 0.00000, 
"""

#
# Initialization
#
i = []
j = []
k = []
l = []
y = []
cn1 = []
cn2 = []

#
# Read data
#
with open(sys.argv[1],'r') as bla:
    for line in bla:
        split = line.split(',')
        for s in split:
            s.strip()
        if int(split[1])>=100:
            i.append(int(split[1].strip()[1:]))
        else:
            i.append(int(split[1].strip()))

        if int(split[2])>=100:
            j.append(int(split[2].strip()[1:]))
        else:
            j.append(int(split[2].strip()))

        if int(split[1])>=100:
            k.append(int(split[1].strip()[0]))
        else:
            k.append(0)

        if int(split[2])>=100:
            l.append(int(split[2].strip()[0]))
        else:
            l.append(0)
        y.append(split[0].strip())
        cn1.append(split[3].strip())
        cn2.append(split[4].strip())

#
# Order data
#
c6 = [ [ [ [ ['-1.0' for i in range(3)] for j in range(5)] for k in range(5)] for l in range (94)] for m in range (94)]
mxc = ['-1.0' for i in range(94)]
for ii,jj,kk,ll,yy,c1,c2 in zip(i,j,k,l,y,cn1,cn2):
    c6[ii-1][jj-1][kk][ll] = [yy,c1,c2]
    c6[jj-1][ii-1][ll][kk] = [yy,c2,c1]
    mxc[ii-1] = max(float(mxc[ii-1]),float(kk))
    mxc[jj-1] = max(float(mxc[jj-1]),float(ll))

#
# Print the code
#
print('static constexpr std::array<unsigned int, 94> maxNC6Raw')
mxcString = '{{ \n '
for i,m in enumerate(mxc): 
    mxcString += '{:1.0f}'.format(m)
    mxcString += ', '
    if ((i+1)%10 == 0):  mxcString += ' \n'
mxcString += '}};'
print(mxcString)
print('')
print('static constexpr std::array<')
print('                    std::array<')
print('                       std::array<')
print('                          std::array<')
print('                             std::array<double,3>,5>,5>,94>,94> c6abRaw')
print('{{')
for i in range(94):
    print('  {{')
    for j in range(94):
        print('    {{')
        for k in range(5):
            print('      {{')
            for l in range(5):
                if l==4: print('        {{'+c6[i][j][k][l][0]+','+c6[i][j][k][l][1]+','+c6[i][j][k][l][2]+'}}' )
                else: print('        {{'+c6[i][j][k][l][0]+','+c6[i][j][k][l][1]+','+c6[i][j][k][l][2]+'}},' )
            if k ==4: print('      }}')
            else: print('      }},')
        if j ==93: print('    }}')
        else: print('    }},')
    if i ==93: print('  }}')
    else: print('  }},')
print('}};')

