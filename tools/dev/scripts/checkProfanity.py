#!/usr/bin/python3
#@file   checkProfanity.py
#
#@date   Nov 9, 2017
#@author Jan Unsleber, David Schnieders
#@copyright \n
# This file is part of the program Serenity.\n\n
# Serenity is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.\n\n
# Serenity is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.\n\n
# You should have received a copy of the GNU Lesser General
# Public License along with Serenity.
# If not, see <http://www.gnu.org/licenses/>.\n

import sys,os


try:
  src = sys.argv[1]
  os.path.isfile(os.path.join(src,"/serenity.cpp"))
except:
  print("Path to src/ folder needed as first argument to this script.")
  sys.exit(1)

profanityList = [
"anal","anus","arse","ass","ballsack","balls","bastard","bitch"
"biatch","bloody","blowjob","blow job","bollock","bollok","boner"
"boob","bugger","bum","butt","buttplug","clitoris","cock","coon"
"crap","cunt","damn","dick","dildo","dyke","fag","feck","fellate"
"fellatio","felching","fuck","f u c k","fudgepacker","fudge packer"
"flange","goddamn","god damn","hell","jerk","jizz","knobend",
"knob end","labia","lmao","lmfao","muff","nigger","nigga","omg","penis"
"piss","poop","prick","pube","pussy","queer","scrotum","sex","shit",
"s hit","sh1t","slut","smegma","spunk","tit","tosser","turd","twat"
"vagina","wank","whore"]

proffiles = []
prof = {}
cwd = os.getcwd()
os.chdir(src)
for root, dirnames, filenames in os.walk('.'):
  for filename in filenames:
    fname = os.path.join(root, filename)
    with open(fname,'r') as f:
      tmp = []
      for l,line in enumerate(f):
        for s in line.split():
          if s.lower() in profanityList:
            tmp.append([l,s])
      if len(tmp) >0: 
        proffiles.append(os.path.join(root, filename))
        prof[os.path.join(root, filename)] = tmp

os.chdir(cwd)

fail = False
if len(proffiles) >0:
  print("Profanity was found in the following files:")
  for f in proffiles:
    for p in prof[f]:
      print("{:s}:{:4d} {:s}".format(f,p[0],p[1]))
  fail= True
  print()

if fail: sys.exit(1);
