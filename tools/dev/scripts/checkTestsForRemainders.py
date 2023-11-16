#!/usr/bin/env python3

######################################################
# To use this script create an empty folder, source 
# the serenity version you want to use and adjust 
# the pathToTest variable below.
#
# When you run this script a list of all available tests 
# is generated and each test is executed in a separate folder. 
# If no files are left behind by the test its folder will
# be deleted. 
#
# When everything is finished the names of the remaining
# folders correspond to the tests that do not clean up
# their files.
######################################################


import sys,os
import subprocess as sp

#adjust to the desired serenity path!
#Example: pathToTest="/home/user/serenity/bin"
pathToTests=""

if pathToTests == "":
	print("Please read the instructions in the script")
	exit()

#generate list of all tests
bashCommand=pathToTests+"/serenity_tests --gtest_list_tests"
output = sp.check_output(bashCommand.split(), stderr=sp.STDOUT)
outString = str(output).split('\\n')

print("Start Selection")

#Outfile for the output if wanted
f = open("out", "w")

testGroup=""
testName=""
#loop over all tests
for line in outString:
	if line[0] != " ":
		testGroup=line.strip()
	if line[0] == " ":
		testName=line.strip()
		print("Running: "+testGroup+testName)

		#mkdir 
		os.mkdir(testGroup+testName)
		os.chdir(testGroup+testName)

		#execute current test
		bashCommand2=pathToTests+"/serenity_tests --gtest_filter="+testGroup+testName
		proc = sp.Popen(bashCommand2.split(), stdout=f)
		proc.communicate()

		#check if the directory is empty and delete if it is
		if len(os.listdir(os.getcwd())) == 0:
			os.chdir("../")
			os.rmdir(testGroup+testName)
		else:
			os.chdir("../")
f.close()
