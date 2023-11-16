#!/usr/bin/env python
# coding: utf-8

#@file   trackRessources.py
#
#@date   Aug 02, 2021
#@author Lars Hellmann
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

#Use this script to track CPU and memory usage for one process
#When executed this script waits for the specified process to 
#start and then records memory and CPU usage every second. When
#the specified proess stops this script dumps all data into a
#csv-file in the specified directory and displays a first plot
#as an overview.

import os
import sys
import getopt
from datetime import datetime
import time
import getpass
import psutil # You might need to download this using 'pip install psutil'
import matplotlib.pyplot as plt
import pandas as pd # You might need to download this using 'pip install pandas'


def checkIfProcessRunning(processName,user):
    #Iterate over the all the running process
    for proc in psutil.process_iter():
        try:
            # Check if process name contains the given name string.
            if processName.lower() in proc.name().lower():
                if user.lower() in proc.username().lower():
                    return True
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass
    return False;

def checkIfProcessRunningPid(pid):
    #Iterate over the all the running process
    for proc in psutil.process_iter():
        try:
            # Check if process name contains the given name string.
            if pid == proc.pid:
                return True
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass
    return False;

def getPid(processName,user):
    for proc in psutil.process_iter():
        try:
            # Check if process name contains the given name string.
            if processName.lower() in proc.name().lower():
                if user.lower() in proc.username().lower():
                    return proc.pid
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            print("Process not found")
            sys.exit(2)
    return 0;

def getProc(processName,user):
    for proc in psutil.process_iter():
        try:
            # Check if process name contains the given name string.
            if processName.lower() in proc.name().lower():
                if user.lower() in proc.username().lower():
                    return proc
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            print("Process not found")
            sys.exit(2)
    return False;

def appendToCsv(path,times,cpu,mem):
    if (os.path.isfile(path+"data.csv") == False):
        #write header
        f = open(path+"data.csv", "w")
        f.write("time cpu mem\n")
        f.close()
    l = zip(times,cpu,mem)
    print(l)
    f = open(path+"data.csv", "a")
    for t in l:
        f.write(str(t[0])+" "+str(t[1])+" "+str(t[2])+"\n")
    f.close()

def trackProcess(proc,path):
    print("PID:"+str(proc.pid))
    starttime=datetime.now()
    times =[]
    cpu = []
    mem = []
    while (proc.pid in psutil.pids()):
        print(proc.pid)
        times.append((datetime.now()-starttime).total_seconds())
        cpu.append(proc.cpu_percent())
        mem.append(proc.memory_full_info().uss / (1024*1024*1024))
        time.sleep(1)
        
    print("times: ", times)
    print("cpu: ",cpu)
    print("mem: ",mem)
    appendToCsv(path,times,cpu,mem)
        
def plotFromCsv(path):
    print("plotting")
    df = pd.read_csv(path+"data.csv", delimiter=' ')
    
    fig,ax = plt.subplots()
    ax.plot(df.time,df.cpu)
    ax.set_xlabel("time / s")
    ax.set_ylabel("cpu / %")
    ax2=ax.twinx()
    ax2.plot(df.time,df.mem,color='orange')
    ax2.set_ylabel("mem / GB")
    plt.show()
    fig.savefig('data.jpg',
            format='jpeg',
            dpi=100,
            bbox_inches='tight')
    
def main(argv):
    user = getpass.getuser()
    process = 'serenity'
    path = os.getcwd()
    try:
        opts, args = getopt.getopt(argv,"hu:p:d:",["user=","process=","path="])
    except getopt.GetoptError:
        print('trackRes.py -u <user> -p <process> -d <directory path>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('trackRes.py -u <user> -p <process name> -d <directory path>')
            sys.exit()
        elif opt in ("-u", "--user"):
            user = arg
        elif opt in ("-p", "--process"):
            process = arg
        elif opt in ("-d", "--dir"):
            path = arg               
    print("User: "+ user)
    print("Proc: "+ process)
    print("Path: "+ path)
    
    #wait until process is started
    while(checkIfProcessRunning(process,user) == False):
        time.sleep(1)
       
    pid = getPid(process,user)
    path += "/"+str(pid)+"_"
    proc = getProc(process,user)
    trackProcess(proc,path)
    plotFromCsv(path)

if __name__ == "__main__":
    main(sys.argv[1:])
