/**
 * @file   FormattedOutput.cpp
 *
 * @date   Apr 14, 2014
 * @author Jan Unsleber
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
/* Include Class Header*/
#include "io/FormattedOutput.h"
/* Include Serenity Internal Headers */
#include "memory/MemoryManager.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <vector>

namespace Serenity {
using namespace std;

void setOutputOptions(const unsigned int digits) {
  cout.precision(digits);
  /**
   * TODO maybe add some cases here
   */
  cout << scientific;
}

void printProgramHead() {
  cout <<endl;
  cout <<"#==============================================================================#"<<endl;
  cout <<"#"<<setw(78)<<" "<<"#"<< endl;
  cout <<"#"<<setw(78)<<"   SSSSS                                                    #########         #"<<endl;
  cout <<"#"<<setw(78)<<"   SS    EEEE RRR  EEEE N   N II TTTTTT YY  YY          ####################  #"<<endl;
  cout <<"#"<<setw(78)<<"   SS    E    R  R E    NN  N II   TT    YYYY        ######################## #"<<endl;
  cout <<"#"<<setw(78)<<"   SSSSS EEEE RRR  EEEE N N N II   TT     YY      #########             ##### #"<<endl;
  cout <<"#"<<setw(78)<<"      SS E    R R  E    N  NN II   TT     YY   #######       OOOOOOOO         #"<<endl;
  cout <<"#"<<setw(78)<<"      SS EEEE R  R EEEE N   N II   TT     YY ########      OOOOOOOOOOOO       #"<<endl;
  cout <<"#"<<setw(78)<<"   SSSSS                                  ########        OOOOOOOOOOOOOO      #"<<endl;
  cout <<"#"<<setw(78)<<"                                     ###########   OOO    OOOOOOOOOOOOOO      #"<<endl;
  cout <<"#"<<setw(78)<<" #########                   ###############      OOOOO   OOOOOOOOOOOOOO      #"<<endl;
  cout <<"#"<<setw(78)<<" ######################################      OO   OOOOO    OOOOOOOOOOOO       #"<<endl;
  cout <<"#"<<setw(78)<<"       #########################         O   OO    OOO       OOOOOOOO         #"<<endl;
  cout <<"#"<<setw(78)<<" "<<"#"<< endl;
  cout <<"#"<<setw(78)<<" "<<"#"<< endl;
  cout <<"#"<< setw(78) << left << centerHeadlines("Serenity")<<"#" << endl;
  cout <<"#"<<setw(78)<<" "<<"#"<< endl;
  cout <<"#"<< setw(78) << left << centerHeadlines("A quantum chemistry code")<<"#" << endl;
  cout <<"#"<< setw(78) << left << centerHeadlines("developed in the workgroup of Johannes Neugebauer")<<"#" << endl;
  cout <<"#"<< setw(78) << left << centerHeadlines("at the WWU Münster.")<<" #" << endl;
  cout <<"#"<<setw(78)<<" "<<"#"<< endl;
  cout <<"#==============================================================================#"<<endl;
  cout <<endl;
}

void printRunStartInfo(){
  // time the program run
  takeTime("the entire run");

  //get time
  time_t     now = time(0);
  struct tm  tstruct;
  char       dateAndTime[50];
  tstruct = *localtime(&now);
  strftime(dateAndTime, sizeof(dateAndTime), "%Y-%m-%d %X", &tstruct);

  /**
   * FIXME works only for Linux
   */
  printSmallCaption("Program started");
  cout <<endl;
  cout <<"    "<<"Time:   " <<dateAndTime<<endl;
  {
    string hostName;
    if (getenv("HOSTNAME")!=NULL) {
      hostName=getenv("HOSTNAME");
    } else {
      hostName="HOSTNAME UNKNOWN";
    }
    cout <<"    "<<"On:                " <<hostName<<endl;
  }

#ifdef _OPENMP
  cout <<"    "<<"Cores:             " <<omp_get_max_threads()<<endl;
#endif
  cout <<"    "<<"Soft Memory Limit: " << MemoryManager::getInstance()->getAvailableSystemMemory() / (1024*1024) <<
          " MB for dynamically stored data." << endl;
  cout <<"    "<<"By:                " <<getenv("USER")<<endl;
  cout <<endl;
}

void printRunEndInfo(){
  //get time
  time_t     now = time(0);
  struct tm  tstruct;
  char       dateAndTime[50];
  tstruct = *localtime(&now);
  strftime(dateAndTime, sizeof(dateAndTime), "%Y-%m-%d %X", &tstruct);

  cout <<endl;
  printSmallCaption("Final Timings");
  Timings::printTimes();
  cout <<endl;
  printSmallCaption("Program ended");
  cout <<"    "<<"Time:  " <<dateAndTime<<endl;
  {
    string hostName;
    if (getenv("HOSTNAME")!=NULL) {
      hostName=getenv("HOSTNAME");
    } else {
      hostName="HOSTNAME UNKNOWN";
    }
    cout <<"    "<<"On:    " <<hostName<<endl;
  }
  timeTaken(0,"the entire run");
  cout <<endl;
}

void printSectionTitle(const string text) {
  cout <<endl;
  cout <<"o------------------------------------------------------------------------------o"<<endl;
  cout << "|" << centerHeadlines(text) << "|" << endl;
  cout <<"o------------------------------------------------------------------------------o"<<endl;
  cout <<endl;
}

void printSubSectionTitle(const string text) {
  cout <<endl;
  cout <<center("------------------------------------------------------------")<<endl;
  cout << center(text) <<endl;
  cout <<center("------------------------------------------------------------")<<endl;
  cout <<endl;
}

void printSmallCaption(const string text) {
  cout <<"  "<<text<<":"<<endl;
  cout <<" ";
    for (unsigned int i=0; i<text.length()+3; ++i){
      cout << "-";
    }
  cout <<endl;
}

void printTableHead(const string text) {
  cout <<"  "<<text<<endl;
  cout <<" ";
    for (unsigned int i=0; i<text.length()+3; ++i){
      cout << "-";
    }
  cout <<endl;
}

void printBigCaption(const string text) {
  cout <<" o";
    for (unsigned int i=0; i<text.length()+3; ++i){
      cout << "-";
    }
  cout << "o"<<endl;
  cout <<" | "<<text<<": |"<<endl;
  cout <<" o";
    for (unsigned int i=0; i<text.length()+3; ++i){
      cout << "-";
    }
  cout << "o"<<endl;
  cout <<endl;
}

void print(const string text) {
  cout <<"    "<< text << endl;
}

string centerHeadlines(const string s) {
   string outpt;
   int l=s.length();
   int pos=(int)((78-l)/2);
   for(int i=0;i<pos;++i) outpt = outpt + " ";
   outpt = outpt +s;
   for(int i=0;i<pos+(78-l)%2;++i) outpt = outpt + " ";
   return outpt;
}

string center(const string s) {
  string outpt;
  int l=s.length();
  int pos=(int)((80-l)/2);
  for(int i=0;i<pos;++i) outpt = outpt + " ";
  outpt = outpt +s;
  for(int i=0;i<pos+(80-l)%2;++i) outpt = outpt + " ";
  return outpt;
}

} /* namespace Serenity */
