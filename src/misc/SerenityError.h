/**
 * @file   SerenityError.h
 * @author Jan Unsleber
 *
 * @date   Oct 25, 2017
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

#ifndef MISC_SERENITYERROR_H_
#define MISC_SERENITYERROR_H_

/* Include Serenity Internal Headers */
#include "io/FormattedOutput.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <iostream>

namespace Serenity {

using namespace std;

/**
 * @class SerenityError SerenityError.h
 * @brief A small class to handle Serenitys error output for thrown errors.
 */
class SerenityError: public std::exception {
public:

    /**
     *  @brief Constructor (C++ STL strings).
     *  @param message The error message.
     */
    explicit SerenityError(const std::string& message):
      _msg(message){
    }

    /**
     * @brief Destructor.
     */
    virtual ~SerenityError() throw (){
    }

    /**
     *  @brief  Returns a pointer to the (constant) error description.
     *  @return A pointer to a const char*. The underlying memory
     *          is in posession of the Exception object. Callers must
     *          not attempt to free the memory.
     */
    virtual const char* what() const throw (){
      //get time
      time_t     now = time(0);
      struct tm  tstruct;
      char       dateAndTime[50];
      tstruct = *localtime(&now);
      strftime(dateAndTime, sizeof(dateAndTime), "%Y-%m-%d %X", &tstruct);

      printSectionTitle("Serenity Crashed");
      cout <<endl;
      printSmallCaption("Error");
      cout <<"    "<<_msg<<std::endl;
      cout <<endl;
      printSmallCaption("Final Timings");
      Timings::printTimes();
      cout <<endl;
      printSmallCaption("Program Crashed");
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

      return _msg.c_str();
    }

protected:
    /// @brief Error message.
    std::string _msg;
};

} /* namespace Serenity */

#endif /* MISC_SERENITYERROR_H_ */
