/**
 * @file Input.h
 *
 * @date Jan 23, 2017
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

#ifndef INPUT_INPUT_H_
#define INPUT_INPUT_H_

/* Include Serenity Internal Headers */
#include "tasks/BasisSetTruncationTask.h"
#include "tasks/CoupledClusterTask.h"
#include "tasks/CubeFileTask.h"
#include "tasks/DispersionCorrectionTask.h"
#include "tasks/DummyTask.h"
#include "tasks/EDATask.h"
#include "tasks/ExportGridTask.h"
#include "tasks/FDETask.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/GeometryOptimizationTask.h"
#include "tasks/GradientTask.h"
#include "tasks/HessianTask.h"
#include "tasks/LRSCFTask.h"
#include "tasks/LocalizationTask.h"
#include "tasks/MP2Task.h"
#include "tasks/MultipoleMomentTask.h"
#include "settings/Options.h"
#include "tasks/PopAnalysisTask.h"
#include "tasks/TDEmbeddingTask.h"
#include "settings/Reflection.h"
#include "tasks/ScfTask.h"
#include "tasks/TSTask.h"
/* Include Std and External Headers */
#include <sstream>


namespace Serenity {

using namespace Serenity::Reflection;
/*
 *
 * @brief Macro for tasks which can be restricted or unrestricted
 * @param TASK task
 * @param MODE scf mode (restricted or unrestricted)
 * @param ACT_SYS Active system(s)
 */
#define createTaskA(TASK,MODE,ACT_SYS) \
  if(MODE==Options::SCF_MODES::RESTRICTED){ \
    auto taskptr = new TASK<Options::SCF_MODES::RESTRICTED>(ACT_SYS); \
    Input::parseTaskSettings(taskptr,settingsStream); \
    task.reset(taskptr); \
  } else { \
    auto taskptr = new TASK<Options::SCF_MODES::UNRESTRICTED>(ACT_SYS); \
    Input::parseTaskSettings(taskptr,settingsStream); \
    task.reset(taskptr); \
  }
/*
 *
 * @brief Macro for tasks which can be restricted or unrestricted
 * @param TASK task
 * @param MODE scf mode (restricted or unrestricted)
 * @param ACT_SYS Active system(s)
 * @param ENV_SYS Environment system(s)
 */
#define createTaskAE(TASK,MODE,ACT_SYS,ENV_SYS) \
    if(MODE==Options::SCF_MODES::UNRESTRICTED){ \
      auto taskptr = new TASK<Options::SCF_MODES::UNRESTRICTED>(ACT_SYS,ENV_SYS); \
      Input::parseTaskSettings(taskptr,settingsStream); \
      task.reset(taskptr); \
    } else { \
      auto taskptr = new TASK<Options::SCF_MODES::RESTRICTED>(ACT_SYS,ENV_SYS); \
      Input::parseTaskSettings(taskptr,settingsStream); \
      task.reset(taskptr); \
    }
/**
 * @class Input Input.h
 * @brief Generating tasks from the input stream
 */
class Input {
public:
  /**
   *
   * @param task task
   * @param stream input stream
   */
  template <class T>
  static void parseTaskSettings(T& task,std::stringstream& stream){
    std::string line;
    std::string word;
    while(getline(stream, line)){
      std::istringstream iss(line);
      word = "";
      iss >> word;

      for (auto& c: word ) c = std::toupper(c);
      // end of settings block
      if (!word.compare("-TASK")){
        break;
      }
      // check for comments and empty lines
      if (word.empty()) continue;
      if (word[0] == '#') continue;
      if (word[0] == '+') {
        std::string blockname = word.erase(0, 1);
        for (auto& c: blockname) c = std::toupper(c);
        if(blockname=="TASK") assert(false);
        while(getline(stream, line)){
          std::istringstream iss2(line);
          word ="";
          iss2 >> word;
          if (word[0] == '-') {
            break;
          }
          // check for comments
          if (word[0] == '#') continue;
          if (word.empty()) continue;
          std::string name = word;
          std::string value;
          if (!(iss2 >> word)){
            throw SerenityError("ERROR: Value missing for keyword '"+name+"'.");
          }
          value = word;
          while(iss2 >> word){
            value+=word;
            value+=" ";
          }

          // set the corresponding fields for the reflection (see Reflection.h)
          bool check = false;

          set_visitor visitor(name,value,check);
          task->visit(task->settings,visitor,blockname);
          if (!check){
            throw SerenityError("ERROR: No keyword '"+name+"' known in this block.");
          }
        }
        continue;
      }

      // unblocked options
      std::string name = word;

      std::string value;

      while(iss >> word){
        value+=word;
        value+=" ";
      }

      if (value.empty()){
        throw SerenityError("ERROR: Value missing for keyword '"+name+"'.");
      }

      // Get rid of last character (whitespace)
      value.pop_back();

      // set the corresponding fields for the reflection (see Reflection.h)
      bool check = false;

      set_visitor visitor(name,value,check);
      task->visit(task->settings,visitor,"");
      if (!check){
        throw SerenityError("ERROR: No keyword '"+name+"' known in this block.");
      }
    }
  }


/**
 *
 * @param systems systems
 * @param name name of the task
 * @param stream input stream
 * @return task or nullpointer if the name couldn't be resolved
 */
  static std::unique_ptr<Task> parseTask(
      std::map<std::string,std::shared_ptr<SystemController>>& systems,
      std::string name,
      std::ifstream& stream){
    if (!name.compare("")) {
      throw SerenityError("ERROR: No task name given.");
    }
    std::unique_ptr<Task> task(nullptr);
    std::string line;
    std::stringstream settingsStream;
    std::vector<std::shared_ptr<SystemController> > activeSystem;
    std::vector<std::shared_ptr<SystemController> > environmentSystem;
    while(getline(stream, line)){
      std::string word;
      std::istringstream iss(line);
      word = "";
      iss >> word;
      bool defaultSettings =false;
      for (auto& c: word) c = std::toupper(c);
      // parse systems
      if (word.empty()){
        continue;
      } else if (word[0] == '#'){
        continue;
      } else if (!word.compare("ACT")
          or !word.compare("ACTIVE")
          or !word.compare("SYSTEM")){
        iss >> word;
        if (!word.compare("")) {
          throw SerenityError("ERROR: No system name given for an active system.");
        }
        try{
          activeSystem.push_back(systems.at(word));
        }catch(...){
          throw SerenityError("ERROR: No system known with name: '"+word+"'.");
        }
        continue;
      } else if (!word.compare("ENV")
          or !word.compare("ENVIRONMENT")){
        iss >> word;
        if (!word.compare("")) {
          throw SerenityError("ERROR: No system name given for an environment system.");
        }
        try{
          environmentSystem.push_back(systems.at(word));
        }catch(...){
          throw SerenityError("ERROR: No system known with name: '"+word+"'.");
        }
        continue;
      } else if(!word.compare("-TASK")){
        // Task block finished...
        // ...check for systems...
        if(activeSystem.size()==0){
          throw SerenityError("ERROR: No active system given for task "+name+".");
        }
        // ...create task using the systems...
        // ...then parse and adjust them ...
        std::string copy = name;
        for (auto& c: copy) c = std::toupper(c);
        if (!copy.compare("BASISSETTASK") or !copy.compare("BST")){
                  createTaskA(BasisSetTruncationTask,activeSystem[0]->getSettings().scfMode,activeSystem[0]);
        } else if (!copy.compare("COUPLEDCLUSTERTASK") or !copy.compare("CC")){
          auto taskptr = new CoupledClusterTask(activeSystem[0]);
          if (!defaultSettings) Input::parseTaskSettings(taskptr,settingsStream);
          task.reset(taskptr);
        } else if (!copy.compare("CUBEFILETASK") or !copy.compare("CUBEFILE") or !copy.compare("CUBE")){
          createTaskAE(CubeFileTask,activeSystem[0]->getSettings().scfMode,activeSystem,environmentSystem);
        } else if (!copy.compare("DISPERSIONCORRECTIONTASK") or !copy.compare("DISPERSION") or !copy.compare("DISP")){
          auto taskptr = new DispersionCorrectionTask(activeSystem[0]);
          if (!defaultSettings) Input::parseTaskSettings(taskptr,settingsStream);
          task.reset(taskptr);
        } else if (!copy.compare("DUMMYTASK") or !copy.compare("DUMMY")){
          auto taskptr = new DummyTask(activeSystem);
          if (!defaultSettings) Input::parseTaskSettings(taskptr,settingsStream);
          task.reset(taskptr);
        } else if (!copy.compare("EXPORTGRIDTASK") or !copy.compare("GRIDTASK") or !copy.compare("GRID")){
          auto taskptr = new ExportGridTask(activeSystem[0]);
          if (!defaultSettings) Input::parseTaskSettings(taskptr,settingsStream);
          task.reset(taskptr);
        } else if (!copy.compare("EDATASK") or !copy.compare("EDA")){
          Options::SCF_MODES scfMode=RESTRICTED;
          for(auto system : activeSystem){
            if(system->getSettings().scfMode==UNRESTRICTED) scfMode=UNRESTRICTED;
          }
          createTaskA(EDATask,scfMode,activeSystem);
        } else if (!copy.compare("FDETASK") or !copy.compare("FDE")){
          createTaskAE(FDETask,activeSystem[0]->getSettings().scfMode,activeSystem[0],environmentSystem);
        } else if (!copy.compare("FREEZEANDTHAWTASK") or !copy.compare("FAT")){
          Options::SCF_MODES scfMode=Options::SCF_MODES::RESTRICTED;
          for(auto system : activeSystem){
            if(system->getSettings().scfMode==Options::SCF_MODES::UNRESTRICTED)
              scfMode=Options::SCF_MODES::UNRESTRICTED;
          }
          createTaskAE(FreezeAndThawTask,scfMode,activeSystem,environmentSystem);
        } else if (!copy.compare("GEOMETRYOPTIMIZATIONTASK") or !copy.compare("GEOOPT") or !copy.compare("OPT")){
          Options::SCF_MODES scfMode=Options::SCF_MODES::RESTRICTED;
          for(auto system : activeSystem){
            if(system->getSettings().scfMode==Options::SCF_MODES::UNRESTRICTED)
              scfMode=Options::SCF_MODES::UNRESTRICTED;
          }
          createTaskAE(GeometryOptimizationTask,scfMode,activeSystem,environmentSystem);
        } else if (!copy.compare("GRADIENTTASK") or !copy.compare("GRADIENT") or !copy.compare("GRAD")){
          Options::SCF_MODES scfMode=Options::SCF_MODES::RESTRICTED;
          for(auto system : activeSystem){
            if(system->getSettings().scfMode==Options::SCF_MODES::UNRESTRICTED)
              scfMode=Options::SCF_MODES::UNRESTRICTED;
          }
          if(scfMode==Options::SCF_MODES::UNRESTRICTED){
            auto taskptr = new GradientTask<Options::SCF_MODES::UNRESTRICTED>(activeSystem,environmentSystem);
            Input::parseTaskSettings(taskptr,settingsStream);
            task.reset(taskptr);
          } else {
            auto taskptr = new GradientTask<Options::SCF_MODES::RESTRICTED>(activeSystem,environmentSystem);
            Input::parseTaskSettings(taskptr,settingsStream);
            task.reset(taskptr);
          }
        } else if (!copy.compare("HESSIANTASK") or !copy.compare("HESSIAN") or !copy.compare("HESS")){
                auto taskptr = new HessianTask(activeSystem,environmentSystem);
                if (!defaultSettings) Input::parseTaskSettings(taskptr,settingsStream);
                task.reset(taskptr);
        } else if (!copy.compare("LOCALIZATIONTASK") or !copy.compare("LOCALIZATION") or !copy.compare("LOC")){
          auto taskptr = new LocalizationTask(activeSystem[0]);
          if (!defaultSettings) Input::parseTaskSettings(taskptr,settingsStream);
          task.reset(taskptr);
        } else if (!copy.compare("LRSCFTASK") or !copy.compare("LRSCF")){
          Options::SCF_MODES scfMode=Options::SCF_MODES::RESTRICTED;
          for(auto system : activeSystem){
             if(system->getSettings().scfMode==Options::SCF_MODES::UNRESTRICTED)
               scfMode=Options::SCF_MODES::UNRESTRICTED;
              }
          createTaskAE(LRSCFTask,scfMode,activeSystem,environmentSystem);
        } else if (!copy.compare("MP2TASK") or !copy.compare("MP2")){
          createTaskA(MP2Task,activeSystem[0]->getSettings().scfMode,activeSystem[0]);
        } else if (!copy.compare("MULTIPOLEMOMENTTASK") or !copy.compare("MULTI")){
          createTaskA(MultipoleMomentTask,activeSystem[0]->getSettings().scfMode,activeSystem[0]);
        } else if (!copy.compare("POPULATIONANALYSISTASK") or !copy.compare("POPULATION") or !copy.compare("POP")){
          createTaskA(PopAnalysisTask,activeSystem[0]->getSettings().scfMode,activeSystem[0]);
        } else if (!copy.compare("PROJECTIONBASEDEMBTASK") or !copy.compare("PBE") or !copy.compare("TDEMBEDDINGTASK") or !copy.compare("TD")){
          Options::SCF_MODES scfMode=Options::SCF_MODES::RESTRICTED;
          for(auto system : activeSystem){
            if(system->getSettings().scfMode==Options::SCF_MODES::UNRESTRICTED)
              scfMode=Options::SCF_MODES::UNRESTRICTED;
          }
          for(auto system : environmentSystem){
            if(system->getSettings().scfMode==Options::SCF_MODES::UNRESTRICTED)
              scfMode=Options::SCF_MODES::UNRESTRICTED;
          }
          createTaskAE(TDEmbeddingTask,scfMode,activeSystem[0],environmentSystem[0]);
        } else if (!copy.compare("SCFTASK") or !copy.compare("SCF")){
          createTaskA(ScfTask,activeSystem[0]->getSettings().scfMode,activeSystem[0]);
        } else if (!copy.compare("TS") or !copy.compare("TSOPT") or !copy.compare("TSTASK")){
          auto taskptr = new TSTask(activeSystem[0],environmentSystem);
          if (!defaultSettings) Input::parseTaskSettings(taskptr,settingsStream);
          task.reset(taskptr);
        } else {
          throw SerenityError("ERROR: No task known with name: '"+name+"'.");
        }
        return task;

      }else{
        // If none of the above:
        // Must be settings keywords...
        // Save in stringstream for parseTaskSettings function
        settingsStream << word << ' ';

        std::string value;
        iss >> value;
        /*
         * If first character is "{", where dealing with a vector:
         * Save everything until "}" in value
         */
        if(value[0]=='{'){
          value.erase(0, 1);
          /*
           * Handle "keyword {val}" case
           */
          if(value.back()=='}'){
            value.pop_back();
          /*
           * Every other case is covered by this loop
           */
          }else{
            settingsStream << value << ' ';
            while(true){
              iss >> value;
              if(value.back()=='}') break;
              settingsStream << value << ' ';
            }
            value.pop_back();
          }
          settingsStream << value << std::endl;
        /*
         * Else, fill single value into settingsStream
         */
        }else{
          settingsStream << value << std::endl;
        }
      }
    }
    throw SerenityError("ERROR: Could not create task "+name+". Please check the input file.");
    return nullptr;
  }
};

} /* namespace Serenity */


#endif /* INPUT_INPUT_H_ */
