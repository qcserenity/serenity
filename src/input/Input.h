/**
 * @file Input.h
 *
 * @date Jan 23, 2017
 * @author Jan Unsleber
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#ifndef INPUT_INPUT_H_
#define INPUT_INPUT_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"
#include "settings/Reflection.h"
#include "system/SystemController.h"
#include "tasks/ActiveSpaceSelectionTask.h"
#include "tasks/BasisSetTruncationTask.h"
#include "tasks/BrokenSymmetryTask.h"
#include "tasks/CoupledClusterTask.h"
#include "tasks/DFTEmbeddedLocalCorrelationTask.h"
#include "tasks/DOSCCTask.h"
#include "tasks/DispersionCorrectionTask.h"
#include "tasks/DummyTask.h"
#include "tasks/EDATask.h"
#include "tasks/ElectronTransferTask.h"
#include "tasks/ElectronicStructureCopyTask.h"
#include "tasks/EvaluateEnergyTask.h"
#include "tasks/ExportCavityTask.h"
#include "tasks/ExportGridTask.h"
#include "tasks/FCIDumpFileWriterTask.h"
#include "tasks/FDETask.h"
#include "tasks/FXDTask.h"
#include "tasks/FiniteFieldTask.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/GWTask.h"
#include "tasks/GeneralizedDOSTask.h"
#include "tasks/GeometryOptimizationTask.h"
#include "tasks/GradientTask.h"
#include "tasks/HessianTask.h"
#include "tasks/ImportCavityTask.h"
#include "tasks/LRSCFTask.h"
#include "tasks/LocalCorrelationTask.h"
#include "tasks/LocalizationTask.h"
#include "tasks/MP2Task.h"
#include "tasks/MultipoleMomentTask.h"
#include "tasks/OrbitalsIOTask.h"
#include "tasks/OrthogonalizationTask.h"
#include "tasks/PCMInteractionEnergyTask.h"
#include "tasks/PlotTask.h"
#include "tasks/PopAnalysisTask.h"
#include "tasks/QuasiRestrictedOrbitalsTask.h"
#include "tasks/ScfTask.h"
#include "tasks/SystemAdditionTask.h"
#include "tasks/SystemSplittingTask.h"
#include "tasks/TDEmbeddingTask.h"
#include "tasks/TSTask.h"
#include "tasks/TopDownStaticEmbeddingTask.h"
#include "tasks/VirtualOrbitalSpaceSelectionTask.h"
#include "tasks/WavefunctionEmbeddingTask.h"
#include "tasks/WriteIntegralsTask.h"

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
#define createTaskA(TASK, MODE, ACT_SYS)                                \
  if (MODE == Options::SCF_MODES::RESTRICTED) {                         \
    auto taskptr = new TASK<Options::SCF_MODES::RESTRICTED>(ACT_SYS);   \
    Input::parseTaskSettings(taskptr, settingsStream);                  \
    task.reset(taskptr);                                                \
  }                                                                     \
  else {                                                                \
    auto taskptr = new TASK<Options::SCF_MODES::UNRESTRICTED>(ACT_SYS); \
    Input::parseTaskSettings(taskptr, settingsStream);                  \
    task.reset(taskptr);                                                \
  }
/*
 *
 * @brief Macro for tasks which can be restricted or unrestricted
 * @param TASK task
 * @param MODE scf mode (restricted or unrestricted)
 * @param ACT_SYS Active system(s)
 * @param ENV_SYS Environment system(s)
 */
#define createTaskAE(TASK, MODE, ACT_SYS, ENV_SYS)                               \
  if (MODE == Options::SCF_MODES::UNRESTRICTED) {                                \
    auto taskptr = new TASK<Options::SCF_MODES::UNRESTRICTED>(ACT_SYS, ENV_SYS); \
    Input::parseTaskSettings(taskptr, settingsStream);                           \
    task.reset(taskptr);                                                         \
  }                                                                              \
  else {                                                                         \
    auto taskptr = new TASK<Options::SCF_MODES::RESTRICTED>(ACT_SYS, ENV_SYS);   \
    Input::parseTaskSettings(taskptr, settingsStream);                           \
    task.reset(taskptr);                                                         \
  }
/*
 *
 * @brief Macro for tasks which can be restricted or unrestricted
 * @param TASK task
 * @param MODE scf mode (restricted or unrestricted)
 * @param SUPER_SYS supersystem(s)
 * @param ACT_SYS Active system(s)
 * @param ENV_SYS Environment system(s)
 */
#define createTaskSAE(TASK, MODE, SUPER_SYS, ACT_SYS, ENV_SYS)                              \
  if (MODE == Options::SCF_MODES::UNRESTRICTED) {                                           \
    auto taskptr = new TASK<Options::SCF_MODES::UNRESTRICTED>(SUPER_SYS, ACT_SYS, ENV_SYS); \
    Input::parseTaskSettings(taskptr, settingsStream);                                      \
    task.reset(taskptr);                                                                    \
  }                                                                                         \
  else {                                                                                    \
    auto taskptr = new TASK<Options::SCF_MODES::RESTRICTED>(SUPER_SYS, ACT_SYS, ENV_SYS);   \
    Input::parseTaskSettings(taskptr, settingsStream);                                      \
    task.reset(taskptr);                                                                    \
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
  template<class T>
  static void parseTaskSettings(T& task, std::stringstream& stream) {
    std::string line;
    std::string word;
    while (getline(stream, line)) {
      std::istringstream iss(line);
      word = "";
      iss >> word;

      for (auto& c : word)
        c = std::toupper(c);
      // end of settings block
      if (!word.compare("-TASK")) {
        break;
      }
      // check for comments and empty lines
      if (word.empty())
        continue;
      if (word[0] == '#')
        continue;
      // sub-block options
      if (word[0] == '+') {
        std::string blockname = word.erase(0, 1);
        for (auto& c : blockname)
          c = std::toupper(c);
        if (blockname == "TASK")
          assert(false);
        while (getline(stream, line)) {
          std::istringstream iss2(line);
          word = "";
          iss2 >> word;
          if (word[0] == '-') {
            break;
          }
          // check for comments
          if (word[0] == '#')
            continue;
          if (word.empty())
            continue;
          std::string name = word;
          std::string value;
          if (!(iss2 >> word)) {
            throw SerenityError("ERROR: Value missing for keyword '" + name + "'.");
          }
          value = word;
          while (iss2 >> word) {
            value += " ";
            value += word;
          }

          // set the corresponding fields for the reflection (see Reflection.h)
          bool check = false;

          set_visitor visitor(name, value, check);
          task->visit(task->settings, visitor, blockname);
          if (!check) {
            throw SerenityError("ERROR: No keyword '" + name + "' known in this block.");
          }
        }
        continue;
      }

      // unblocked options
      std::string name = word;

      std::string value;

      while (iss >> word) {
        value += word;
        value += " ";
      }

      if (value.empty()) {
        throw SerenityError("ERROR: Value missing for keyword '" + name + "'.");
      }

      // Get rid of last character (whitespace)
      value.pop_back();

      // set the corresponding fields for the reflection (see Reflection.h)
      bool check = false;

      set_visitor visitor(name, value, check);
      task->visit(task->settings, visitor, "");
      if (!check) {
        set_visitor visitorGeneral(name, value, check);
        visit_each(task->generalSettings, visitorGeneral);
        if (!check)
          throw SerenityError("ERROR: No keyword '" + name + "' known in this block.");
      }
    }
  }

  /**
   * @brief Determines which SCFMode a task should be run in.
   * @return The new SCFMode.
   */
  static Options::SCF_MODES determineSCFMode(std::vector<std::shared_ptr<SystemController>> activeSystems = {},
                                             std::vector<std::shared_ptr<SystemController>> environmentSystems = {},
                                             std::vector<std::shared_ptr<SystemController>> superSystems = {}) {
    Options::SCF_MODES runMode = Options::SCF_MODES::RESTRICTED;
    for (auto system : activeSystems)
      runMode = (system->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) ? Options::SCF_MODES::UNRESTRICTED : runMode;
    for (auto system : environmentSystems)
      runMode = (system->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) ? Options::SCF_MODES::UNRESTRICTED : runMode;
    for (auto system : superSystems)
      runMode = (system->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) ? Options::SCF_MODES::UNRESTRICTED : runMode;
    return runMode;
  }

  /**
   *
   * @param systems systems
   * @param name name of the task
   * @param stream input stream
   * @return task or nullpointer if the name couldn't be resolved
   */
  static std::unique_ptr<Task> parseTask(std::map<std::string, std::shared_ptr<SystemController>>& systems,
                                         std::string name, std::ifstream& stream) {
    if (!name.compare("")) {
      throw SerenityError("ERROR: No task name given.");
    }
    std::unique_ptr<Task> task(nullptr);
    std::string line;
    std::stringstream settingsStream;
    std::vector<std::shared_ptr<SystemController>> activeSystem;
    std::vector<std::shared_ptr<SystemController>> environmentSystem;
    std::vector<std::shared_ptr<SystemController>> supersystem;
    while (getline(stream, line)) {
      std::string word;
      std::istringstream iss(line);
      word = "";
      iss >> word;
      for (auto& c : word)
        c = std::toupper(c);
      // parse systems
      if (word.empty()) {
        continue;
      }
      else if (word[0] == '#') {
        continue;
      }
      else if (!word.compare("ACT") or !word.compare("ACTIVE") or !word.compare("SYSTEM")) {
        iss >> word;
        if (!word.compare("")) {
          throw SerenityError("ERROR: No system name given for an active system.");
        }
        try {
          activeSystem.push_back(systems.at(word));
        }
        catch (...) {
          throw SerenityError("ERROR: No system known with name: '" + word + "'.");
        }
        continue;
      }
      else if (!word.compare("ENV") or !word.compare("ENVIRONMENT")) {
        iss >> word;
        if (!word.compare("")) {
          throw SerenityError("ERROR: No system name given for an environment system.");
        }
        try {
          environmentSystem.push_back(systems.at(word));
        }
        catch (...) {
          throw SerenityError("ERROR: No system known with name: '" + word + "'.");
        }
        continue;
      }
      else if (!word.compare("SUP") or !word.compare("SUPER")) {
        iss >> word;
        if (!word.compare("")) {
          throw SerenityError("ERROR: No system name given for a supersystem.");
        }
        try {
          supersystem.push_back(systems.at(word));
        }
        catch (...) {
          throw SerenityError("ERROR: No system known with name: '" + word + "'.");
        }
        continue;
      }
      else if (!word.compare("-TASK")) {
        std::string copy = name;
        for (auto& c : copy)
          c = std::toupper(c);
        // Task block finished...
        // ...check for systems...
        if (activeSystem.size() == 0) {
          if ((copy.compare("ACTIVESPACETASK") && copy.compare("AS") && copy.compare("ACTIVESPACESELECTIONTASK")) or
              supersystem.size() == 0)
            throw SerenityError("ERROR: No active system given for task " + name + ".");
        }
        // ...create task using the systems...
        // ...then parse and adjust them ...
        if (!copy.compare("BROKENSYMMETRYTASK") or !copy.compare("BS")) {
          auto taskptr = new BrokenSymmetryTask(activeSystem, environmentSystem);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("BASISSETTASK") or !copy.compare("BST") or !copy.compare("BASISSETTRUNCATIONTASK")) {
          createTaskA(BasisSetTruncationTask, activeSystem[0]->getSCFMode(), activeSystem[0]);
        }
        else if (!copy.compare("COUPLEDCLUSTERTASK") or !copy.compare("CC")) {
          auto taskptr = new CoupledClusterTask(activeSystem[0], environmentSystem);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("PLOTTASK") or !copy.compare("PLOT") or !copy.compare("CUBEFILETASK") or
                 !copy.compare("CUBEFILE") or !copy.compare("CUBE")) {
          createTaskAE(PlotTask, determineSCFMode(activeSystem, environmentSystem), activeSystem, environmentSystem);
        }
        else if (!copy.compare("DISPERSIONCORRECTIONTASK") or !copy.compare("DISPERSION") or !copy.compare("DISP")) {
          auto taskptr = new DispersionCorrectionTask(activeSystem[0]);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("DFTEMB") or !copy.compare("DFTEMBEDDING") or !copy.compare("DFTEMBLC") or
                 !copy.compare("DFTEMBEDDEDLOCALCORRELATIONTASK")) {
          auto taskptr = new DFTEmbeddedLocalCorrelationTask(activeSystem[0], environmentSystem,
                                                             (supersystem.size() < 1) ? nullptr : supersystem[0]);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("DOSCC") or !copy.compare("DOSCCTASK")) {
          auto taskptr = new DOSCCTask(activeSystem);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("DUMMYTASK") or !copy.compare("DUMMY")) {
          auto taskptr = new DummyTask(activeSystem);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("ElectronTransferTask") or !copy.compare("ET") or !copy.compare("FDEET")) {
          auto taskptr = new ElectronTransferTask(activeSystem);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("ENERGYTASK") or !copy.compare("ENERGY") or !copy.compare("E") or
                 !copy.compare("EVALUATEENERGYTASK")) {
          createTaskA(EvaluateEnergyTask, determineSCFMode(activeSystem), activeSystem);
        }
        else if (!copy.compare("EXPORTGRIDTASK") or !copy.compare("GRIDTASK") or !copy.compare("GRID")) {
          auto taskptr = new ExportGridTask(activeSystem[0]);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("EDATASK") or !copy.compare("EDA")) {
          createTaskA(EDATask, determineSCFMode(activeSystem), activeSystem);
        }
        else if (!copy.compare("FCIDUMP") or !copy.compare("FCIDUMPFILEWRITERTASK")) {
          createTaskAE(FCIDumpFileWriterTask, activeSystem[0]->getSCFMode(), activeSystem[0], environmentSystem);
        }
        else if (!copy.compare("FDETASK") or !copy.compare("FDE")) {
          createTaskAE(FDETask, activeSystem[0]->getSCFMode(), activeSystem[0], environmentSystem);
        }
        else if (!copy.compare("FXDTASK") or !copy.compare("FXD")) {
          createTaskA(FXDTask, activeSystem[0]->getSCFMode(), activeSystem[0]);
        }
        else if (!copy.compare("FINITEFIELD") or !copy.compare("FF") or !copy.compare("FINITEFIELDTASK")) {
          createTaskAE(FiniteFieldTask, determineSCFMode({activeSystem[0]}, environmentSystem), activeSystem[0],
                       environmentSystem);
        }
        else if (!copy.compare("FREEZEANDTHAWTASK") or !copy.compare("FAT")) {
          createTaskAE(FreezeAndThawTask, determineSCFMode(activeSystem, environmentSystem), activeSystem, environmentSystem);
        }
        else if (!copy.compare("TDSTATIC") or !copy.compare("STATIC") or !copy.compare("TOPDOWNSTATICEMBEDDINGTASK")) {
          Options::SCF_MODES scfMode = Options::SCF_MODES::RESTRICTED;
          for (auto system : activeSystem) {
            if (system->getSCFMode() == Options::SCF_MODES::UNRESTRICTED)
              scfMode = Options::SCF_MODES::UNRESTRICTED;
          }
          createTaskAE(TopDownStaticEmbeddingTask, scfMode, activeSystem, environmentSystem);
        }
        else if (!copy.compare("GDOS") or !copy.compare("GENERALIZEDDOS") or !copy.compare("GENERALIZEDDOSTASK")) {
          createTaskAE(GeneralizedDOSTask, determineSCFMode(activeSystem, environmentSystem), activeSystem, environmentSystem);
        }
        else if (!copy.compare("GEOMETRYOPTIMIZATIONTASK") or !copy.compare("GEOOPT") or !copy.compare("OPT")) {
          createTaskAE(GeometryOptimizationTask, determineSCFMode(activeSystem, environmentSystem), activeSystem,
                       environmentSystem);
        }
        else if (!copy.compare("GRADIENTTASK") or !copy.compare("GRADIENT") or !copy.compare("GRAD")) {
          createTaskAE(GradientTask, determineSCFMode(activeSystem, environmentSystem), activeSystem, environmentSystem);
        }
        else if (!copy.compare("HESSIANTASK") or !copy.compare("HESSIAN") or !copy.compare("HESS")) {
          createTaskAE(HessianTask, determineSCFMode(activeSystem, environmentSystem), activeSystem, environmentSystem);
        }
        else if (!copy.compare("LOCALCORRELATIONTASK") or !copy.compare("LC")) {
          auto taskptr = new LocalCorrelationTask(activeSystem[0], environmentSystem);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("LOCALIZATIONTASK") or !copy.compare("LOCALIZATION") or !copy.compare("LOC")) {
          auto taskptr = new LocalizationTask(activeSystem[0], environmentSystem);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("LRSCFTASK") or !copy.compare("LRSCF") or !copy.compare("RICC2")) {
          createTaskAE(LRSCFTask, determineSCFMode(activeSystem, environmentSystem), activeSystem, environmentSystem);
        }
        else if (!copy.compare("GWTASK") or !copy.compare("GW")) {
          createTaskAE(GWTask, determineSCFMode(activeSystem, environmentSystem), activeSystem, environmentSystem);
        }
        else if (!copy.compare("MP2TASK") or !copy.compare("MP2")) {
          createTaskAE(MP2Task, determineSCFMode({activeSystem[0]}, environmentSystem), activeSystem[0], environmentSystem);
        }
        else if (!copy.compare("MULTIPOLEMOMENTTASK") or !copy.compare("MULTI")) {
          auto taskptr = new MultipoleMomentTask(activeSystem);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("POPULATIONANALYSISTASK") or !copy.compare("POPULATION") or !copy.compare("POP") or
                 !copy.compare("POPANALYSISTASK")) {
          createTaskA(PopAnalysisTask, activeSystem[0]->getSCFMode(), activeSystem[0]);
        }
        else if (!copy.compare("PCMINTERACTIONENERGYTASK") or !copy.compare("PCM")) {
          auto taskptr = new PCMInteractionEnergyTask(activeSystem[0]);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("PROJECTIONBASEDEMBTASK") or !copy.compare("PBE") or !copy.compare("TDEMBEDDINGTASK") or
                 !copy.compare("TD")) {
          if (environmentSystem.size() < 1)
            throw SerenityError("ERROR: At least one environment system is needed!");
          createTaskAE(TDEmbeddingTask, determineSCFMode({activeSystem[0]}, {environmentSystem[0]}), activeSystem[0],
                       environmentSystem[0]);
        }
        else if (!copy.compare("QRO") or !copy.compare("QUASIRESTRICTEDORBITALS") or
                 !copy.compare("QUASIRESTRICTEDORBITALSTASK")) {
          createTaskAE(QuasiRestrictedOrbitalsTask, activeSystem[0]->getSCFMode(), activeSystem[0], environmentSystem);
        }
        else if (!copy.compare("READORBS") or !copy.compare("READ") or !copy.compare("WRITEORBS") or
                 !copy.compare("WRITE") or !copy.compare("IO") or !copy.compare("ORBITALSIOTASK")) {
          createTaskA(OrbitalsIOTask, activeSystem[0]->getSCFMode(), activeSystem[0]);
        }
        else if (!copy.compare("SCFTASK") or !copy.compare("SCF")) {
          createTaskA(ScfTask, activeSystem[0]->getSCFMode(), activeSystem[0]);
        }
        else if (!copy.compare("ORTHOGONALIZATIONTASK") or !copy.compare("ORTHO")) {
          if (supersystem.size() > 0) {
            createTaskAE(OrthogonalizationTask, determineSCFMode(activeSystem, {supersystem[0]}), activeSystem,
                         supersystem[0]);
          }
          else {
            createTaskAE(OrthogonalizationTask, determineSCFMode(activeSystem), activeSystem, nullptr);
          }
        }
        else if (!copy.compare("ACTIVESPACETASK") or !copy.compare("AS") or !copy.compare("ACTIVESPACESELECTIONTASK")) {
          createTaskSAE(ActiveSpaceSelectionTask, determineSCFMode(activeSystem, environmentSystem, supersystem),
                        supersystem, activeSystem, environmentSystem);
        }
        else if (!copy.compare("VIRTUALORBITALSPACESELECTIONTASK") or !copy.compare("VOSSTASK") or !copy.compare("VOSS")) {
          createTaskAE(VirtualOrbitalSpaceSelectionTask, determineSCFMode(activeSystem, environmentSystem),
                       activeSystem, environmentSystem);
        }
        else if (!copy.compare("SPLIT") or !copy.compare("SPLITTINGTASK") or !copy.compare("SYSTEMSPLITTINGTASK")) {
          createTaskAE(SystemSplittingTask, determineSCFMode({activeSystem[0]}, environmentSystem), activeSystem[0],
                       environmentSystem);
        }
        else if (!copy.compare("ADD") or !copy.compare("ADDITIONTASK") or !copy.compare("SYSTEMADDITIONTASK")) {
          createTaskAE(SystemAdditionTask, determineSCFMode({activeSystem[0]}, environmentSystem), activeSystem[0],
                       environmentSystem);
        }
        else if (!copy.compare("COPY") or !copy.compare("ESC") or !copy.compare("ELECTRONICSTRUCTURECOPYTASK")) {
          createTaskAE(ElectronicStructureCopyTask, determineSCFMode({activeSystem[0]}), activeSystem[0], environmentSystem);
        }
        else if (!copy.compare("EXPORTCAVITY") or !copy.compare("EXPORTCAVITYTASK")) {
          auto taskptr = new ExportCavityTask(activeSystem[0]);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("IMPORTCAVITY") or !copy.compare("IMPORTCAVITYTASK")) {
          auto taskptr = new ImportCavityTask(activeSystem[0]);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("TS") or !copy.compare("TSOPT") or !copy.compare("TSTASK")) {
          auto taskptr = new TSTask(activeSystem[0], environmentSystem);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("WFEMB") or !copy.compare("WAVEFUNCTIONEMBEDDING") or
                 !copy.compare("WAVEFUNCTIONEMBEDDINGTASK")) {
          auto taskptr = new WavefunctionEmbeddingTask(activeSystem[0], environmentSystem);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else if (!copy.compare("INT") or !copy.compare("WRITEINTS") or !copy.compare("INTEGRALS") or
                 !copy.compare("WRITEINTEGRALSTASK")) {
          auto taskptr = new WriteIntegralsTask(activeSystem[0]);
          Input::parseTaskSettings(taskptr, settingsStream);
          task.reset(taskptr);
        }
        else {
          throw SerenityError("ERROR: No task known with name: '" + name + "'.");
        }
        return task;
      }
      else {
        // If none of the above:
        // Must be settings keywords...
        // Save in stringstream for parseTaskSettings function
        settingsStream << word << ' ';

        std::string value, valueCopy;
        iss >> value;
        /*
         * If first character is "{", where dealing with a vector:
         * Save everything until "}" in value
         */
        if (value[0] == '{') { // This checks and reads input in the form of '{val1 val2 val3}'
          std::string vector;  // Create temporary string to validate the vector input
          vector.append(value);
          while (iss >> value) {
            vector.append(" " + value); // Get the remaining part of the line
          }
          if (vector.back() != '}') // Check for closing bracket
            throw SerenityError("Missing closing bracket in the input. The line is: " + line);
          else
            vector.pop_back();                                                                // Remove closing bracket
          vector.erase(0, 1);                                                                 // Remove leading bracket
          if (vector.find('}') != std::string::npos or vector.find('{') != std::string::npos) // Check for additional
                                                                                              // brackets in the line
            throw SerenityError("Wrong bracket in the input. The line is: " + line);
          std::istringstream local_iss(vector); // Create istreamstring to pass the vector to the settingsstream
          while (local_iss >> value) {
            settingsStream << value << ' ';
          }
          settingsStream << std::endl;
        }
        else {
          settingsStream << value << std::endl;
        }
      }
    }
    throw SerenityError("ERROR: Could not create task " + name + ". Please check the input file.");
    return nullptr;
  }
};

} /* namespace Serenity */

#endif /* INPUT_INPUT_H_ */
