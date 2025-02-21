/**
 * @file    MemoryManager.h
 * @author  Thomas Dresselhaus, Jan Unsleber
 *
 * @date    14. Juli 2016
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
#ifndef MEMORYMANAGER_H
#define MEMORYMANAGER_H
/* Include Std and External Headers */
#include <memory>
#include <mutex>

namespace Serenity {
/**
 * @class MemoryManager MemoryManager.h
 * @brief Manages the distribution of available 'dynamic' memory.
 *
 * Some routines can gain a significant speed-up when granted huge chunks of memory. In larger
 * calculations, not all of that data can be kept in memory anymore and the corresponding objects
 * will fall back to alternative more CPU intense but less memory-consuming ways of performing
 * their jobs. This class gives access to the stats on available memory for these purposes.
 *
 * The manager does now only give information for the classes that are about to generate huge chunks of data,
 * the idea to control the generation form here was abandoned not because it does not necessarily make sense,
 * but because the manager would have to manage every object created to have an actual feeling for the memory.
 * E.g. 100 density matrices were not be managed but ERIs for one of them were, it is easy to see that the
 * matrices themselves can easily be the bigger chunck of memory.
 *
 * At the moment the idea is to have a class that can tell if big chunks can still be added, and it does
 * so using the system info or a set limit for the memory.
 * The memory set, however does not account for overhead and memory of small objects, there is de-facto
 * no memory management. The program will need a certain amount of memory.
 * Each object should recognize its size on its own and take appropriate measures
 */
class MemoryManager {
 public:
  virtual ~MemoryManager() = default;
  /**
   * @returns the only ever existing instance of the MemoryManager (singleton)
   */
  static std::shared_ptr<MemoryManager> getInstance();

  /**
   *
   * @return The amount of memory used by the current process in total right now, value in bytes.
   */
  long long getSerenityMemoryUsage();
  /**
   *
   * @return The amount of memory used by the system right now, value in bytes.
   */
  long long getSystemMemoryUsage();
  /**
   * @brief Checks for free dynamic memory
   * @return The amount of memory available on the system for Serenity right now, value in bytes.
   */
  long long getAvailableSystemMemory();
  /**
   * @brief The maximum system RAM
   * @returns The maximum system memory available.
   */
  long long getSystemPhysicalMemorySize();

 private:
  MemoryManager();
  // Stop the compiler generating methods to copy the object
  MemoryManager(MemoryManager const& copy) = delete;            // Not Implemented
  MemoryManager& operator=(MemoryManager const& copy) = delete; // Not Implemented

  long long _softMaxiumMemory;
  std::mutex _lock;
};

} /* namespace Serenity */

#endif /* MEMORYMANAGER_H */
