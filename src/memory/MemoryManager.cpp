/**
 * @file MemoryManager.cpp
 * @author Thomas Dresselhaus
 *
 * @date Juli 14, 2014
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

/* Include Class Header*/
#include "memory/MemoryManager.h"
/* Include Std and External Headers */
#include <stdio.h>
#include <string.h>
#include <cassert>
#include <cstdlib>
#include <sstream>
#if __APPLE__ || __MACH__
#include <mach/mach_host.h>
#include <mach/mach_vm.h>
#include <sys/sysctl.h>
#include <unistd.h>
#elif __linux__ || __unix__ || __unix
#include <linux/sysinfo.h>
#include <sys/sysinfo.h>
#endif

namespace Serenity {

MemoryManager::MemoryManager() : _softMaxiumMemory(-1) {
  if (getenv("SERENITY_MEMORY") != NULL) {
    std::stringstream convert;
    convert << getenv("SERENITY_MEMORY");
    convert >> _softMaxiumMemory;
    _softMaxiumMemory *= 1024 * 1024;
    assert(_softMaxiumMemory < getSystemPhysicalMemorySize());
    if (_softMaxiumMemory == 0)
      _softMaxiumMemory = -1;
  }
}

std::shared_ptr<MemoryManager> MemoryManager::getInstance() {
  // The only instance
  // Guaranteed to be lazy initialized
  // Guaranteed that it will be destroyed correctly
  static std::shared_ptr<MemoryManager> instance(new MemoryManager());
  return instance;
}

#if __APPLE__ || __MACH__
long long MemoryManager::getSerenityMemoryUsage() {
  assert(false && "Deprecated on Mac OS");
  return 0;
}
long long MemoryManager::getSystemMemoryUsage() {
  assert(false && "Deprecated on Mac OS");
  return 0;
}
long long MemoryManager::getSystemPhysicalMemorySize() {
  assert(false && "Deprecated on Mac OS");
  return 0;
}
long long MemoryManager::getAvailableSystemMemory() {
  mach_msg_type_number_t count = HOST_VM_INFO_COUNT;
  vm_statistics_data_t vmstat;
  if (KERN_SUCCESS != host_statistics(mach_host_self(), HOST_VM_INFO, (host_info_t)&vmstat, &count)) {
    assert(false && "Error while reading RAM infomation!");
  }
  long long freePages = vmstat.free_count + vmstat.inactive_count; // also accounts for cached files
  return freePages * getpagesize();
}

#elif __linux__ || __unix__ || __unix
long long MemoryManager::getSerenityMemoryUsage() {
  // Will lock here and unlock on return
  std::lock_guard<std::mutex> lock(_lock);
  /*
   * This works for linux only
   */
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];

  while (fgets(line, 128, file) != NULL) {
    if (strncmp(line, "VmRSS:", 6) == 0) {
      // This assumes that a digit will be found and the line ends in " Kb".
      int i = strlen(line);
      const char* p = line;
      while (*p < '0' || *p > '9')
        p++;
      line[i - 3] = '\0';
      result = atoi(p);
      break;
    }
  }
  fclose(file);
  return result * 1024;
}
long long MemoryManager::getSystemMemoryUsage() {
  struct sysinfo memInfo;

  char buf[60]; /* actual lines we expect are ~30 chars or less */
  unsigned long cached = 0;
  FILE* file = fopen("/proc/meminfo", "r");
  while (fgets(buf, sizeof(buf), file) != NULL) {
    if (sscanf(buf, "Cached: %lu %*s\n", &cached) == 1)
      break;
  }
  fclose(file);

  sysinfo(&memInfo);
  long long physMemUsed = memInfo.totalram - memInfo.freeram - memInfo.bufferram - cached * 1024;
  // Multiply in next statement to avoid int overflow on right hand side...
  physMemUsed *= memInfo.mem_unit;
  return physMemUsed;
}

long long MemoryManager::getSystemPhysicalMemorySize() {
  struct sysinfo memInfo;

  sysinfo(&memInfo);
  long long totalPhysMem = memInfo.totalram;
  // Multiply in next statement to avoid int overflow on right hand side...
  totalPhysMem *= memInfo.mem_unit;
  return totalPhysMem;
}

long long MemoryManager::getAvailableSystemMemory() {
  long long freeMem;
  if (_softMaxiumMemory == -1) {
    freeMem = getSystemPhysicalMemorySize() - getSystemMemoryUsage();
  }
  else {
    freeMem = _softMaxiumMemory - getSerenityMemoryUsage();
  }
  return freeMem;
}
#endif
} /* namespace Serenity */
