/*
 * Memory usage info
 */

#include "pymol/memory.h"

#ifdef __linux__
#include <sys/sysinfo.h>
#include <sys/types.h>
#include <unistd.h>
#elif defined(__APPLE__)
#import <mach/mach.h>
#import <mach/mach_host.h>
#elif defined(_WIN32)
#include <windows.h>
// must come *after* windows.h
# include <psapi.h>
#endif

#include <cstdio>
#include <iostream>

namespace pymol
{

/**
 * Get the memory usage of the current process in bytes.
 *
 * May be inaccurate (e.g. on Linux this is the resident set size VmRSS)
 */
size_t memory_usage()
{
#ifdef _WEBGL
  return 0;
#elif defined(__linux__)
  size_t vmRSS = 0;
  if (auto fp = std::fopen("/proc/self/statm", "rb")) {
    std::fscanf(fp, "%*zu%zu", &vmRSS);
    std::fclose(fp);
  }
  return vmRSS * sysconf(_SC_PAGESIZE);
#elif defined(__APPLE__)
  task_vm_info_data_t data{};
  mach_msg_type_number_t data_count = TASK_VM_INFO_COUNT;
  task_info(mach_task_self(), TASK_VM_INFO, task_info_t(&data), &data_count);
  return data.phys_footprint;
#elif defined(_WIN32)
  PROCESS_MEMORY_COUNTERS pmc{};
  if (auto hProcess =
          unique_ptr_take_ownership(GetCurrentProcess(), CloseHandle)) {
    GetProcessMemoryInfo(hProcess.get(), &pmc, sizeof(pmc));
  }
  return pmc.WorkingSetSize;
#else
#error "platform not handled"
#endif
}

/**
 * Get the "available" memory in bytes.
 *
 * TODO: On macOS, this value is often very small and doesn't reflect the
 * memory pressure at all.
 */
size_t memory_available()
{
#ifdef _WEBGL
  return 0;
#elif defined(__linux__)
  size_t memAvail = 0;
  if (auto fp = unique_ptr_take_ownership(
          std::fopen("/proc/meminfo", "rb"), std::fclose)) {
    char buffer[80];
    while (std::fgets(buffer, sizeof(buffer), fp.get())) {
      if (std::sscanf(buffer, "MemAvailable: %zu", &memAvail) > 0) {
        return memAvail * 1000 /* kB */;
      }
    }
  }
  return 0;
#elif defined(__APPLE__)
  // TODO is this meaningful?
  mach_port_t const port = mach_host_self();
  vm_size_t pagesize = 0;
  host_page_size(port, &pagesize);
  vm_statistics_data_t data{};
  mach_msg_type_number_t data_count = sizeof(data) / sizeof(integer_t);
  host_statistics(port, HOST_VM_INFO, host_info_t(&data), &data_count);
  return data.free_count * pagesize;
#elif defined(_WIN32)
  MEMORYSTATUSEX statex{};
  statex.dwLength = sizeof(statex);
  GlobalMemoryStatusEx(&statex);
  return statex.ullAvailPhys;
#else
#error "platform not handled"
#endif
}

} // namespace pymol
