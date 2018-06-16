/*
 * CumulativeTimeMeasurement.cpp
 *
 *  Created on: 24 Mar 2017
 *      Author: mikolajr
 */

#include "CumulativeTimeMeasurement.h"

#include <sys/time.h> // timeval
#include <sys/resource.h> // rusage

namespace SignificantPattern {

double CumulativeTimeMeasurement::measureTime()
{
  struct rusage t;
  struct timeval tv;
  getrusage(RUSAGE_SELF, &t);
  tv = t.ru_utime;
  return tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

} /* namespace SignificantPattern */
