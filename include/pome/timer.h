#ifndef INCLUDE_TIMER_H
#define INCLUDE_TIMER_H

#include <chrono>
#include <iostream>

namespace pome {

using AwesomeClock = std::chrono::high_resolution_clock;

class Timer {
 public:
  Timer();
  ~Timer();

 protected:
  std::chrono::time_point<AwesomeClock> m_start;
  std::chrono::time_point<AwesomeClock> m_end;
  std::chrono::duration<double> m_duration;
};

}  // namespace pome

#endif