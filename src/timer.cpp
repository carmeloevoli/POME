#include "pome/timer.h"

namespace pome {

Timer::Timer() { m_start = AwesomeClock::now(); }

Timer::~Timer() {
  m_end = std::chrono::high_resolution_clock::now();
  m_duration = m_end - m_start;
  std::cerr << "Timer took : " << m_duration.count() << " s." << std::endl;
}

}  // namespace pome
