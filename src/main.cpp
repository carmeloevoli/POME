#include "git_revision.h"
#include "logging.h"
#include "pome.h"
#include "timer.h"

int main(int argc, char *argv[]) {
  log_startup_information();
  try {
    POME::Timer timer;

    POME::ModelState state;

    POME::Pome pome(state);
    pome.buildSource();
    pome.dumpSource();

    pome.evolve();

  } catch (const std::exception &e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return 0;
}