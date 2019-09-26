   
#ifndef HANDY_TIMER
#define HANDY_TIMER

#include <chrono>

// A simple timer class to measure time spent
// Lekhraj 

namespace measure {
typedef std::chrono::nanoseconds TimeT;
struct timer
{
    std::chrono::steady_clock::time_point start_t;
    timer() { start(); }
    void start() {
         start_t = std::chrono::steady_clock::now();
       
    }
    double elapsed() {
      TimeT duration = std::chrono::duration_cast< TimeT> 
                            (std::chrono::steady_clock::now() - start_t);
        return duration.count();
    }
};
}

#endif
