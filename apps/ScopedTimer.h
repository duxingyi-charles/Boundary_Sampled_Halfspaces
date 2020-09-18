#pragma once

#include <chrono>
#include <iostream>
#include <string>

template <typename Timer = std::chrono::high_resolution_clock>
class ScopedTimer {
   public:
    ScopedTimer(const std::string& name = "default") : m_name(name) {
        m_tic = Timer::now();
    }

    ~ScopedTimer() {
        auto toc = Timer::now();
        std::chrono::duration<double> t = toc - m_tic;
        std::cout << "[" << m_name << "]: " << t.count() << "s"
                  << std::endl;
    }

   private:
    std::string m_name;
    std::chrono::time_point<Timer> m_tic;
};
