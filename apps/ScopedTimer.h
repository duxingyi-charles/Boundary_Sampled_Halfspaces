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
        std::cout << "[" << m_name << "]: " << toc() << "s"
                  << std::endl;
    }

    double toc() const {
        const auto t = Timer::now();
        std::chrono::duration<double> d = t - m_tic;
        return d.count();
    }

   private:
    std::string m_name;
    std::chrono::time_point<Timer> m_tic;
};
