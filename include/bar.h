#ifndef PROGRESS_BAR_H
#define PROGRESS_BAR_H

#include <iostream>

class ProgressBar {
public:
  ProgressBar() {
  }
  virtual ~ProgressBar() {
    std::cout << std::endl;
  }
  void print(float progress) const {
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
  }
private:
  int barWidth = 70;
};

#endif
