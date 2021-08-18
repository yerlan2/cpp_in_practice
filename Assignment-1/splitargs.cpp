// Copyright year Negmetulla_Yerlan 180107219@stu.sdu.edu.kz

#include <iostream>

int main(int argc, char **argv) {
  argv++;
  int c = 0;
  while (*argv) {
    if (c < 4)
      std::cout << *argv << std::endl;
    else 
      std::cerr << *argv << std::endl;
    argv++;
    c++;
  }
  return 0;
}
