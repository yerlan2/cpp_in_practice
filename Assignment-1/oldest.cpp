// Copyright year Negmetulla_Yerlan 180107219@stu.sdu.edu.kz

#include <iostream>
#include <string>
#include <stdlib.h>

using namespace std;
int main(int argc, char **argv) {
  string str = "ls -tr | head -";
  if (*(++argv))
    str += *argv;
  const char *cmd = str.c_str();
  system(cmd);
  return 0;
}
