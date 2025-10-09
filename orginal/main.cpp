#include "../inc/common.h"

int main(int argc, const char *argv[]) {
  if (argc != 2) {
    cout << "Server wrong input parameters!" << endl;
    exit(1);
  }

  string filepath = argv[1];

  return 0;
}