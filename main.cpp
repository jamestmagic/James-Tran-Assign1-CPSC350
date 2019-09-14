#define _USE_MATH_DEFINES
#include "Assignment1.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <cmath>
#include <time.h>

using namespace std;

int main(int argc, char** argv){
  string dnaFile = argv[1];
  Assignment1 a;
  a.allOperationsAndGenerations(dnaFile);
  char answer = 'y';
  while(answer=='y'){
    cout << "Would you like to proccess another list?(y/n)" << endl;
    cin >> answer;
    if(towlower(answer) == 'y'){
      cout << "Please enter the name of your new file." << endl;
      cin >> dnaFile;
      a.allOperationsAndGenerations(dnaFile);
    }
  }
  return 0;


}
