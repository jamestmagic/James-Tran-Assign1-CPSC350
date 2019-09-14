#include "Assignment1.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

Assignment1::Assignment1()
{

}
Assignment1::~Assignment1()
{
  cout << "object deleted" << endl;
}

void readStrings(string t){

}

double Assignment1::findMean(double sum, double lines){ //lines is the amount of lines of DNA strings in the file
  return sum/lines;
}

double Assignment1::findVariance(double meanDiff, double lines){ //meanDiff is the sum of all the mean differences squared
  return meanDiff/lines;
}

double Assignment1::findSD(double variance){
  return sqrt(variance);
}

double Assignment1::findRelativeProb(double numInstances, double total){
  return numInstances/total; //divides number of instances of a given nucleotide by the total amount of nucleotides
}
