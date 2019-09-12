#include "Assignment1.h"
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

int main(int argc, char** argv){
  string dnaFile = argv[1];
  Assignment1 a;

  ifstream readFile(dnaFile);

  double numLines = 0; //number of lines in the file or number of dna strings
  double numChars = 0; //total number of nucleotides in all lines
  double meanDiff = 0; //end sum of the differences of the means and number of characters per lines
                      //used to calculate the variance
  string dnaString; //


  while(getline(readFile, dnaString)){
    ++numLines;
    numChars = numChars + dnaString.length();
  }

  double mean = a.findMean(numChars, numLines);

  while(getline(readFile, dnaString)){
    meanDiff += pow((dnaString.length() - mean),2); //takes the difference of the mean and squares it
  }

  double variance = a.findVariance(meanDiff, numLines);
  //cout << "Mean = " << mean << endl;
  return 0;


}
