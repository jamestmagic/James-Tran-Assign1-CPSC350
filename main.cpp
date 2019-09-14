#include "Assignment1.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

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
    ++numLines; //counts number of lines
    numChars = numChars + dnaString.length();
  }

  double mean = a.findMean(numChars, numLines);
  cout << "mean = " << mean << endl;

  ifstream readFile1(dnaFile);
  while(getline(readFile1, dnaString)){
    meanDiff += pow((dnaString.length() - mean),2); //takes the difference of the mean and squares it
  }

  double variance = a.findVariance(meanDiff, numLines);
  double standardDeviation = a.findSD(variance);
  cout<< "variance = " << variance << endl;
  cout<< "sd = " << standardDeviation << endl;
  //nucleotides for relative probability
  char nucleotide;
  int nA = 0; //how many instances of this particular nucleotide
  int nT = 0;
  int nC = 0;
  int nG = 0;
  int total = 0; //total amount of nucleotides in the file
  ifstream readFile2(dnaFile);
  while(getline(readFile2, dnaString)){
    for(int i = 0; i < dnaString.length();++i){ //iterate through dnaString
      if(tolower(dnaString[i]) == 'a'){//checks if the character matches the nucleotide
        ++nA; //adds to number of given nucleotide
        ++total; //adds to total amount of nucleotides in file
      }
      else if(tolower(dnaString[i]) == 'c'){//checks if the character matches the nucleotide
        ++nC; //adds to number of given nucleotide
        ++total; //adds to total amount of nucleotides in file
      }
      else if(tolower(dnaString[i]) == 't'){//checks if the character matches the nucleotide
        ++nT; //adds to number of given nucleotide
        ++total; //adds to total amount of nucleotides in file
      }
      else if(tolower(dnaString[i]) == 'g'){//checks if the character matches the nucleotide
        ++nG; //adds to number of given nucleotide
        ++total; //adds to total amount of nucleotides in file
      }
    }
  }
  
  cout << "a: " << nA << endl;
  //cout << "Mean = " << mean << endl;
  return 0;


}
