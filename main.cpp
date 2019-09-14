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
  double nA,nT,nC,nG = 0; //how many instances of this particular nucleotide
  //bigram instances
  double aa,ac,at,ag,ca,cc,ct,cg,ta,tc,tt,tg,ga,gc,gt,gg = 0; //instances of bigrams
  double total, bigramTotal = 0; //total amount of nucleotides in the file, total amount of bigrams
  ifstream readFile2(dnaFile);
  while(getline(readFile2, dnaString)){
    for(int i = 0; i < dnaString.length();++i){ //iterate through dnaString
      ++total; //adds to total number of nucleotides
      if(tolower(dnaString[i]) == 'a'){//checks if the character matches the nucleotide
        ++nA; //adds to number of given nucleotide
      }
      else if(tolower(dnaString[i]) == 'c'){//checks if the character matches the nucleotide
        ++nC; //adds to number of given nucleotide
      }
      else if(tolower(dnaString[i]) == 't'){//checks if the character matches the nucleotide
        ++nT; //adds to number of given nucleotide
      }
      else if(tolower(dnaString[i]) == 'g'){//checks if the character matches the nucleotide
        ++nG; //adds to number of given nucleotide
      }
    }

  }

  ifstream readFile3(dnaFile);

  while(getline(readFile3, dnaString)){
    for(int j = 0; j<dnaString.length()-1;++j){
      ++bigramTotal; //tracks how many bigrams there are in total throughout the file
      if(tolower(dnaString[j]) == 'a'){
        if(tolower(dnaString[j+1]) == 'a'){
          ++aa;
        }
        else if(tolower(dnaString[j+1]) == 'c'){
          ++ac;
        }
        else if(tolower(dnaString[j+1]) == 't'){
          ++at;
        }
        else if(tolower(dnaString[j+1]) == 'g'){
          ++ag;
        }
      }
      else if(tolower(dnaString[j]) == 'c'){
        if(tolower(dnaString[j+1]) == 'a'){
          ++ca;
        }
        else if(tolower(dnaString[j+1]) == 'c'){
          ++cc;
        }
        else if(tolower(dnaString[j+1]) == 't'){
          ++ct;
        }
        else if(tolower(dnaString[j+1]) == 'g'){
          ++cg;
        }
      }
      else if(tolower(dnaString[j]) == 't'){
        if(tolower(dnaString[j+1]) == 'a'){
          ++ta;
        }
        else if(tolower(dnaString[j+1]) == 'c'){
          ++tc;
        }
        else if(tolower(dnaString[j+1]) == 't'){
          ++tt;
        }
        else if(tolower(dnaString[j+1]) == 'g'){
          ++tg;
        }
      }
      else if(tolower(dnaString[j]) == 'g'){
        if(tolower(dnaString[j+1]) == 'a'){
          ++ga;
        }
        else if(tolower(dnaString[j+1]) == 'c'){
          ++gc;
        }
        else if(tolower(dnaString[j+1]) == 't'){
          ++gt;
        }
        else if(tolower(dnaString[j+1]) == 'g'){
          ++gg;
        }
      }
    }
  }

  double probA = a.findRelativeProb(nA, total);
  double probT = a.findRelativeProb(nT, total);
  double probC = a.findRelativeProb(nC, total);
  double probG = a.findRelativeProb(nG, total);
  cout << "a: " << nA << endl;
  cout << "t: " << nT << endl;
  cout << "c: " << nC << endl;
  cout << "g: " << nG << endl;
  cout << "Total: " << total << endl;
  cout << "Probability of A: " << probA <<endl;
  cout << "Probability of T: " << probT <<endl;
  cout << "Probability of C: " << probC <<endl;
  cout << "Probability of G: " << probG <<endl;
  cout << "\n" << endl;

  cout << "Bigram total: " << bigramTotal <<endl;
  cout << "AA:" << aa << endl;
  double probAA = a.findRelativeProb(aa, bigramTotal);
  double probAC = a.findRelativeProb(ac, bigramTotal);
  double probAT = a.findRelativeProb(at, bigramTotal);
  double probAG = a.findRelativeProb(ag, bigramTotal);
  double probCA = a.findRelativeProb(ca, bigramTotal);
  double probCC = a.findRelativeProb(cc, bigramTotal);
  double probCT = a.findRelativeProb(ct, bigramTotal);
  double probCG = a.findRelativeProb(cg, bigramTotal);
  double probTA = a.findRelativeProb(ta, bigramTotal);
  double probTC = a.findRelativeProb(tc, bigramTotal);
  double probTT = a.findRelativeProb(tt, bigramTotal);
  double probTG = a.findRelativeProb(tg, bigramTotal);
  double probGA = a.findRelativeProb(ga, bigramTotal);
  double probGC = a.findRelativeProb(gc, bigramTotal);
  double probGT = a.findRelativeProb(gt, bigramTotal);
  double probGG = a.findRelativeProb(gg, bigramTotal);
  cout << "Probability of AA: " << probAA << endl;
  cout << "Probability of AC: " << probAC << endl;
  cout << "Probability of AT: " << probAT << endl;
  cout << "Probability of AG: " << probAG << endl;
  cout << "Probability of CA: " << probCA << endl;
  cout << "Probability of CC: " << probCC << endl;
  cout << "Probability of CT: " << probCT << endl;
  cout << "Probability of CG: " << probCG << endl;
  cout << "Probability of TA: " << probTA << endl;
  cout << "Probability of TC: " << probTC << endl;
  cout << "Probability of TT: " << probTT << endl;
  cout << "Probability of TG: " << probTG << endl;
  cout << "Probability of GA: " << probGA << endl;
  cout << "Probability of GC: " << probGC << endl;
  cout << "Probability of GT: " << probGT << endl;
  cout << "Probability of GG: " << probGG << endl;
  //cout << "Mean = " << mean << endl;

  ofstream writeFile;
  writeFile.open("jamestran.out");
  writeFile << "mean = " << mean << endl;
  writeFile<< "variance = " << variance << endl;
  writeFile<< "sd = " << standardDeviation << endl;
  writeFile << "a: " << nA << endl;
  writeFile << "t: " << nT << endl;
  writeFile << "c: " << nC << endl;
  writeFile << "g: " << nG << endl;
  writeFile << "Total: " << total << endl;
  writeFile << "Probability of A: " << probA <<endl;
  writeFile << "Probability of T: " << probT <<endl;
  writeFile << "Probability of C: " << probC <<endl;
  writeFile << "Probability of G: " << probG <<endl;
  writeFile << "\n" << endl;
  writeFile << "Bigram total: " << bigramTotal <<endl;
  writeFile << "Probability of AA: " << probAA << endl;
  writeFile << "Probability of AC: " << probAC << endl;
  writeFile << "Probability of AT: " << probAT << endl;
  writeFile << "Probability of AG: " << probAG << endl;
  writeFile << "Probability of CA: " << probCA << endl;
  writeFile << "Probability of CC: " << probCC << endl;
  writeFile << "Probability of CT: " << probCT << endl;
  writeFile << "Probability of CG: " << probCG << endl;
  writeFile << "Probability of TA: " << probTA << endl;
  writeFile << "Probability of TC: " << probTC << endl;
  writeFile << "Probability of TT: " << probTT << endl;
  writeFile << "Probability of TG: " << probTG << endl;
  writeFile << "Probability of GA: " << probGA << endl;
  writeFile << "Probability of GC: " << probGC << endl;
  writeFile << "Probability of GT: " << probGT << endl;
  writeFile << "Probability of GG: " << probGG << endl;

  writeFile.close();
  return 0;


}
