#define _USE_MATH_DEFINES
#include "Assignment1.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <cmath>
#include <time.h>

using namespace std;

Assignment1::Assignment1()//constructor
{

}
Assignment1::~Assignment1()//destructor
{
  cout << "object deleted" << endl;
}

double Assignment1::findMean(double sum, double lines){ //lines is the amount of lines of DNA strings in the file
  return sum/lines;
}

double Assignment1::findVariance(double meanDiff, double lines){ //meanDiff is the sum of all the mean differences squared
  return meanDiff/lines;
}

double Assignment1::findSD(double variance){//calculates standardDeviation from variance
  return sqrt(variance);
}

double Assignment1::findRelativeProb(double numInstances, double total){
  return numInstances/total; //divides number of instances of a given nucleotide by the total amount of nucleotides
}

void Assignment1::allOperationsAndGenerations(string fileName){
  srand(time(NULL)); //time seed for rand()

  ifstream readFile(fileName); //starts stream

  double numLines = 0; //number of lines in the file or number of dna strings
  double numChars = 0; //total number of nucleotides in all lines
  double meanDiff = 0; //end sum of the differences of the means and number of characters per lines
                      //used to calculate the variance
  string dnaString; //placeholder string for iterating through file streams


  while(getline(readFile, dnaString)){
    ++numLines; //counts number of lines
    numChars = numChars + dnaString.length(); //keeps track of total number of characters in file
  }

  double mean = findMean(numChars, numLines);//calculates mean based on number of characters and number of lines in file

  ifstream readFile1(fileName); //starts stream
  while(getline(readFile1, dnaString)){
    meanDiff += pow((dnaString.length() - mean),2); //takes the difference of the mean and squares it
  }

  //calculates variance and standardDeviation based on data collected from file
  double variance = findVariance(meanDiff, numLines);
  double standardDeviation = findSD(variance);
  //nucleotides for relative probability
  double nA,nT,nC,nG = 0; //how many instances of this particular nucleotide
  //bigram instances
  double aa,ac,at,ag,ca,cc,ct,cg,ta,tc,tt,tg,ga,gc,gt,gg = 0; //instances of bigrams
  double total, bigramTotal = 0; //total amount of nucleotides in the file, total amount of bigrams
  ifstream readFile2(fileName); //starts stream again for new iteration
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

  ifstream readFile3(fileName); //starts stream

  while(getline(readFile3, dnaString)){
    for(int j = 0; j<dnaString.length()-1;++j){
      ++bigramTotal; //tracks how many bigrams there are in total throughout the file
      //if statements check for specific bigram permutations
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
  //relative probability of individual nucleotides
  double probA = findRelativeProb(nA, total);
  double probT = findRelativeProb(nT, total);
  double probC = findRelativeProb(nC, total);
  double probG = findRelativeProb(nG, total);

  //relative probability of nucleotide bigrams
  double probAA = findRelativeProb(aa, bigramTotal);
  double probAC = findRelativeProb(ac, bigramTotal);
  double probAT = findRelativeProb(at, bigramTotal);
  double probAG = findRelativeProb(ag, bigramTotal);
  double probCA = findRelativeProb(ca, bigramTotal);
  double probCC = findRelativeProb(cc, bigramTotal);
  double probCT = findRelativeProb(ct, bigramTotal);
  double probCG = findRelativeProb(cg, bigramTotal);
  double probTA = findRelativeProb(ta, bigramTotal);
  double probTC = findRelativeProb(tc, bigramTotal);
  double probTT = findRelativeProb(tt, bigramTotal);
  double probTG = findRelativeProb(tg, bigramTotal);
  double probGA = findRelativeProb(ga, bigramTotal);
  double probGC = findRelativeProb(gc, bigramTotal);
  double probGT = findRelativeProb(gt, bigramTotal);
  double probGG = findRelativeProb(gg, bigramTotal);

  //Ranges for the probability of each nucleotide bigram
  double aaRange = probAA;
  double acRange = probAA+probAC;
  double atRange = probAA+probAC+probAT;
  double agRange = probAA+probAC+probAT+probAG;
  double caRange = probAA+probAC+probAT+probAG+probCA;
  double ccRange = probAA+probAC+probAT+probAG+probCA+probCC;
  double ctRange = probAA+probAC+probAT+probAG+probCA+probCC+probCT;
  double cgRange = probAA+probAC+probAT+probAG+probCA+probCC+probCT+probCG;
  double taRange = probAA+probAC+probAT+probAG+probCA+probCC+probCT+probCG+probTA;
  double tcRange = probAA+probAC+probAT+probAG+probCA+probCC+probCT+probCG+probTA+probTC;
  double ttRange = probAA+probAC+probAT+probAG+probCA+probCC+probCT+probCG+probTA+probTC+probTT;
  double tgRange = probAA+probAC+probAT+probAG+probCA+probCC+probCT+probCG+probTA+probTC+probTT+probTG;
  double gaRange = probAA+probAC+probAT+probAG+probCA+probCC+probCT+probCG+probTA+probTC+probTT+probTG+probGA;
  double gcRange = probAA+probAC+probAT+probAG+probCA+probCC+probCT+probCG+probTA+probTC+probTT+probTG+probGA+probGC;
  double gtRange = probAA+probAC+probAT+probAG+probCA+probCC+probCT+probCG+probTA+probTC+probTT+probTG+probGA+probGC+probGT;
  double ggRange = probAA+probAC+probAT+probAG+probCA+probCC+probCT+probCG+probTA+probTC+probTT+probTG+probGA+probGC+probGT+probGG;


  ofstream writeFile; //starts an out stream to write to a file
  writeFile.open("jamestran.out"); //opens a new .out file

  //printing to the file
  writeFile << "James Tran" << endl;
  writeFile << "2318908" << endl;
  writeFile << "CPSC 350" << endl;
  writeFile << "Assignment 1 \n \n" << endl;
  writeFile << "Mean: " << mean << endl;
  writeFile<< "Variance: " << variance << endl;
  writeFile<< "SD: " << standardDeviation << endl;
  writeFile << "A: " << nA << endl;
  writeFile << "T: " << nT << endl;
  writeFile << "C: " << nC << endl;
  writeFile << "G: " << nG << endl;
  writeFile << "Probability of A: " << probA <<endl;
  writeFile << "Probability of T: " << probT <<endl;
  writeFile << "Probability of C: " << probC <<endl;
  writeFile << "Probability of G: " << probG <<endl;
  writeFile << "Total: " << total << endl;
  writeFile << "\n" << endl;
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
  writeFile << "Bigram total: " << bigramTotal << "\n\n" <<endl;

  //Probability ranges for the individual nucleotides
  double aRange = probA;
  double cRange = probA + probC;
  double tRange = probA + probC + probT;
  double gRange = probA + probC + probT + probG;

  //start of dna string generation from statistics above
  for(int i = 0;i<1000;++i){ //creates 1000 strings
    double a1 = ((double)rand()/RAND_MAX);//generates random number between 0 and 1
    double b1 = ((double)rand()/RAND_MAX);//generates random number between 0 and 1
    //cout << "a1 = " << a1 <<endl;
    //cout << "b1 = " << b1 <<endl;
    double c1 = sqrt(-2*log(a1)) * cos(2*M_PI*b1);
    //cout << "c1 = " << c1 <<endl;
    int d1 = (standardDeviation*c1) + mean; //how many nucleotides in this string
    //cout << "length of a string: " << d1 << endl;
    for(int j = 0;j<d1/2;++j){//creates d1 pairs of nucleotides
      double a2 = ((double)rand()/RAND_MAX);//generate random number from 0 to 1
      //if statements check if the random number falls within the range of the nucleotides' probability
      //generates first nucleotide and then completes bigram with same random instance based on bigram probability
      if(a2 >= 0 && a2 < aRange){
        writeFile << "A";
        if(a2 >= 0 && a2 < aaRange){
          writeFile << "A";
        }else if(a2 >= aaRange && a2 < acRange){
          writeFile << "C";
        }
        else if(a2 >= acRange && a2 < atRange){
          writeFile << "T";
        }
        else if(a2 >= atRange && a2 < agRange){
          writeFile << "G";
        }
      }
      else if(a2 >= aRange && a2 < cRange){
        writeFile << "C";
        if(a2 >= 0 && a2 < caRange){
          writeFile << "A";
        }else if(a2 >= caRange && a2 < ccRange){
          writeFile << "C";
        }
        else if(a2 >= ccRange && a2 < ctRange){
          writeFile << "T";
        }
        else if(a2 >= ctRange && a2 < cgRange){
          writeFile << "G";
        }
      }
      else if(a2 >= cRange && a2 < tRange){
        writeFile << "T";
        if(a2 >= 0 && a2 < taRange){
          writeFile << "A";
        }else if(a2 >= taRange && a2 < tcRange){
          writeFile << "C";
        }
        else if(a2 >= tcRange && a2 < ttRange){
          writeFile << "T";
        }
        else if(a2 >= ttRange && a2 < tgRange){
          writeFile << "G";
        }
      }
      else if(a2 >= tRange && a2 < gRange){
        writeFile << "G";
        if(a2 >= 0 && a2 < gaRange){
          writeFile << "A";
        }else if(a2 >= gaRange && a2 < gcRange){
          writeFile << "C";
        }
        else if(a2 >= gcRange && a2 < gtRange){
          writeFile << "T";
        }
        else if(a2 >= gtRange && a2 < ggRange){
          writeFile << "G";
        }
      }
    }
    writeFile << "\n";//starts new string
  }


  writeFile.close();//close out stream


}
