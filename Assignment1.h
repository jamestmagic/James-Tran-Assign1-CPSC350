#include <iostream> //preprocessor directive

using namespace std;

class Assignment1{
  public:
    Assignment1(); //constructor
    ~Assignment1(); //destructor
    void readStrings(string t);
    double findMean(double total, double lines);
    double findVariance(double meanDiff, double lines);
    double findSD(double variance);
    double findRelativeProb(double numInstances, double total);
};
