#ifndef XCSMACROS_H
#define XCSMACROS_H


/**
 * Most parameter-names are chosen similar to the 'An Algorithmic Description of XCS' ( Butz&Wilson, IlliGAL report 2000017).
 */

const double alpha=0.1; //The fall of rate in the fitness evaluation.
const double beta=0.2; //The learning rate for updating fitness, prediction, prediction error, and action set size estimate in XCS's classifiers.
const double gama=0.95; //The discount rate in multi-step problems.
const double delta=0.1; //The fraction of the mean fitness of the population below which the fitness of a classifier may be considered in its vote for deletion.
const double nu=5.0; //Specifies the exponent in the power function for the fitness evaluation.
const double theta_GA=25.0; //The threshold for the GA application in an action set.
const double epsilon_0= 10.0; //The error threshold under which the accuracy of a classifier is set to one.
const int theta_del = 20; //Specified the threshold over which the fitness of a classifier may be considered in its deletion probability.
const int crossoverType = 2; // 0 uniform, 1 onePoint, and 2 twoPoint crossover.
const double pX = 0.5; // 0.8; // 0.04; //0.04; //The probability of applying crossover in an offspring classifier.
const double pM = 0.5; //0.04; //0.8; //The probability of mutating one allele and the action in an offspring classifier.
const float p_ol = 0.5;  // probability of using a filter from observed list
const int mutationType = 0; // 0 niche, and 1 general mutation.
const double P_dontcare = 0; // 0.33; //The probability of using a don't care symbol in an allele when covering.
const double pDontcare_filter = .33; // this is for don't care within filter
const double predictionErrorReduction = 1.0; //0.25; //The reduction of the prediction error when generating an offspring classifier.
const double fitnessReduction=0.1; //The reduction of the fitness when generating an offspring classifier.
const int theta_sub = 20; //The experience of a classifier required to be a subsumer.
const int theta_promising = 10; // A classifier is "promising" if its experience > theta_promising and error < epsilon_0
const double predictionIni=10.0; //The initial prediction value when generating a new classifier (e.g in covering).
const double predictionErrorIni=0.0; //The initial prediction error value when generating a new classifier (e.g in covering).
const double fitnessIni=0.01; //The initial fitness value when generating a new classifier (e.g in covering).
const bool doGASubsumption = true; //Specifies if GA subsumption should be executed.
const bool doActSetSubsumption = false;//true; //Specifies if action set subsumption should be executed.
const double tournamentSize = 0.4; //The fraction of classifiers participating in a tournament from an action set.
const double forceDifferentInTournament = 0.0;
const bool doGAErrorBasedSelect = false;
const double selectTolerance = 0.0;
const char dontcare='#'; //The don't care symbol (normally '#')
const float m_0 = 0.5; //to be used in mutation
const float m = 0.1; // how much to mutuate one single ellel
const float epsilon = 0.5; // the probability of exploration for epsilon greedy strategy
const int N_filter_ol = 2500;  // maximum limit of filters in the managed observed list


const long _M = 2147483647; //Constant for the random number generator (modulus of PMMLCG = 2^31 -1).
const long _A = 16807; //Constant for the random number generator (default = 16807).


void setSeed(long newSeed);
long getSeed();
double drand();

int irand(int n); // returns a number from 0 to n-1 inclusive


#endif // XCSMACROS_H
