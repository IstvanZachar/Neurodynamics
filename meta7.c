// Deme-besed metapopulational search
// by István Zachar
// patterns are represented as (-1.0, 1.0) floats

// Version history
//   3: 
//   4: added changing environment; removed Nk landscape entirely
//   5: bugfixes at testing stored values (does not affect experiments)
//   6: implemented TYPE 7 simulated AANN-s returning closest stored pattern
//   7: Added ENVPOP == 4: ENVPOPMUT input is only active when after MAXTRAINT


#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <string.h>

#include "randomGenerator.c" // own RNG library (same as András's, with minor bugfix and extra functionality) 
#include "colors.c"          // for monitoring purposes only


#define TOSTRING(x)     (#x)  // macros to work with definitions like `FITNESS (hammingFitness)`
#define FUNCTIONNAME(x) (TOSTRING(x))

#define SEED (393052)

// 0 = unchanged input
// 1 = mutated output (P_mut = 1/N)
// 2 = correlated output 
// 3 = random output
// 4 = autoassociative attractor neural network
// 5 = autoassociative attractor neural network, with no spurious pattern
// 6 = simulated AANN, with 1 memory slot each (last learned)
// 7 = simulated AANN, returning closest stored (with some noise?)
#define TYPE 4 // type of operator to turn input to output

// Fitness landscape & problem size
#define N           100  // number of neurons (200) // if Nk landscape is used, the maximal value of N is 63!
#define P           10   // partition size; must be an integer that divides N (only if blockFitness is used)
#define T           2   // number of optima per block (only if blockFitness is used)
#define B       	  (N/P) // number of blocks per sequence (only if blockFitness is used)
// NOTE: pearsonFitness cannot deal with sequences of uniform bits!!!
#define FITNESS 		hammingFitness // fitness function (1.0 is best, minimum need not be 0.0): pearsonFitness, hammingFitness, blockFitness
#define FITNESSNAME	FUNCTIONNAME(FITNESS)

// Populations
#define SPATIAL  0       // use spatial lattice? (only matters if (DX*DY) > 1)
#define DX       1       // 1st dimension of spatial deme lattice
#define DY       1       // 2nd dimension of spatial deme lattice
#define DN       (DX*DY) // total deme number
#define NN       100      // population size of a deme (number of networks)

// Selection & evolution
#define SHUFFLE       1 // shuffle output population of previous generation for input to the next generation
#define ELITESEL      0    // select best from population OR select one randomly for replication
#define REPLACE       0 // how to replace patterns in the output population: 0=replace one of the worst with the selected; 1=replace all of the output with the selected; 2=replace one randomly
#define RECOMBINATION 0.0 //(0.1)    // probability of recombining two output sequences; if not recombination, migration
#define MIGRATION     0.0 //(0.004)   // probability of migration (only if DN > 1)
#define MUTATION      (1.0/(float)N) // per-bit mutation rate of selected re-entered output of previous generation
#define RECALLNOISE   (0.0) // (0.002) // when simulating AANN (TYPE == 7), mutate stored pattern to yield actual output 

// Environment & simulation time
#define MAXT       12000   // number of generations
#define ENVNUM     2       // number of different environments to cycle through (use ENVNUM = 1 and ENVLEN > MAXT for stable environment; ENVNUM = 0 won't work; Note, that if ENVNUM > 1, global optimum is set accordingly!)
#define ENVLEN     2000    // length of a stable period
#define ENVPOP	   1       // 0 = keep previous population; 1 = random new population; 2 = uniform orthogonal-to-actual-global-optimum population
#define ENVPOPMUT  0.5  // only if ENVPOP = 3
#define ENVWRES    0 // reset weights matrix to random at changes
#define ENDT       0       //(RESOLUTION) // if termination condition is satisfied for this many steps, terminate simulation; use 0 for no termination condition
#define RESOLUTION 1 //(ENVLEN/2)
#define VERBOSE    0

// Autoassociative attractor neural networks
#define UPDATENUM  5    // number of update steps when recall
#define THRESHOLD  0.0  // firing threshold
#define TRAINNUM   40    // number of networks to be trained with the same output (only if TYPE >= 4)
#define TRAINNOISE 0.0  // per-bit mutation rate of re-learned pattern (only if TYPE >= 4)
#define MAXTRAINT  10000 // train up to this time only

// Attractor network testing (only if TYPE >= 5)
#define STORETEST     1    // test storedness of output (only used if TYPE == 4; otherwise sometimes it is automatically on)
#define STOREMAX      N    // maximal number of stored patterns per network; only used for tracking learned patterns; theoretical capacity is: N/sqrt(2*log(N));
#define TESTBEST      0 //(ENVLEN/2) // test recall of `globalOptimumAll` at every `TESTBEST` timestep; TESTBEST == 0 does not test at all 
#define TESTNOISE     0.1  // noise applied for iterative recall of the same input
#define TESTTOLERANCE 0.0  // accept equivalence if output is closer than (or equal to) tolerance (in relative Hamming distance)



// General landscape variables
float globalOptimumS[N];
float globalOptimumW;
float globalOptimumAll[ENVNUM][N]; // store optimum for each environment (for Pearson or Hamming fitnesses)
// Block landscape
float blockOptimaSequence[B][T][P]; // store local optima sequences for each block (for block fitness only)
float blockOptimaFitness [B][T];    // store local optima fitnesses for each block (for block fitness only)

// Auxiliary
int matrixNeighbors[DN][8]; // spatial deme neighbors
int STORE = (TYPE == 5) || (TYPE == 6) || (TYPE == 7) || ((TYPE == 4) && STORETEST); // store learned patterns in `stored` array
int PN; // if `STORE`, store this many learned patterns
int learnCounter[DN][NN];

typedef struct DEME {
  float w[NN], minW, maxW, avgW, input[NN][N], output[NN][N], test[NN][N], maxOutput[N], weight[NN][N][N];
	int   minP, maxP, storedQ[NN], storedP[NN];
	float storedD[NN], stored[NN][STOREMAX][N];
} DEME;
DEME deme[DN];






// mathematical functions

float pearsonCorrelation(float *u, float *v) { // Pearson product-moment correlation coefficient of vectors `u` and `v`
  int i;
  float uAvg = 0.0, vAvg = 0.0, uSqd = 0.0, vSqd = 0.0, cov = 0.0;
  
  for(i = 0; i < N; i++) {
		uAvg += u[i];
		vAvg += v[i];
  }
  uAvg = uAvg/(float)N;
  vAvg = vAvg/(float)N;
  for(i = 0; i < N; i++) {
		cov  += (u[i]-uAvg)*(v[i]-vAvg);
		uSqd += (u[i]-uAvg)*(u[i]-uAvg);
		vSqd += (v[i]-vAvg)*(v[i]-vAvg);
  }
	//printf("... %f %f %f\n", cov, uSqd, vSqd);
	if(cov*uSqd*vSqd == 0) {printf("Uniform sequence encountered in `pearsonCorrelation`. Aborting.\n"); exit(1);};
  return(cov/sqrt(uSqd)/sqrt(vSqd)); 
}

float pearsonFitness(float *v) {
	return(pearsonCorrelation(v, globalOptimumS));
}

int hammingDistanceN(float *v, float *u, int n) { // standard HD up to length `n`
  int i, d = 0;
  for(i = 0; i < n; i++) if(v[i] != u[i]) d++;
  return(d);
}

int hammingDistance(float *v, float *u) { // standard HD of length N
  int i, d = 0;
  for(i = 0; i < N; i++) if(v[i] != u[i]) d++;
  return(d);
}

float relativeHammingDistance(float *v, float *u) { // HD/N
  return((float)hammingDistance(v, u)/(float)N);
}

float hammingFitness(float *v) { // 1 - (HD/N)
  return(1.0 - relativeHammingDistance(v, globalOptimumS));
}

int randomMinPosition(float *v, int n) { // selects the position of the smallest value in vector `v` up to length `n`; if there are multiple instances, selects one randomly.
	float found;
	int i, count = 0, r, pos[n];
	for(i = 0; i < n; i++) pos[i] = -1;
	for(i = 0; i < n; i++) {
		if(i == 0 || v[i] < found) {
			count = 0;
			found = v[i];
		} 
		if(v[i] == found) {
			pos[count] = i;
			count++;
		}
	}
	if(count==1) r = 0;
	else         r = randl(count);
	// for(i = 0; i < n; i++) printf("%f ", v[i]);   printf("\n");
	// for(i = 0; i < n; i++) printf("%d ", pos[i]); printf("\n");
	// printf("MIN:%f   COUNT:%d   REFPOS:%d   POS:%d\n", found, count, r, pos[r]);
	return(pos[r]);
}

int randomMaxPosition(float *v, int n) { // selects the position of the largest value in vector `v` up to length `n`; if there are multiple instances, selects one randomly.
	float found;
	int i, count = 0, r, pos[n];
	for(i = 0; i < n; i++) pos[i] = -1;
	for(i = 0; i < n; i++) {
		if(i == 0 || v[i] > found) {
			count = 0;
			found = v[i];
		} 
		if(v[i] == found) {
			pos[count] = i;
			count++;
		}
	}
	if(count==1) r = 0;
	else         r = randl(count);
	// for(i = 0; i < n; i++) printf("%f ", v[i]);   printf("\n");
	// for(i = 0; i < n; i++) printf("%d ", pos[i]); printf("\n");
	// printf("MAX:%f   COUNT:%d   REFPOS:%d   POS:%d\n", found, count, r, pos[r]);
	return(pos[r]);
}

int firstMaxPosition(float *v, int n) { // selects the position of the largest value in vector `v` up to length `n`; if there are multiple instances, selects last (to the right).
	float max = -999.;
	int i, pos = -1;
	for(i = 0; i < n; i++) if(v[i] > max) {
		max = v[i]; 
		pos = i;
	}
	return(pos);
}


// pattern functions

int samePatternQ(float *u, float *v) { // boolean test of pattern identity
	int i = 0, q = 1;
  while((i < N) && q) {
		if(u[i] != v[i]) q = 0;
		i++;
	}
	return(q);
}

void copyPattern(float *to, float *from) {
	int i;
	for(i = 0; i < N; i++) to[i] = from[i];
}

void invertPattern(float *v) { // flip each bit in vector `v`
	int i;
	for(i = 0; i < N; i++) v[i] = -1.0 * v[i];	
}

void printPattern(float *v) {
	int i;
	for(i = 0; i < N; i++) {
		if      (v[i] == -1.) printf("-");
		else if (v[i] ==  1.) printf("+");
		else                  printf(".");
	}
	printf("\n");
	fflush(stdout);
}

const char *vectorToString(float *v, int n) { // cannot be called multiple times in e.g. printf!
	int i;
	static char s[N+1]; // maximum length of vector to be printed is N
	s[0] = '\0';
	s[n] = '\0';
  for(i = 0; i < n; i++) if(v[i] == -1.0) s[i] = '-'; else {
		if(v[i] == 1.0) s[i] = '+'; else s[i] = '.';
	}
	return(s);
}

void mutatePattern(float *v, float mut) { // mutate pattern per-digit mutation rate `mut`
	int i;
	if(mut > 0.0) for(i = 0; i < N; i++) if(randd() < mut) v[i] = -1. * v[i];
}

void exactlyMutatePattern(float *v, float mut) { // mutates exactly `mut*N` bits (maximally N)
	int i = 0, m = mut*N, r, ref[N] = {0};	
	if(mut > 0) {
		while(i < m && i < N) {
			do r = randl(N); while(ref[r]);
			v[r] = -1. * v[r];
			ref[r] = 1;
			i++;
		}
	}
}

void mutateSinglebitPattern(float *v, float mut) {
	int i;
	i = randl(N/2);
	v[N/2 + i] = -1. * v[N/2 + i];
}

void randomPattern(float *v, float mut) { // generates a random pattern with `mut` probability that a bit is +1
	int i;
	for(i = 0; i < N; i++) v[i] = (randd() < mut) ? 1.0 : -1.0;
}

void alternatingPattern(float *v, int l, int start) { // alternating +1.0 and -1.0 with partition length `l`, starting with `start` (1 -> +1.0, else -> -1.0) 
	int i, b = 1;
	for(i = 0; i < N; i++) {
		v[i] = (b)?(+1.0):(-1.0);
		if(((i+1) % l) == 0) b = 1 - b;
	}
	if(start != 1) invertPattern(v);
}

void continuousPattern(float *v, int from, int to) { // [---..+++..---] where + subpattern runs from `from` to `to-1`
	int i, t;
	// TODO: Is not handled when `from` > `to`
	for(i = 0; (i < from) && (i < N); i++) v[i] = -1.0;
	for(     ; (i < to)   && (i < N); i++) v[i] = +1.0;
	for(     ; (i < N);               i++) v[i] = -1.0;
}

void correlatePattern(float *v, float corr) { // generate a partially uncorrelated version of vector `v` (saving back to `v`)
	int i;
  float x[N];
	if(corr < 1.0) {
		if(corr == -1.0) {
			for(i = 0; i < N; i++) v[i] *= -1.0;
		}	else {
			for(i = 0; i < N; i++) x[i] = v[i];
			while(pearsonCorrelation(x, v) >= corr) {
				i = randl(N);     
				v[i] *= -1.0;
			}
		}
	}
}

void recombinePattern(float * u, float * v) { // two-point recombination
	int i, p1 = randl( N ), p2 = randl( N );
	float temp;
	if(p1 > p2) {
		p1 = p1+p2;
		p2 = p1-p2;
		p1 = p1-p2;
	}
	for(i = p1; i < p2; i++) {
		temp = u[i];
		u[i] = v[i];
		v[i] = temp;
	}
}


// network functions

void setWeights(void) {
  int d, i, j, k;
  for(d = 0; d < DN; d++) for(i = 0; i < NN; i++) for(j = 0; j < N; j++) for(k = 0; k < N; k++) deme[d].weight[i][j][k] = 0.0;
}

void randomWeights(void) {
  int d, i, j, k;
  for(d = 0; d < DN; d++) for(i = 0; i < NN; i++) for(j = 0; j < N; j++) for(k = 0; k < N; k++) deme[d].weight[i][j][k] = randd();
}

void storePattern(float *v, int d, int n) { // tests if `v` is stored in network `n` at deme `d` or not, and pushes/stores it to most recent position (on the right)
  // If there is an identical in `stored`, remove it, shift those right of it one step to the left and store the new one at the last position
	int j = PN-1;
	while(!samePatternQ(v, deme[d].stored[n][j]) && j >= 0) j--;
	if(j < 0) j = 0; // if already present, start shift at position `j`, if not present, start shift at position 0
	for( ; j < PN-1; j++) copyPattern(deme[d].stored[n][j], deme[d].stored[n][j+1]);
	copyPattern(deme[d].stored[n][PN-1], v);
}

void trainNetwork(float *v, int d, int n) { // training network `n` in deme `d` with vector `v`, updating `weight`
  int i, j;
	float h[N] = {0}, f = 1.0/(float)N;
	for(i = 0; i < N; i++) for(j = 0; j < N; j++) h[i] += deme[d].weight[n][i][j] * v[j];
	for(i = 0; i < N; i++) for(j = 0; j < N; j++) {
		if(i == j) deme[d].weight[n][i][j]  = 0.0;
		else       deme[d].weight[n][i][j] += f*v[i]*v[j] - f*v[i]*h[j] - f*v[j]*h[i];
	}
}

void trainNetworksRandom(int n) { // train each network on `n` random patterns (and store these patterns if `STORE == 1`)
	int d, i, j;
	float test[N];
	for(d = 0; d < DN; d++) for(i = 0; i < NN; i++) for(j = 0; j < n; j++) {
		randomPattern(test, 0.5);
		if(STORE) storePattern(test, d, i);
		trainNetwork(test, d, i);
	}
}

void updateOutput(float *v, int d, int n, int i) { // update of neuron `i` in network `n` in deme `d` with threshold neuron model
   int j;
   float h = 0.0;
   for(j = 0; j < N; j++) if(j != i) h += (deme[d].weight[n][i][j]) * v[j];
   if(h >= THRESHOLD) v[i] =  1.0;
   else               v[i] = -1.0;
}

void updateNetwork(float *v, int d, int n) { // updating network 'n' in deme `d` with input `v`; output is written into `v`
   int i, u;
   for(u = 0; u < UPDATENUM; u++) for(i = 0; i < N; i++) updateOutput(v, d, n, randl(N));
}

void shuffleOutputs(int d) { // shuffle `output` population for deme `d`
  // Fisher-Yates shuffle from: http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
	int i, j;
	float temp[N];
	for (i = NN - 1; i > 0; i--) {
			j = randl(i+1);
			copyPattern(temp,              deme[d].output[j]);
			copyPattern(deme[d].output[j], deme[d].output[i]);
			copyPattern(deme[d].output[i], temp);
	}
}

float testAgainst(int e, int iter) { // tests all networks of all demes `iter` times against the global optimum of environment `e`, returning average per-bit relative distance
	int d, n, i, sum = 0;
	float v[N];
	for(d = 0; d < DN; d++) {
		for(n = 0; n < NN; n++) {
			for(i = 0; i < iter; i++) {
				copyPattern(v, globalOptimumAll[e]);
				mutatePattern(v, TESTNOISE);
				updateNetwork(v, d, n);
				sum += hammingDistance(v, globalOptimumAll[e]);
			}
		}
	}
	return((float)sum/(float)(N*NN*DN*iter));
}

int closestStored(float *v, int d, int n) { // returns the index of the stored pattern that is closest to `v` in deme `d`, network `n`
	int j;
	float hd[PN];
	for(j = 0; j < PN; j++) hd[j] = relativeHammingDistance(v, deme[d].stored[n][j]);
	return(randomMinPosition(hd, PN));
}


// Landscaping

void setBlockOptimaDefault(int env) { // sets the same T optima (sequence and fitness) for the B blocks in environment `env` {1, 2, ...}
  // E1: first [++++...], rest [-+-+...]
	// E2: first [----...], rest [+-+-...]
	int i, j, k;
	for(i = 0; i < B; i++) {
		for(j = 0; j < T; j++) {
			if(j == 0) {
				for(k = 0; k < P; k++) blockOptimaSequence[i][j][k] = (env == 1)?(+1.0):(-1.0);
				blockOptimaFitness[i][j] = 10.0;
			} else {
				for(k = 0; k < P; k++) blockOptimaSequence[i][j][k] = (((env == 1) && (k % 2)) || ((env == 2) && ((k+1) % 2)))?(+1.0):(-1.0);
				blockOptimaFitness[i][j] = 1.0;
			}
		}
	}
}

float blockGlobalOptimumFitness(void) { // calculates global optimum using `blockOptimaSequence` and `blockOptimaFitness`
	int i, j, bestP;
	float bestW, sum = 0.0;
	for(i = 0; i < B; i++) {
		bestW = 0.0; // best fitness
		bestP = -99; // position of best fitness
		for(j = 0; j < T; j++) { // find best fitnessed target for given block; assuming there are no to identical best fitnesses
			if(blockOptimaFitness[i][j] > bestW) {
				bestP = j;
				bestW = blockOptimaFitness[i][j];
			}
		}
		for(j = 0; j < T; j++) if(j == bestP) sum += bestW; else sum += 1.0/(1.0 + hammingDistanceN(blockOptimaSequence[i][bestP], blockOptimaSequence[i][j], P));
	}
	return(sum);
}

void blockGlobalOptimumSequence(float *v) { // calculates global optimum using `blockOptimaSequence` and `blockOptimaFitness`
	int i, j, bestP;
	float bestW;
	for(i = 0; i < B; i++) {
		bestW = 0.0; // best fitness
		bestP = -99; // position of best fitness
		for(j = 0; j < T; j++) { // find best fitnessed target for given block; assuming there are no to identical best fitnesses
			if(blockOptimaFitness[i][j] > bestW) {
				bestP = j;
				bestW = blockOptimaFitness[i][j];
			}
		}
		for(j = 0; j < P; j++) v[i*P + j] = blockOptimaSequence[i][bestP][j];
	}
}

float blockFitness(float *v) { // building block fitness; rescales range (min, max) to (min/max, 1.0)
	int i, j, d;
	float c, block[P], sum = 0.0;
	for(i = 0; i < B; i++) {
		for(j = 0; j < P; j++) block[j] = v[i*P + j];
		for(j = 0; j < T; j++) {
			d = hammingDistanceN(block, blockOptimaSequence[i][j], P);
			if(d == 0) c = blockOptimaFitness[i][j]; else c = 1.0/(1.0+(float)d);
			sum += c;
			//printf("\t%d %d %d %f %f\n", i, j, d, c, sum);
		}
	}
	return(sum/globalOptimumW);
}

void setLandscape(int ls, int env) { // sets `globalOptimumS` and `globalOptimumW`
	int i, j;
	
	if((ls == 1) || (ls == 2)) {	// Pearson & Hamming fitnesses
		copyPattern(globalOptimumS, globalOptimumAll[env]);
		globalOptimumW = 1.0;
	} else if(ls == 3) { // Block fitness		
		setBlockOptimaDefault(env); // generate optima for environment
		blockGlobalOptimumSequence(globalOptimumS);
		globalOptimumW = blockGlobalOptimumFitness();
		if(VERBOSE) {
			printf("\tBlock target optima ({N=%d P=%d B=%d T=%d):\n", N, P, B, T);
			for(i = 0; i < B; i++) {
				printf("\tBLOCK #%d\n", i);
				for(j = 0; j < T; j++) printf("\t\tT%d %s %f\n", j, vectorToString(blockOptimaSequence[i][j], P), blockOptimaFitness[i][j]);
			}
		}
	}
	
	if(VERBOSE) {
		printf("\tnew environment: E%d\n", env);
		printf("\tglobal optimum:  %s %f\n\n", vectorToString(globalOptimumS, N), globalOptimumW);
	}
}


// Deme functions

void setNeighborDemes(int n, int m) { // sets up a lookup table for neughboring positions for each matrix element; 8 neighbours, no periodic boundary
	int i, j, k, pos[n][m], c = 0, p;
  for(i = 0; i < n; i++) for(j = 0; j < m; j++) {
		pos[i][j] = c++;
		for(k = 0; k < 8; k++) matrixNeighbors[i*m + j][k] = -99;
	}
	// for(i = 0; i < n; i++) for(j = 0; j < m; j++) printf("%2d%s", pos[i][j], (j == m-1)?"\n":" ");
	// for(i = 0; i < n*m; i++) for(j = 0; j < 8; j++) printf("%d%s", matrixNeighbors[p][j], (j == 7)?"\n":" ");
	for(i = 0; i < n; i++) for(j = 0; j < m; j++) {
		c = 0;
		p = i*m + j;
		if(i > 0) {
			if(j > 0)   matrixNeighbors[p][c++] = pos[i-1][j-1];
			            matrixNeighbors[p][c++] = pos[i-1][j];
			if(j < m-1) matrixNeighbors[p][c++] = pos[i-1][j+1];
		}
	  if(j > 0)   matrixNeighbors[p][c++] = pos[i][j-1];
	  if(j < m-1) matrixNeighbors[p][c++] = pos[i][j+1];
		if(i < n-1) {
			if(j > 0)   matrixNeighbors[p][c++] = pos[i+1][j-1];
			            matrixNeighbors[p][c++] = pos[i+1][j];
			if(j < m-1) matrixNeighbors[p][c++] = pos[i+1][j+1];
		}
	}
}






int main(int argc, char** argv) {
	int i, j, k, d, e, t = 0, terminateSum = 0, environmentC = 0, LANDSCAPE = 0, learnC = 0;
	float bestTest[ENVNUM] = {0};

	
	seed(SEED); // seed RNG of randomGenerator.c with `unsigned long`
	
	
	
	// Landscape setup
	if(strcmp("pearsonFitness", FITNESSNAME) == 0) LANDSCAPE = 1; else
	if(strcmp("hammingFitness", FITNESSNAME) == 0) LANDSCAPE = 2; else
	if(strcmp("blockFitness",   FITNESSNAME) == 0) LANDSCAPE = 3; else {
		printf("Invalid fitness function %s. Aborting.\n", FITNESSNAME);
		exit(1);
	}
	if(LANDSCAPE == 1) for(e = 0; e < ENVNUM; e++) alternatingPattern(globalOptimumAll[e], N/2, (e==0));   else // pearsonFitness; [++..--..] for E0, [--..++..] for any other Ei
	if(LANDSCAPE == 2) for(e = 0; e < ENVNUM; e++) randomPattern(globalOptimumAll[e], (e==0)?(1.0):(0.0));      // hammingFitness; [+++++..] for E0, [-----..] for any other Ei
	//if(LANDSCAPE == 2) for(e = 0; e < ENVNUM; e++) continuousPattern(globalOptimumAll[e], e*N/ENVNUM, (e+1)*N/ENVNUM);      // hammingFitness; [+++..---] for E(first), [--..++..--] for intermediate E, [---..+++] for E(last)
	if(VERBOSE) {
		printf("LANDSCAPE: %s %d\n", FITNESSNAME, LANDSCAPE);
		for(i = 0; i < ENVNUM; i++) printPattern(globalOptimumAll[i]);
	}

		
	// Spatial setup
	if(SPATIAL) setNeighborDemes(DX, DY);
	
	
	// Network setup
	if(((TYPE == 4) && (STORETEST)) || (TYPE == 5) || (TYPE == 7)) PN = STOREMAX; else PN = 1;
	for(d = 0; d < DN; d++) for(i = 0; i < NN; i++) randomPattern(deme[d].output[i], 0.5);
	if((TYPE == 4) || (TYPE == 5) || (TYPE == 7)) {
		setWeights();
		trainNetworksRandom(STOREMAX); // NOTE: MUST be set to STOREMAX so that `stored` is correctly initialized and contains patterns at all positions. DO NOT set it to PN, as it can be zero which would result in faulty AANN-s
	} else if(TYPE == 6) {
		for(d = 0; d < DN; d++) for(i = 0; i < NN; i++) for(j = 0; j < PN; j++) randomPattern(deme[d].stored[i][j], 0.5);
	}
	
	
	if(VERBOSE) {
		printf("Method: %d\n", TYPE);
		printf("\tDN=%d (%d x %d)   spatial=%d   NN=%d \n", DN, DX, DY, SPATIAL, NN);
		printf("\tstore=%d   N=%d   C(exp)=%f   C(act)=%d\n", STORE, N, N/sqrt(2.0*log((float)N)), PN);
		printf("\tN(env)=%d   L(env)=%d    pop.reset? %d\n", ENVNUM, ENVLEN, ENVPOP);
		printf("\n");
	}
	
	
	
	// Generations
	while((t <= MAXT) && (((ENDT > 0) && (terminateSum < ENDT)) || (ENDT == 0))) {
		int bestD = -99, bestP = -99, recombineC = 0, migrateC = 0;
		float avgW = 0, bestW = -99., storedC = 0.0;
		
		
		// test recall properties at middle and end of environment period
		if((TESTBEST) > 0 && ((t % (TESTBEST)) == 0)) for(i = 0; i < ENVNUM; i++) bestTest[i] = testAgainst(i, 100);
		
		
		// Environment & landscaping
		if(((ENVNUM > 1) && ((t % ENVLEN) == 0)) || ((ENVNUM == 1) && (t == 0))) {
					
			setLandscape(LANDSCAPE, environmentC);
			if((ENVWRES) && (t < MAXTRAINT)) {
				setWeights();
				trainNetworksRandom(STOREMAX);
			}
			// printf("CHANGED ENVIRONMENT (T=%d)\n", t);
			// printf("\tnew optimum: E%d %s %f\n", environmentC, vectorToString(globalOptimumS, N), globalOptimumW);
			// printf("\tmean learnC: %f (/network/deme/%d steps)\n", (float)(learnC-learnCOld)/(float)(DN*NN), ENVLEN);
			if(environmentC == (ENVNUM-1)) environmentC = 0; else environmentC++;
			// new output population
			if(ENVPOP == 0) {} else 
			if(ENVPOP == 1) for(d = 0; d < DN; d++) for(i = 0; i < NN; i++) randomPattern(deme[d].output[i], 0.5); else
		  if(ENVPOP == 2) for(d = 0; d < DN; d++) for(i = 0; i < NN; i++) {
				copyPattern(deme[d].output[i], globalOptimumS);
				invertPattern(deme[d].output[i]);
				//printPattern(deme[d].output[i]);
			}
			if(ENVPOP == 3) for(d = 0; d < DN; d++) for(i = 0; i < NN; i++) { // use mutated optimum cue
				copyPattern(deme[d].output[i], globalOptimumS);
				mutatePattern(deme[d].output[i], ENVPOPMUT);
			}
			if(ENVPOP == 4) { // use random input when learning AND use noisy optimum when not learning anymore
				if(t < MAXTRAINT) {
					for(d = 0; d < DN; d++) for(i = 0; i < NN; i++) randomPattern(deme[d].output[i], 0.5);
				} else {
					for(d = 0; d < DN; d++) for(i = 0; i < NN; i++) { // when still learning, use previous population; after learning, use mutated optimum cue
						copyPattern(deme[d].output[i], globalOptimumS);
						mutatePattern(deme[d].output[i], ENVPOPMUT);
					}
				}
			}
			learnC = 0;
		}

		
		for(d = 0; d < DN; d++) {
			
			/* Update output based on new input */
			if(t > 0) {
				if(SHUFFLE) shuffleOutputs(d);
				for(i = 0; i < NN; i++) copyPattern(deme[d].input[i], deme[d].output[i]); // save input for comparison (for TYPE == 5 or for STORETEST == 1)
				
				if(TYPE == 0) {                                                                       } else  // 0 = unchanged input
				if(TYPE == 1) {for(i = 0; i < NN; i++) mutatePattern(deme[d].output[i], 1.0/(float)N);} else  // 1 = mutated output
				if(TYPE == 2) {for(i = 0; i < NN; i++) correlatePattern(deme[d].output[i], 1.0 );     } else  // 2 = correlated output
				if(TYPE == 3) {for(i = 0; i < NN; i++) randomPattern(deme[d].output[i], 0.5 );        } else  // 3 = random output
				if(TYPE == 4) {for(i = 0; i < NN; i++) updateNetwork(deme[d].output[i], d, i);        } else  // 4 = AANN
				if(TYPE == 5) {for(i = 0; i < NN; i++) updateNetwork(deme[d].output[i], d, i);        } else  // 5 = AANN with no spurious output (iterate recall until stored is returned)
				if(TYPE == 6) {for(i = 0; i < NN; i++) copyPattern(deme[d].output[i], deme[d].stored[i][0]);} else // 6 = simulated AANN with 1-slot memory, always returning last learned
			  if(TYPE == 7) {for(i = 0; i < NN; i++) {
					int j;
					j = closestStored(deme[d].output[i], d, i);
					copyPattern(deme[d].output[i], deme[d].stored[i][j]);
					mutatePattern(deme[d].output[i], RECALLNOISE);
					} // 7 = simulated AANN returning closest learned pattern
				}
			}
				
			for(i = 0; i < NN; i++) {
				
				/* Test storedness */
				//if(((TYPE == 4) && (STORETEST)) || (TYPE == 5)) { 
				if((STORE && STORETEST) || (TYPE == 5)) { 
					float hd[PN], minhd = 1.0;
					
					for(j = 0; j < PN; j++) {
						hd[j] = relativeHammingDistance(deme[d].output[i], deme[d].stored[i][j]);
						if(hd[j] < minhd) minhd = hd[j];
					}
					// printf("D%d N%d: ", d, i);
					// for(j = 0; j < PN; j++) printf("%1.2f ", hd[j]);
					// printf("\n");
					
					/* Only for AANN without spurious patterns */
					/* If output `test` is not close to any stored pattern, repeat input->output */
					/* This practically eliminates the possibility to encounter non-stored "spurious" (or random) output patterns */
					while ((TYPE == 5) && (minhd > TESTTOLERANCE)) {
						copyPattern(deme[d].test[i], deme[d].input[i]);
						mutatePattern(deme[d].test[i], TESTNOISE);
						updateNetwork(deme[d].test[i], d, i);
						for(j = 0; j < PN; j++) {
							hd[j] = relativeHammingDistance(deme[d].test[i], deme[d].stored[i][j]);
							if(hd[j] < minhd) minhd = hd[j];
						}
					}
					
					if(minhd <= TESTTOLERANCE) {
						deme[d].storedP[i] = randomMinPosition(hd, PN);
						deme[d].storedD[i] = hd[deme[d].storedP[i]];
						deme[d].storedQ[i] = 1;
						// if(deme[d].storedD[i] != minhd) {
							// printf("ERROR: randomMinPosition could not find the appropriate minimum! Aborting.\n");
							// exit(1);
						// };
					} else {
						deme[d].storedP[i] = -9;
						deme[d].storedD[i] = minhd;
						deme[d].storedQ[i] = 0;
					}
					storedC += deme[d].storedQ[i];						
				}
				
				
				deme[d].w[i] = FITNESS(deme[d].output[i]);
			}
			
			/* Find best solution (chose one randomly, if there are multiple with identical HD-s) */
			deme[d].avgW = 0.0;
			for(i = 0; i < NN; i ++) deme[d].avgW += deme[d].w[i]; 
			deme[d].maxP = randomMaxPosition(deme[d].w, NN); // NOTE: this might chose a position that is not an attractor, however there might be another that is.
			deme[d].maxW = deme[d].w[deme[d].maxP];
			deme[d].minP = randomMinPosition(deme[d].w, NN); // NOTE: this might chose a position that is not an attractor, however there might be another that is.
			deme[d].minW = deme[d].w[deme[d].minP];
			deme[d].avgW = deme[d].avgW/(float)NN;
			copyPattern(deme[d].maxOutput, deme[d].output[deme[d].maxP]);
			if(deme[d].maxW > bestW) {
				bestD = d;
				bestW = deme[d].maxW;
				bestP = deme[d].maxP;
			}
			avgW += deme[d].avgW;
		}
		
		
		
		if(VERBOSE) {
		//if((t > 10980) && (t < 10983)) {
			printf("T%d\n", t);
			for(d = 0; d < DN; d++) {
				printf("  D%d\n", d);
				for(i = 0; i < NN; i++) {
					if((d == bestD) && (i == bestP)) {
						printCyan ("\tN%d %s %f", i, vectorToString(deme[d].output[i], N), deme[d].w[i]);
					} else if(i == deme[d].maxP) {
						printGreen("\tN%d %s %f", i, vectorToString(deme[d].output[i], N), deme[d].w[i]);
					} else if(i == deme[d].minP) {
						printRed  ("\tN%d %s %f", i, vectorToString(deme[d].output[i], N), deme[d].w[i]);
					} else {
						printf    ("\tN%d %s %f", i, vectorToString(deme[d].output[i], N), deme[d].w[i]);
					}
					if(TYPE >= 4) printf("\t%d #%d: %f", deme[d].storedQ[i], deme[d].storedP[i], deme[d].storedD[i]);
					printf("\n");
					
					// for(j = 0; j < NN; j++) {
						// printf("\tS%d ", j);
						// printPattern(deme[d].stored[i][j]);
						// printf("\n");
					// }
					
					fflush(stdout);
				}
			}
		}
		

		if(t > 0) { // skip this in the first turn
		
			// Selection (recombination or mutation)
			for(d = 0; d < DN; d++) {
				int p1;
				float w1, v1[N];
				
				if(ELITESEL) p1 = deme[d].maxP; else p1 = randl(NN);
				
				if(randd() < RECOMBINATION) {
					int p2;
					float w2, v2[N];
					
					/* Recombine */
					recombineC++;
					copyPattern(v1, deme[d].output[p1]);
					do p2 = randl(NN); while (p2 == p1);
					if((DN > 1) && (randd() < MIGRATION)) {
						/* Migrate */
						int e;
						
						if(SPATIAL) {
							do e = matrixNeighbors[d][randl(8)]; while (e < 0); // select neighboring deme on lattice
						} else {
							do e = randl(DN); while (e == d); // select random other deme
						}
						
						copyPattern(v2, deme[e].output[p2]);
						migrateC++;
						//printf("\t\t%d -> %d\n", e, d);
					} else {
						copyPattern(v2, deme[d].output[p2]);
					}
					recombinePattern(v1, v2);
					w1 = FITNESS(v1);
					w2 = FITNESS(v2);
					if(w1 < w2) {
						copyPattern(v1, v2);
						w1 = w2;
					}
				} else {
					/* Mutate */
					copyPattern(v1, deme[d].output[p1]);
					mutatePattern(v1, (MUTATION));
					w1 = FITNESS(v1);
				}
				
				/* Replication & learn */
				if(w1 > deme[d].minW) {

					if(REPLACE == 0) {copyPattern(deme[d].output[deme[d].minP], v1);} else // replace the worst patterns in the output (choosing one randomly, if there are multiple)
					if(REPLACE == 1) {for(i = 0; i < NN; i++); copyPattern(deme[d].output[i], v1);} else // replace ALL of the output population with the new output
					if(REPLACE == 2) {copyPattern(deme[d].output[randl(NN)], v1);} // replace ONE pattern randomly
					
					if((TYPE >= 4) && (t < MAXTRAINT)) { // Learning
						int pos[NN] = {0};
						float toLearn[N];
						for(k = 0; k < ((TRAINNUM>NN)?NN:TRAINNUM); k++) { // the best pattern is assured to be trained to TRAINNUM *different* networks (faster convergence)
							do i = randl(NN); while(pos[i]);
							pos[i] = 1;
							copyPattern(toLearn, v1);
							mutatePattern(toLearn, TRAINNOISE);
							trainNetwork(toLearn, d, i);

							/* Store if learned */
							learnC++;
							if(STORE) storePattern(toLearn, d, i); // Store learnt patterns in `stored` only if there is no identical already
						}
					}
				}			
			} // deme `d`
		}
	
	
		//if(VERBOSE && (bestW == 1.0)) terminateSum++; else terminateSum = 0;
		if((bestW == 1.0)) terminateSum++; else terminateSum = 0;

		
		
		if(!(t % (RESOLUTION))) {

			printf("%d\t%f\t%f\t%d\t%d\t%f\t%d\t%f\t", t, bestW, avgW/(float)DN, deme[bestD].storedQ[bestP], deme[bestD].storedP[bestP], deme[bestD].storedD[bestP], learnC, (float)storedC/(float)(DN*NN));
			if(STORE) for(d = 0; d < DN; d++) {
				for(i = 0; i < NN; i++) printf("%d\t", deme[d].storedP[i]);
				//for(i = 0; i < NN; i++) printf("%f\t", deme[d].storedD[i]);
				if(DN > 1) printf("\n");
			}
			if((TESTBEST) > 0) for(e = 0; e < ENVNUM; e++) printf("%f\t", bestTest[e]);
			printf("\n");
			fflush(stdout);
		}
		
		
		t++;
	}
	
	
	
	
  return(0);
}
