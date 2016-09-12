// Storkey's attractor network
// by András Szilágyi, István Zachar
// patterns are represented as (-1, 1) floats
// version 13


#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <inttypes.h>
#include <time.h>
#include <string.h>
#include <windows.h>

#include "randomGenerator.c" // own RNG library (same as András's, with minor bugfix and extra functionality) 
#include "ran4.c" // for Nk's pseudorandom generator; ran4 should only be used for Nk, not for any other random number, as it would mess up the seed's state
#include "colors.c" // for monitoring purposes only



#define SEED (223)

// Autoassociative attractor neural networks
#define AANN      (0)   // use autoassociative attractor neural networks or simply pass input to output without change
#define UPDATENUM (100) // number of update steps (200)
#define LEARNNUM  (1)   // number of training steps (1)
#define THRESHOLD (0.0) // firing threshold

// Fitness landscape & problem size
#define N       (40)  // number of neurons (200) // if Nk landscape is used, the maximal value of N is 63!
#define P       (20)   // partition size; must be an integer that divides N
#define T       (2)   // number of optima per block
#define B       (N/P) // number of blocks per sequence
#define FITNESS (blockFitness) // fitness function (1.0 is best, minimum need not be 0.0): pearsonFitness, hammingFitness, NkFitness, blockFitness
#define NKP     (N)   // Nk landscape phenotype length (generally safe to have N = P
#define K       (6)   // Nk landscape number of epistatic neighbours (uniform for each gene; subsequence length is K+1)

// Multinetwork population
#define DN          (5)  // deme number
#define NN          (10) // number of networks (20)
#define INPUTNOISE  (1.0/(float)N) // (0.01)  // per-bit mutation rate of re-entered output of previous generation
#define LEARNNOISE  (0.005) // per-bit mutation rate of re-learned output of previous generation (only if AANN > 0)
#define TRAINNUM    (1)    // number of networks to be trained with the same output (though differently noised) (only if AANN > 0)
#define REC         (0.5)   // probability of recombining two output sequences
#define MIGRATION   (0.0004) // migration rate
#define MAXT        (5000000)  // number of generations
#define ENDFRAME    (10)    // if the global distance is zero for the last ENDFRAME steps, terminate run
#define GIVEUP      (5)     // when best fitness does not change for GIVEUP time, introduce new variation (recombination, increased mutation, etc.)

// Monitor
#define MONITOR    (0)
#define RESOLUTION (10000)
#define VERBOSE    (1)

// Testing
#define PN               (30)    // maximal number of stored patterns (N/sqrt(2*ln(N)) * 1.1); only used for tracking learned patterns in some tests
#define TESTTOLERANCE    (0.02)  // register attractorness only if output is closer than threshold (in relative Hamming distance)
#define TESTNOISE        (0.1)   //(0.1)  // test attractorness by feeding in noisy variants with this noise-per-bit rate
#define TESTITERATION    (10)    //(10)   // test attractorness by feeding in this many samples iteratively 
#define TESTGOOD         (5)     //(5)  // if this many iterations are reconstructed correctly from TESTITERATION, attractorness is attributed
#define STOREDCORRELATE  (0.0)   // correlation of stored patterns
#define RECALLINPUTNOISE (0.1)   // per-bit mutation over input when recall

#define TESTITER    (1000) // for testing purposes




float storedChanged   [NN][PN][N];
float storedOriginal  [NN][PN][N];
float storedRandom50  [NN][PN][N];
float storedRandom75  [NN][PN][N];
float storedRandom90  [NN][PN][N];
float storedRectangle [NN][PN][N];
float storedMutated   [NN][PN][N];
float storedCorrelated[NN][PN][N];
float storedSingleBit [NN][PN][N];

float output[DN*NN][N];
float weight[NN][N][N];
float globalBest[N];

// Nk landscape
const uint64_t MAX  = (1ULL << (N  )); // size of sequence space (2^N), max decimal value of a sequence is MAX-1
const uint64_t MAXK = (1ULL << (K+1)); // size of subsequence space (2^(K+1)), max decimal value of a subsequence of length K+1 is MAXK-1
const uint64_t MAXB = (1ULL << (N-1)); // largest sequence with a single 1 bit (2^(N-1), represents the sequence [10000..0]; for printing purposes)
uint64_t MAP[NKP] = {0};    // epistatic map, uniform for each sequence
unsigned long RAN4SEED;

// block landscape
float blockOptimaSequence[B][T][P];
float blockOptimaFitness [B][T];
float blockGlobalMax;

// temporary
int counter = 0;
int correct = 0;





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
  return(cov/sqrt(uSqd)/sqrt(vSqd)); 
}

float pearsonFitness(float *v) {
	return(pearsonCorrelation(v, globalBest));
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
  return(1.0 - relativeHammingDistance(v, globalBest));
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

void randomPattern(float *v, float mut) {
	int i;
	for(i = 0; i < N; i++) v[i] = (randd() < mut) ? 1.0 : -1.0;
}

void randomCorrelatePattern(float *v, float corr) { // generate a partially uncorrelated version of vector `v` (saving back to `v`)
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

void recombinePattern(float * u, float * v) { // two point recombination
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

void resetWeights(void) {
  int i, j, k;
  for(k = 0; k < NN; k++) for(i = 0; i < N; i++) for(j = 0; j < N; j++) weight[k][i][j] = 0.0;
}

void resetStoredPatterns(void) { // set various stored patterns
  int i, j, k;
	float master[N];
	for(i = 0; i < NN; i++) {
		randomPattern(master, 0.5);
		for(j = 0; j < PN; j++) {
			randomPattern(storedChanged[i][j], 0.5);
			copyPattern(storedOriginal[i][j], storedChanged[i][j]);
			copyPattern(storedMutated[i][j], master);
			copyPattern(storedCorrelated[i][j], master);
			mutatePattern(storedMutated[i][j], 0.02);
			randomCorrelatePattern(storedCorrelated[i][j], STOREDCORRELATE);
			randomPattern(storedRandom50[i][j], 0.5);
			randomPattern(storedRandom75[i][j], 0.75);
			randomPattern(storedRandom90[i][j], 0.9);
			
			for(k = 0; k < N; k++) {
				storedRectangle[i][j][k] = (k >= j*(N/PN) && k < (j+1)*(N/PN)) ? 1.0 : -1.0;
				storedSingleBit[i][j][k] = (j == k) ? 1.0 : -1.0;
			}
			//printPattern(storedRectangle[i][j]);
			//printPattern(storedOriginal[i][j]);
			//printPattern(storedMutated[i][j]);
			//printPattern(storedCorrelated[i][j]);
		}
	}	
}

void train_network(float *v, int m) { // training network 'm' with vector `v` (palimpsest updating: "replacing" older patterns in weight matrix)
  int i, j;
  float h[N], f = (1.0)/(float)N;
	
  for(i = 0; i < N; i++) {
    h[i] = 0.0;
    for(j = 0; j < N; j++) h[i] += weight[m][i][j] * v[j];
  }
  for(i = 0; i < N; i++) for(j = 0; j < N; j++) {
      if(i == j) weight[m][i][j]  = 0.0;
      else       weight[m][i][j] += f*v[i]*v[j] - f*v[i]*h[j] - f*v[j]*h[i];
    }
}

void update_output(int m, int n) { // update of neuron 'n' in network 'm' with threshold neuron model
   int t;
   float h = 0.0;
   for(t = 0; t < N; t++) if(t != n) h += weight[m][n][t] * output[m][t];
   if(h >= THRESHOLD) output[m][n] =  1.0; // linear threshold
   else               output[m][n] = -1.0;
}

void update_network(int m) { // updating network 'm'
   int t;
   for(t = 0; t < N; t++) update_output(m, randl(N));
}


void recallPatternAgainst(float *v, float vv[NN][PN][N], float mut) { // recall pattern `v` against patterns of `vv`
	int i, j, u, closest = -1, expectedQ = -1, attractorQ = -1;
	float d[PN], dOC, dIE, dOE, input[N];
	
	for(i = 0; i < NN; i++) {
		copyPattern(input, v);
		exactlyMutatePattern(input, mut);
		copyPattern(output[i], input);
		//printPattern(output[i]);
	  for(u = 0; u < UPDATENUM; u++) update_network(i);
		//printPattern(output[i]);
		for(j = 0; j < PN; j++) d[j] = relativeHammingDistance(output[i], vv[i][j]);
	  closest = randomMinPosition(d, PN);
		dOC = d[closest];                                  // Output against Closest pattern
		dOE = relativeHammingDistance(output[i], v);        // Output against Expected pattern
		dIE = relativeHammingDistance(input, v);     // Input  against Expected pattern
		attractorQ = samePatternQ(output[i], vv[i][closest]);
		expectedQ  = samePatternQ(output[i], v);
		if(VERBOSE) {
			printf("\tN#%2d   C=%2d   d(IE)=%.3f   d(OC)=%.3f   d(OE)=%.3f \n", i, closest, dIE, dOC, dOE);
			fflush(stdout);
		} else {
			printf("\tN:%d\tC:%d\td(OC):%f\tOK? %d\n", counter, closest, dOC, dOC < (3.0/(float)N));
			fflush(stdout);
			counter++;
			correct += dOC < (3.0/(float)N);
		}
	}
}

void recallAgainst(float vv[NN][PN][N], float mut) { // recall patterns of `vv` against patterns of `vv`
	int i, j, k, u, closest = -1, expectedQ = -1, attractorQ = -1, expectedSum = 0, attractorSum = 0, count = 0;
	float d[PN], dOC, dIE, dOE, dOT, input[N];
	
	for(i = 0; i < NN; i++) for(j = 0; j < PN; j++) {
		count++;
		copyPattern(input, vv[i][j]);
		mutatePattern(input, mut);
		//if(mut > 0.0) mutateSinglebitPattern(input, 99.0); // TODO REMOVE
		// printPattern(vv[i][j]);
		// printPattern(input);
		copyPattern(output[i], input);
	  for(u = 0; u < UPDATENUM; u++) update_network(i);
		for(k = 0; k < PN; k++) d[k] = relativeHammingDistance(output[i], vv[i][k]);
	  closest = randomMinPosition(d, PN);
		dOC = d[closest];                                  // Output against Closest attractor
		dOE = relativeHammingDistance(output[i], vv[i][j]); // Output against Expected attractor
		dIE = relativeHammingDistance(input, vv[i][j]);     // Input  against Expected attractor
		dOT = relativeHammingDistance(output[i], globalBest);     // Output against globally optimal Target
		attractorQ = samePatternQ(output[i], vv[i][closest]);
		expectedQ  = samePatternQ(output[i], vv[i][j]);
		attractorSum += attractorQ;
		expectedSum  += expectedQ;
		if(VERBOSE) {
			printf("N/E:%2d/%2d   C=%2d   d(IE)=%.3f   d(OC)=%.3f   d(OE)=%.3f   d(OT)=%.3f   O?C:%d    O?E:%d\n", i, j, closest, dIE, dOC, dOE, dOT, attractorQ, expectedQ);
			fflush(stdout);
		}
	}
	if(VERBOSE) {
		printf("RECALL expectedQ  (IE): failed=%d\tok=%d \n",   count-expectedSum, expectedSum);
		printf("RECALL attractorQ (IC): failed=%d\tok=%d \n\n", count-attractorSum, attractorSum);
		fflush(stdout);
	} else {
		printf("%d\t%d\n", count - expectedQ, expectedQ);
		fflush(stdout);
	}
}

int storedPatternQTest(float *v, float vv[NN][PN][N], int nn, float mut) {
	int i, j, u, s = -99, stored = -99, stableQ = 1, storedQ = 1, iter = 10;
	float out[N];
	
	for(i = 0; i < iter; i++) {
		s = -99;
		copyPattern(output[nn], v);
		mutatePattern(output[nn], mut);
		for(u = 0; u < UPDATENUM; u++) update_network(nn);
		for(j = 0; j < PN; j++) if(samePatternQ(output[nn], vv[nn][j])) s = j;
		if(i == 0) {
			copyPattern(out, output[nn]);
			stored = s;
		}
		stableQ = stableQ && samePatternQ(output[nn], out);
		storedQ = storedQ && s >= 0 && s == stored;
		//printf("\t\tC:%d\n", s);
	}
	printf("\tstableQ: %d   storedQ: %d   stored#: %d\t (i:%d, m:%f)\n", stableQ, storedQ, stored, iter, mut);
	fflush(stdout);
	return(stableQ && storedQ);
}

int storedAttractorQ(float *v, float vv[NN][PN][N], int nn) { // tests if the recalled output `v` can be reproduced AND it is a stored pattern in `vv[nn]`, of network `nn`
	int i = 0, j, pos = -99, sum = 0;
	float d[PN], minD = 1.0;
	
	for(j = 0; j < PN; j++) {
		d[j] = relativeHammingDistance(v, vv[nn][j]);
		if(d[j] < minD) minD = d[j];
	}
	pos = randomMinPosition(d, PN); // NOTE: this returns one randomly, but if there are multiple identical distances, those are the same sequences - which is impossible inside a single network)
	
	do {
		copyPattern(output[nn], v);
		exactlyMutatePattern(output[nn], TESTNOISE);
		for(j = 0; j < UPDATENUM; j++) update_network(nn);
		sum += (relativeHammingDistance(output[nn], v) <= TESTTOLERANCE);
		i++;
	} while (i < TESTITERATION && sum < TESTGOOD);
	
	printf("\tN#%d   C:%2d   d:%.3f   s:%2d   A?:%d   %s\n", nn, pos, minD, sum, (sum >= TESTGOOD) && (minD <= TESTTOLERANCE), (sum >= TESTGOOD && minD > TESTTOLERANCE)?"!":"");
	fflush(stdout);
	
	return(sum >= TESTGOOD) && (minD <= TESTTOLERANCE);
}

int testAttractorQ(float *v, int nn) { // tests if the recalled outputs of `v` can be reproduced (for at max TESTTOLERANCE relative distance)
	int i = 0, j, c = 0;
	float d;
	
	do {
		copyPattern(output[nn], v);
		exactlyMutatePattern(output[nn], TESTNOISE);
		//mutatePattern(output[nn], TESTNOISE);
		for(j = 0; j < UPDATENUM; j++) {
			update_network(nn);
			//if(t == 41 && nn == 2 && j >= UPDATENUM-2) printPattern(output[nn]);
		}
		d = relativeHammingDistance(output[nn], v);
		c += (d <= TESTTOLERANCE);
		// if(t == 2 || t == 41) {
			// if(i == 0)  printf("\tTEST #%d\n", nn);
			// printf("\t\t%d %f %d\n", i, d, c);
		// }
		i++;
	} while (i < TESTITERATION && c < TESTGOOD);
	
	//if(t == 2 || t == 41) printf("\t\tEND: %d\n", c >= TESTGOOD);
	return(c >= TESTGOOD);
}


float averageDistanceStored(int *ref, float vv[NN][PN][N], int nn) { // calculates average distance of each possible pair of patterns in `vv`, according to the reference vector `ref` indicating whether `vv[nn][i]` is filled or not
	int i, j, c = 0;
	float sum = 0.0;
	
	for(i = 0; i < PN; i++) {
		if(ref[i] > 0) {
			for(j = i+1; j < PN; j++) {
				if(ref[j] > 0) {
					c++;
					sum += relativeHammingDistance(vv[nn][i], vv[nn][j]);
				}
			}
		}
	}
	return(sum/(float)c);	
}

float averageDistance(int *ref, float vv[TESTITER][N]) { // calculates average distance of each possible pair of patterns in `vv`, according to the reference vector `ref` indicating whether `vv[i]` is filled or not
	int i, j, c = 0;
	float sum = 0.0;
	
	for(i = 0; i < TESTITER; i++) {
		if(ref[i] > 0) {
			for(j = 0; j < i; j++) {
				if(ref[j] > 0) {
					c++;
					sum += relativeHammingDistance(vv[i], vv[j]);
				}
			}
		}
	}
	return(sum/(float)c);	
}

void attractorStatistics(int nn) { // generates and inputs TESTITER random patterns to the network `nn` to find its attractors
	int nSum = 0, i, j, k, attQ, closest, c = 0;
	int sSum = 0, sDiffSum = 0, sRef[PN] = {0};
	int eSum = 0, eDiffSum = 0, eRef[TESTITER] = {0};
	float v[TESTITER][N], d[PN];

	for(i = 0; i < TESTITER; i++) {
		k = -1;
		randomPattern(output[nn], 0.5);
		for(j = 0; j < UPDATENUM; j++) update_network(nn);
		copyPattern(v[i], output[nn]);
		attQ = testAttractorQ(v[i], nn);	
		for(j = 0; j < PN; j++) d[j] = relativeHammingDistance(storedChanged[nn][j], v[i]);
		closest = randomMinPosition(d, PN);
		if(attQ) {
			if(d[closest] == 0) {
				sSum++;
				sRef[closest]++;
			} else {
				k = 0;
				eSum++;
				while(!samePatternQ(v[k], v[i]) && (k < i)) k++;
				eRef[k]++;
			}
		} else nSum++;
		// printf("\t#%3d A:%d P:%2d d:%.3f    %s   [%d]\n", i, attQ, closest, d[closest], (attQ ? ((d[closest] == 0) ? "S" : "X") : "-"), k); // "S" = stored attractor; "E" = not stored (extra) attractor
		// fflush(stdout);
	}
	
	for(i = 0; i < PN;       i++) sDiffSum += (sRef[i] ? 1 : 0);
	for(i = 0; i < TESTITER; i++) eDiffSum += (eRef[i] ? 1 : 0);
	if(VERBOSE) {
		printf("ATTRACTOR STATISTICS of network #%d: \n", nn);
		printf("\tstored:       %d, %d different (of %d), avg.dist.int:%f \n", sSum, sDiffSum, PN, averageDistanceStored(sRef, storedChanged, nn));
		printf("\tdistribution: ["); for(i = 0; i < PN; i++)       printf("%d ", sRef[i]); printf("\b]\n");
		printf("\tnot stored:   %d, %d different, avg.dist.int:%f\n", eSum, eDiffSum, averageDistance(eRef, v));
		//printf("\tdistribution: ["); for(i = 0; i < TESTITER; i++) printf("%d ", eRef[i]); printf("\b]\n");
		printf("\tno attractor: %d\n", nSum);
		printf("\ttotal:        %d\n", sSum+eSum+nSum);
		fflush(stdout);
	}
}



// Nk landscape

uint64_t rotateRight(uint64_t num, int step) { // bitwise rotation to the right (circular)
	int s = step % N;
	return((num >> s) | ((num << (N-s)) % MAX));
}

uint64_t rotateLeft(uint64_t num, int step) { // bitwise rotation to the left (circular)
	int s = step % N;
	return(((num << s) % MAX) | (num >> (N-s)));
}

uint64_t fromBinaryVector(float *v) { // converts binary vector `v` of {-1, +1} to a binary integer
	int i;
	uint64_t num = 0;
	for(i = 0; i < N; i++) num += rotateRight(MAXB, i)*(v[i] > 0); // NOTE: Consider how `v` is represented! {0, 1} or {-1, +1}?
	return(num);
}

double NkFitness(float *v, float *u) { // Calculates the Nk fitness value for `v`, with [111..1] having w=1.0 and [000..0] having w=0.0 fitness; expecting a binary float vector `v`; ignores `u`
	int i;
	uint64_t seq, nth, num;
	double sum = 0.;
	
	num = fromBinaryVector(v);
	for(i = 0; i < NKP; i++) {
		seq = rotateLeft(num & MAP[i], K+1+i); // extract epistatic neighborhood of size `K+1` at position `i` from sequence `num`, according to the epistatic map `MAP`
		if((0 < seq) && (seq < MAXK-1)) {
			nth = seq + i*MAXK;
			sum += seededNthRan4(RAN4SEED, (unsigned long)nth);
		} else if(seq == MAXK-1) {
			sum += 1.0; // subsequence is of best sequence [111..1], has fitness contribution 1.0;
		} // else do not increment `sum`
		//printf("\t#%d\t%"PRIu64"\t%"PRIu64"\t%f\n", i, seq, i*MAXK, sum);
	}
	//return(sum/(double)P);
	return(0.1 + 0.9*(sum/(double)NKP)); // NOTE
}


// Block landscape

void setBlockOptimaDefault( ) { // sets the same T optima (sequence and fitness) for the B blocks: first is [1111...], rest is [0101...]
	int i, j, k;
	for(i = 0; i < B; i++) {
		for(j = 0; j < T; j++) {
			if(j == 0) {
				for(k = 0; k < P; k++) blockOptimaSequence[i][j][k] = 1.0; // [1111...]
				blockOptimaFitness[i][j] = 10.0;
			} else {
				for(k = 0; k < P; k++) blockOptimaSequence[i][j][k] = (k % 2)?(1.0):(-1.0); // [0101...]
				blockOptimaFitness[i][j] = 1.0;
			}
		}
	}
}

float blockGlobalOptimumFitness( ) { // calculates global optimum using `blockOptimaSequence` and `blockOptimaFitness`
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
	return(sum/blockGlobalMax);
}








int main(int argc, char** argv) {
	int i;


	// if(TRAINNUM > NN) { // error checking
		// printf("TRAINNUM > NN. Aborting.\n");
		// exit(1);
	// }


	seed(SEED); // seed RNG of randomGenerator.c with `unsigned long`
	RAN4SEED = randl(LONG_MAX); // seed for Nk algorithm (`ran4` seed, required so that `ran4` depends on the session, i.e. on the overall randomizer state set by `SEED`)

	
	// reset wheights and training patterns
	resetWeights();
	resetStoredPatterns();

	
	if(0) { // Simple recall
		int u, i, j;
		
		for(u = 0; u < LEARNNUM; u++) for(i = 0; i < NN; i++) for(j = 0; j < PN; j++) train_network(storedOriginal[i][j], i);
		recallAgainst(storedOriginal, RECALLINPUTNOISE);
		recallAgainst(storedOriginal, RECALLINPUTNOISE);
		
		exit(0);
	}

	if(0) { // Recall with many random patterns to asses basin sizes
		int nn = 0, u, j, dummy = 100, hd[PN];
		uint64_t i, stored[PN] = {0}, iter = 100000;
		float v[N];
		
		//train with dummy patterns to initialize weight matrix
		for(i = 0; i < 1000; i++) {
			randomPattern(v, 0.5);
			for(u = 0; u < LEARNNUM; u++) train_network(v, nn);
		}	
		for(u = 0; u < LEARNNUM; u++) for(j = 0; j < PN; j++) train_network(storedOriginal[nn][j], nn);
	
		for(i = 0; i < iter; i++) {			
			randomPattern(output[nn], 0.5);
			for(u = 0; u < UPDATENUM; u++) update_network(nn);
			for(j = 0; j < PN; j++) {
				hd[j] = hammingDistance(output[nn], storedOriginal[nn][j]);
				stored[j] += (hd[j] == 0);
				//printf("%d ", hd[j]);
				//fflush(stdout);
			}
		}
		
		printf("\n");
		for(j = 0; j < PN; j++) printf("%"PRIu64" ", stored[j]);
		printf("   (of %"PRIu64") \n", iter);
		
		exit(0);
	}
	
	if(0) { // Recall with increasing load (Fig. 1 at Storkey XXX)
		int nn = 0, u, i, j, b = 0, c;
	  float d;
		
		for(i = 0; i < PN; i++) for(u = 0; u < LEARNNUM; u++) train_network(storedCorrelated[nn][i], nn);
		
		for(i = 0; i < PN; i++) {
			c = 0;
			for(j = 0; j < TESTITERATION; j++) {
				copyPattern(output[nn], storedCorrelated[nn][i]);
				exactlyMutatePattern(output[nn], TESTNOISE);
				for(u = 0; u < UPDATENUM; u++) update_network(nn);
				d = relativeHammingDistance(storedCorrelated[nn][i], output[nn]);
				c += (d <= TESTTOLERANCE);
				//printf("\t#%d %f %d\n", i, d, c);
			}
			b += (c >= TESTGOOD);
		}
		
		printf("%d %d\n", PN, b);
		
		exit(0);
	}
	
	if(0) { // Retraining
		int u, i, j, k;
	
		for(u = 0; u < LEARNNUM; u++) for(i = 0; i < NN; i++) for(j = 0; j < PN; j++) train_network(storedOriginal[i][j], i);
		recallAgainst(storedOriginal, 0.0);
		recallAgainst(storedOriginal, RECALLINPUTNOISE);
		
		if(VERBOSE) printf("\nDUMMY TRAINING...\n\n\n");
		for(k = 0; k < 100; k++) {
			resetStoredPatterns();
			for(u = 0; u < LEARNNUM; u++) for(i = 0; i < NN; i++) for(j = 0; j < PN; j++) train_network(storedRandom90[i][j], i);
		}
		
		resetStoredPatterns();
		for(u = 0; u < LEARNNUM; u++) for(i = 0; i < NN; i++) for(j = 0; j < PN; j++) train_network(storedRandom90[i][j], i);
		recallAgainst(storedRandom90, RECALLINPUTNOISE);
		exit(0);
	}
	
	if(0) { // Testing attractor properties of outputs of a single network (exact test, samePatternQ)
		float in[N], mut = 0.1;
		int u, i, j, nn = 0, c = 0, iter = 100;
		
		printf("TEST\n");
		for(u = 0; u < LEARNNUM; u++) for(i = 0; i < NN; i++) for(j = 0; j < PN; j++) train_network(storedOriginal[i][j], i);

		for(i = 0; i < iter; i++) {
			//randomPattern(in, 0.5);
			copyPattern(in, storedOriginal[randl(NN)][randl(PN)]);
			c += storedPatternQTest(in, storedOriginal, nn, mut);
		}
		printf("M:%f   STORED && STABLE: %d/%d\n", mut, c, iter);
		fflush(stdout);
		exit(0);
	}
	
	if(0) { // Testing attractor properties of outputs of a single network (loose test, relativeHammingDistance < 0.05)
		float v[N];
		int e, u, i, j, nn = 0, c = 0, iter = 100, dummy = 500, epochs = 2;
		
		for(e = 0; e < epochs; e++) {
			
			if(e > 0) {
				c = 0;
				for(u = 0; u < LEARNNUM; u++) for(j = 0; j < PN; j++) train_network(storedOriginal[nn][j], nn);
				for(i = 0; i < iter; i++) {
					//randomPattern(output[nn], 0.5);
					copyPattern(output[nn], storedOriginal[nn][randl(PN)]);
					//copyPattern(output[nn], storedOriginal[nn][i]);
					for(j = 0; j < UPDATENUM; j++) update_network(nn);
					copyPattern(v, output[nn]);
					//c += testAttractorQ(v, nn);
					c += storedAttractorQ(v, storedOriginal, nn);
				}
				printf("\nSTORED ATTRACTORS: %d/%d   (i: %d/%d, m: %.3f, t: %.3f)\n", c, iter, TESTGOOD, TESTITERATION, TESTNOISE, TESTTOLERANCE);
				fflush(stdout);
			}

			if(epochs > 1) for(i = 0; i < dummy; i++) {
				randomPattern(v, 0.5);
				for(u = 0; u < LEARNNUM; u++) train_network(v, nn);
			}
		}
			
		exit(0);
	}

	if(0) { // find attractors of a network: test network with random input and each found to be an attractor is tested against 1) stored ones if it is a stored attractor or 2) against all previous non-stored ("extra") attractors, if it is a non-stored one.
		int i, u, nn = 0, dummy = 1000;
		float v[N];
		
		for(i = 0; i < dummy; i++) {
			randomPattern(v, 0.5);
			for(u = 0; u < LEARNNUM; u++) train_network(v, nn);
		}
			
		for(u = 0; u < LEARNNUM; u++) for(i = 0; i < PN; i++) train_network(storedOriginal[nn][i], nn);
		attractorStatistics(nn);
		
		exit(0);
	}
	
	if(1) { // Multiple networks searching for global optimum
		int i, j, k, u, t = 0, terminateSum = 0, bestNet, worstNet, dummy = 500, giveupC = 0;
		float bestOutput[N], input[NN][N], v[N], maxW = -99., bestW, lastW, worstW, w[NN];
		
		//printf("Multiple networks searching for global optimum\n");
		
		

		
		
		// Landscaping
		// Pearson correlation
	  randomPattern(globalBest, 0.5); 
		// Nk landscape
		randomPattern(globalBest, 1.0); // [111..1]
		for(i = 0; i < K+1; i++) MAP[0] = MAP[0] + pow(2, N-i-1); // generate epistatic map; here a circular sequence is used with K neighbours to the right
		for(i = 0; i < NKP; i++)   MAP[i] = rotateRight(MAP[0], i);
		// Block landscape
		randomPattern(globalBest, 1.0); // [111..1]
		setBlockOptimaDefault(); // generate optima
		blockGlobalOptimumSequence(v);
		blockGlobalMax = blockGlobalOptimumFitness();
		if(1) {
			printf("Block target optima ({N=%d P=%d B=%d):\n", N, P, B);
			for(i = 0; i < B; i++) {
				printf("BLOCK #%d\n", i);
				for(j = 0; j < T; j++) {
					printf("\t");
					for(k = 0; k < P; k++) printf("%s", (blockOptimaSequence[i][j][k] > 0)?"+":"-");
					printf("\t%f\n", blockOptimaFitness[i][j]);
				}
			}
			printf("Global optimum:\n");
			printPattern(v);
			printf("\t%f\n\n", blockGlobalMax);
		}
		
		
		if(AANN) { // train with dummy patterns to initialize weight matrix
			for(i = 0; i < NN; i++) for(j = 0; j < dummy; j++) {
				randomPattern(v, 0.5);
				for(u = 0; u < LEARNNUM; u++) train_network(v, i);
			}
		}
	
	
		// start with the same random pattern as input for all networks
		//randomPattern(v, 0.5);
		//randomPattern(v, 0.0);// [000..0] as worst for initial input
		//for(i = 0; i < NN; i++) copyPattern(input[i], v);
		for(i = 0; i < NN; i++) randomPattern(input[i], 0.5);
		
		
		// Initial setup
		for(i = 0; i < NN; i++) w[i] = FITNESS(input[i]);
		bestNet  = randomMaxPosition(w, NN); // NOTE: this might chose a position that is not an attractor, however there might be another that is.
		bestW    = w[bestNet];
	  worstNet = randomMinPosition(w, NN); // NOTE: this might chose a position that is not an attractor, however there might be another that is.
		worstW   = w[worstNet];
		lastW    = bestW;
		printf("%d %f %d %d %d\n", t, bestW, -1, -1, -1);
		fflush(stdout);
		t++;
		
		
		// generations
		while(t < MAXT && terminateSum < ENDFRAME) {
			int attractorQ = 0, storedQ = 0, storedLastQ = 0, recombineQ = 0, learnQ = 0;
			

		// if((t > 0) && !(t % 1000)) { // test attractor statistics
			// printf("T%d ", t);
			// for(i = 0; i < NN; i++) attractorStatistics(i);
		// }
			
			// update network
			for(i = 0; i < NN; i++) {
				copyPattern(output[i], input[i]);
				if(AANN) for(u = 0; u < UPDATENUM; u++) update_network(i); // if no AANN is used, output will be simply the input
			}
						
			//find best solution (chose one randomly, if there are multiple with identical HD-s)
			for(i = 0; i < NN; i++) {
				w[i] = FITNESS(output[i]);
				//printf("\t#%d\t%f\t", i, w[i]); printPattern(output[i]); if(i == (NN-1)) printf("\n");
			}
			bestNet  = randomMaxPosition(w, NN); // NOTE: this might chose a position that is not an attractor, however there might be another that is.
			bestW    = w[bestNet];
			worstNet = randomMinPosition(w, NN);
			worstW   = w[worstNet];
			copyPattern(bestOutput, output[bestNet]);
			

			if(AANN && MONITOR) { // monitor all outputs for attractorness and other properties
				int hd, storedQSum = 0, storedNotLastQSum = 0;
				int storedMinHD[NN], attractorQList[NN] = {0}, storedIList[NN] = {0}, storedNList[NN] = {0}, storedQList[NN] = {0}, storedDist[NN];
				float test[N], orig[N];
				
				for(i = 0; i < NN; i++) {
					storedIList[i] = -1;
					storedMinHD[i] = N;
					storedDist[i]  = N;
					
					// test whether output is an attractor
					copyPattern(orig, output[i]);
					copyPattern(test, output[i]);
					attractorQList[i] = testAttractorQ(test, i);
					
					// test whether output is a stored (learnt) pattern; checking all and counting matches
					for(j = 0; j < PN; j++) {
						hd = hammingDistance(orig, storedChanged[i][j]);
						if(hd < storedMinHD[i]) storedMinHD[i] = hd;
						if(samePatternQ(orig, storedChanged[i][j])) {
							storedNList[i]++;
							storedIList[i] = j;
						}
					}
					storedQList[i]    =  storedNList[i] > 0;
					storedQSum        += storedQList[i];
					storedNotLastQSum += storedQList[i] && (storedIList[i] != (PN-1));
					storedDist[i]     =  hammingDistance(orig, globalBest);
				}
				
				attractorQ  = attractorQList[bestNet];
				storedQ     = storedQList[bestNet]; // independent of `attractorQ`: could equal to stored pattern but unable to recall it multiple times correctly
				storedLastQ = (storedIList[bestNet] == (PN-1));
				
				if(1) {
					printf("\tT = %d\n", t);
					printf("\tINDEX      [");
					for(i = 0; i < NN; i++) if(i == bestNet) printGreen("%2d ", i); else printf("%2d ", i);
					printf("\b]   best: \t#%2d\n", bestNet);
					printf("\tATTRACTOR? [");
					for(i = 0; i < NN; i++) if(i == bestNet) printGreen("%2d ", attractorQList[i]); else printf("%2d ", attractorQList[i]);
					printf("\b]   attractor?\t%d\n", attractorQList[bestNet]);
					printf("\tMIN HD     [");
					for(i = 0; i < NN; i++) if(i == bestNet) printGreen("%2d ", storedMinHD[i]); else printf("%2d ", storedMinHD[i]);
					printf("\b]   min HD\t%d\n", storedMinHD[bestNet]);
					printf("\tSTORED?    [");
					for(i = 0; i < NN; i++) if(i == bestNet) printGreen("%2d ", storedQList[i]); else printf("%2d ", storedQList[i]);
					printf("\b]   stored?\t%d   (sum: %d)\n", storedQList[bestNet], storedQSum);
					printf("\tSTORED IND [");
					for(i = 0; i < NN; i++) if(i == bestNet) printGreen("%2d ", storedIList[i]); else printf("%2d ", storedIList[i]);
					printf("\b]   index: \t%d\n", storedIList[bestNet]);
					printf("\tSTORED CNT [");
					for(i = 0; i < NN; i++) if(i == bestNet) printGreen("%2d ", storedNList[i]); else printf("%2d ", storedNList[i]);
					printf("\b]   count: \t%d\n", storedNList[bestNet]);
					printf("\tGLOBAL HD  [");
					for(i = 0; i < NN; i++) if(i == bestNet) printGreen("%2d ", storedDist[i]); else printf("%2d ", storedDist[i]);
					printf("\b]   HD:    \t%d   w: %f\n", storedDist[bestNet], bestW);
				}
			} else if(AANN) { // monitor only the output of the best net
			  int storedI = -1;
				float orig[N];
				
				copyPattern(orig, output[bestNet]);
				storedQ = 0;
				for(j = 0; j < PN; j++) if(samePatternQ(orig, storedChanged[bestNet][j])) {
						storedQ = 1;
						storedI = j;
					}
				storedLastQ = (storedI == (PN-1));
				attractorQ = testAttractorQ(orig, bestNet);
			}

			
			if(AANN && randd() < 0.9) { // Retrain network(s) with mutated output
				// NOTE: the strict condition `bestW > maxW` usually makes the system stuck in an attractor
				// if(randd() < .5) { // learn only half of the time
				// if(1) { // learn always
				// if(bestW > maxW) { // learn only when better
				int pos[NN] = {0};
				learnQ = 1;
				for(k = 0; k < TRAINNUM; k++) { // the best pattern is assured to be trained to *different* networks! (converges faster than retraining the same network multiple times)
					do i = randl(NN); while(pos[i]);
					pos[i] = 1;
					copyPattern(v, bestOutput);
					mutatePattern(v, LEARNNOISE);
					for(u = 0; u < LEARNNUM; u++) train_network(v, i);
					
					// Store learnt patterns in `storedChanged`; only store if there is no identical already stored.
					// If there is an identical, remove it, shift those behind it one step to the left and store the new one at the last position
					j = PN-1;
					while(!samePatternQ(v, storedChanged[i][j]) && j >= 0) j--;
					if(j < 0) j = 0; // if already present,  start shift at position `j`, if not present, start shift at position 0.
					for( ; j < PN-1; j++) copyPattern(storedChanged[i][j], storedChanged[i][j+1]);
					copyPattern(storedChanged[i][PN-1], v);
				}
			}
			
						
			/*
			// Monitor best fitness
			if(bestW > maxW) maxW = bestW;
			giveupC += (bestW == lastW);
			lastW = bestW;
			
		
		  // Mutate output population
			// Recombine the best of the outputs with a random other
			// at an inner point chosen randomly from uniform distribution
			// if `givup` reaches 5, the input of the next generation will be the recombinant
			// only carries on `bestOutput` and discard partner
			if(giveupC >= GIVEUP || (randd() < REC)) {
				giveupC = 0;
				//randomPattern(otherOutput, 0.5);
				//copyPattern(bestOutput, other2Output);
				recombinePattern(bestOutput, otherOutput);
				recombineQ = 1;
			}
			for(i = 0; i < NN; i++) {
				copyPattern(input[i], bestOutput);
				mutatePattern(input[i], INPUTNOISE);
			}
			*/			
			
			if(!AANN) { // Selection (recombination or mutation)
				float newW, new[N];
				
				for(i = 0; i < NN; i++) copyPattern(input[i], output[i]);
				if(randd() < REC) {
					int p1 = randl(NN), p2 = randl(NN);
					float w1, w2;
					while(p2 == p1) p2 = randl(NN);
					recombinePattern(output[p1], output[p2]);
					w1 = FITNESS(output[p1]);
					w2 = FITNESS(output[p2]);
					if(w1 >= w2) {
						copyPattern(new, output[p1]);
						newW = w1;
					}	else {
						copyPattern(new, output[p2]);
						newW = w2;
					}
					recombineQ = 1;
				} else {
					int p = randl(NN);
					copyPattern(new, output[p]);
					mutatePattern(new, INPUTNOISE);
					newW = FITNESS(new);
				}
				if(newW > worstW) {
					copyPattern(input[worstNet], new);
					//printf("\tREP\t%f -> %d\t", FITNESS(new), worstNet); printPattern(input[worstNet]);
					//printf("\n");
				}
			}
		
		
		
			if(bestW == 1.0) terminateSum++; else terminateSum = 0;
	
			
			if(!(t % RESOLUTION)) {
				printf("%d %f %d %d %d\t", t, bestW, attractorQ, storedQ, storedLastQ);
				printPattern(output[bestNet]);
				fflush(stdout);
			}
			
			t++;
		}
	}
	
	
	
  return(0);
}
