// Deme-besed metapopulational search
// by István Zachar
// patterns are represented as (-1.0, 1.0) floats
// version 1


#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>

#include "randomGenerator.c" // own RNG library (same as András's, with minor bugfix and extra functionality) 
#include "ran4.c"            // for Nk's pseudorandom generator; ran4 should only be used for Nk, not for any other random number, as it would mess up the seed's state
#include "colors.c"          // for monitoring purposes only


#define SEED 227063

// 0 = unchanged input
// 1 = mutated output (P_mut = 1/N)
// 2 = correlated output 
// 3 = random output
// 4 = autoassociative attractor neural network
// 5 = autoassociative attractor neural network, with no spurious pattern
#define TYPE 4 // type of operator to turn input to output

// Fitness landscape & problem size
#define N       20  // number of neurons (200) // if Nk landscape is used, the maximal value of N is 63!
#define P       10   // partition size; must be an integer that divides N
#define T       2   // number of optima per block
#define B       N/P // number of blocks per sequence
#define FITNESS blockFitness //(blockFitness) // fitness function (1.0 is best, minimum need not be 0.0): pearsonFitness, hammingFitness, NkFitness, blockFitness
#define F       N   // Nk landscape phenotype length (generally safe to have N = P
#define K       6   // Nk landscape number of epistatic neighbours (uniform for each gene; subsequence length is K+1)

// Demes
#define SPATIAL 1 // use spatial lattice?
#define DX 10     // 1st dimension of spatial deme lattice
#define DY 10     // 2nd dimension of spatial deme lattice
#define DN DX*DY  // total deme number

// Multinetwork population
#define NN            10    // population size of a deme (number of networks)
#define ELITESEL      0    // select best from population OR select one randomly for replication
#define RECOMBINATION 0.1 //(0.1)    // probability of recombining two output sequences; if not recombination, migration
#define MIGRATION     0.004 //(0.004)   // probability of migration; only if DN > 1
#define MUTATION      1.0/(float)N // per-bit mutation rate of re-entered output of previous generation
#define LEARNING      1.    // probability of learning

// Simulation time
#define MAXT       10000  // number of generations
#define RESOLUTION 1
#define VERBOSE    0
#define ENDFRAME   500 //(RESOLUTION)

// Autoassociative attractor neural networks
#define UPDATENUM  5 // (100) // number of update steps (100)
#define THRESHOLD  0. // firing threshold
#define TRAINNUM   5    // number of networks to be trained with the same output (though differently noised) (only if TYPE >= 4)
#define LEARNNOISE 0. // per-bit mutation rate of re-learned output of previous generation (only if TYPE >= 1)

// Attractor network testing (only if TYPE >= 5)
#define TESTSTORED    0   // test storedness of output (always on if TYPE >= 5)
#define PN            N   // maximal number of stored patterns (N/sqrt(2*ln(N)) * 1.1); only used for tracking learned patterns in some tests
#define TESTTOLERANCE 0. // accept equivalence if output is closer than (or equal to) tolerance (in relative Hamming distance)
#define TESTNOISE     0.05// noise applied for iterative reentry of the same input



// Landscapes
float globalBest[N];

// Nk landscape
// const uint64_t MAXN = (1ULL << (N  )); // size of sequence space (2^N), max decimal value of a sequence is MAXN-1
// const uint64_t MAXK = (1ULL << (K+1)); // size of subsequence space (2^(K+1)), max decimal value of a subsequence of length K+1 is MAXK-1
// const uint64_t MAXB = (1ULL << (N-1)); // largest sequence with a single 1 bit (2^(N-1), represents the sequence [10000..0]; for printing purposes)
// uint64_t MAP[F] = {0};    // epistatic map, uniform for each sequence
unsigned long RAN4SEED;

// block landscape
float blockOptimaSequence[B][T][P];
float blockOptimaFitness [B][T];
float blockGlobalMax;

typedef struct DEME {
  float w[NN], minW, maxW, input[NN][N], output[NN][N], maxOutput[N], weight[NN][N][N];
	int   minP, maxP, storedQ[NN], storedP[NN];
	float storedD[NN], stored[NN][PN][N];
} DEME;

DEME deme[DN];

// spatial deme neighbors
int matrixNeighbors[DN][8];
int STORE = (TYPE >= 5) || ((TYPE == 4) && TESTSTORED); // store learned patterns in `stored` array






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

void setWeights(void) {
  int d, i, j, k;
  for(d = 0; d < DN; d++) for(i = 0; i < NN; i++) for(j = 0; j < N; j++) for(k = 0; k < N; k++) deme[d].weight[i][j][k] = 0.0;
}

void storePattern(float *v, int d, int i) { // tests if `v` is stored in network `i` at deme `d` or not, and pushes/stores it to most recent position (on the right)
  // If there is an identical in `stored`, remove it, shift those behind it one step to the left and store the new one at the last position
	int j = PN-1;
	while(!samePatternQ(v, deme[d].stored[i][j]) && j >= 0) j--;
	if(j < 0) j = 0; // if already present, start shift at position `j`, if not present, start shift at position 0
	for( ; j < PN-1; j++) copyPattern(deme[d].stored[i][j], deme[d].stored[i][j+1]);
	copyPattern(deme[d].stored[i][PN-1], v);
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

void trainNetworksRandom(int n) { // train each network on `n` random patterns and store these patterns
	int d, i, j;
	//float v[N];
	for(d = 0; d < DN; d++) for(i = 0; i < NN; i++) for(j = 0; j < n; j++) {
		// randomPattern(v, 0.5);
		// trainNetwork(v, d, i);
		randomPattern(deme[d].stored[i][j], 0.5);
		if(STORE) storePattern(deme[d].stored[i][j], d, i);
		trainNetwork(deme[d].stored[i][j], d, i);
	}
}

void updateOutput(float *v, int d, int n, int i) { // update of neuron `i` in network `n` in deme `d` with threshold neuron model
   int j;
   float h = 0.0;
   for(j = 0; j < N; j++) if(j != i) h += (deme[d].weight[n][i][j]) * v[j];
   if(h >= THRESHOLD) v[i] =  1.0;
   else               v[i] = -1.0;
}

void updateNetwork(float *v, int d, int n) { // updating network 'n' in deme `d`
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


// Nk landscape
/*
uint64_t rotateRight(uint64_t num, int step) { // bitwise rotation to the right (circular)
	int s = step % N;
	return((num >> s) | ((num << (N-s)) % MAXN));
}

uint64_t rotateLeft(uint64_t num, int step) { // bitwise rotation to the left (circular)
	int s = step % N;
	return(((num << s) % MAXN) | (num >> (N-s)));
}

uint64_t fromBinaryVector(float *v) { // converts binary vector `v` of {-1, +1} to a binary integer
	int i;
	uint64_t num = 0;
	for(i = 0; i < N; i++) num += rotateRight(MAXB, i)*(v[i] > 0); // NOTE: Consider how `v` is represented! {0, 1} or {-1, +1}?
	return(num);
}

double NkFitness(float *v, float *u) { // Nk fitness of `v`, with [111..1] having w=1.0 and [000..0] having w=0.0 fitness; expecting a binary float vector `v`; ignores `u`
	int i;
	uint64_t seq, nth, num;
	double sum = 0.;
	
	num = fromBinaryVector(v);
	for(i = 0; i < F; i++) {
		seq = rotateLeft(num & MAP[i], K+1+i); // extract epistatic neighborhood of size `K+1` at position `i` from sequence `num`, according to the epistatic map `MAP`
		if((0 < seq) && (seq < MAXK-1)) {
			nth = seq + i*MAXK;
			sum += seededNthRan4(RAN4SEED, (unsigned long)nth);
		} else if(seq == MAXK-1) {
			sum += 1.0; // subsequence is of best sequence [111..1], has fitness contribution 1.0;
		} // else do not increment `sum`
		//printf("\t#%d\t%"PRIu64"\t%"PRIu64"\t%f\n", i, seq, i*MAXK, sum);
	}
	return(0.1 + 0.9*(sum/(double)F)); // NOTE
}
*/

// Block landscape

void setBlockOptimaDefault(void) { // sets the same T optima (sequence and fitness) for the B blocks: first is [1111...], rest is [0101...]
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
	return(sum/blockGlobalMax);
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

	seed(SEED); // seed RNG of randomGenerator.c with `unsigned long`
	RAN4SEED = randl(LONG_MAX); // seed for Nk algorithm (`ran4` seed, required so that `ran4` depends on the session, i.e. on the overall randomizer state set by `SEED`)

	
	if(1) { // Multiple networks searching for global optimum
		int i, j, k, d, t = 0, terminateSum = 0;
		

		// Landscaping
		/* Pearson correlation */
	  //randomPattern(globalBest, 0.5); 
		/* Nk landscape */
		// randomPattern(globalBest, 1.0); // [111..1]
		// for(i = 0; i < K+1; i++) MAP[0] = MAP[0] + pow(2, N-i-1); // generate epistatic map; here a circular sequence is used with K neighbours to the right
		// for(i = 0; i < F;   i++) MAP[i] = rotateRight(MAP[0], i);
		/* Block landscape */		
		setBlockOptimaDefault(); // generate optima
		blockGlobalOptimumSequence(globalBest);
		blockGlobalMax = blockGlobalOptimumFitness();
		if(VERBOSE) {
			printf("Block target optima ({N=%d P=%d B=%d T=%d):\n", N, P, B, T);
			for(i = 0; i < B; i++) {
				printf("BLOCK #%d\n", i);
				for(j = 0; j < T; j++) printf("\tT%d %s %f\n", j, vectorToString(blockOptimaSequence[i][j], P), blockOptimaFitness[i][j]);
			}
			printf("Global:\t%s %f\n\n", vectorToString(globalBest, N), blockGlobalMax);
		}
		
		
		// Initial setup
		if(SPATIAL) setNeighborDemes(DX, DY);
		for(d = 0; d < DN; d++) for(i = 0; i < NN; i++) randomPattern(deme[d].output[i], 0.5);
		
		// Network setup
		if(TYPE >= 4) {
			setWeights();
			trainNetworksRandom(PN);
		}
		
		
		// Generations
		while(t < MAXT && terminateSum < ENDFRAME) {
			int bestD = -99, bestP = -99, recombineQ = 0, migrateQ = 0, learnQ = 0;
			float bestW = -99., storedC = 0.0;
			
			
			for(d = 0; d < DN; d++) {
				
				/* Update output based in new input */
				if(t > 0) {
					shuffleOutputs(d);
					if(TYPE >= 5) for(i = 0; i < NN; i++) copyPattern(deme[d].input[i], deme[d].output[i]);
					
					if(TYPE == 0) {                                                                  } else // 0 = unchanged input
					if(TYPE == 1) {for(i = 0; i < NN; i++) mutatePattern(   deme[d].output[i], 1.0/(float)N);} else // 1 = mutated output
					if(TYPE == 2) {for(i = 0; i < NN; i++) correlatePattern(deme[d].output[i], 1.0 );} else // 2 = correlated output
					if(TYPE == 3) {for(i = 0; i < NN; i++) randomPattern(   deme[d].output[i], 0.5 );} else // 3 = random output
					if(TYPE >= 4) {for(i = 0; i < NN; i++) updateNetwork(   deme[d].output[i], d, i);}      // 4 = AANN
																																																	// 5 = AANN with no spurious output (iterate until stored)
				}
					
				for(i = 0; i < NN; i++) {
					
					/* Test storedness */
					if(STORE) { 
						float hd[PN], minhd = 1.0;
						
						for(j = 0; j < PN; j++) {
							hd[j] = relativeHammingDistance(deme[d].output[i], deme[d].stored[i][j]);
							if(hd[j] < minhd) minhd = hd[j];
						}
						
						while ((TYPE >= 5) && (minhd > TESTTOLERANCE)) { // repeat input->output if output is not stored
							copyPattern(deme[d].output[i], deme[d].input[i]);
							mutatePattern(deme[d].output[i], TESTNOISE);
							updateNetwork(deme[d].output[i], d, i);
							for(j = 0; j < PN; j++) {
								hd[j] = relativeHammingDistance(deme[d].output[i], deme[d].stored[i][j]);
								if(hd[j] < minhd) minhd = hd[j];
							}
						}
						
						deme[d].storedP[i] = randomMinPosition(hd, PN);
						deme[d].storedD[i] = hd[deme[d].storedP[i]];
						deme[d].storedQ[i] = (deme[d].storedD[i] <= TESTTOLERANCE);
						storedC += deme[d].storedQ[i];						
					}
					
					deme[d].w[i] = FITNESS(deme[d].output[i]);
				}
				
				/* Find best solution (chose one randomly, if there are multiple with identical HD-s) */
				deme[d].maxP = randomMaxPosition(deme[d].w, NN); // NOTE: this might chose a position that is not an attractor, however there might be another that is.
				deme[d].maxW = deme[d].w[deme[d].maxP];
				deme[d].minP = randomMinPosition(deme[d].w, NN); // NOTE: this might chose a position that is not an attractor, however there might be another that is.
				deme[d].minW = deme[d].w[deme[d].minP];
				copyPattern(deme[d].maxOutput, deme[d].output[deme[d].maxP]);
				if(deme[d].maxW > bestW) {
					bestD = d;
					bestW = deme[d].maxW;
					bestP = deme[d].maxP;
				}
			}
			
			
			
			if(VERBOSE) {
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
						recombineQ++;
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
							migrateQ++;
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
						mutatePattern(v1, MUTATION);
						w1 = FITNESS(v1);
					}
					
					/* Replication & learn */
					if((w1 > deme[d].minW)) {

						copyPattern(deme[d].output[deme[d].minP], v1); // replace ONE of the worst patterns in the output
						//for(i = 0; i < NN; i++); copyPattern(deme[d].output[i], v1); // replace ALL of the output population with the new output
						
						if((TYPE >= 4) && (randd() < LEARNING)) { // Learning
							int pos[NN] = {0};
							float toLearn[N];
							for(k = 0; k < ((TRAINNUM>NN)?NN:TRAINNUM); k++) { // the best pattern is assured to be trained to TRAINNUM *different* networks (faster convergence)
								learnQ++;
								do i = randl(NN); while(pos[i]);
								pos[i] = 1;
								copyPattern(toLearn, v1);
								mutatePattern(toLearn, LEARNNOISE);
								trainNetwork(toLearn, d, i);
								if(STORE) storePattern(toLearn, d, i); // Store learnt patterns in `stored` only if there is no identical already
							}
						}
					}			
				} // deme `d`
			}
		
		
			//if(VERBOSE && (bestW == 1.0)) terminateSum++; else terminateSum = 0;
			if((bestW == 1.0)) terminateSum++; else terminateSum = 0;
	
			
			if(!(t % RESOLUTION)) {
				//printf("%d\t%f\tR%d\tM%d\tL%d\t%d.%d.%s\n", t, bestW, recombineQ, migrateQ, learnQ, bestD, bestP, vectorToString(deme[bestD].output[bestP], N));
				printf("%d\t%f\t%d\t%d\t%f\n", t, bestW, deme[bestD].storedQ[bestP], learnQ, (float)storedC/(float)(DN*NN));
				fflush(stdout);
			}
			
			t++;
		}
	}
	
	
	
  return(0);
}