#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

/******* RANDOM NUMBER GENERATOR **********/
#define AA 471
#define B 1586
#define CC 6988
#define DD 9689
#define M 16383
#define RIMAX 2147483648.0        /* = 2^31 */
#define RandomInteger (++nd, ra[nd & M] = ra[(nd-AA) & M] ^ ra[(nd-B) & M] ^ ra[(nd-CC) & M] ^ ra[(nd-DD) & M])
void seed(long seed);
static long ra[M+1], nd;


void seed(long seed)
{
 int  i;

 if(seed<0) { puts("SEED error."); exit(1); }
 ra[0]= (long) fmod(16807.0*(double)seed, 2147483647.0);
 for(i=1; i<=M; i++)
 {
  ra[i] = (long)fmod( 16807.0 * (double) ra[i-1], 2147483647.0);
 }
}

long randl(long num)      /* random number between 0 and num-1 */
{
 return(RandomInteger % num);
}

double randd(void)
{
 return((double) RandomInteger / RIMAX);
}
/****************************************/



#define N (200)  // number of neurons
#define NN (20) // number of networks
#define PN (50) // depth of memory
#define Th (0.0) // firing threshold
#define NOUP (100) // number of update steps
#define mu_l (0.005) // mutation rate of best output (for input) 

#define GENERATIONS (10) // number of generations

float pattern[NN][PN][N];
float spattern[NN][N];
float output[NN][N];
float weight[NN][N][N];
float optvect[N];
int optnum;
unsigned long int zeed = 0;

void init(void)
{
  int i,j,k;
  
  for(k=0;k<NN;k++)
    for(i=0;i<N;i++)
      for(j=0;j<N;j++)
        weight[k][i][j]=0.0;
}

void piksr2(int n,float *arr, float *brr)
{
  int i,j;
  float a,b;

  for (j=2;j<=n;j++)
  {
    a=arr[j];
    b=brr[j];
    i=j-1;
    while (i > 0 && arr[i] > a)
    {
      arr[i+1]=arr[i];
      brr[i+1]=brr[i];
      i--;
    }
    arr[i+1]=a;
    brr[i+1]=b;
  }
}


// set the stored (random) patterns
void set_stored_pattern(void)
{
  int m,n,k;

  for(k=0;k<NN;k++)
  {
     for(m = 0; m < PN; m++)
     {
       for(n = 0; n < N; n++)
       {
         if(randd()>0.5)
	   pattern[k][m][n]=-1.0;
	 else
	   pattern[k][m][n]=1.0;

       }
    }   
  }
}


// set the special patterns

void set_special_patterns(void)
{
  int n,k;

  for(k=0;k<NN;k++)
  {
     for(n = 0; n < N; n++)
     {
       if(n<k*N/(NN-1))
         spattern[k][n]=1.0;
       else
         spattern[k][n]=-1.0;
       }
  }
}


// training network 'm' with vector vec -- palimpsest
void train_network_pali(int m, float *vec)
{
  int i,j;
  float h[N];
  
  for(i=0;i<N;i++)
  {
    h[i]=0.0;
    for(j=0;j<N;j++)    
      h[i]+=weight[m][i][j]*vec[j];
  }

  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      if(i==j)
	weight[m][i][j]=0.0;
      else
        weight[m][i][j]+=1.0/((float)N)*vec[i]*vec[j]-1.0/((float)N)*vec[i]*h[j]-1.0/((float)N)*vec[j]*h[i];
    }
  }
}

float relat_Hamm(float* in, float* ou)
{
  int t;
  float rhd=0;
  
  for(t=0;t<N;t++)
    if(in[t]!=ou[t])
      rhd=rhd+1.0;
    
  rhd=rhd/(float)N;
  
  return(1.0-rhd);
  
}

// update of neuron 'n' in network 'm'
void update_output(int m, int n)
{
   int t;
   float h=0.0;
   
   for(t = 0; t < N; t++)
     if(t != n )
       h += weight[m][n][t] * output[m][t];

   if(h>=Th)
      output[m][n] =  1.0; // linear threshold
   else
      output[m][n] = -1.0;

}

// updating network 'm'
void update_net(int m) 
{
   int t;

   for(t=0;t<N;t++)
     update_output(m, randl(N));
}


int main(int argc, char** argv)
{
 
   int k,m,n,up,rou,rr,pro;
   float vect[N],pc[NN+1],pcopt,ref[NN+1];

   seed(13788);

   init();

   set_stored_pattern();
   set_special_patterns();

// train all networks with PN-1 random patterns
   for(m=0;m<NN;m++)
     for(k=0;k<PN-1;k++)
       train_network_pali(m,pattern[m][k]);

// train networks with special patterns
   for(m=1;m<NN;m++)
     train_network_pali(m,spattern[m]);


//start with the same random pattern (similar to the worst) for all
  for(n=0;n<N;n++)
  {
    if(randd()<0.9)
      vect[n]=-1.0;
    else
      vect[n]=1.0;
  }

  for(m=0;m<NN;m++)
    for(n=0;n<N;n++)
      output[m][n]=vect[n];    

  printf("Round\t# of best net\tfitness\n");
  for(rou=1;rou<=GENERATIONS;rou++)
  {
    for(m=0;m<NN;m++)
      for(up=0;up<NOUP;up++)
        update_net(m);
    
      
//sort and find the best solution
    optnum=-1;    
    for(m=0;m<NN;m++)
    {
      ref[m+1]=m;
      for(n=0;n<N;n++)
	vect[n]=output[m][n];
      pc[m+1]=relat_Hamm(vect,spattern[NN-1]);
    } 

    piksr2(NN,pc,ref);
   
    for(m=NN;m>0;m--)
      if(pc[NN]>pc[m])
        break;
      
    rr=randl(NN-m);
    pro=(int)ref[NN-rr];
   
    pcopt=pc[NN];
    optnum=pro;

    
  
// elitism
   for(n=0;n<N;n++)
   {
     vect[n]=output[optnum][n];
     optvect[n]=vect[n]; 
   }


// set inputs for next round
    for(m=0;m<NN;m++)
    {
      for(n=0;n<N;n++)
      {
        output[m][n]=vect[n];
        if(randd()<mu_l)
          output[m][n]=-1.0*output[m][n];
      }
    }
    printf("%d\t%d\t%f\n",rou,optnum,pcopt);
    fflush(stdout);

  }

  return(0);
}
