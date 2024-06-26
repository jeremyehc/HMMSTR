/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   15 December 1997
**      File:   hmm.h
**      Purpose: datastructures used for HMM.
**      Organization: University of Maryland
**
**	Update:
**	Author: Tapas Kanungo
**	Purpose: include <math.h>. Not including this was
**		creating a problem with forward.c
**      $Id: hmm.h,v 1.9 1999/05/02 18:38:11 kanungo Exp kanungo $
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <glib.h> //added

typedef struct {
	int N;		/* number of states;  Q={1,2,...,N} */
	int M; 		/* number of observation symbols; V={1,2,...,M}*/
	double	**A;	/* A[1..N][1..N]. a[i][j] is the transition prob
			   of going from state i at time t to state j
			   at time t+1 */
	double	**B;	/* B[1..N][1..M]. b[j][k] is the probability of
			   of observing symbol k in state j */
	double	*pi;	/* pi[1..N] pi[i] is the initial state distribution. */
} HMM;
void ReadHMM(FILE *fp, HMM *phmm);
void PrintHMM(FILE *fp, HMM *phmm);
void InitHMM(HMM *phmm, int N, int M, int seed);
void CopyHMM(HMM *phmm1, HMM *phmm2);
void FreeHMM(HMM *phmm);
void PrintDelta(FILE *fp, double **delta, HMM *phmm, int T);
void PrintTraceBack(FILE *fp, int **delta, HMM *phmm, int T);
void ReadSequence(FILE *fp, int *pT, int **pO);
void ReadSequences(FILE *fp, int *pN, int **pT, int ***pO);
void PrintSequences(int *pN, int *size, int ***O);
void PrintSequence(FILE *fp, int T, int *O);
void PrintPath(FILE *fp, int T_new,GList *list);
void PrintLikelihoods(FILE *fp, int T_new, GList *likelihoods);
void GenSequenceArray(HMM *phmm, int seed, int T, int *O, int *q);
int GenInitalState(HMM *phmm);
int GenNextState(HMM *phmm, int q_t);
int GenSymbol(HMM *phmm, int q_t);
void PrintAllOuts(FILE *fp, int T_new, GList *list,double ** delta, int N, int T);


void Forward(HMM *phmm, int T, int *O, double **alpha, double *pprob);
void ForwardWithScale(HMM *phmm, int T, int *O, double **alpha,
        double *scale, double *pprob);
void Backward(HMM *phmm, int T, int *O, double **beta, double *pprob);
void BackwardWithScale(HMM *phmm, int T, int *O, double **beta,
        double *scale, double *pprob);
void BaumWelch(HMM *phmm,HMM *phmm_orig,FILE *fp, int *niter,
	double *plogprobinit, double *plogprobfinal,char *log_out,char *hmm_fp);

double *** AllocXi(int T, int N);
void FreeXi(double *** xi, int T, int N);
void ComputeGamma(HMM *phmm, int T, double **alpha, double **beta,
        double **gamma);
void ComputeXi(HMM* phmm, int T, int *O, double **alpha, double **beta,
        double ***xi);
void Viterbi(HMM *phmm, int T, int *O, double **delta, int **psi,int **traceback_dir,
        GList **list, int *T_new,int *q, double *pprob);
void ViterbiLog(HMM *phmm, int T, int *O, double **delta, int **psi,int **traceback_dir,
        GList **list,GList **likelihoods_ints, GList **likelihoods_dec, int *T_new, int *q, double *pprob,int repeat_len,int prefix_idx);

/* random number generator related functions*/

int hmmgetseed(void);
void hmmsetseed(int seed);
double hmmgetrand(void);

#define MAX(x,y)        ((x) > (y) ? (x) : (y))
#define MIN(x,y)        ((x) < (y) ? (x) : (y))
