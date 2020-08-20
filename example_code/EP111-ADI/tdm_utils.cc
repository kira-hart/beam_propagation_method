#include<tdm_utils.h>


void TDM_Apply(const comple *a, const comple *b, const comple *c, const comple *src, comple *tar, int stride, int n) {
  /**
   * colculates tar = M*src where M is a tri-diagonal matrix specified as follows:
   * n - number of equations
   * a - suB-diagonal (means it is the diagonal below the main diagonal) -- indexed from 1..n-1
   * b - the main diagonal
   * c - suP-diagonal (means it is the diagonal above the main diagonal) -- indexed from 0..n-2
   **/

  tar[0]   = b[0]*src[0] + c[0]*src[ stride ];
  for(int i=1;i<n-1;i++) tar[stride*i] = a[i]*src[stride*(i-1)] + b[i]*src[stride*i] + c[i]*src[stride*(i+1)];
  tar[stride*(n-1)] = a[n-1]*src[stride*(n-2)] + b[n-1]*src[stride*(n-1)];
}

void TDM_Solve(const comple *a, comple *b, const comple *c, comple *src, comple *tar, int stride, int n)
{
        /**
         * n - number of equations
         * a - suB-diagonal (means it is the diagonal below the main diagonal) -- indexed from 1..n-1
         * b - the main diagonal
         * c - suP-diagonal (means it is the diagonal above the main diagonal) -- indexed from 0..n-2
         * v - right hand side
         * x - solution
         * NOTE: array b will be DESTROYED!!!
         */
        for (int i = 1; i < n; i++)
        {
                comple m = a[i]/b[i-1];
                b[i] = b[i] - m*c[i-1];
                src[stride*i] = src[stride*i] - m*src[stride*(i-1)];
        }
 
        tar[stride*(n-1)] = src[stride*(n-1)]/b[n-1];
	
        for (int i = n - 2; i >= 0; i--)
	  tar[stride*i]=(src[stride*i]-c[i]*tar[stride*(i+1)])/b[i];
}

