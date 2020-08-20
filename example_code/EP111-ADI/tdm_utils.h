#ifndef __tdm_utils_h__
#define __tdm_utils_h__

#include<complex>
typedef std::complex<double> comple;


void TDM_Apply(const comple *a, const comple *b, const comple *c, const comple *src, comple *tar, int stride, int n);
void TDM_Solve(const comple *a,       comple *b, const comple *c,       comple *src, comple *tar, int stride, int n);



#endif // __tdm_utils_h__
