#include <fftw3.h>

#define ODD 0
#define EVEN 1
#define PARITY(a) (a)%2 ? ODD : EVEN
#define FFTLog_SWAP(a,b) do { typeof(a) temp = a; a = b; b = temp; } while (0)

typedef struct FFTLog_complex {
  double re;
  double im;
  double amp;
  double arg;
} FFTLog_complex;

typedef struct {
  int N;
  fftw_plan p_forward;
  fftw_plan p_backward;
  fftw_complex *an;
  fftw_complex *cm;
  fftw_complex *um;
  double min;
  double max;
  double q;
  double mu;
  double kr;
} FFTLog_config;

void FFTLogPT(FFTLog_config *fc, double (*func)(double));
FFTLog_config *FFTLogPT_init(int N, double min, double max, double nu);
void FFTLogPT_free(FFTLog_config *fc);
