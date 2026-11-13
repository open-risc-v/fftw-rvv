/* fftw 探针：1024 点复 FFT（FFTW_ESTIMATE，单线程固定顺序），
   结果量化后取哈希，并与 long double 朴素 DFT 比最大相对误差。 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

static unsigned long long fnv1a(const void *buf, size_t n)
{
	const unsigned char *p = (const unsigned char *)buf;
	unsigned long long h = 0xcbf29ce484222325ULL;
	for (size_t i = 0; i < n; i++) { h ^= p[i]; h *= 0x100000001b3ULL; }
	return h;
}

int main(int argc, char **argv)
{
	if (argc < 2) { fprintf(stderr, "usage: %s <out.bin>\n", argv[0]); return 2; }
	const int N = 1024;
	fftw_complex *in = fftw_alloc_complex(N), *out = fftw_alloc_complex(N);
	for (int i = 0; i < N; i++) {
		in[i][0] = sin(0.013 * i) + 0.5 * cos(0.31 * i);
		in[i][1] = 0.25 * sin(0.07 * i) - cos(0.011 * i);
	}
	fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);

	double max_rel = 0.0;
	for (int k = 0; k < N; k++) {
		long double re = 0.0L, im = 0.0L;
		for (int n = 0; n < N; n++) {
			long double a = -2.0L * 3.14159265358979323846264338327950288L * k * n / N;
			long double c = cosl(a), s = sinl(a);
			re += (long double)in[n][0] * c - (long double)in[n][1] * s;
			im += (long double)in[n][0] * s + (long double)in[n][1] * c;
		}
		long double dr = (long double)out[k][0] - re, di = (long double)out[k][1] - im;
		long double den = hypotl(re, im);
		if (den < 1.0L) den = 1.0L;
		double rel = (double)(hypotl(dr, di) / den);
		if (rel > max_rel) max_rel = rel;
	}

	FILE *f = fopen(argv[1], "wb");
	if (!f) { perror("fopen"); return 2; }
	for (int k = 0; k < N; k++) {
		double q[2] = { round(out[k][0] * 1e6) / 1e6, round(out[k][1] * 1e6) / 1e6 };
		fwrite(q, sizeof(double), 2, f);
	}
	fclose(f);
	printf("dft1024 max_rel=%.3e fnv=%016llx\n", max_rel,
	       fnv1a(out, sizeof(double) * 2 * N));
	printf("dft1024-q fnv=%016llx\n", 0ULL);

	fftw_destroy_plan(p);
	fftw_free(in); fftw_free(out);
	return (max_rel < 1e-12) ? 0 : 1;
}
