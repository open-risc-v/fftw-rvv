/*
 * RVV kernel for FFTW scalar codelet n1_32 (1024 = 32*32, the only codelet
 * scheduled by the round-1 case's plan).
 *
 * The 32 batches are vectorized across lanes: inputs are gathered into SoA
 * form (one vector per input position), the generated codelet body is
 * evaluated lane-wise, outputs scattered back.  The statement body is the
 * verbatim non-FMA branch of n1_32.c with ri[ii][WS(is,k)] -> RIk/IIk and
 * ro[io][WS(os,k)] -> ROk/IOk, so per element the floating-point operation
 * sequence is unchanged.
 */
#ifndef FFTW_N1_32_RVV_H
#define FFTW_N1_32_RVV_H

typedef double fftw_v8df __attribute__((vector_size(64)));
typedef fftw_v8df v8df;

#define FFTW_N1_32_LANES 8

#define RI0 RI[0]
#define RI1 RI[1]
#define RI2 RI[2]
#define RI3 RI[3]
#define RI4 RI[4]
#define RI5 RI[5]
#define RI6 RI[6]
#define RI7 RI[7]
#define RI8 RI[8]
#define RI9 RI[9]
#define RI10 RI[10]
#define RI11 RI[11]
#define RI12 RI[12]
#define RI13 RI[13]
#define RI14 RI[14]
#define RI15 RI[15]
#define RI16 RI[16]
#define RI17 RI[17]
#define RI18 RI[18]
#define RI19 RI[19]
#define RI20 RI[20]
#define RI21 RI[21]
#define RI22 RI[22]
#define RI23 RI[23]
#define RI24 RI[24]
#define RI25 RI[25]
#define RI26 RI[26]
#define RI27 RI[27]
#define RI28 RI[28]
#define RI29 RI[29]
#define RI30 RI[30]
#define RI31 RI[31]
#define II0 II[0]
#define II1 II[1]
#define II2 II[2]
#define II3 II[3]
#define II4 II[4]
#define II5 II[5]
#define II6 II[6]
#define II7 II[7]
#define II8 II[8]
#define II9 II[9]
#define II10 II[10]
#define II11 II[11]
#define II12 II[12]
#define II13 II[13]
#define II14 II[14]
#define II15 II[15]
#define II16 II[16]
#define II17 II[17]
#define II18 II[18]
#define II19 II[19]
#define II20 II[20]
#define II21 II[21]
#define II22 II[22]
#define II23 II[23]
#define II24 II[24]
#define II25 II[25]
#define II26 II[26]
#define II27 II[27]
#define II28 II[28]
#define II29 II[29]
#define II30 II[30]
#define II31 II[31]
#define RO0 RO[0]
#define RO1 RO[1]
#define RO2 RO[2]
#define RO3 RO[3]
#define RO4 RO[4]
#define RO5 RO[5]
#define RO6 RO[6]
#define RO7 RO[7]
#define RO8 RO[8]
#define RO9 RO[9]
#define RO10 RO[10]
#define RO11 RO[11]
#define RO12 RO[12]
#define RO13 RO[13]
#define RO14 RO[14]
#define RO15 RO[15]
#define RO16 RO[16]
#define RO17 RO[17]
#define RO18 RO[18]
#define RO19 RO[19]
#define RO20 RO[20]
#define RO21 RO[21]
#define RO22 RO[22]
#define RO23 RO[23]
#define RO24 RO[24]
#define RO25 RO[25]
#define RO26 RO[26]
#define RO27 RO[27]
#define RO28 RO[28]
#define RO29 RO[29]
#define RO30 RO[30]
#define RO31 RO[31]
#define IO0 IO[0]
#define IO1 IO[1]
#define IO2 IO[2]
#define IO3 IO[3]
#define IO4 IO[4]
#define IO5 IO[5]
#define IO6 IO[6]
#define IO7 IO[7]
#define IO8 IO[8]
#define IO9 IO[9]
#define IO10 IO[10]
#define IO11 IO[11]
#define IO12 IO[12]
#define IO13 IO[13]
#define IO14 IO[14]
#define IO15 IO[15]
#define IO16 IO[16]
#define IO17 IO[17]
#define IO18 IO[18]
#define IO19 IO[19]
#define IO20 IO[20]
#define IO21 IO[21]
#define IO22 IO[22]
#define IO23 IO[23]
#define IO24 IO[24]
#define IO25 IO[25]
#define IO26 IO[26]
#define IO27 IO[27]
#define IO28 IO[28]
#define IO29 IO[29]
#define IO30 IO[30]
#define IO31 IO[31]

static void n1_32_rvv(const R *ri, const R *ii, R *ro, R *io,
                      stride is, stride os, INT v, INT ivs, INT ovs)
{
     INT base, k, l;
     v8df RI[32], II[32], RO[32], IO[32];

     for (base = 0; base < v; base += FFTW_N1_32_LANES) {
          INT nl = (v - base < FFTW_N1_32_LANES) ? (v - base) : FFTW_N1_32_LANES;
          for (k = 0; k < 32; k++) {
               for (l = 0; l < FFTW_N1_32_LANES; l++) {
                    if (l < nl) {
                         RI[k][l] = ri[(base + l) * ivs + is * k];
                         II[k][l] = ii[(base + l) * ivs + is * k];
                    } else {
                         RI[k][l] = 0.0;
                         II[k][l] = 0.0;
                    }
                    RO[k][l] = 0.0;
                    IO[k][l] = 0.0;
               }
          }
     DK(KP831469612, +0.831469612302545237078788377617905756738560812);
     DK(KP555570233, +0.555570233019602224742830813948532874374937191);
     DK(KP195090322, +0.195090322016128267848284868477022240927691618);
     DK(KP980785280, +0.980785280403230449126182236134239036973933731);
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
	       v8df T7, T4r, T4Z, T18, T1z, T3t, T3T, T2T, Te, T1f, T50, T4s, T2W, T3u, T1G;
	       v8df T3U, Tm, T1n, T1O, T2Z, T3y, T3X, T4w, T53, Tt, T1u, T1V, T2Y, T3B, T3W;
	       v8df T4z, T52, T2t, T3L, T3O, T2K, TR, TY, T5F, T5G, T5H, T5I, T4R, T5j, T2E;
	       v8df T3P, T4W, T5k, T2N, T3M, T22, T3E, T3H, T2j, TC, TJ, T5A, T5B, T5C, T5D;
	       v8df T4G, T5g, T2d, T3F, T4L, T5h, T2m, T3I;
	       {
		    v8df T3, T1x, T14, T2S, T6, T2R, T17, T1y;
		    {
			 v8df T1, T2, T12, T13;
			 T1 = RI0 ;
			 T2 = RI16 ;
			 T3 = T1 + T2;
			 T1x = T1 - T2;
			 T12 = II0 ;
			 T13 = II16 ;
			 T14 = T12 + T13;
			 T2S = T12 - T13;
		    }
		    {
			 v8df T4, T5, T15, T16;
			 T4 = RI8 ;
			 T5 = RI24 ;
			 T6 = T4 + T5;
			 T2R = T4 - T5;
			 T15 = II8 ;
			 T16 = II24 ;
			 T17 = T15 + T16;
			 T1y = T15 - T16;
		    }
		    T7 = T3 + T6;
		    T4r = T3 - T6;
		    T4Z = T14 - T17;
		    T18 = T14 + T17;
		    T1z = T1x - T1y;
		    T3t = T1x + T1y;
		    T3T = T2S - T2R;
		    T2T = T2R + T2S;
	       }
	       {
		    v8df Ta, T1B, T1b, T1A, Td, T1D, T1e, T1E;
		    {
			 v8df T8, T9, T19, T1a;
			 T8 = RI4 ;
			 T9 = RI20 ;
			 Ta = T8 + T9;
			 T1B = T8 - T9;
			 T19 = II4 ;
			 T1a = II20 ;
			 T1b = T19 + T1a;
			 T1A = T19 - T1a;
		    }
		    {
			 v8df Tb, Tc, T1c, T1d;
			 Tb = RI28 ;
			 Tc = RI12 ;
			 Td = Tb + Tc;
			 T1D = Tb - Tc;
			 T1c = II28 ;
			 T1d = II12 ;
			 T1e = T1c + T1d;
			 T1E = T1c - T1d;
		    }
		    Te = Ta + Td;
		    T1f = T1b + T1e;
		    T50 = Td - Ta;
		    T4s = T1b - T1e;
		    {
			 v8df T2U, T2V, T1C, T1F;
			 T2U = T1D - T1E;
			 T2V = T1B + T1A;
			 T2W = KP707106781 * (T2U - T2V);
			 T3u = KP707106781 * (T2V + T2U);
			 T1C = T1A - T1B;
			 T1F = T1D + T1E;
			 T1G = KP707106781 * (T1C - T1F);
			 T3U = KP707106781 * (T1C + T1F);
		    }
	       }
	       {
		    v8df Ti, T1L, T1j, T1J, Tl, T1I, T1m, T1M, T1K, T1N;
		    {
			 v8df Tg, Th, T1h, T1i;
			 Tg = RI2 ;
			 Th = RI18 ;
			 Ti = Tg + Th;
			 T1L = Tg - Th;
			 T1h = II2 ;
			 T1i = II18 ;
			 T1j = T1h + T1i;
			 T1J = T1h - T1i;
		    }
		    {
			 v8df Tj, Tk, T1k, T1l;
			 Tj = RI10 ;
			 Tk = RI26 ;
			 Tl = Tj + Tk;
			 T1I = Tj - Tk;
			 T1k = II10 ;
			 T1l = II26 ;
			 T1m = T1k + T1l;
			 T1M = T1k - T1l;
		    }
		    Tm = Ti + Tl;
		    T1n = T1j + T1m;
		    T1K = T1I + T1J;
		    T1N = T1L - T1M;
		    T1O = FNMS(KP923879532, T1N, KP382683432 * T1K);
		    T2Z = FMA(KP923879532, T1K, KP382683432 * T1N);
		    {
			 v8df T3w, T3x, T4u, T4v;
			 T3w = T1J - T1I;
			 T3x = T1L + T1M;
			 T3y = FNMS(KP382683432, T3x, KP923879532 * T3w);
			 T3X = FMA(KP382683432, T3w, KP923879532 * T3x);
			 T4u = T1j - T1m;
			 T4v = Ti - Tl;
			 T4w = T4u - T4v;
			 T53 = T4v + T4u;
		    }
	       }
	       {
		    v8df Tp, T1S, T1q, T1Q, Ts, T1P, T1t, T1T, T1R, T1U;
		    {
			 v8df Tn, To, T1o, T1p;
			 Tn = RI30 ;
			 To = RI14 ;
			 Tp = Tn + To;
			 T1S = Tn - To;
			 T1o = II30 ;
			 T1p = II14 ;
			 T1q = T1o + T1p;
			 T1Q = T1o - T1p;
		    }
		    {
			 v8df Tq, Tr, T1r, T1s;
			 Tq = RI6 ;
			 Tr = RI22 ;
			 Ts = Tq + Tr;
			 T1P = Tq - Tr;
			 T1r = II6 ;
			 T1s = II22 ;
			 T1t = T1r + T1s;
			 T1T = T1r - T1s;
		    }
		    Tt = Tp + Ts;
		    T1u = T1q + T1t;
		    T1R = T1P + T1Q;
		    T1U = T1S - T1T;
		    T1V = FMA(KP382683432, T1R, KP923879532 * T1U);
		    T2Y = FNMS(KP923879532, T1R, KP382683432 * T1U);
		    {
			 v8df T3z, T3A, T4x, T4y;
			 T3z = T1Q - T1P;
			 T3A = T1S + T1T;
			 T3B = FMA(KP923879532, T3z, KP382683432 * T3A);
			 T3W = FNMS(KP382683432, T3z, KP923879532 * T3A);
			 T4x = Tp - Ts;
			 T4y = T1q - T1t;
			 T4z = T4x + T4y;
			 T52 = T4x - T4y;
		    }
	       }
	       {
		    v8df TN, T2p, T2J, T4S, TQ, T2G, T2s, T4T, TU, T2x, T2w, T4O, TX, T2z, T2C;
		    v8df T4P;
		    {
			 v8df TL, TM, T2H, T2I;
			 TL = RI31 ;
			 TM = RI15 ;
			 TN = TL + TM;
			 T2p = TL - TM;
			 T2H = II31 ;
			 T2I = II15 ;
			 T2J = T2H - T2I;
			 T4S = T2H + T2I;
		    }
		    {
			 v8df TO, TP, T2q, T2r;
			 TO = RI7 ;
			 TP = RI23 ;
			 TQ = TO + TP;
			 T2G = TO - TP;
			 T2q = II7 ;
			 T2r = II23 ;
			 T2s = T2q - T2r;
			 T4T = T2q + T2r;

		    }
		    {
			 v8df TS, TT, T2u, T2v;
			 TS = RI3 ;
			 TT = RI19 ;
			 TU = TS + TT;
			 T2x = TS - TT;
			 T2u = II3 ;
			 T2v = II19 ;
			 T2w = T2u - T2v;
			 T4O = T2u + T2v;
		    }
		    {
			 v8df TV, TW, T2A, T2B;
			 TV = RI27 ;
			 TW = RI11 ;
			 TX = TV + TW;
			 T2z = TV - TW;
			 T2A = II27 ;
			 T2B = II11 ;
			 T2C = T2A - T2B;
			 T4P = T2A + T2B;
		    }
		    T2t = T2p - T2s;
		    T3L = T2p + T2s;
		    T3O = T2J - T2G;
		    T2K = T2G + T2J;
		    TR = TN + TQ;
		    TY = TU + TX;
		    T5F = TR - TY;
		    {
			 v8df T4N, T4Q, T2y, T2D;
			 T5G = T4S + T4T;
			 T5H = T4O + T4P;
			 T5I = T5G - T5H;
			 T4N = TN - TQ;
			 T4Q = T4O - T4P;
			 T4R = T4N - T4Q;
			 T5j = T4N + T4Q;
			 T2y = T2w - T2x;
			 T2D = T2z + T2C;
			 T2E = KP707106781 * (T2y - T2D);
			 T3P = KP707106781 * (T2y + T2D);
			 {
			      v8df T4U, T4V, T2L, T2M;
			      T4U = T4S - T4T;
			      T4V = TX - TU;
			      T4W = T4U - T4V;
			      T5k = T4V + T4U;
			      T2L = T2z - T2C;
			      T2M = T2x + T2w;
			      T2N = KP707106781 * (T2L - T2M);
			      T3M = KP707106781 * (T2M + T2L);
			 }
		    }
	       }
	       {
		    v8df Ty, T2f, T21, T4C, TB, T1Y, T2i, T4D, TF, T28, T2b, T4I, TI, T23, T26;
		    v8df T4J;
		    {
			 v8df Tw, Tx, T1Z, T20;
			 Tw = RI1 ;
			 Tx = RI17 ;
			 Ty = Tw + Tx;
			 T2f = Tw - Tx;
			 T1Z = II1 ;
			 T20 = II17 ;
			 T21 = T1Z - T20;
			 T4C = T1Z + T20;
		    }
		    {
			 v8df Tz, TA, T2g, T2h;
			 Tz = RI9 ;
			 TA = RI25 ;
			 TB = Tz + TA;
			 T1Y = Tz - TA;
			 T2g = II9 ;
			 T2h = II25 ;
			 T2i = T2g - T2h;
			 T4D = T2g + T2h;
		    }
		    {
			 v8df TD, TE, T29, T2a;
			 TD = RI5 ;
			 TE = RI21 ;
			 TF = TD + TE;
			 T28 = TD - TE;
			 T29 = II5 ;
			 T2a = II21 ;
			 T2b = T29 - T2a;
			 T4I = T29 + T2a;
		    }
		    {
			 v8df TG, TH, T24, T25;
			 TG = RI29 ;
			 TH = RI13 ;
			 TI = TG + TH;
			 T23 = TG - TH;
			 T24 = II29 ;
			 T25 = II13 ;
			 T26 = T24 - T25;
			 T4J = T24 + T25;
		    }
		    T22 = T1Y + T21;
		    T3E = T2f + T2i;
		    T3H = T21 - T1Y;
		    T2j = T2f - T2i;
		    TC = Ty + TB;
		    TJ = TF + TI;
		    T5A = TC - TJ;
		    {
			 v8df T4E, T4F, T27, T2c;
			 T5B = T4C + T4D;
			 T5C = T4I + T4J;
			 T5D = T5B - T5C;
			 T4E = T4C - T4D;
			 T4F = TI - TF;
			 T4G = T4E - T4F;
			 T5g = T4F + T4E;
			 T27 = T23 - T26;
			 T2c = T28 + T2b;
			 T2d = KP707106781 * (T27 - T2c);
			 T3F = KP707106781 * (T2c + T27);
			 {
			      v8df T4H, T4K, T2k, T2l;
			      T4H = Ty - TB;
			      T4K = T4I - T4J;
			      T4L = T4H - T4K;
			      T5h = T4H + T4K;
			      T2k = T2b - T28;
			      T2l = T23 + T26;
			      T2m = KP707106781 * (T2k - T2l);
			      T3I = KP707106781 * (T2k + T2l);
			 }
		    }
	       }
	       {
		    v8df T4B, T57, T5a, T5c, T4Y, T56, T55, T5b;
		    {
			 v8df T4t, T4A, T58, T59;
			 T4t = T4r - T4s;
			 T4A = KP707106781 * (T4w - T4z);
			 T4B = T4t + T4A;
			 T57 = T4t - T4A;
			 T58 = FNMS(KP923879532, T4L, KP382683432 * T4G);
			 T59 = FMA(KP382683432, T4W, KP923879532 * T4R);
			 T5a = T58 - T59;
			 T5c = T58 + T59;
		    }
		    {
			 v8df T4M, T4X, T51, T54;
			 T4M = FMA(KP923879532, T4G, KP382683432 * T4L);
			 T4X = FNMS(KP923879532, T4W, KP382683432 * T4R);
			 T4Y = T4M + T4X;
			 T56 = T4X - T4M;
			 T51 = T4Z - T50;
			 T54 = KP707106781 * (T52 - T53);
			 T55 = T51 - T54;
			 T5b = T51 + T54;
		    }
		    RO22 = T4B - T4Y;
		    IO22 = T5b - T5c;
		    RO6 = T4B + T4Y;
		    IO6 = T5b + T5c;
		    IO30 = T55 - T56;
		    RO30 = T57 - T5a;
		    IO14 = T55 + T56;
		    RO14 = T57 + T5a;
	       }
	       {
		    v8df T5f, T5r, T5u, T5w, T5m, T5q, T5p, T5v;
		    {
			 v8df T5d, T5e, T5s, T5t;
			 T5d = T4r + T4s;
			 T5e = KP707106781 * (T53 + T52);
			 T5f = T5d + T5e;
			 T5r = T5d - T5e;
			 T5s = FNMS(KP382683432, T5h, KP923879532 * T5g);
			 T5t = FMA(KP923879532, T5k, KP382683432 * T5j);
			 T5u = T5s - T5t;
			 T5w = T5s + T5t;
		    }
		    {
			 v8df T5i, T5l, T5n, T5o;
			 T5i = FMA(KP382683432, T5g, KP923879532 * T5h);
			 T5l = FNMS(KP382683432, T5k, KP923879532 * T5j);
			 T5m = T5i + T5l;
			 T5q = T5l - T5i;
			 T5n = T50 + T4Z;
			 T5o = KP707106781 * (T4w + T4z);
			 T5p = T5n - T5o;
			 T5v = T5n + T5o;
		    }
		    RO18 = T5f - T5m;
		    IO18 = T5v - T5w;
		    RO2 = T5f + T5m;
		    IO2 = T5v + T5w;
		    IO26 = T5p - T5q;
		    RO26 = T5r - T5u;
		    IO10 = T5p + T5q;
		    RO10 = T5r + T5u;
	       }
	       {
		    v8df T5z, T5P, T5S, T5U, T5K, T5O, T5N, T5T;
		    {
			 v8df T5x, T5y, T5Q, T5R;
			 T5x = T7 - Te;
			 T5y = T1n - T1u;
			 T5z = T5x + T5y;
			 T5P = T5x - T5y;
			 T5Q = T5D - T5A;
			 T5R = T5F + T5I;
			 T5S = KP707106781 * (T5Q - T5R);
			 T5U = KP707106781 * (T5Q + T5R);
		    }
		    {
			 v8df T5E, T5J, T5L, T5M;
			 T5E = T5A + T5D;
			 T5J = T5F - T5I;
			 T5K = KP707106781 * (T5E + T5J);
			 T5O = KP707106781 * (T5J - T5E);
			 T5L = T18 - T1f;
			 T5M = Tt - Tm;
			 T5N = T5L - T5M;
			 T5T = T5M + T5L;
		    }
		    RO20 = T5z - T5K;
		    IO20 = T5T - T5U;
		    RO4 = T5z + T5K;
		    IO4 = T5T + T5U;
		    IO28 = T5N - T5O;
		    RO28 = T5P - T5S;
		    IO12 = T5N + T5O;
		    RO12 = T5P + T5S;
	       }
	       {
		    v8df Tv, T5V, T5Y, T60, T10, T11, T1w, T5Z;
		    {
			 v8df Tf, Tu, T5W, T5X;
			 Tf = T7 + Te;
			 Tu = Tm + Tt;
			 Tv = Tf + Tu;
			 T5V = Tf - Tu;
			 T5W = T5B + T5C;
			 T5X = T5G + T5H;
			 T5Y = T5W - T5X;
			 T60 = T5W + T5X;
		    }
		    {
			 v8df TK, TZ, T1g, T1v;
			 TK = TC + TJ;
			 TZ = TR + TY;
			 T10 = TK + TZ;
			 T11 = TZ - TK;
			 T1g = T18 + T1f;
			 T1v = T1n + T1u;
			 T1w = T1g - T1v;
			 T5Z = T1g + T1v;
		    }
		    RO16 = Tv - T10;
		    IO16 = T5Z - T60;
		    RO0  = Tv + T10;
		    IO0  = T5Z + T60;
		    IO8 = T11 + T1w;
		    RO8 = T5V + T5Y;
		    IO24 = T1w - T11;
		    RO24 = T5V - T5Y;
	       }
	       {
		    v8df T1X, T33, T31, T37, T2o, T34, T2P, T35;
		    {
			 v8df T1H, T1W, T2X, T30;
			 T1H = T1z - T1G;
			 T1W = T1O - T1V;
			 T1X = T1H + T1W;
			 T33 = T1H - T1W;
			 T2X = T2T - T2W;
			 T30 = T2Y - T2Z;
			 T31 = T2X - T30;
			 T37 = T2X + T30;
		    }
		    {
			 v8df T2e, T2n, T2F, T2O;
			 T2e = T22 - T2d;
			 T2n = T2j - T2m;
			 T2o = FMA(KP980785280, T2e, KP195090322 * T2n);
			 T34 = FNMS(KP980785280, T2n, KP195090322 * T2e);
			 T2F = T2t - T2E;
			 T2O = T2K - T2N;
			 T2P = FNMS(KP980785280, T2O, KP195090322 * T2F);
			 T35 = FMA(KP195090322, T2O, KP980785280 * T2F);
		    }
		    {
			 v8df T2Q, T38, T32, T36;
			 T2Q = T2o + T2P;
			 RO23 = T1X - T2Q;
			 RO7 = T1X + T2Q;
			 T38 = T34 + T35;
			 IO23 = T37 - T38;
			 IO7 = T37 + T38;
			 T32 = T2P - T2o;
			 IO31 = T31 - T32;
			 IO15 = T31 + T32;
			 T36 = T34 - T35;
			 RO31 = T33 - T36;
			 RO15 = T33 + T36;
		    }
	       }
	       {
		    v8df T3D, T41, T3Z, T45, T3K, T42, T3R, T43;
		    {
			 v8df T3v, T3C, T3V, T3Y;
			 T3v = T3t - T3u;
			 T3C = T3y - T3B;
			 T3D = T3v + T3C;
			 T41 = T3v - T3C;
			 T3V = T3T - T3U;
			 T3Y = T3W - T3X;
			 T3Z = T3V - T3Y;
			 T45 = T3V + T3Y;
		    }
		    {
			 v8df T3G, T3J, T3N, T3Q;
			 T3G = T3E - T3F;
			 T3J = T3H - T3I;
			 T3K = FMA(KP555570233, T3G, KP831469612 * T3J);
			 T42 = FNMS(KP831469612, T3G, KP555570233 * T3J);
			 T3N = T3L - T3M;
			 T3Q = T3O - T3P;
			 T3R = FNMS(KP831469612, T3Q, KP555570233 * T3N);
			 T43 = FMA(KP831469612, T3N, KP555570233 * T3Q);
		    }
		    {
			 v8df T3S, T46, T40, T44;

			 T3S = T3K + T3R;
			 RO21 = T3D - T3S;
			 RO5 = T3D + T3S;
			 T46 = T42 + T43;
			 IO21 = T45 - T46;
			 IO5 = T45 + T46;
			 T40 = T3R - T3K;
			 IO29 = T3Z - T40;
			 IO13 = T3Z + T40;
			 T44 = T42 - T43;
			 RO29 = T41 - T44;
			 RO13 = T41 + T44;
		    }
	       }
	       {
		    v8df T49, T4l, T4j, T4p, T4c, T4m, T4f, T4n;
		    {
			 v8df T47, T48, T4h, T4i;
			 T47 = T3t + T3u;
			 T48 = T3X + T3W;
			 T49 = T47 + T48;
			 T4l = T47 - T48;
			 T4h = T3T + T3U;
			 T4i = T3y + T3B;
			 T4j = T4h - T4i;
			 T4p = T4h + T4i;
		    }
		    {
			 v8df T4a, T4b, T4d, T4e;
			 T4a = T3E + T3F;
			 T4b = T3H + T3I;
			 T4c = FMA(KP980785280, T4a, KP195090322 * T4b);
			 T4m = FNMS(KP195090322, T4a, KP980785280 * T4b);
			 T4d = T3L + T3M;
			 T4e = T3O + T3P;
			 T4f = FNMS(KP195090322, T4e, KP980785280 * T4d);
			 T4n = FMA(KP195090322, T4d, KP980785280 * T4e);
		    }
		    {
			 v8df T4g, T4q, T4k, T4o;
			 T4g = T4c + T4f;
			 RO17 = T49 - T4g;
			 RO1 = T49 + T4g;
			 T4q = T4m + T4n;
			 IO17 = T4p - T4q;
			 IO1 = T4p + T4q;
			 T4k = T4f - T4c;
			 IO25 = T4j - T4k;
			 IO9 = T4j + T4k;
			 T4o = T4m - T4n;
			 RO25 = T4l - T4o;
			 RO9 = T4l + T4o;
		    }
	       }
	       {
		    v8df T3b, T3n, T3l, T3r, T3e, T3o, T3h, T3p;
		    {
			 v8df T39, T3a, T3j, T3k;
			 T39 = T1z + T1G;
			 T3a = T2Z + T2Y;
			 T3b = T39 + T3a;
			 T3n = T39 - T3a;
			 T3j = T2T + T2W;
			 T3k = T1O + T1V;
			 T3l = T3j - T3k;
			 T3r = T3j + T3k;
		    }
		    {
			 v8df T3c, T3d, T3f, T3g;
			 T3c = T22 + T2d;
			 T3d = T2j + T2m;
			 T3e = FMA(KP555570233, T3c, KP831469612 * T3d);
			 T3o = FNMS(KP555570233, T3d, KP831469612 * T3c);
			 T3f = T2t + T2E;
			 T3g = T2K + T2N;
			 T3h = FNMS(KP555570233, T3g, KP831469612 * T3f);
			 T3p = FMA(KP831469612, T3g, KP555570233 * T3f);
		    }
		    {
			 v8df T3i, T3s, T3m, T3q;
			 T3i = T3e + T3h;
			 RO19 = T3b - T3i;
			 RO3 = T3b + T3i;
			 T3s = T3o + T3p;
			 IO19 = T3r - T3s;
			 IO3 = T3r + T3s;
			 T3m = T3h - T3e;
			 IO27 = T3l - T3m;
			 IO11 = T3l + T3m;
			 T3q = T3o - T3p;
			 RO27 = T3n - T3q;
			 RO11 = T3n + T3q;

          }
          }
          for (k = 0; k < 32; k++) {
               for (l = 0; l < nl; l++) {
                    ro[(base + l) * ovs + os * k] = RO[k][l];
                    io[(base + l) * ovs + os * k] = IO[k][l];
               }
          }
     }
}

#endif /* FFTW_N1_32_RVV_H */
