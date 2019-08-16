#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

#include "typedefs.h"
#include "typedefs1.h"
#include "rnnoise.h"


/*----------------------------------------------------------------------------*/
/*						          瀹忓畾涔?                                     */
/*----------------------------------------------------------------------------*/
#define Word16 short
#define Word32 long
#define Float32 float
#define PI           3.1415926535897932
#define FS           16000           //sampling rate
#define FRM_LEN      160             //length of a frame(16k)
#define NUM          6               //number of microphone
#define D            FRM_LEN         //half of the number of subband
#define ORD2         (6*2*D)         //order of prototype low pass filter
#define T60          200             //reverberation time
#define ORD3         (T60*FS/1000/D) //order of subband filter
#define ORD4         6         
#define ORD6         10

#define DIAMETER     0.08f           //diameter of microphone array(cm)
#define SPEED        343             //velocity of sound
#define AZ_STEP      18              //step size in search
#define AZ_NUM       (20*AZ_STEP)    //number of azimuth
#define INIT_LEN     64              //Initialization number of noise spectrum
#define MED_NUM      7               //order of median filter 
#define AA           0.95f 
#define BB           11    

#define FRM_LEN2     480             //length of a frame(48k)
#define FRM_LEN3     160
#define ORD5         72              //order of anti-alias filter
#define DELAY_LEN    (ORD5-1)    

#define FRM_LEN1    480
#define ORD48       73

#define LAMBDA      0.9

///VAD
#define TRUE			1
#define FALSE			0
//#define	NUM_CHAN		16
#define NUM_CHAN        20
#define	LO_CHAN			0
#define	MID_CHAN		5
#define	HI_CHAN			15

#define	UPDATE_THLD		35
#define	METRIC_THLD		45
#define	INDEX_THLD		12
#define	SETBACK_THLD	12
#define	SNR_THLD		6
#define	INDEX_CNT_THLD	5
#define	UPDATE_CNT_THLD	50
#define	NORM_ENRG		1.0f	     // use (32768.0 * 32768.0) for fractional 
#define	NOISE_FLOOR		1.0f
#define	MIN_CHAN_ENRG	0.0625f
#define	INE			    16.0f
#define	MIN_GAIN		-13.0f
#define	HYSTER_CNT_THLD	6		     // forced update constants... 
#define	HIGH_TCE_DB		50.0f	     // 50 - 10log10(NORM_ENRG) 
#define	LOW_TCE_DB		30.0f	     // 30 - 10log10(NORM_ENRG) 
#define	TCE_RANGE		(HIGH_TCE_DB-LOW_TCE_DB)
#define	DEV_THLD		28.0f

const Word16 lpf3[ORD5] = { -63, -112, -119, -22, 179, 399, 497, 373, 63, -248, -336, -109, 290, 548, 410, -93, -612, -709, -215, 572, 1039, 714, -300, -1287, -1387, -294, 1352, 2263, 1415, -1031, -3492, -3746, -339, 6262, 13598, 18408, 18408, 13598, 6262, -339, -3746, -3492, -1031, 1415, 2263, 1352, -294, -1387, -1287, -300, 714, 1039, 572, -215, -709, -612, -93, 410, 548, 290, -109, -336, -248, 63, 373, 497, 399, 179, -22, -119, -112, -63 };
const Word16 lpf48[ORD48] = { -2, -2, 0, 22, 18, 1, -136, -104, 0, 512,448,-32,-1408,-1280,256,3840,3328,-1024,-8192,-6144,3584,14336,10240,-4096,-22528,-16384,8192,32767,28672,-8192,-32768,-28672,16384,32767,32767,20480,22528,20480,32767,32767,16384,-28672,-32768,-8192,28672,32767,8192,-16384,-22528,-4096,10240,14336,3584,-6144,-8192,-1024,3328,3840,256,-1280,-1408,-32,448,512,0,-104,-136,1,18,22,0,-2,-2 };

typedef struct
{
	Word16 x[DELAY_LEN];
} DOWNSAMPLE_STR;

void down_sample_init(DOWNSAMPLE_STR *st)
{
	Word32 i;

	for (i = 0; i < DELAY_LEN; i++) st->x[i] = 0;
}

void down_sample_run(
	DOWNSAMPLE_STR *st,
	Word16 in_sp[FRM_LEN2],
	Word16 out_sp[FRM_LEN3]
)
{
	Word16 delay[FRM_LEN2 + DELAY_LEN];
	Word32 acc2, i, j;

	for (j = 0; j < DELAY_LEN; j++) delay[j] = st->x[j];
	for (j = 0; j < FRM_LEN2; j++) delay[DELAY_LEN + j] = in_sp[j];
	for (j = 0; j < DELAY_LEN; j++) st->x[j] = in_sp[FRM_LEN2 - DELAY_LEN + j];

	for (i = 0; i < FRM_LEN3; i++)
	{
		for (acc2 = j = 0; j < ORD5; j++) acc2 += delay[3 * i + j] * lpf3[j];
		out_sp[i] = (Word16)(acc2 >> 16);
	}
}

typedef struct
{
	Word16 x[ORD48];
} INTERP_STR;

void interp_init(INTERP_STR* st)
{
	Word32 i;
	for (i = 0; i < ORD48; i++) st->x[i] = 0;
}

void interp_run(
	INTERP_STR* st,
	Word16 in_sp[FRM_LEN],
	Word16 out_sp[FRM_LEN1]
)
{
	Word16 delay[ORD48 + FRM_LEN1];
	Word32 acc2, i, j;

	for (j = 0; j < ORD48; j++) delay[j] = st->x[j];
	for (j = 0; j < FRM_LEN; j++)
	{
		delay[ORD48 + 3 * j] = in_sp[j];
		delay[ORD48 + 3 * j + 1] = 0;
		delay[ORD48 + 3 * j + 2] = 0;
		//delay[ORD48 + 3 * j] = delay[ORD48 + 3 * j + 1] = 0;
		//delay[ORD48 + 3 * j + 2] = in_sp[j];
	}
	for (j = 0; j < ORD48; j++) st->x[j] = delay[FRM_LEN1 + j];

	for (i = 0; i < FRM_LEN1; i++)
	{
		for (acc2 = j = 0; j < ORD48; j++) acc2 += delay[1 + i + j] * lpf48[j];
		out_sp[i] = (Word16)(acc2 >> 16);
	}
}

enum SOUND_SOURCE_CATEGORY
{
	HUMAN_FACE,
	NOISE_SOURCE
};

//function of agc initialization
void AGC_init(Agc_t *st)
{
	int minLevel, maxLevel, agcMode, fs;
	AGC_config_t agcConfig;

	minLevel = 0;
	maxLevel = 255;
	agcMode = kAgcModeFixedDigital;
	fs = 16000;
	st->lastError = 0;

	AGC_Init(st, minLevel, maxLevel, agcMode, fs);

	agcConfig.compressionGaindB = 6;
	agcConfig.limiterEnable = 1;
	agcConfig.targetLevelDbfs = 3;
	AGC_set_config(st, agcConfig);
}

void INV_3(double A[3][3], int n) // three-order matrix inversion  //blockinv1
{
	long i, j;
	double tp, tmp, temp[2];
	double a[2][2];

	temp[0] = A[0][2] / A[2][2]; temp[1] = A[1][2] / A[2][2];
	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 2; j++)
		{
			a[i][j] = A[i][j] - temp[i] * A[2][j];
		}
	}
	//INV_2(a, 2);
	tp = a[0][1] / a[1][1];
	a[0][0] = 1.0 / (a[0][0] - tp * a[1][0]);
	a[0][1] = -a[0][0] * tp;
	a[1][1] = (1.0 - a[1][0] * a[0][1]) / a[1][1];
	a[1][0] = a[0][1];
	for (tmp = i = 0; i < 2; i++)
	{
		for (A[i][2] = j = 0; j < 2; j++)
		{
			A[i][2] -= a[i][j] * temp[j];//负号变+=为-=
		}
		tmp += A[2][i] * A[i][2];
		A[2][i] = A[i][2];
	}
	A[2][2] = (1 - tmp) / A[2][2];
	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 2; j++)
		{
			A[i][j] = a[i][j];
		}
	}
}
void Matrix_Mul_3(double a[3][3], double b[3][3], double c[3][3], int n) //three-order matrix multiplication,c=a*b
{
	int i, j, k;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (c[i][j] = k = 0; k < n; k++) c[i][j] += a[i][k] * b[k][j];
		}
	}
}
void INV(double A[NUM][NUM], int n)
{
	int i, j;
	double a[3][3], b[3][3], c[3][3], d[3][3], temp1[3][3], temp2[3][3], a1[3][3], b1[3][3], d1[3][3];
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			a[i][j] = A[i][j];
			b[i][j] = A[i][j + 3];
			c[i][j] = A[i + 3][j];
			d[i][j] = A[i + 3][j + 3];
		}
	}
	INV_3(d, 3);
	Matrix_Mul_3(b, d, temp1, 3);
	Matrix_Mul_3(temp1, c, a1, 3);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			a1[i][j] = a[i][j] - a1[i][j];
		}
	}
	INV_3(a1, 3);
	Matrix_Mul_3(a1, temp1, b1, 3);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			b1[i][j] = -b1[i][j];
		}
	}
	Matrix_Mul_3(c, b1, temp2, 3);
	for (i = 0; i < 3; i++)
	{
		temp2[i][i] = 1 - temp2[i][i];
		for (j = 0; j < 3; j++)
		{
			if (j != i) temp2[i][j] = -temp2[i][j];
		}
	}
	Matrix_Mul_3(d, temp2, d1, 3);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			A[i][j] = a1[i][j];
			A[i][j + 3] = b1[i][j];
			A[i + 3][j] = b1[j][i];
			A[i + 3][j + 3] = d1[i][j];
		}
	}
}

/*void INV(double A[NUM][NUM], int n) //实数矩阵求逆的函数，对矩阵A求逆，求逆的结果还放在矩阵A中  //blockinv2
{
	long i, j, k;
	double tp, temp1[2][2], temp2[2][2], tmp[2][2], temp3[2][4], temp4[2][2], temp5[2][4], temp6[4][4], temp7[4][4];

	//INV_2(a, 2);
	tp=A[4][5]/A[5][5];
	A[4][4]=1.0/(A[4][4]-tp*A[5][4]);
	A[4][5]=-A[4][4]*tp;
	A[5][5]=(1.0-A[5][4]*A[4][5])/A[5][5];
	A[5][4]=A[4][5];
	for (i=2; i<4; i++)
	{
		for (j=4; j<n; j++)
		{
			temp1[i-2][j-4] = A[i][4]*A[4][j] + A[i][5]*A[5][j];
		}
	}
	for (i=0; i<2; i++)
	{
		for (j=2; j<4; j++)
		{
			A[i+2][j] -= temp1[i][0]*A[4][j] + temp1[i][1]*A[5][j];
		}
	}
	tp=A[2][3]/A[3][3];
	A[2][2]=1.0/(A[2][2]-tp*A[3][2]);
	A[2][3]=-A[2][2]*tp;
	A[3][3]=(1.0-A[3][2]*A[2][3])/A[3][3];
	A[3][2]=A[2][3];
	for (i=2; i<4; i++)
	{
		for (j=0; j<2; j++)
		{
			A[i][j+4] = -A[i][2]*temp1[0][j] - A[i][3]*temp1[1][j];
		}
	}
	for (i=4; i<n; i++)
	{
		for (j=4; j<n; j++)
		{
			temp2[i-4][j-4] = -A[i][2]*A[2][j] - A[i][3]*A[3][j];
		}
	}
	temp2[0][0] += 1.0;
	temp2[1][1] += 1.0;
	for (i=4; i<n; i++)
	{
		for (j=0; j<2; j++)
		{
			tmp[i-4][j] = A[i][4]*temp2[0][j] + A[i][5]*temp2[1][j];
		}
	}
	for (i=0; i<2; i++)
	{
		for (j=0; j<2; j++)
		{
			A[i+4][j+4] = tmp[i][j];
		}
	}
	for (i=4; i<n; i++)
	{
		for (j=2; j<4; j++)
		{
			A[i][j] = A[j][i];
		}
	}////4*4的A22的逆√
	for (i=0; i<2; i++)
	{
		for (j=0; j<4; j++)
		{
			for (temp3[i][j]=0,k=2; k<n; k++) temp3[i][j] += A[i][k]*A[k][j+2];
		}
	}
	for (i=0; i<2; i++)
	{
		for (j=0; j<2; j++)
		{
			for (temp4[i][j]=k=0; k<4; k++) temp4[i][j] += temp3[i][k]*A[k+2][j];
		}
	}
	for (i=0; i<2; i++)
	{
		for (j=0; j<2; j++)
		{
			A[i][j] -= temp4[i][j];
		}
	}
	tp=A[0][1]/A[1][1];
	A[0][0]=1.0/(A[0][0]-tp*A[1][0]);
	A[0][1]=-A[0][0]*tp;
	A[1][1]=(1.0-A[1][0]*A[0][1])/A[1][1];
	A[1][0]=A[0][1];
	for (i=0; i<2; i++)
	{
		for (j=0; j<4; j++)
		{
			for (temp5[i][j]=k=0; k<2; k++) temp5[i][j] += A[i][k]*temp3[k][j];
		}
	}
	for (i=0; i<2; i++)
	{
		for (j=0; j<4; j++)
		{
			A[i][j+2] = -temp5[i][j];
		}
	}
	for (i=0; i<4; i++)
	{
		for (j=0; j<4; j++)
		{
			for (temp6[i][j]=k=0; k<2; k++) temp6[i][j] -= A[i+2][k]*A[k][j+2];
		}
		temp6[i][i] += 1.0;
	}
	for (i=0; i<4; i++)
	{
		for (j=0; j<4; j++)
		{
			for (temp7[i][j]=k=0; k<4; k++) temp7[i][j] += A[i+2][k+2]*temp6[k][j];
		}
	}
	for (i=0; i<4; i++)
	{
		for (j=0; j<4; j++)
		{
			A[i+2][j+2] = temp7[i][j];
		}
	}
	for (i=2; i<n; i++)
	{
		for (j=0; j<2; j++)
		{
			A[i][j] = A[j][i];
		}
	}
}*/

/*void INV(double A[NUM][NUM], int n)   //normal
{
	Word16 i, j, k, JS[NUM], IS[NUM];
	double d, temp;

	for (k=0; k<n; k++)
	{
		for (d=0, i=k; i<n; i++)
		{
			for (j=k; j<n; j++)
			{
				if (fabs(A[i][j])>d)
				{
					d = fabs(A[i][j]);
					IS[k] = i;
					JS[k] = j;
				}
			}
		}
		if (d+1.0 == 1.0) return;

		if (IS[k] != k)
		{
			for (j = 0;j<n;j++)
			{
				temp = A[k][j];
				A[k][j] = A[IS[k]][j];
				A[IS[k]][j] = temp;
			}
		}
		if (JS[k] != k)
		{
			for (i = 0;i<n;i++)
			{
				temp = A[i][k];
				A[i][k] = A[i][JS[k]];
				A[i][JS[k]] = temp;
			}
		}
		A[k][k] = 1.0 / A[k][k];
		for (j=0; j<n; j++)
		{
			if (j != k) A[k][j] = A[k][j] * A[k][k];
		}
		for (i=0; i<n; i++)
		{
			if (i != k)
			{
				for (j = 0;j<n;j++)
				{
					if (j != k) A[i][j] = A[i][j] - A[i][k] * A[k][j];
				}
			}
		}
		for (i = 0;i<n;i++)
		{
			if (i != k) A[i][k] = -A[i][k] * A[k][k];
		}
	}
	for (k = n - 1;k >= 0;k--)
	{
		for (j = 0;j<n;j++)
		{
			if (JS[k] != k)
			{
				temp = A[k][j];
				A[k][j] = A[JS[k]][j];
				A[JS[k]][j] = temp;
			}
		}
		for (i = 0;i<n;i++)
		{
			if (IS[k] != k)
			{
				temp = A[i][k];
				A[i][k] = A[i][IS[k]];
				A[i][IS[k]] = temp;
			}
		}
	}
}*/

void Matrix_Mul(double a[NUM][NUM], double b[NUM][NUM], double c[NUM][NUM], int n)
{
	int i, j, k;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (c[i][j] = k = 0; k < n; k++) c[i][j] += a[i][k] * b[k][j];
		}
	}
}

void inverse(double A[NUM][NUM], double B[NUM][NUM], double C_re[NUM][NUM], double C_im[NUM][NUM], int n)
{
	int i, j;
	double  temp_A[NUM][NUM], temp_AB[NUM][NUM], temp_BAB[NUM][NUM];

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++) temp_A[i][j] = A[i][j];
	}
	INV(A, NUM);
	Matrix_Mul(A, B, temp_AB, NUM);
	Matrix_Mul(B, temp_AB, temp_BAB, NUM);
	for (i = 0; i < NUM; i++)
	{
		for (j = 0; j < NUM; j++)	C_re[i][j] = temp_A[i][j] + temp_BAB[i][j];
	}
	INV(C_re, NUM);
	Matrix_Mul(temp_AB, C_re, C_im, NUM);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++) C_im[i][j] = -C_im[i][j];
	}
}


typedef struct
{
	Word16  buf_mic[NUM][ORD2];
	Word16  buf_spk[ORD2];
	Float32 spk_ana_re[D][ORD3];
	Float32 spk_ana_im[D][ORD3];
	Float32 spk_spec[D][ORD3];
	Float32 spk_ener[ORD3];
	Float32 fbf_re[D][ORD6];
	Float32 fbf_im[D][ORD6];
	Float32 bmIN_re[NUM][D][ORD6];
	Float32 bmIN_im[NUM][D][ORD6];
	Float32 pn;

	Float32 mic0_power[D];
	Float32 spk_power[D];
	Float32 e0_power[D];
	Float32 mic0_e0_power_re[D];
	Float32 mic0_e0_power_im[D];
	Float32 mic0_spk_power_re[D];
	Float32 mic0_spk_power_im[D];

	Word16 divergeState;
	Word16 micState;
	Word16 echoState;

	Float32 hNlMic0SpkAvgMin;
	Float32 hNlFbMin, hNlFbLocalMin;
	//Word16 hNlNewMin, hNlMinCtr;
	//Float32 overDrive, overDriveSm;
	Word16 nlp_mode;
	Word16 mult;

	Float32 h_r[NUM][D][ORD3];
	Float32 h_i[NUM][D][ORD3];

	Float32 hb_r[NUM][D][ORD3];
	Float32 hb_i[NUM][D][ORD3];

	Float32 echo_r[D][ORD4];
	Float32 echo_i[D][ORD4];
	Float32 noise_spectrum[NUM][D];
	Float32 gr[NUM][D];
	Float32 gi[NUM][D];
	Float32 hr[NUM][D];
	Float32 hi[NUM][D];
	Float32 sr[NUM][D][ORD6];
	Float32 si[NUM][D][ORD6];
	Float32 syn[6][ORD2 * 2];
	Float32 syn_all_around[ORD2 * 2];

	Float32 sub_az[D];
	Float32 audio_pit;
	Float32 audio_az[MED_NUM];
	Float32 video_pit;
	Float32 video_az;
	Word32  source_category;
	Word16  init_farme_cnt;
	Word16  voice_frame_cnt;

	NS_STRUCT my_str_anc;
	NS_STRUCT anc_all_around;
	Agc_t my_str_agc;
	Agc_t agc_all_around;

	Float32 pearson;
	Float32 alpha;

	Word16 frame_success;
	Word16 frame_noise;
	Word16 frame_voice;
	Word16 frame_locate;
	///VAD
	Word16 vad_first;
	Word16 hyster_cnt;
	Word16 last_update_cnt;
	Word16 update_cnt;
	Word32 frame_cnt;
	Float32 ch_enrg[NUM_CHAN];
	Float32 ch_noise[NUM_CHAN];
	Float32 ch_enrg_long_db[NUM_CHAN];

	Float32 qxx_re[D][NUM][NUM];/// 
	Float32 qxx_im[D][NUM][NUM];///
	Float32 qnn_re[D][NUM][NUM];/// 
	Float32 qnn_im[D][NUM][NUM];///

	double aec_time;
	double srp_time;
	double mvdr_time;

	DOWNSAMPLE_STR down_str;
	INTERP_STR interp_str;
	DenoiseState *rnn_str;

	//Word32 rnn_vad;
	Word32 bf_first;
	Float32 bf_re[NUM][D], bf_im[NUM][D];

	//RLS
	Float32 p_rls_re[D][ORD3][ORD3];
	Float32 p_rls_im[D][ORD3][ORD3];
	Float32 h_rls_re[D][ORD3];
	Float32 h_rls_im[D][ORD3];
	Float32 k_rls_re[D][ORD3];
	Float32 k_rls_im[D][ORD3];
} AEC_SRP_ST;



Word32 frame;
Float32 prototype_filter[ORD2];
Float32 cos_tab0[D][2 * D], sin_tab0[D][2 * D];
Float32 cos_tab1[5][5], sin_tab1[5][5];
Float32 cos_tab2[D], sin_tab2[D];
Float32 wr[AZ_NUM][NUM][D], wi[AZ_NUM][NUM][D]; //AZ_NUM*NUM*FFT_LEN/2
int inMicLevel = 0, outMicLevel = 0;
uint8_t saturationWarning;
const Word16 rev[320] = {
0,64,128,192,256,32,96,160,224,288,16,80,144,208,272,48,112,176,240,304,8,72,136,200,264,40,104,168,232,296,
24,88,152,216,280,56,120,184,248,312,4,68,132,196,260,36,100,164,228,292,20,84,148,212,276,52,116,180,244,308,
12,76,140,204,268,44,108,172,236,300,28,92,156,220,284,60,124,188,252,316,2,66,130,194,258,34,98,162,226,290,
18,82,146,210,274,50,114,178,242,306,10,74,138,202,266,42,106,170,234,298,26,90,154,218,282,58,122,186,250,314,
6,70,134,198,262,38,102,166,230,294,22,86,150,214,278,54,118,182,246,310,14,78,142,206,270,46,110,174,238,302,
30,94,158,222,286,62,126,190,254,318,1,65,129,193,257,33,97,161,225,289,17,81,145,209,273,49,113,177,241,305,
9,73,137,201,265,41,105,169,233,297,25,89,153,217,281,57,121,185,249,313,5,69,133,197,261,37,101,165,229,293,
21,85,149,213,277,53,117,181,245,309,13,77,141,205,269,45,109,173,237,301,29,93,157,221,285,61,125,189,253,317,
3,67,131,195,259,35,99,163,227,291,19,83,147,211,275,51,115,179,243,307,11,75,139,203,267,43,107,171,235,299,
27,91,155,219,283,59,123,187,251,315,7,71,135,199,263,39,103,167,231,295,23,87,151,215,279,55,119,183,247,311,
15,79,143,207,271,47,111,175,239,303,31,95,159,223,287,63,127,191,255,319 };
const Float32 ptrGCoh[2] = { 0.8f, 0.2f };
const Float32 kTargetSupp[3] = { -6.9f, -11.5f, -18.4f };
const Float32 min_overdrive[3] = { 1.0f, 2.0f, 5.0f };
const Float32 weightCurve[160] = {
	0.0000f, 0.1000f, 0.1239f, 0.1338f, 0.1413f, 0.1477f,
	0.1534f, 0.1585f, 0.1631f, 0.1675f, 0.1716f, 0.1755f,
	0.1792f, 0.1827f, 0.1861f, 0.1893f, 0.1924f, 0.1955f,
	0.1984f, 0.2013f, 0.2040f, 0.2067f, 0.2094f, 0.2119f,
	0.2145f, 0.2169f, 0.2193f, 0.2217f, 0.2240f, 0.2263f,
	0.2285f, 0.2307f, 0.2329f, 0.2350f, 0.2371f, 0.2392f,
	0.2412f, 0.2432f, 0.2452f, 0.2471f, 0.2490f, 0.2509f,
	0.2528f, 0.2547f, 0.2565f, 0.2583f, 0.2601f, 0.2619f,
	0.2636f, 0.2654f, 0.2671f, 0.2688f, 0.2704f, 0.2721f,
	0.2738f, 0.2754f, 0.2770f, 0.2786f, 0.2802f, 0.2818f,
	0.2833f, 0.2849f, 0.2864f, 0.2879f, 0.2894f, 0.2909f,
	0.2924f, 0.2939f, 0.2954f, 0.2968f, 0.2983f, 0.2997f,
	0.3011f, 0.3025f, 0.3039f, 0.3053f, 0.3067f, 0.3081f,
	0.3094f, 0.3108f, 0.3121f, 0.3135f, 0.3148f, 0.3161f,
	0.3174f, 0.3187f, 0.3200f, 0.3213f, 0.3226f, 0.3239f,
	0.3252f, 0.3264f, 0.3277f, 0.3289f, 0.3302f, 0.3314f,
	0.3326f, 0.3338f, 0.3351f, 0.3363f, 0.3375f, 0.3387f,
	0.3399f, 0.3410f, 0.3422f, 0.3434f, 0.3446f, 0.3457f,
	0.3469f, 0.3480f, 0.3492f, 0.3503f, 0.3515f, 0.3526f,
	0.3537f, 0.3548f, 0.3559f, 0.3571f, 0.3582f, 0.3593f,
	0.3604f, 0.3614f, 0.3625f, 0.3636f, 0.3647f, 0.3658f,
	0.3668f, 0.3679f, 0.3690f, 0.3700f, 0.3711f, 0.3721f,
	0.3732f, 0.3742f, 0.3752f, 0.3763f, 0.3773f, 0.3783f,
	0.3794f, 0.3804f, 0.3814f, 0.3824f, 0.3834f, 0.3844f,
	0.3854f, 0.3864f, 0.3874f, 0.3884f, 0.3894f, 0.3904f,
	0.3913f, 0.3923f, 0.3933f, 0.3942f, 0.3952f, 0.3962f,
	0.3971f, 0.3981f, 0.3990f, 0.4000f
};

///VAD
//const Word16 ch_tbl[NUM_CHAN][2] = {{ 2,  3}, { 4,  5}, { 6,  7}, { 8,  9}, {10, 11}, {12, 13}, {14, 16}, {17, 19}, {20, 22}, {23, 26}, {27, 30}, {31, 35}, {36, 41}, {42, 48}, {49, 55}, {56, 63}};
//const Word16 ch_tbl[NUM_CHAN][2] = {{ 2,  3}, { 4,  5}, { 6,  8}, { 9,  11}, {12, 14}, {15, 17}, {18, 20}, {21, 24}, {25, 28}, {29, 33}, {34, 38}, {39, 44}, {45, 51}, {52, 60}, {61, 69}, {70, 79}};
const Word16 ch_tbl[NUM_CHAN][2] = { {2, 3}, {4, 5}, {6, 7}, {8, 10}, {11, 13}, {13, 15}, {16, 18}, {19, 22}, {22, 25}, {26, 30}, {30, 34}, {35, 40}, {41, 46}, {47, 54}, {55, 63}, {64, 74}, {75, 88}, {89, 106}, {107, 128}, {129, 154} };
//const Word16 ch_tbl[NUM_CHAN][2] = {{ 2,  3}, { 4,  5}, { 6,  7}, { 8,  10}, {11, 13}, {14, 16}, {17, 20}, {21, 24}, {25, 28}, {29, 33}, {34, 38}, {39, 44}, {45, 51}, {52, 60}, {61, 69}, {70, 79}};

const Word16 vm_tbl[90] = { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 8, 8, 9, 9, 10, 10, 11, 12, 12, 13, 13, 14, 15, 15, 16, 17, 17, 18, 19, 20, 20, 21, 22, 23, 24, 24, 25, 26, 27, 28, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50 };

/*const Float32 overDriveCurve[160]={
	1.0000, 1.0793, 1.1122, 1.1374, 1.1586, 1.1773,
	1.1943, 1.2098, 1.2243, 1.2379, 1.2508, 1.2630,
	1.2747, 1.2859, 1.2967, 1.3071, 1.3172, 1.3270,
	1.3365, 1.3457, 1.3547, 1.3634, 1.3720, 1.3803,
	1.3885, 1.3965, 1.4044, 1.4121, 1.4196, 1.4271,
	1.4344, 1.4416, 1.4486, 1.4556, 1.4624, 1.4692,
	1.4758, 1.4824, 1.4889, 1.4953, 1.5016, 1.5078,
	1.5140, 1.5200, 1.5261, 1.5320, 1.5379, 1.5437,
	1.5494, 1.5551, 1.5608, 1.5664, 1.5719, 1.5774,
	1.5828, 1.5881, 1.5935, 1.5987, 1.6040, 1.6092,
	1.6143, 1.6194, 1.6244, 1.6295, 1.6344, 1.6394,
	1.6443, 1.6491, 1.6540, 1.6588, 1.6635, 1.6682,
	1.6729, 1.6776, 1.6822, 1.6868, 1.6914, 1.6959,
	1.7004, 1.7049, 1.7093, 1.7137, 1.7181, 1.7225,
	1.7268, 1.7312, 1.7354, 1.7397, 1.7439, 1.7482,
	1.7524, 1.7565, 1.7607, 1.7648, 1.7689, 1.7730,
	1.7770, 1.7811, 1.7851, 1.7891, 1.7931, 1.7970,
	1.8009, 1.8049, 1.8088, 1.8126, 1.8165, 1.8203,
	1.8242, 1.8280, 1.8318, 1.8355, 1.8393, 1.8430,
	1.8467, 1.8505, 1.8541, 1.8578, 1.8615, 1.8651,
	1.8687, 1.8724, 1.8760, 1.8795, 1.8831, 1.8867,
	1.8902, 1.8937, 1.8972, 1.9007, 1.9042, 1.9077,
	1.9111, 1.9146, 1.9180, 1.9214, 1.9248, 1.9282,
	1.9316, 1.9350, 1.9384, 1.9417, 1.9450, 1.9484,
	1.9517, 1.9550, 1.9582, 1.9615, 1.9648, 1.9680,
	1.9713, 1.9745, 1.9777, 1.9810, 1.9842, 1.9873,
	1.9905, 1.9937, 1.9969, 2.0000
};*/

//clock_t start1, start2, start3, finish1, finish2, finish3;

//fft of 320
void fft_320(Float32 *tmp1, Float32 *tmp2)
{
	Word16 i, k, j, L, b, r, m;
	Float32 acc2, acc3, x_re[320], x_im[320], x_re1[320], x_im1[320], tr, ti, sr, si, rr, ri;

	for (i = 0; i < 320; i++) x_re1[i] = tmp1[rev[i]];
	for (i = 0; i < 320; i++) x_im1[i] = tmp2[rev[i]];

	for (i = 0; i < 64; i++)
	{
		for (k = 0; k < 5; k++)
		{
			for (acc2 = 0.0f, j = 0; j < 5; j++) acc2 += x_re1[j + i * 5] * cos_tab1[k][j] - x_im1[j + i * 5] * sin_tab1[k][j];
			x_re[k + i * 5] = acc2;
			for (acc3 = 0.0f, j = 0; j < 5; j++) acc3 += x_re1[j + i * 5] * sin_tab1[k][j] + x_im1[j + i * 5] * cos_tab1[k][j];
			x_im[k + i * 5] = acc3;
		}
	}

	for (b = L = 1; L < 7; L++)
	{
		b = 2 * b;//b=2^L
		m = 5 * (b / 2);
		for (j = 0; j < (64 / b); j++)
		{
			for (i = 0; i < m; i++)
			{
				r = i * 64 / b;

				rr = cos_tab2[r];
				ri = sin_tab2[r];

				tr = x_re[j * 2 * m + i];
				ti = x_im[j * 2 * m + i];
				sr = x_re[j * 2 * m + i + m];
				si = x_im[j * 2 * m + i + m];

				x_re[j * 2 * m + i] = tr + sr * rr - si * ri;
				x_im[j * 2 * m + i] = ti + sr * ri + si * rr;
				x_re[j * 2 * m + i + m] = tr - sr * rr + si * ri;
				x_im[j * 2 * m + i + m] = ti - sr * ri - si * rr;
			}
		}
	}

	memcpy(tmp1, x_re, 320 * 4);
	memcpy(tmp2, x_im, 320 * 4);
}

static int CmpFloat(const void* a, const void* b) {
	const float* da = (const float*)a;
	const float* db = (const float*)b;

	return (*da > *db) - (*da < *db);
}




void aec_srp_gsc_init(AEC_SRP_ST *st,            //褰撳墠寰呭鐞嗛�氶亾鍙橀噺缁撴瀯浣撶殑鎸囬拡
	Float32 appointed_pitch   //鎸囧畾淇话瑙? 鍙栧�艰寖鍥? [-PI/2, PI/2], 鍗曚綅:寮у害
)
{
	Word32 i, j, k, m, n;
	Float32 theta, tao;

	for (j = 0; j < ORD2; j++) prototype_filter[j] = (Float32)(0.54 - 0.46*cos((2.0*j + 1) / ORD2)) * (Float32)sin(PI*(2 * j - ORD2 + 1) / (4 * D)) / (Float32)(PI*(2 * j - ORD2 + 1) / 2.0);

	//甯告暟鏁扮粍鍒濆鍖?
	for (k = 0; k < D; k++)
	{
		cos_tab2[k] = (Float32)cos(2 * PI*k / (2 * D));
		sin_tab2[k] = (Float32)sin(2 * PI*k / (2 * D));
		for (j = 0; j < 2 * D; j++) cos_tab0[k][j] = (Float32)cos(2 * PI*k*j / (2 * D));
		for (j = 0; j < 2 * D; j++) sin_tab0[k][j] = (Float32)sin(2 * PI*k*j / (2 * D));
	}

	for (k = 0; k < 5; k++)
	{
		for (j = 0; j < 5; j++) cos_tab1[k][j] = (Float32)cos(2 * PI*k*j / 5);
		for (j = 0; j < 5; j++) sin_tab1[k][j] = (Float32)sin(2 * PI*k*j / 5);
	}

	//鏍规嵁鎸囧畾鐨勪刊浠拌, 鍒濆鍖栧悇涓柟浣嶈, 棰戠巼鐐圭殑鏉冮噸鍥犲瓙
	for (n = 0; n < AZ_NUM; n++)
	{
		theta = (Float32)(n * 2 * PI / AZ_NUM);
		for (i = 0; i < NUM; i++)
		{
			tao = (Float32)(PI / D)*(Float32)(FS*DIAMETER / 2 / SPEED)*(Float32)cos(appointed_pitch)*(Float32)cos(theta - 2 * PI*i / NUM);
			for (k = 0; k < D; k++)
			{
				wr[n][i][k] = (Float32)cos(k*tao);
				wi[n][i][k] = (Float32)sin(k*tao);
			}
		}
	}

	//閫氶亾鍙橀噺鍒濆鍖?
	for (j = 0; j < ORD2; j++) st->buf_spk[j] = 0;
	for (m = 0; m < NUM; m++) { for (j = 0; j < ORD2; j++) st->buf_mic[m][j] = 0; }
	for (k = 0; k < D; k++) { for (j = 0; j < ORD3; j++) st->spk_ana_re[k][j] = st->spk_ana_im[k][j] = 0; }

	for (m = 0; m < NUM; m++)
	{
		for (k = 0; k < D; k++) { for (j = 0; j < ORD3; j++) st->hb_r[m][k][j] = st->hb_i[m][k][j] = st->h_r[m][k][j] = st->h_i[m][k][j] = 0; }
	}
	for (k = 0; k < D; k++) { for (j = 0; j < ORD4; j++) st->echo_r[k][j] = st->echo_i[k][j] = 0; }
	for (k = 0; k < D; k++) { for (j = 0; j < ORD3; j++) st->spk_spec[k][j] = 0; }
	for (j = 0; j < ORD3; j++) st->spk_ener[j] = 0;

	for (k = 1; k < D; k++)
	{
		for (j = 0; j < ORD6; j++)
		{
			st->fbf_re[k][j] = 0;
			st->fbf_im[k][j] = 0;
		}
	}

	for (m = 0; m < NUM; m++)
	{
		for (j = 0; j < D; j++)
		{
			st->gr[m][j] = 0;
			st->gi[m][j] = 0;
			st->hr[m][j] = 0;
			st->hi[m][j] = 0;
		}
	}

	for (m = 0; m < NUM; m++)
	{
		for (k = 1; k < D; k++)
		{
			for (j = 0; j < ORD6; j++)
			{
				st->sr[m][k][j] = 0;
				st->si[m][k][j] = 0;
			}
		}
	}

	for (m = 0; m < NUM; m++)
	{
		for (k = 0; k < D; k++)
		{
			for (j = 0; j < ORD6; j++)
			{
				st->bmIN_re[m][k][j] = 0;
				st->bmIN_im[m][k][j] = 0;
			}
		}
	}

	for (k = 0; k < D; k++) st->sub_az[k] = 0;
	st->audio_pit = 0;
	for (j = 0; j < MED_NUM; j++) st->audio_az[j] = 0;
	st->video_pit = -1;
	st->video_az = -1;
	st->source_category = HUMAN_FACE;

	for (i = 0; i < 6; i++)
	{
		for (j = 0; j < ORD2 * 2; j++) st->syn[i][j] = 0;
	}
	for (j = 0; j < ORD2 * 2; j++)
	{
		st->syn_all_around[j] = 0;
	}

	st->init_farme_cnt = 0;
	st->voice_frame_cnt = 0;
	st->alpha = 1;

	NS_Init(&(st->my_str_anc), 16000, 2);
	NS_Init(&(st->anc_all_around), 16000, 2);
	AGC_init(&(st->my_str_agc));
	AGC_init(&(st->agc_all_around));

	for (k = 0; k < D; k++)
	{
		st->mic0_power[k] = 0;
		st->spk_power[k] = 0;
		st->e0_power[k] = 0;
		st->mic0_e0_power_re[k] = 0;
		st->mic0_e0_power_im[k] = 0;
		st->mic0_spk_power_re[k] = 0;
		st->mic0_spk_power_im[k] = 0;
	}
	st->divergeState = 0;
	st->micState = 0;
	st->echoState = 0;
	st->hNlMic0SpkAvgMin = 1;
	st->nlp_mode = 1;
	st->hNlFbMin = 1;
	st->hNlFbLocalMin = 1;
	//st->hNlNewMin = 0;
	//st->hNlMinCtr = 0;
	//st->overDrive = 2;
	//st->overDriveSm = 2;
	st->mult = 2;

	st->pn = 0.0f;

	st->frame_success = 0;
	st->frame_noise = 0;
	st->frame_voice = 0;
	st->frame_locate = 0;
	///VAD
	st->vad_first = TRUE;
	st->hyster_cnt = 0;
	st->last_update_cnt = 0;
	st->update_cnt = 0;
	st->frame_cnt = 0;
	for (i = 0; i < NUM_CHAN; i++) st->ch_enrg_long_db[i] = st->ch_enrg[NUM_CHAN] = st->ch_noise[NUM_CHAN] = 0;

	for (i = 0; i < D; i++)///
	{
		for (m = 0; m < NUM; m++)
		{
			for (n = 0; n < NUM; n++) st->qxx_re[i][m][n] = st->qxx_im[i][m][n] = st->qnn_re[i][m][n] = st->qnn_im[i][m][n] = 0.0;
			st->qnn_re[i][m][m] = 1.0;
		}
	}

	st->aec_time = 0;
	st->srp_time = 0;
	st->mvdr_time = 0;

	down_sample_init(&(st->down_str));
	interp_init(&(st->interp_str));
	st->rnn_str = rnnoise_create(NULL);

	//st->rnn_vad = 0;
	st->bf_first = 1;
	for (m = 0; m < NUM; m++)
	{
		for (k = 1; k < D; k++)
		{
			st->bf_re[m][k] = 0.0f;
			st->bf_im[m][k] = 0.0f;
		}
	}

	//RLS
	for (k = 0; k < D; k++)
	{
		for (i = 0; i < ORD3; i++)
		{
			for (j = 0; j < ORD3; j++)
			{
				st->p_rls_re[k][i][j] = 0.0f;
				st->p_rls_im[k][i][j] = 0.0f;
			}
		}
	}
	for (k = 0; k < D; k++)
	{
		for (j = 0; j < ORD3; j++)
		{
			st->h_rls_re[k][j] = 0.0f;
			st->h_rls_im[k][j] = 0.0f;
			st->p_rls_re[k][j][j] = 1000.0f;
		}
	}
}

void aec_srp_gsc_free(AEC_SRP_ST* st)
{
	rnnoise_destroy(st->rnn_str);
}


void set_position(AEC_SRP_ST *st,          //褰撳墠寰呭鐞嗛�氶亾鍙橀噺缁撴瀯浣撶殑鎸囬拡
	Float32 video_pit,       //鍥惧儚璇嗗埆鍑虹殑澹版簮淇话瑙? 鍙栧�艰寖鍥? [-PI/2, PI/2], 鍗曚綅:寮у害
	Float32 video_az,        //鍥惧儚璇嗗埆鍑虹殑澹版簮鏂逛綅瑙? 鍙栧�艰寖鍥? [0, 2*PI), 鍗曚綅:寮у害
	Word32 source_category  //鍥惧儚璇嗗埆鍑虹殑澹版簮绫诲埆, 鍙栧�艰寖鍥? 瑙丼OUND_SOURCE_CATEGORY, 鍗曚綅:寮у害
)
{
	st->video_pit = video_pit;
	st->video_az = video_az;
	st->source_category = source_category;
}

/////////////////////////////////////////////////////////////////////////
//鑾峰彇闊抽瀹氫綅缁撴灉
//杩斿洖鍊? 闊抽瀹氫綅鐨勬柟浣嶈
/////////////////////////////////////////////////////////////////////////
Float32 get_position(AEC_SRP_ST *st        //褰撳墠寰呭鐞嗛�氶亾鍙橀噺缁撴瀯浣撶殑鎸囬拡
)
{
	Word32 i, j;
	Float32 tmp[MED_NUM], acc2;

	for (j = 0; j < MED_NUM; j++) tmp[j] = st->audio_az[j];

	for (i = 0; i < MED_NUM - 1; i++)
	{
		for (j = MED_NUM - 1; j > i; j--)
		{
			if (tmp[j] < tmp[j - 1])
			{
				acc2 = tmp[j];
				tmp[j] = tmp[j - 1];
				tmp[j - 1] = acc2;
			}
		}
	}

	return(tmp[MED_NUM / 2]);
}

void aec_srp_gsc(AEC_SRP_ST *st,        //褰撳墠寰呭鐞嗛�氶亾鍙橀噺缁撴瀯浣撶殑鎸囬拡
	Word16 mic_sp[NUM][FRM_LEN],
	Word16 spk_sp[FRM_LEN],      //褰撳墠寰呭鐞嗗枃鍙�氶亾杈撳叆淇″彿
	Word16 out_sp[6][FRM_LEN],      //褰撳墠澧炲己鍚庣殑淇″彿
	Word16 all_around_out[FRM_LEN],
	Word16 vad_out[FRM_LEN]
)
{
	Float32 mic_ana_re[NUM][D], mic_ana_im[NUM][D];    //楹﹀厠淇″彿鍒嗘瀽婊ゆ尝杈撳嚭
	Float32 tmp[2 * D], tmp2[2 * D], e_r[NUM][D], e_i[NUM][D], ener[D], mic0_spec[D], spectrum[NUM][D], spectrum2[NUM][D], mic0_spec1[D];
	Float32 re[D], im[D], re2[NUM][D], im2[NUM][D], re3[NUM][D], im3[NUM][D], re4[NUM][D], im4[NUM][D], re5[2 * D], im5[2 * D], re6[NUM][D], im6[NUM][D], re_all_around[2 * D], im_all_around[2 * D];
	Float32 acc2, acc3, pearson, aec_step, mic0_ener, cxx[NUM], cxy, cyy, peak, att[D], power_sub[D], peak_sub[D];
	Word16 temp1[FRM_LEN], temp2[FRM_LEN], tmp1_all_around[FRM_LEN], tmp2_all_around[FRM_LEN];
	Word32  i, j, k, m, n, vad_cnt, tot_az, sub_azimuth, num_sub[D];
	Float32 mic0_power_sum, e0_power_sum;
	Float32 coh_mic0_e0[D], coh_mic0_spk[D];
	Float32 hNl_mic0_e0_Avg, hNl_mic0_spk_Avg;
	Float32 hNlFb = 0, hNlFblow = 0;
	Float32 hNl[D]; //鏉冮噸婊ゆ尝鍣?
	Word16 minPrefBand = 2, prefBandSize = 60;
	Float32 hNlPref[60];
	Float32 prefBandQuant = 0.75f, prefBandQuantLow = 0.5f;
	Float32 energy[D];
	Float32 pn, ps;
	Float32 snr;
	Word16 az_num;
	Float32 temp;
	Float32 weight_re[NUM], weight_im[NUM];
	Float32 den_re, den_im;
	Float32 tao[NUM], phrase_re[NUM][D], phrase_im[NUM][D];
	///VAD
	Word16 update_flag, vm_sum, ch_snr[NUM_CHAN];
	Float32 ch_enrg_db[NUM_CHAN], tne, tce, vv, ch_enrg_dev, alpha;

	Float32 tr, ti, matrix_re[NUM], matrix_im[NUM], Fmvdr_re[NUM], Fmvdr_im[NUM];///
	double inv1_re[NUM][NUM], inv1_im[NUM][NUM], inv_re[NUM][NUM], inv_im[NUM][NUM];///

	clock_t start1, start2, start3, finish1, finish2, finish3;

	Word16 tmp16[FRM_LEN], tmp48[FRM_LEN1];
	Float32 rnn_tmp[FRM_LEN1];
	Float32 vad_prob;

	//RLS
	Float32 rls_re[D], rls_im[D];
	Float32 rls_tmp_re[D][ORD3], rls_tmp_im[D][ORD3];
	Float32 klman_den;
	Float32 p_tmp_re[ORD3][ORD3], p_tmp_im[ORD3][ORD3];


	start1 = clock();
	//鏇存柊鍠囧彮淇″彿缂撳瓨
	for (j = 0; j < ORD2 - FRM_LEN; j++) st->buf_spk[j] = st->buf_spk[j + FRM_LEN];
	for (j = 0; j < FRM_LEN; j++) st->buf_spk[ORD2 - FRM_LEN + j] = spk_sp[j];

	//鏇存柊楹﹀厠淇″彿缂撳瓨
	for (m = 0; m < NUM; m++)
	{
		for (j = 0; j < ORD2 - FRM_LEN; j++) st->buf_mic[m][j] = st->buf_mic[m][j + FRM_LEN];
		for (j = 0; j < FRM_LEN; j++) st->buf_mic[m][ORD2 - FRM_LEN + j] = mic_sp[m][j];
	}

	for (k = 1; k < D; k++)
	{
		//鏇存柊鍠囧彮淇″彿鍒嗘瀽婊ゆ尝鍣ㄨ緭鍑虹紦瀛?
		for (j = 0; j < ORD3 - 1; j++) st->spk_ana_re[k][j] = st->spk_ana_re[k][j + 1];
		for (j = 0; j < ORD3 - 1; j++) st->spk_ana_im[k][j] = st->spk_ana_im[k][j + 1];
		//鏇存柊鍠囧彮淇″彿甯у箙搴﹁氨
		for (j = 0; j < ORD3 - 1; j++) st->spk_spec[k][j] = st->spk_spec[k][j + 1];
	}

	//鏇存柊鍠囧彮淇″彿甯ц兘閲忕紦瀛?
	for (j = 0; j < ORD3 - 1; j++) st->spk_ener[j] = st->spk_ener[j + 1];

	//鍠囧彮淇″彿瀛愬甫鍒嗘瀽
	for (i = 0; i < 2 * D; i++)
	{
		for (acc2 = 0.0f, j = 0; j < ORD2 / (2 * D); j++)	acc2 += st->buf_spk[j * 2 * D + i] * prototype_filter[j * 2 * D + i];
		tmp[2 * D - 1 - i] = acc2;
		tmp2[i] = 0;
	}

	fft_320(tmp, tmp2);

	for (k = 0; k < D; k++)
	{
		st->spk_ana_re[k][ORD3 - 1] = tmp[k];
		st->spk_ana_im[k][ORD3 - 1] = tmp2[k];
	}

	//楹﹀厠淇″彿瀛愬甫鍒嗘瀽
	for (m = 0; m < NUM; m++)
	{
		for (i = 0; i < 2 * D; i++)
		{
			for (acc2 = 0.0f, j = 0; j < ORD2 / (2 * D); j++)	acc2 += st->buf_mic[m][j * 2 * D + i] * prototype_filter[j * 2 * D + i];
			tmp[2 * D - 1 - i] = acc2;
			tmp2[i] = 0;
		}
		fft_320(tmp, tmp2);
		for (k = 1; k < D; k++)
		{
			//for(acc2=0.0f, j=0; j<2*D; j++) acc2 += tmp[j]*cos_tab0[k][j];
			mic_ana_re[m][k] = tmp[k];
			//for(acc3=0.0f, j=0; j<2*D; j++) acc3 += tmp[j]*sin_tab0[k][j];
			mic_ana_im[m][k] = tmp2[k];
		}
	}

	//浼拌鍠囧彮淇″彿骞呭害璋? 鎬昏兘閲?
	for (acc3 = 0.0f, k = 1; k < D; k++)
	{
		acc2 = st->spk_ana_re[k][ORD3 - 1] * st->spk_ana_re[k][ORD3 - 1] + st->spk_ana_im[k][ORD3 - 1] * st->spk_ana_im[k][ORD3 - 1];
		acc3 += acc2;
		st->spk_spec[k][ORD3 - 1] = (Float32)sqrt(acc2);
	}
	st->spk_ener[ORD3 - 1] = acc3;

	//浼拌绗?涓害鍏嬩俊鍙峰箙搴﹁氨, 鎬昏兘閲?
	for (mic0_ener = 0.0f, k = 1; k < D; k++)
	{
		acc2 = mic_ana_re[0][k] * mic_ana_re[0][k] + mic_ana_im[0][k] * mic_ana_im[0][k];
		mic0_ener += acc2;
		mic0_spec[k] = (Float32)sqrt(acc2);
	}

	//double talk妫�娴
	for (pearson = 0.0f, j = 0; j < ORD3; j++)
	{
		for (acc2 = 0.0f, k = 1; k < D; k++) acc2 += st->spk_spec[k][j] * mic0_spec[k];
		acc3 = (Float32)sqrt(mic0_ener*st->spk_ener[j]);
		acc2 = (acc3 > 1024.0f) ? acc2 / acc3 : 0;
		if (pearson < acc2)
			pearson = acc2;
	}

	pearson = (Float32)pow(pearson, 2);
	//st->pearson = pearson * pearson;

	//瀛愬甫鑷�傚簲婊ゆ尝鍣ㄥ彉姝ラ暱
	//if(pearson>0.9)      aec_step = 0.1f;
	//else if(pearson>0.7) aec_step = 0.05f;
	//else if(pearson>0.5) aec_step = 0.01f;
	//else if(pearson>0.3) aec_step = 0.006f;
	//else                 aec_step = 0.0001f;
	aec_step = (Float32)(10 * pearson + sqrt(10 * pearson))*(10 * pearson + sqrt(10 * pearson)) / 100.0f;
	//aec_step = (Float32)(7*pearson+sqrt(7*pearson))*(7*pearson+sqrt(7*pearson))/100.0f;

	//RLS
	for (k = 1; k < D; k++)
	{
		//滤波
		for (acc2 = 0.0f, j = 0; j < ORD3; j++)
		{
			acc2 += st->spk_ana_re[k][j] * st->h_rls_re[k][j] + st->spk_ana_im[k][j] * st->h_rls_im[k][j];
		}
		for (acc3 = 0.0f, j = 0; j < ORD3; j++)
		{
			acc3 += st->spk_ana_im[k][j] * st->h_rls_re[k][j] - st->spk_ana_re[k][j] * st->h_rls_im[k][j];
		}
		//计算误差
		rls_re[k] = mic_ana_re[0][k] - acc2;
		rls_im[k] = mic_ana_im[0][k] - acc3;
		//计算Kalman增益
		for (acc2 = acc3 = 0.0f, i = 0; i < ORD3; i++)
		{
			rls_tmp_re[k][i] = rls_tmp_im[k][i] = 0.0f;
			for (j = 0; j < ORD3; j++)
			{
				rls_tmp_re[k][i] += st->p_rls_re[k][i][j] * st->spk_ana_re[k][j] - st->p_rls_im[k][i][j] * st->spk_ana_im[k][j];
				rls_tmp_im[k][i] += st->p_rls_re[k][i][j] * st->spk_ana_im[k][j] + st->p_rls_im[k][i][j] * st->spk_ana_re[k][j];
			}
			acc2 += st->spk_ana_re[k][i] * rls_tmp_re[k][i] + st->spk_ana_im[k][i] * rls_tmp_im[k][i];
			acc3 += st->spk_ana_re[k][i] * rls_tmp_im[k][i] - st->spk_ana_im[k][i] * rls_tmp_re[k][i];
		}
		klman_den = sqrt(acc2 * acc2 + acc3 * acc3);
		for (j = 0; j < ORD3; j++)
		{
			st->k_rls_re[k][j] = rls_tmp_re[k][j] / (LAMBDA + klman_den);
			st->k_rls_im[k][j] = rls_tmp_im[k][j] / (LAMBDA + klman_den);
		}
		//更新滤波器系数
		for (j = 0; j < ORD3; j++)
		{
			st->h_rls_re[k][j] += st->k_rls_re[k][j] * rls_re[k] + st->k_rls_im[k][j] * rls_im[k];
			st->h_rls_im[k][j] += -st->k_rls_re[k][j] * rls_im[k] + st->k_rls_im[k][j] * rls_re[k];
		}
		//更新自相关的逆矩阵
		for (i = 0; i < ORD3; i++)
		{
			for (j = 0; j < ORD3; j++)
			{
				p_tmp_re[i][j] = st->k_rls_re[k][i] * st->spk_ana_re[k][j] + st->k_rls_im[k][i] * st->spk_ana_im[k][j];
				p_tmp_im[i][j] = -st->k_rls_re[k][i] * st->spk_ana_im[k][j] + st->k_rls_im[k][i] * st->spk_ana_re[k][j];
			}
		}
		for (i = 0; i < ORD3; i++)
		{
			for (j = 0; j < ORD3; j++)
			{
				for (acc2 = acc3 = 0.0f, n = 0; n < ORD3; n++)
				{
					acc2 += p_tmp_re[i][n] * st->p_rls_re[k][n][j] - p_tmp_im[i][n] * st->p_rls_im[k][n][j];
					acc3 += p_tmp_re[i][n] * st->p_rls_im[k][n][j] + p_tmp_im[i][n] * st->p_rls_re[k][n][j];
				}
				st->p_rls_re[k][i][j] = (st->p_rls_re[k][i][j] - acc2) / LAMBDA;
				st->p_rls_im[k][i][j] = (st->p_rls_im[k][i][j] - acc3) / LAMBDA;
			}
		}
	}

	//鏇存柊妯℃嫙鍥炴尝淇″彿瀛愬甫缂撳瓨
	for (k = 1; k < D; k++)
	{
		for (j = 0; j < ORD4 - 1; j++) st->echo_r[k][j] = st->echo_r[k][j + 1];
		for (j = 0; j < ORD4 - 1; j++) st->echo_i[k][j] = st->echo_i[k][j + 1];
	}

	//鍥炴尝鎶垫秷涓婚�氳矾mic 0
	for (k = 1; k < D; k++)
	{
		//浼拌鍠囧彮淇″彿瀛愬甫鑳介噺
		for (acc2 = 0.0f, j = 0; j < ORD3; j++) acc2 += st->spk_ana_re[k][j] * st->spk_ana_re[k][j] + st->spk_ana_im[k][j] * st->spk_ana_im[k][j];
		ener[k] = acc2;

		//浼拌妯℃嫙鍥炴尝瀛愬甫淇″彿
		for (acc2 = 0.0f, j = 0; j < ORD3; j++) acc2 += st->spk_ana_re[k][j] * st->h_r[0][k][j] - st->spk_ana_im[k][j] * st->h_i[0][k][j];
		for (acc3 = 0.0f, j = 0; j < ORD3; j++) acc3 += st->spk_ana_im[k][j] * st->h_r[0][k][j] + st->spk_ana_re[k][j] * st->h_i[0][k][j];
		st->echo_r[k][ORD4 - 1] = acc2;
		st->echo_i[k][ORD4 - 1] = acc3;

		//浼拌璇樊瀛愬甫淇″彿
		e_r[0][k] = mic_ana_re[0][k] - acc2;
		e_i[0][k] = mic_ana_im[0][k] - acc3;
	}

	for (k = 1; k < D; k++)
	{
		//鍠囧彮淇″彿瀛愬甫鑳介噺瓒冲澶?
		if (ener[k] > 4.0f)
		{
			acc2 = aec_step / ener[k]; //归一化步长
			for (j = 0; j < ORD3; j++) st->h_r[0][k][j] += acc2 * (st->spk_ana_re[k][j] * e_r[0][k] + st->spk_ana_im[k][j] * e_i[0][k]);
			for (j = 0; j < ORD3; j++) st->h_i[0][k][j] += acc2 * (st->spk_ana_re[k][j] * e_i[0][k] - st->spk_ana_im[k][j] * e_r[0][k]);
		}
	}

	//鍥炴尝鎶垫秷娆￠�氳矾mic 1<->(NUM-1)

	//浼拌妯℃嫙鍥炴尝淇″彿瀛愬甫鑳介噺
	for (k = 1; k < D; k++)
	{

		for (acc2 = 0.0f, j = 0; j < ORD4; j++) acc2 += st->echo_r[k][j] * st->echo_r[k][j] + st->echo_i[k][j] * st->echo_i[k][j];
		ener[k] = acc2;
	}

	for (m = 1; m < NUM; m++)
	{
		for (k = 1; k < D; k++)
		{
			//浼拌妯℃嫙鍥炴尝瀛愬甫淇″彿
			for (acc2 = 0.0f, j = 0; j < ORD4; j++) acc2 += st->echo_r[k][j] * st->h_r[m][k][j] - st->echo_i[k][j] * st->h_i[m][k][j];
			for (acc3 = 0.0f, j = 0; j < ORD4; j++) acc3 += st->echo_i[k][j] * st->h_r[m][k][j] + st->echo_r[k][j] * st->h_i[m][k][j];

			//浼拌璇樊瀛愬甫淇″彿
			e_r[m][k] = mic_ana_re[m][k] - acc2;
			e_i[m][k] = mic_ana_im[m][k] - acc3;
		}
		for (k = 1; k < D; k++)
		{
			//鍠囧彮淇″彿瀛愬甫鑳介噺瓒冲澶?
			if (ener[k] > 4.0f)
			{
				acc2 = aec_step / ener[k];
				for (j = 0; j < ORD4; j++) st->h_r[m][k][j] += acc2 * (st->echo_r[k][j] * e_r[m][k] + st->echo_i[k][j] * e_i[m][k]);
				for (j = 0; j < ORD4; j++) st->h_i[m][k][j] += acc2 * (st->echo_r[k][j] * e_i[m][k] - st->echo_i[k][j] * e_r[m][k]);
			}
		}
	}

	//for (i = ORD2 * 2 - 1; i >= 2 * D; i--) st->syn[0][i] = st->syn[0][i - 2 * D];
	//for (i = 0; i < 2 * D; i++)
	//{
	//	for(acc2=0.0f, k=1; k<D; k++) acc2 += e_r[0][k]*cos_tab0[k][i] - e_i[0][k]*sin_tab0[k][i];
	//	st->syn[0][i]=acc2;	    
	//}
	//for (i = 0; i < FRM_LEN; i++)
	//{
	//	for(acc2=0.0f, j=0; j<ORD2/D; j++) acc2 += prototype_filter[j*D+i]*st->syn[0][j*2*D+(j&1)*D+i];
	//	out_sp[0][i]=(short)(acc2*D*32);
	//}


	//璁＄畻鍔熺巼璋卞拰浜掑姛鐜囪氨
	for (mic0_power_sum = e0_power_sum = 0, k = 1; k < D; k++)
	{
		st->mic0_power[k] = ptrGCoh[0] * st->mic0_power[k] + ptrGCoh[1] * (mic_ana_re[0][k] * mic_ana_re[0][k] + mic_ana_im[0][k] * mic_ana_im[0][k]);
		st->spk_power[k] = ptrGCoh[0] * st->spk_power[k] + ptrGCoh[1] * (st->spk_ana_re[k][ORD3 - 1] * st->spk_ana_re[k][ORD3 - 1] + st->spk_ana_im[k][ORD3 - 1] * st->spk_ana_im[k][ORD3 - 1]);
		st->e0_power[k] = ptrGCoh[0] * st->e0_power[k] + ptrGCoh[1] * (e_r[0][k] * e_r[0][k] + e_i[0][k] * e_i[0][k]);
		st->mic0_e0_power_re[k] = ptrGCoh[0] * st->mic0_e0_power_re[k] + ptrGCoh[1] * (mic_ana_re[0][k] * e_r[0][k] + mic_ana_im[0][k] * e_i[0][k]);
		st->mic0_e0_power_im[k] = ptrGCoh[0] * st->mic0_e0_power_im[k] + ptrGCoh[1] * (mic_ana_re[0][k] * e_i[0][k] - mic_ana_im[0][k] * e_r[0][k]);
		st->mic0_spk_power_re[k] = ptrGCoh[0] * st->mic0_spk_power_re[k] + ptrGCoh[1] * (mic_ana_re[0][k] * st->spk_ana_re[k][ORD3 - 1] + mic_ana_im[0][k] * st->spk_ana_im[k][ORD3 - 1]);
		st->mic0_spk_power_im[k] = ptrGCoh[0] * st->mic0_spk_power_im[k] + ptrGCoh[1] * (mic_ana_re[0][k] * st->spk_ana_im[k][ORD3 - 1] - mic_ana_im[0][k] * st->spk_ana_re[k][ORD3 - 1]);

		mic0_power_sum += st->mic0_power[k];
		e0_power_sum += st->e0_power[k];
	}

	//鍒ゅ埆鏄惁鍙戞暎
	if (st->divergeState == 0)
	{
		if (e0_power_sum > mic0_power_sum) { st->divergeState = 1; }
	}
	else
	{
		if (e0_power_sum * 1.05f < mic0_power_sum) { st->divergeState = 0; }
	}

	if (st->divergeState == 1)
	{
		memcpy(e_r[0], mic_ana_re[0], D * sizeof(float));
		memcpy(e_i[0], mic_ana_im[0], D * sizeof(float));
	}

	//璁＄畻鍚勪釜瀛愬甫鐨勪簰鐩稿叧
	for (k = 1; k < D; k++)
	{
		coh_mic0_e0[k] = (st->mic0_e0_power_re[k] * st->mic0_e0_power_re[k] + st->mic0_e0_power_im[k] * st->mic0_e0_power_im[k]) / (st->mic0_power[k] * st->e0_power[k] + 1e-10f);
		coh_mic0_spk[k] = (st->mic0_spk_power_re[k] * st->mic0_spk_power_re[k] + st->mic0_spk_power_im[k] * st->mic0_spk_power_im[k]) / (st->mic0_power[k] * st->spk_power[k] + 1e-10f);
	}

	//璁＄畻浜掔浉鍏冲钩鍧囧�?
	hNl_mic0_spk_Avg = 0;
	for (k = minPrefBand; k < minPrefBand + prefBandSize; k++)
	{
		hNl_mic0_spk_Avg += coh_mic0_spk[k];
	}
	hNl_mic0_spk_Avg /= prefBandSize;
	hNl_mic0_spk_Avg = 1 - hNl_mic0_spk_Avg;

	hNl_mic0_e0_Avg = 0;
	for (k = minPrefBand; k < minPrefBand + prefBandSize; k++)
	{
		hNl_mic0_e0_Avg += coh_mic0_e0[k];
	}
	hNl_mic0_e0_Avg /= prefBandSize;

	if (hNl_mic0_spk_Avg < 0.75f && hNl_mic0_spk_Avg < st->hNlMic0SpkAvgMin)
	{
		st->hNlMic0SpkAvgMin = hNl_mic0_spk_Avg;
	}

	//鍒ゆ柇鏄惁鏈変汉澹帮紝鏄惁鏈夊洖澹?
	if (hNl_mic0_e0_Avg > 0.18f && hNl_mic0_spk_Avg > 0.28f)
	{
		st->micState = 1;
	}
	else if (hNl_mic0_e0_Avg < 0.18f || hNl_mic0_spk_Avg < 0.28f)
	{
		st->micState = 0;
	}

	if (st->hNlMic0SpkAvgMin == 1)
	{
		st->echoState = 0;
		//st->overDrive = min_overdrive[st->nlp_mode];
		if (st->micState == 1)
		{
			memcpy(hNl, coh_mic0_e0, D * sizeof(float));
			hNlFb = hNl_mic0_e0_Avg;
			hNlFblow = hNl_mic0_e0_Avg;
		}
		else
		{
			for (k = 1; k < D; k++) { hNl[k] = (1 - coh_mic0_spk[k]); }
			hNlFb = hNl_mic0_spk_Avg;
			hNlFblow = hNl_mic0_spk_Avg;
		}
	}
	else
	{
		if (st->micState == 1)
		{
			st->echoState = 0;
			memcpy(hNl, coh_mic0_e0, D * sizeof(float));
			hNlFb = hNl_mic0_e0_Avg;
			hNlFblow = hNl_mic0_e0_Avg;
		}
		else
		{
			st->echoState = 1;
			for (k = 1; k < D; k++) { hNl[k] = 0.05*min(coh_mic0_e0[k], 1 - coh_mic0_spk[k]); }
			memcpy(hNlPref, &hNl[minPrefBand], sizeof(float)*prefBandSize);
			qsort(hNlPref, prefBandSize, sizeof(float), CmpFloat);
			hNlFb = hNlPref[(int)floor(prefBandQuant*(prefBandSize - 1))];
			hNlFblow = hNlPref[(int)floor(prefBandQuantLow*(prefBandSize - 1))];
		}
	}

	if (hNlFblow < 0.6f && hNlFblow < st->hNlFbLocalMin)
	{
		st->hNlFbLocalMin = hNlFblow;
		st->hNlFbMin = hNlFblow;
		//st->hNlNewMin = 1;
		//st->hNlMinCtr = 0;
	}

	st->hNlFbLocalMin = min(st->hNlFbLocalMin + 0.0008f / st->mult, 1);
	st->hNlMic0SpkAvgMin = min(st->hNlMic0SpkAvgMin + 0.0006f / st->mult, 1);

	//if(st->hNlNewMin == 1)
   // {
	//    st->hNlMinCtr++;
	//}

   // if(st->hNlMinCtr == 2)
   // {
	  //  st->hNlNewMin = 0;
	  //  st->hNlMinCtr = 0;
		//st->overDrive = max(kTargetSupp[st->nlp_mode]/((float)log(st->hNlFbMin + 1e-10f) + 1e-10f), min_overdrive[st->nlp_mode]);
	//}

	//Smooth the overdrive
	//if(st->overDrive < st->overDriveSm)
	//{
	//    st->overDriveSm = 0.99f * st->overDriveSm + 0.01f * st->overDrive;
	//}
	//else
	//{
	//    st->overDriveSm = 0.9f * st->overDriveSm + 0.1f * st->overDrive;
	//}

	for (k = 1; k < D; k++)
	{
		if (hNl[k] > hNlFb)
		{
			hNl[k] = weightCurve[k] * hNlFb + (1 - weightCurve[k]) * hNl[k];
		}

		hNl[k] *= 1.0f;
		//hNl[k] = pow(hNl[k], st->overDriveSm * overDriveCurve[k]);

		for (j = 0; j < NUM; j++)
		{
			e_r[j][k] *= (hNl[k]);
			e_i[j][k] *= (hNl[k]);
		}

		//e_r[0][k] *= (hNl[k]*1.5);
		//e_i[0][k] *= (hNl[k]*1.5);
	}

	//全向麦克风输出
	for (i = ORD2 * 2 - 1; i >= 2 * D; i--) st->syn_all_around[i] = st->syn_all_around[i - 2 * D];

	for (i = 1; i < D; i++)
	{
		re_all_around[2 * D - i] = re_all_around[i] = e_r[0][i];
		im_all_around[2 * D - i] = -e_i[0][i];
		im_all_around[i] = e_i[0][i];
	}
	re_all_around[D] = re_all_around[0] = im_all_around[D] = im_all_around[0] = 0;

	fft_320(re_all_around, im_all_around);

	for (i = 0; i < 2 * D; i++)
	{
		st->syn_all_around[i] = re_all_around[i];
	}

	for (i = 0; i < FRM_LEN; i++)
	{
		for (acc2 = 0.0f, j = 0; j < ORD2 / D; j++) acc2 += prototype_filter[j*D + i] * st->syn_all_around[j * 2 * D + (j & 1)*D + i];
		//out_sp[i]=(short)(acc2*D*32);
		tmp16[i] = tmp1_all_around[i] = (short)(acc2*D * 10);
		//temp1[i]=(short)(acc2*D*32);
	}
	interp_run(&(st->interp_str), tmp16, tmp48); //3倍插值
	for (i = 0; i < FRM_LEN1; i++) rnn_tmp[i] = tmp48[i];
	vad_prob = rnnoise_process_frame(st->rnn_str, rnn_tmp, rnn_tmp);//执行rnn降噪
	for (i = 0; i < FRM_LEN; i++) vad_out[i] = vad_prob > 0.3 ? 4000 : 0; //输出vad检测结果
	for (i = 0; i < FRM_LEN1; i++) tmp48[i] = rnn_tmp[i];
	down_sample_run(&(st->down_str), tmp48, tmp16);//3倍抽取
	for (i = 0; i < FRM_LEN; i++) all_around_out[i] = tmp16[i];

	/*NS_run(&(st->anc_all_around), tmp1_all_around, tmp2_all_around);
	AGC
	AGC_Process(&(st->agc_all_around), tmp2_all_around, NULL, FRM_LEN, all_around_out, NULL, inMicLevel, &outMicLevel, 0, (uint8_t*)(&saturationWarning));
	*/

	finish1 = clock();
	st->aec_time += (double)(finish1 - start1) / CLOCKS_PER_SEC;
	//for(m=0; m<NUM; m++)
	//{
	//	for(k=1; k<D; k++)
	//	{
	//		e_r[m][k] = mic_ana_re[m][k];
	//		e_i[m][k] = mic_ana_im[m][k];
	//	}
	//}

	for (i = ORD2 * 2 - 1; i >= 2 * D; i--) st->syn[0][i] = st->syn[0][i - 2 * D];
	for (i = 0; i < 2 * D; i++)
	{
		for(acc2=0.0f, k=1; k<D; k++) acc2 += e_r[0][k]*cos_tab0[k][i] - e_i[0][k]*sin_tab0[k][i];
		//for (acc2 = 0.0f, k = 1; k < D; k++) acc2 += rls_re[k] * cos_tab0[k][i] - rls_im[k] * sin_tab0[k][i];
		st->syn[0][i] = acc2;
	}
	for (i = 0; i < FRM_LEN; i++)
	{
		for (acc2 = 0.0f, j = 0; j < ORD2 / D; j++) acc2 += prototype_filter[j*D + i] * st->syn[0][j * 2 * D + (j & 1)*D + i];
		out_sp[0][i] = (short)(acc2*D * 20);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////
	//  澹版簮瀹氫綅
	//////////////////////////////////////////////////////////////////////////////////////////////////
	//浼拌骞呭害璋?

	start2 = clock();
	for (m = 0; m < NUM; m++)
	{
		for (cxx[m] = 0, k = 1; k < D; k++)
		{
			acc2 = e_r[m][k] * e_r[m][k] + e_i[m][k] * e_i[m][k];
			cxx[m] += acc2;
			spectrum[m][k] = (Float32)sqrt(acc2);
		}
	}

	//自动电平对齐
	for (acc2 = 0, m = 0; m < NUM; m++) acc2 += cxx[m];
	for (m = 0; m < NUM; m++)
	{
		if (cxx[m] > 0.0f)
		{
			acc3 = (Float32)sqrt(acc2 / NUM / cxx[m]);
			for (k = 1; k < D; k++) re6[m][k] = e_r[m][k] * acc3;
			for (k = 1; k < D; k++) im6[m][k] = e_i[m][k] * acc3;
		}
		else
		{
			for (k = 1; k < D; k++) re6[m][k] = 0;
			for (k = 1; k < D; k++) im6[m][k] = 0;
		}
	}

	//PHAT权值
	for (m = 0; m < NUM; m++)
	{
		for (k = 1; k < D; k++)
		{
			acc2 = (Float32)sqrt(re6[m][k] * re6[m][k] + im6[m][k] * im6[m][k]);
			if (acc2 > 0.0f)
			{
				re2[m][k] = re6[m][k] / acc2;
				im2[m][k] = im6[m][k] / acc2;
			}
		}
	}


	//计算功率权值
	for (peak = 0.0f, k = 1; k < D; k++)
	{
		acc3 = re6[0][k] * re6[0][k] + im6[0][k] * im6[0][k];
		mic0_spec1[k] = (Float32)sqrt(acc3);
		if (peak < mic0_spec1[k]) peak = mic0_spec1[k];
	}

	if (peak > 0.1f)
	{
		//for (k = 1; k < D; k++) mic0_spec[k] = max(2.0f + 1.0f*(Float32)log10(mic0_spec1[k] / peak), 0);
		//for (peak = 0.0f, k = 1; k < D; k++) { if (peak < mic0_spec[k]) peak = mic0_spec[k]; }
		//for (k = 1; k < D; k++) mic0_spec[k] = mic0_spec[k] / peak;
		for (k = 1; k < D; k++)
		{
			acc3 = mic0_spec1[k] / peak;
			if (acc3 >= 0.01f) mic0_spec[k] = 1;
			else mic0_spec[k] = 0;
		}
	}



	//for (i = ORD2 * 2 - 1; i >= 2 * D; i--) st->syn[0][i] = st->syn[0][i - 2 * D];
	//for (i = 0; i < 2 * D; i++)
	//{
	//	for(acc2=0.0f, k=1; k<D; k++) acc2 += e_r[0][k]*cos_tab0[k][i] - e_i[0][k]*sin_tab0[k][i];
	//	st->syn[0][i]=acc2;	    
	//}
	//for (i = 0; i < FRM_LEN; i++)
	//{
	//	for(acc2=0.0f, j=0; j<ORD2/D; j++) acc2 += prototype_filter[j*D+i]*st->syn[0][j*2*D+(j&1)*D+i];
	//	out_sp[0][i]=(short)(acc2*D*32);
	//}

	//old  VAD
	/*for(vad_cnt=m=0; m<NUM; m++)
	{
		if(st->init_farme_cnt<INIT_LEN)
		{
			for(k=1; k<D; k++) st->noise_spectrum[m][k] += spectrum[m][k]/INIT_LEN;
		}
		else
		{
			for(cyy=0, k=1; k<D; k++) cyy += st->noise_spectrum[m][k]*st->noise_spectrum[m][k];
			for(cxy=0, k=1; k<D; k++) cxy += spectrum[m][k]*st->noise_spectrum[m][k];
			if( (cxx[m]>4096.0f)&&(cyy>0.0f) )
			{
				//璁＄畻璋辩浉鍏崇郴鏁?
				pearson = cxy*cxy/cxx[m]/cyy;

				if( (pearson>0.6)&&(0.25<cxx[m]/cyy)&&(cxx[m]/cyy<4) )
				{
					//蹇�熸洿鏂拌儗鏅櫔澹?
					acc2 = min((Float32)pow(pearson, 6), 0.05f);
					for(k=1; k<D; k++) st->noise_spectrum[m][k] = (1.0f-acc2)*st->noise_spectrum[m][k] + acc2*spectrum[m][k];
				}
				else
				{
					vad_cnt++;
					//鎱㈤�熸洿鏂拌儗鏅櫔澹?
					acc2 = (Float32)pow(pearson, 14);
					for(k=1; k<D; k++) st->noise_spectrum[m][k] = (1.0f-acc2)*st->noise_spectrum[m][k] + acc2*spectrum[m][k];
				}
			}
		}
	}*/

	//start2 = clock();
	//new VAD
	alpha = (st->vad_first == TRUE) ? 1.0f : 0.55f;
	for (i = LO_CHAN; i <= HI_CHAN; i++)
	{
		for (vv = 0.0f, j = ch_tbl[i][0]; j <= ch_tbl[i][1]; j++) vv += re6[0][j] * re6[0][j] + im6[0][j] * im6[0][j];
		st->ch_enrg[i] = max((1.0f - alpha)*st->ch_enrg[i] + alpha * vv / (ch_tbl[i][1] - ch_tbl[i][0] + 1), MIN_CHAN_ENRG);
	}

	st->frame_cnt++;
	/* Initialize channel noise estimate to channel energy of vad_first four frames */
	if (st->frame_cnt < 6) { for (i = LO_CHAN; i <= HI_CHAN; i++) st->ch_noise[i] = max(st->ch_enrg[i], INE); }

	/* Compute the channel SNR indices */
	for (i = LO_CHAN; i <= HI_CHAN; i++) ch_snr[i] = (int)((max(10.0f*log10(st->ch_enrg[i] / st->ch_noise[i]), 0.0f) + 0.1875f) / 0.375f);

	/* Compute the sum of voice metrics */
	for (vm_sum = 0, i = LO_CHAN; i <= HI_CHAN; i++) vm_sum += vm_tbl[min(ch_snr[i], 89)];
	/* Compute the total noise estimate (tne) and total channel energy estimate (tce) */
	for (tne = 0.0f, i = LO_CHAN; i <= HI_CHAN; i++) tne += st->ch_noise[i];
	for (tce = 0.0f, i = LO_CHAN; i <= HI_CHAN; i++) tce += st->ch_enrg[i];
	/* Calculate log spectral deviation */
	for (i = LO_CHAN; i <= HI_CHAN; i++) ch_enrg_db[i] = 10.0f*(Float32)log10(st->ch_enrg[i]);

	if (st->vad_first == TRUE) { for (i = LO_CHAN; i <= HI_CHAN; i++) st->ch_enrg_long_db[i] = ch_enrg_db[i]; }

	for (ch_enrg_dev = 0.0f, i = LO_CHAN; i <= HI_CHAN; i++) ch_enrg_dev += (Float32)fabs(st->ch_enrg_long_db[i] - ch_enrg_db[i]);
	/* Calculate long term integration constant as a function of total channel energy (tce) *//* (i.e., high tce (-40 dB) -> slow integration (alpha = 0.99), low tce (-60 dB) -> fast integration (alpha = 0.50) */
	alpha = max(min(0.99f - (0.49f / TCE_RANGE)*(HIGH_TCE_DB - 10.0f*(Float32)log10(tce)), 0.9f), 0.5f);
	/* Calc long term log spectral energy */
	for (i = LO_CHAN; i <= HI_CHAN; i++) st->ch_enrg_long_db[i] = alpha * st->ch_enrg_long_db[i] + (1.0f - alpha)*ch_enrg_db[i];

	/* Set or reset the update flag */
	update_flag = FALSE;
	if (vm_sum <= UPDATE_THLD)
	{
		st->update_cnt = 0;
		update_flag = TRUE;
	}
	else if ((tce > NOISE_FLOOR) && (ch_enrg_dev < DEV_THLD))
	{
		st->update_cnt++;
		if (st->update_cnt >= UPDATE_CNT_THLD) update_flag = TRUE;
	}

	if (st->update_cnt == st->last_update_cnt) st->hyster_cnt++;
	else                                     st->hyster_cnt = 0;

	st->last_update_cnt = st->update_cnt;
	if (st->hyster_cnt > HYSTER_CNT_THLD) st->update_cnt = 0;

	/* Update the channel noise estimates */
	if (update_flag == TRUE) { for (i = LO_CHAN; i <= HI_CHAN; i++) st->ch_noise[i] = max(0.9f*st->ch_noise[i] + 0.1f*st->ch_enrg[i], MIN_CHAN_ENRG); }
	st->vad_first = FALSE;
	//alpha = (update_flag == TRUE) ? AA : 1.0f;


	if (update_flag == 1) for (i = 0; i < FRM_LEN; i++) out_sp[1][i] = 0;
	else         for (i = 0; i < FRM_LEN; i++) out_sp[1][i] = 4000;


	//if(vad_cnt == 6) for(i=0; i<FRM_LEN; i++) out_sp[2][i] = 4000;
	//else         for(i=0; i<FRM_LEN; i++) out_sp[2][i] = 0;

	if (st->frame_cnt >= 6)
	{
		if (update_flag == 1)
		{
			snr = 0.0f;
			st->frame_noise++;
			for (pn = 0.0f, k = 1; k < D; k++)
			{
				pn += re6[0][k] * re6[0][k] + im6[0][k] * im6[0][k];
			}

			if (st->frame_noise > 30)
			{
				if ((pn / st->pn <= 1.2f) && (st->pn / pn <= 1.2f))
				{
					st->pn = 0.9*st->pn + 0.1*pn;
				}
			}
			else
			{
				st->pn = 0.9*st->pn + 0.1*pn;
			}

			//st->pn = 0.9*st->pn + 0.1*pn;
			//st->pn = pn;
		}
		else
		{
			st->frame_voice++;
			for (ps = 0.0f, k = 1; k < D; k++)
			{
				ps += re6[0][k] * re6[0][k] + im6[0][k] * im6[0][k];
			}
			if (st->pn > 0.0f)
			{
				snr = 10.0f * (Float32)log10(ps / st->pn);
			}
			else
			{
				snr = 0.0f;
			}


			/*if(snr > 15.0f)
			{
				if(frame>=178)
				{
					i=i;
				}
				st->frame_success++;
				for(k=1; k<D; k++) { peak_sub[k] = 0; }
				for(peak=0.0f, j=n=0; n<AZ_NUM; n+=AZ_STEP)
				{
					for(acc2=0.0f, k=1; k<D; k++)
					{
						//鐩镐綅鍔犳潈
						for(re[k]=0, m=0; m<NUM; m++) re[k] += re2[m][k]*wr[n][m][k] + im2[m][k]*wi[n][m][k];
						for(im[k]=0, m=0; m<NUM; m++) im[k] += im2[m][k]*wr[n][m][k] - re2[m][k]*wi[n][m][k];
						//璁＄畻姝ゆ柟鍚戠殑淇″彿鍔熺巼
						power_sub[k] = re[k]*re[k]+im[k]*im[k];
						if(peak_sub[k]<power_sub[k])
						{
							peak_sub[k] = power_sub[k];
							num_sub[k] = n;
						}
						acc2 += power_sub[k]*mic0_spec[k];
					}
					if(peak<acc2)
					{
						peak=acc2;
						j=n;
					}
				}


				//瀛愬甫娉㈣揪鏂瑰悜绮炬悳绱?
				for(k=1; k<D; k++)
				{
					for(sub_azimuth=num_sub[k], i=num_sub[k]-AZ_STEP/2; i<num_sub[k]+AZ_STEP/2; i++)
					{
						n=i;
						if(n<0) n += AZ_NUM;
						if(n>AZ_NUM) n -= AZ_NUM;
						for(re[k]=0, m=0; m<NUM; m++) re[k] += re2[m][k]*wr[n][m][k] + im2[m][k]*wi[n][m][k];
						for(im[k]=0, m=0; m<NUM; m++) im[k] += im2[m][k]*wr[n][m][k] - re2[m][k]*wi[n][m][k];

						acc2 = re[k]*re[k]+im[k]*im[k];
						if(peak_sub[k]<acc2) {peak_sub[k]=acc2; sub_azimuth=n;}
					}
					st->sub_az[k] = sub_azimuth*2.0f*(Float32)PI/AZ_NUM;
				}

				//鍏ㄥ甫娉㈣揪鏂瑰悜绮炬悳绱?
				for(tot_az=j, i=j-AZ_STEP/2; i<=j+AZ_STEP/2; i++)
				{
					//瑙掑害鏄犲皠鍒癧0, 2*PI)鍖洪棿
					n=i;
					if(n<0) n += AZ_NUM;
					if(n>AZ_NUM) n -= AZ_NUM;

					for(acc2=0.0f, k=1; k<D; k++)
					{
						//鐩镐綅鍔犳潈
						for(re[k]=0, m=0; m<NUM; m++) re[k] += re2[m][k]*wr[n][m][k] + im2[m][k]*wi[n][m][k];
						for(im[k]=0, m=0; m<NUM; m++) im[k] += im2[m][k]*wr[n][m][k] - re2[m][k]*wi[n][m][k];
						//璁＄畻姝ゆ柟鍚戠殑淇″彿鍔熺巼
						acc2 += (re[k]*re[k]+im[k]*im[k])*mic0_spec[k];
					}
					if(peak<acc2)
					{
						peak=acc2;
						tot_az=n;
					}
				}

				for(j=0; j<MED_NUM-1; j++) st->audio_az[j]=st->audio_az[j+1];
				st->audio_az[MED_NUM-1] = tot_az*2.0f*(Float32)PI/AZ_NUM;
				if(fabs(get_position(st)-PI/2)<=PI/18)
				//if(fabs(st->audio_az[MED_NUM-1]-PI/2)<=PI/18)
				//if((fabs(get_position(st)-2*PI)<=PI/18)||(fabs(get_position(st))<=PI/18))
				//if((fabs(st->audio_az[MED_NUM-1]-2*PI)<=PI/36)||(fabs(st->audio_az[MED_NUM-1])<=PI/36))
				{
					st->frame_locate++;
				}
			}*/
		}
	}

	else
	{
		snr = 0.0f;
	}


	finish2 = clock();
	st->srp_time += (double)(finish2 - start2) / CLOCKS_PER_SEC;


	//娉㈣揪鏂瑰悜绮楁悳绱?
	for (k = 1; k < D; k++) { peak_sub[k] = 0; }
	for (peak = 0.0f, j = n = 0; n < AZ_NUM; n += AZ_STEP)
	{
		for (acc2 = 0.0f, k = 1; k < D; k++)
		{
			//鐩镐綅鍔犳潈    			
			for (re[k] = 0, m = 0; m < NUM; m++) re[k] += re2[m][k] * wr[n][m][k] + im2[m][k] * wi[n][m][k];
			for (im[k] = 0, m = 0; m < NUM; m++) im[k] += im2[m][k] * wr[n][m][k] - re2[m][k] * wi[n][m][k];
			//璁＄畻姝ゆ柟鍚戠殑淇″彿鍔熺巼
			power_sub[k] = re[k] * re[k] + im[k] * im[k];
			if (peak_sub[k] < power_sub[k])
			{
				peak_sub[k] = power_sub[k];
				num_sub[k] = n;
			}
			acc2 += power_sub[k] * mic0_spec[k];
		}
		if (peak < acc2)
		{
			peak = acc2;
			j = n;
		}
	}


	//瀛愬甫娉㈣揪鏂瑰悜绮炬悳绱?
	for (k = 1; k < D; k++)
	{
		for (sub_azimuth = num_sub[k], i = num_sub[k] - AZ_STEP / 2; i < num_sub[k] + AZ_STEP / 2; i++)
		{
			n = i;
			if (n < 0) n += AZ_NUM;
			if (n > AZ_NUM) n -= AZ_NUM;
			for (re[k] = 0, m = 0; m < NUM; m++) re[k] += re2[m][k] * wr[n][m][k] + im2[m][k] * wi[n][m][k];
			for (im[k] = 0, m = 0; m < NUM; m++) im[k] += im2[m][k] * wr[n][m][k] - re2[m][k] * wi[n][m][k];

			acc2 = re[k] * re[k] + im[k] * im[k];
			if (peak_sub[k] < acc2) { peak_sub[k] = acc2; sub_azimuth = n; }
		}
		st->sub_az[k] = sub_azimuth * 2.0f*(Float32)PI / AZ_NUM;
	}

	//鍏ㄥ甫娉㈣揪鏂瑰悜绮炬悳绱?
	for (tot_az = j, i = j - AZ_STEP / 2; i <= j + AZ_STEP / 2; i++)
	{
		//瑙掑害鏄犲皠鍒癧0, 2*PI)鍖洪棿
		n = i;
		if (n < 0) n += AZ_NUM;
		if (n > AZ_NUM) n -= AZ_NUM;

		for (acc2 = 0.0f, k = 1; k < D; k++)
		{
			//鐩镐綅鍔犳潈    			
			for (re[k] = 0, m = 0; m < NUM; m++) re[k] += re2[m][k] * wr[n][m][k] + im2[m][k] * wi[n][m][k];
			for (im[k] = 0, m = 0; m < NUM; m++) im[k] += im2[m][k] * wr[n][m][k] - re2[m][k] * wi[n][m][k];
			//璁＄畻姝ゆ柟鍚戠殑淇″彿鍔熺巼
			acc2 += (re[k] * re[k] + im[k] * im[k])*mic0_spec[k];
		}
		if (peak < acc2)
		{
			peak = acc2;
			tot_az = n;
		}
	}

	if (snr > 15.0f)
	{
		st->frame_success++;
		for (j = 0; j < MED_NUM - 1; j++) st->audio_az[j] = st->audio_az[j + 1];
		st->audio_az[MED_NUM - 1] = tot_az * 2.0f*(Float32)PI / AZ_NUM;
		if (fabs(get_position(st) - PI / 2) <= PI / 18)
			//if(fabs(st->audio_az[MED_NUM-1]-PI/2)<=PI/18)
			//if((fabs(get_position(st)-2*PI)<=PI/18)||(fabs(get_position(st))<=PI/18))
			//if((fabs(st->audio_az[MED_NUM-1]-2*PI)<=PI/36)||(fabs(st->audio_az[MED_NUM-1])<=PI/36))
		{
			st->frame_locate++;
		}
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////
	//  GSC娉㈡潫鎴愬舰+鍑犱綍鐩插垎绂?
	//////////////////////////////////////////////////////////////////////////////////////////////////
	if ((0 <= st->video_az) && (st->video_az < 2 * PI) && (-PI / 2.0 < st->video_pit) && (st->video_pit < PI / 2))
	{
		peak = st->video_az;
	}
	else
	{
		peak = get_position(st);
	}

	az_num = (Word32)floor(peak*AZ_NUM / 2.0f / PI + 0.5);
	if (az_num < 0) az_num += AZ_NUM;
	if (az_num > AZ_NUM) az_num -= AZ_NUM;

	for (i = 0; i < NUM; i++)
	{
		tao[i] = (Float32)(PI / D)*(Float32)(FS*DIAMETER / 2 / SPEED)*(Float32)cos(peak - 2 * PI*i / NUM);
		for (k = 0; k < D; k++)
		{
			phrase_re[i][k] = (Float32)cos(k*(tao[i] - tao[0]));
			phrase_im[i][k] = -(Float32)sin(k*(tao[i] - tao[0]));
		}
	}


	start3 = clock();
	//MVDR
	if (frame >= 2150)
	{
		i = i;
	}
	//alpha = (update_flag == TRUE) ? AA : 1.0f;
	alpha = vad_prob > 0.3 ? 1.0f : AA;
	if (st->bf_first != 1)
	{
		for (k = 1; k < D; k++)
		{
			for (m = 0; m < NUM; m++)
			{
				for (n = 0; n <= m; n++)
				{
					//tr = re6[m][k] * re6[n][k] + im6[m][k] * im6[n][k];
					//ti = im6[m][k] * re6[n][k] - re6[m][k] * im6[n][k];
					tr = st->bf_re[m][k] * st->bf_re[n][k] + st->bf_im[m][k] * st->bf_im[n][k];
					ti = st->bf_im[m][k] * st->bf_re[n][k] - st->bf_re[m][k] * st->bf_im[n][k];
					st->qxx_re[k][m][n] = 0.2f*tr + 0.8f*st->qxx_re[k][m][n];
					st->qxx_im[k][m][n] = 0.2f*ti + 0.8f*st->qxx_im[k][m][n];
					st->qnn_re[k][m][n] = (1 - alpha)*tr + alpha * st->qnn_re[k][m][n];
					st->qnn_im[k][m][n] = (1 - alpha)*ti + alpha * st->qnn_im[k][m][n];
					inv1_re[m][n] = st->qxx_re[k][m][n] + (BB - 1)*st->qnn_re[k][m][n];
					inv1_im[m][n] = st->qxx_im[k][m][n] + (BB - 1)*st->qnn_im[k][m][n];
					//inv1_re[m][n] = st->qxx_re[k][m][n] ;
					//inv1_im[m][n] = st->qxx_im[k][m][n] ;
				}
			}
			for (m = 0; m < NUM; m++)
			{
				for (n = m + 1; n < NUM; n++)
				{
					st->qxx_re[k][m][n] = st->qxx_re[k][n][m];
					st->qxx_im[k][m][n] = -st->qxx_im[k][n][m];
					st->qnn_re[k][m][n] = st->qnn_re[k][n][m];
					st->qnn_im[k][m][n] = -st->qnn_im[k][n][m];
					inv1_re[m][n] = inv1_re[n][m];
					inv1_im[m][n] = -inv1_im[n][m];
				}
			}
			inverse(inv1_re, inv1_im, inv_re, inv_im, NUM);

			for (m = 0; m < NUM; m++)
			{
				matrix_re[m] = st->qxx_re[k][m][0] - st->qnn_re[k][m][0];
				matrix_im[m] = st->qxx_im[k][m][0] - st->qnn_im[k][m][0];
			}

			//for(m=0; m<NUM; m++)
			//{
			//	acc3 = sqrt(matrix_re[m]*matrix_re[m] + matrix_im[m]*matrix_im[m]);
			//	if(acc3>0.0f)
			//	{
			//		matrix_re[m] /= acc3;
			//		matrix_im[m] /= acc3;
			//	}
			//	else
			//	{
			//		matrix_re[m] = 0.0f;
			//		matrix_im[m] = 0.0f;
			//	}
			//}
			e_r[0][k] = e_i[0][k] = 0.0f;
			for (m = 0; m < NUM; m++)
			{
				tr = ti = 0.0f;
				for (n = 0; n < NUM; n++)
				{
					tr += (Float32)(inv_re[m][n] * matrix_re[n] - inv_im[m][n] * matrix_im[n]);
					ti += (Float32)(inv_re[m][n] * matrix_im[n] + inv_im[m][n] * matrix_re[n]);
					//tr += (Float32)(inv_re[m][n] * re2[n][k] - inv_im[m][n] * im2[n][k]);
					//ti += (Float32)(inv_re[m][n] * im2[n][k] + inv_im[m][n] * re2[n][k]);
				}
				//e_r[0][k] += tr * re6[m][k] + ti * im6[m][k];
				//e_i[0][k] += tr * im6[m][k] - ti * re6[m][k];
				e_r[0][k] += tr * st->bf_re[m][k] + ti * st->bf_im[m][k];
				e_i[0][k] += tr * st->bf_im[m][k] - ti * st->bf_re[m][k];
			}
		}

		/*alpha = (update_flag == TRUE) ? 0.7f : 1.0f;
		for(k=1; k<D; k++)
		{
			for(m=0; m<NUM; m++)
			{
				for(n=0; n<=m; n++)
				{
					tr = re6[m][k]*re6[n][k] + im6[m][k]*im6[n][k];
					ti = im6[m][k]*re6[n][k] - re6[m][k]*im6[n][k];
					//st->qxx_re[k][m][n] = 0.2f*tr + 0.8f*st->qxx_re[k][m][n];
					//st->qxx_im[k][m][n] = 0.2f*ti + 0.8f*st->qxx_im[k][m][n];
					st->qnn_re[k][m][n] = (1-alpha)*tr + alpha*st->qnn_re[k][m][n];
					st->qnn_im[k][m][n] = (1-alpha)*ti + alpha*st->qnn_im[k][m][n];
					inv1_re[m][n] = st->qnn_re[k][m][n];
					inv1_im[m][n] = st->qnn_im[k][m][n];
					//inv1_re[m][n] = st->qxx_re[k][m][n] ;
					//inv1_im[m][n] = st->qxx_im[k][m][n] ;
				}
			}
			for(m=0; m<NUM; m++)
			{
				for(n=m+1; n<NUM; n++)
				{
					//st->qxx_re[k][m][n] = st->qxx_re[k][n][m];
					//st->qxx_im[k][m][n] = -st->qxx_im[k][n][m];
					//st->qnn_re[k][m][n] = st->qnn_re[k][n][m];
					//st->qnn_im[k][m][n] = -st->qnn_im[k][n][m];
					inv1_re[m][n] = inv1_re[n][m];
					inv1_im[m][n] = -inv1_im[n][m];
				}
			}
			inverse(inv1_re, inv1_im, inv_re, inv_im, NUM);

			//for(m=0; m<NUM; m++)
			//{
			//    matrix_re[m] = st->qxx_re[k][m][0] - st->qnn_re[k][m][0];
			//    matrix_im[m] = st->qxx_im[k][m][0] - st->qnn_im[k][m][0];
			//}

			if(k==1)
			{
				k=k;
			}
			for(den_re=den_im=0.0f,m=0; m<NUM; m++)
			{
				weight_re[m] = weight_im[m] = 0.0f;
				for(n=0; n<NUM; n++)
				{
					//weight_re[m] += inv_re[m][n]*wr[az_num][n][k] + inv_im[m][n]*wi[az_num][n][k];
					//weight_im[m] += inv_im[m][n]*wr[az_num][n][k] - inv_re[m][n]*wi[az_num][n][k];
					weight_re[m] += inv_re[m][n]*phrase_re[n][k] + inv_im[m][n]*phrase_im[n][k];
					weight_im[m] += inv_im[m][n]*phrase_re[n][k] - inv_re[m][n]*phrase_im[n][k];
				}
				den_re += phrase_re[m][k]*weight_re[m] - phrase_im[m][k]*weight_im[m];
				den_im += phrase_re[m][k]*weight_im[m] + phrase_im[m][k]*weight_re[m];
			}

			acc2 = den_re*den_re + den_im*den_im;

			for(e_r[0][k]=e_i[0][k]=0.0f,m=0; m<NUM; m++)
			{
				weight_re[m] = (weight_re[m]*den_re+weight_im[m]*den_im)/acc2;
				weight_im[m] = (weight_im[m]*den_re-weight_re[m]*den_im)/acc2;

				e_r[0][k] += weight_re[m]*re6[m][k] + weight_im[m]*im6[m][k];
				e_i[0][k] += weight_re[m]*im6[m][k] - weight_im[m]*re6[m][k];
			}

		}*/

		//finish3 = clock();
		//st->mvdr_time += (double)(finish3 - start3) / CLOCKS_PER_SEC;

		for (i = ORD2 * 2 - 1; i >= 2 * D; i--) st->syn[2][i] = st->syn[2][i - 2 * D];
		for (i = 0; i < 2 * D; i++)
		{
			for (acc2 = 0.0f, k = 1; k < D; k++) acc2 += e_r[0][k] * cos_tab0[k][i] - e_i[0][k] * sin_tab0[k][i];
			st->syn[2][i] = acc2;
		}
		for (i = 0; i < FRM_LEN; i++)
		{
			for (acc2 = 0.0f, j = 0; j < ORD2 / D; j++) acc2 += prototype_filter[j*D + i] * st->syn[2][j * 2 * D + (j & 1)*D + i];
			out_sp[2][i] = (short)(acc2*D * 20);
		}
	}
	else
	{
		for (i = 0; i < FRM_LEN; i++)
		{
			out_sp[2][i] = 0;
		}
	}

	finish3 = clock();
	st->mvdr_time += (double)(finish3 - start3) / CLOCKS_PER_SEC;

	st->bf_first = 0;

	for (m = 0; m < NUM; m++)
	{
		for (k = 1; k < D; k++)
		{
			st->bf_re[m][k] = re6[m][k];
			st->bf_im[m][k] = im6[m][k];
		}
	}



	/*for(k=1; k<D; k++)
	{
		for(j=0; j<ORD6-1; j++) st->fbf_re[k][j] = st->fbf_re[k][j+1];
		for(j=0; j<ORD6-1; j++) st->fbf_im[k][j] = st->fbf_im[k][j+1];
	}

	for(m=0; m<NUM; m++)
	{
		for(k=0; k<D; k++)
		{
			for(j=0; j<ORD6-1; j++) st->bmIN_re[m][k][j] = st->bmIN_re[m][k][j+1];
			for(j=0; j<ORD6-1; j++) st->bmIN_im[m][k][j] = st->bmIN_im[m][k][j+1];
		}
	}

	//璁＄畻淇″彿瀵归綈
	for(k=1; k<D; k++)
	{
		for(m=0; m<NUM; m++)
		{
			re3[m][k] = re6[m][k]*wr[n][m][k] + im6[m][k]*wi[n][m][k];
			st->bmIN_re[m][k][ORD6-1] = re6[m][k];

			im3[m][k] = im6[m][k]*wr[n][m][k] - re6[m][k]*wi[n][m][k];
			st->bmIN_im[m][k][ORD6-1] = im6[m][k];
		}
		//for(m=0; m<NUM; m++) st->bmIN_im[m][k][ORD6-1] = im3[m][k] = im6[m][k]*wr[n][m][k] - re6[m][k]*wi[n][m][k];
		//for(m=0; m<NUM; m++) re3[m][k] = e_r[m][k]*wr[n][m][k] + e_i[m][k]*wi[n][m][k];
		//for(m=0; m<NUM; m++) im3[m][k] = e_i[m][k]*wr[n][m][k] - e_r[m][k]*wi[n][m][k];
	}

	//浼犵粺鍔犳潈姹傚拰娉㈡潫褰㈡垚
	for(k=1; k<D; k++)
	{
		for(acc2=0.0f, m=0; m<NUM; m++) acc2 += re3[m][k];
		for(acc3=0.0f, m=0; m<NUM; m++) acc3 += im3[m][k];
		st->fbf_re[k][ORD6-1] = re[k] = acc2/NUM;
		st->fbf_im[k][ORD6-1] = im[k] = acc3/NUM;
	}

	for(k=1; k<D; k++)
	{
		for(acc2=0.0f, j=0; j<ORD6; j++) acc2 += st->fbf_re[k][j]*st->fbf_re[k][j] + st->fbf_im[k][j]*st->fbf_im[k][j];
		energy[k] = acc2;
	}

	//if(vad_cnt == NUM) acc3 = 1.0f;
	//else               acc3 = 0.001f;

	for(m=0; m<NUM; m++)
	{
		for(k=1; k<D; k++)
		{
			for(acc2=0.0f, j=0; j<ORD6; j++) acc2 += st->fbf_re[k][j]*st->sr[m][k][j] - st->fbf_im[k][j]*st->si[m][k][j];
			for(acc3=0.0f, j=0; j<ORD6; j++) acc3 += st->fbf_im[k][j]*st->sr[m][k][j] + st->fbf_re[k][j]*st->si[m][k][j];

			//re4[m][k] = re3[m][k] - acc2;
			//im4[m][k] = im3[m][k] - acc3;

			re4[m][k] = st->bmIN_re[m][k][ORD6-1] - acc2;
			im4[m][k] = st->bmIN_im[m][k][ORD6-1] - acc3;
		}

		//if(vad_cnt == NUM) acc3 = 0.5f;
		//if(snr >= 20.0) acc3 = 0.8f;
		//else               acc3 = 0.00001f;
		//for(k=1; k<D; k++)
		//{
		//	if(energy[k]>4.0f)
		//	{
		//		for(j=0; j<ORD6; j++) st->sr[m][k][j] += acc3/energy[k] * (st->fbf_re[k][j]*re4[m][k] + st->fbf_im[k][j]*im4[m][k]);
		//		for(j=0; j<ORD6; j++) st->si[m][k][j] += acc3/energy[k] * (st->fbf_re[k][j]*im4[m][k] - st->fbf_im[k][j]*re4[m][k]);
		//	}
		//}
		if(snr >= 20.0f)
		{
			acc3=0.8f;
			for(k=1; k<D; k++)
			{
				for(j=0; j<ORD6; j++) st->sr[m][k][j] += acc3/energy[k] * (st->fbf_re[k][j]*re4[m][k] + st->fbf_im[k][j]*im4[m][k]);
				for(j=0; j<ORD6; j++) st->si[m][k][j] += acc3/energy[k] * (st->fbf_re[k][j]*im4[m][k] - st->fbf_im[k][j]*re4[m][k]);
			}
		}
	}



	//鍔犳潈姹傚拰娉㈡潫鎴愬舰杈撳嚭
	for (i = ORD2 * 2 - 1; i >= 2 * D; i--) st->syn[1][i] = st->syn[1][i - 2 * D];
	for (i = 0; i < 2 * D; i++)
	{
		for (acc2 = 0.0f, k = 1; k < D; k++) acc2 += re[k] * cos_tab0[k][i] - im[k] * sin_tab0[k][i];
		st->syn[1][i] = acc2;
	}
	for (i = 0; i < FRM_LEN; i++)
	{
		for (acc2 = 0.0f, j = 0; j < ORD2 / D; j++) acc2 += prototype_filter[j*D + i] * st->syn[1][j * 2 * D + (j & 1)*D + i];
		out_sp[1][i] = (short)(acc2*D * 32);
	}

	//閫氳繃Block鐭╅樀
	//for(k=1; k<D; k++)
	//{
	//	for(m=0; m<NUM-1; m++) re4[m][k] = re3[m+1][k] - re3[m][k];
	//	for(m=0; m<NUM-1; m++) im4[m][k] = im3[m+1][k] - im3[m][k];
	//}

	//浼拌瀛愬甫鑳介噺
	for(k=1; k<D; k++)
	{
		for(acc2=0.0f, m=0; m<NUM; m++) acc2 += re4[m][k]*re4[m][k] + im4[m][k]*im4[m][k];
		ener[k]=acc2;
	}

	//multichannel ANC
	for(k=1; k<D; k++)
	{
		//浼拌瀛愬甫妯℃嫙淇″彿
		for(acc2=0.0f, m=0; m<NUM; m++) acc2 += re4[m][k]*st->gr[m][k] - im4[m][k]*st->gi[m][k];
		for(acc3=0.0f, m=0; m<NUM; m++) acc3 += im4[m][k]*st->gr[m][k] + re4[m][k]*st->gi[m][k];
		//浼拌璇樊瀛愬甫淇″彿
		e_r[0][k] = re[k] - acc2;
		e_i[0][k] = im[k] - acc3;
	}

	//GSC娉㈡潫鎴愬舰杈撳嚭
	for (i = ORD2 * 2 - 1; i >= 2 * D; i--) st->syn[2][i] = st->syn[2][i - 2 * D];
	for (i = 0; i < 2 * D; i++)
	{
		for (acc2 = 0.0f, k = 1; k < D; k++) acc2 += e_r[0][k] * cos_tab0[k][i] - e_i[0][k] * sin_tab0[k][i];
		st->syn[2][i] = acc2;
	}
	for (i = 0; i < FRM_LEN; i++)
	{
		for (acc2 = 0.0f, j = 0; j < ORD2 / D; j++) acc2 += prototype_filter[j*D + i] * st->syn[2][j * 2 * D + (j & 1)*D + i];
		out_sp[2][i] = (short)(acc2*D * 32);
	}

	//鑷�傚簲鎶垫秷
	//if(() acc3=0.5f;
	//else             acc3=0.0001f;
	//for(k=1; k<D; k++)
	//{
	//   if(ener[k]>4.0f)
	//	{
	//	    acc2 = acc3/ener[k];
	 //	    for(m=0; m<NUM; m++) st->gr[m][k] += acc2*(re4[m][k]*e_r[0][k] + im4[m][k]*e_i[0][k]);
	//		for(m=0; m<NUM; m++) st->gi[m][k] += acc2*(re4[m][k]*e_i[0][k] - im4[m][k]*e_r[0][k]);
	//	}
	//}
	if(snr < 1.0f)
	{
		acc3=0.8f;
		for(k=1; k<D; k++)
		{
			if(ener[k]>0)
			{
				for(m=0; m<NUM; m++) st->gr[m][k] += acc3/ener[k]*(re4[m][k]*e_r[0][k] + im4[m][k]*e_i[0][k]);
				for(m=0; m<NUM; m++) st->gi[m][k] += acc3/ener[k]*(re4[m][k]*e_i[0][k] - im4[m][k]*e_r[0][k]);
			}
		}
	}*/

	//几何盲分离
	for (k = 1; k < D; k++)
	{
		//att[k] = (Float32)pow(0.5 + 0.5*cos(st->sub_az[k] - peak), 6);
		att[k] = (Float32)pow(0.5 + 0.5*cos(tot_az*2.0f*(Float32)PI / AZ_NUM - (PI / 2 - PI / 36)), 6);
		e_r[0][k] *= att[k];
		e_i[0][k] *= att[k];
	}

	//瀛愬甫鍚堟垚
	for (i = ORD2 * 2 - 1; i >= 2 * D; i--) st->syn[3][i] = st->syn[3][i - 2 * D];

	for (i = 1; i < D; i++)
	{
		re5[2 * D - i] = re5[i] = e_r[0][i];
		im5[2 * D - i] = -e_i[0][i];
		im5[i] = e_i[0][i];
	}
	re5[D] = re5[0] = im5[D] = im5[0] = 0;

	fft_320(re5, im5);

	for (i = 0; i < 2 * D; i++)
	{
		st->syn[3][i] = re5[i];
	}

	//鍑犱綍鐩插垎绂昏緭鍑?
	for (i = 0; i < FRM_LEN; i++)
	{
		for (acc2 = 0.0f, j = 0; j < ORD2 / D; j++) acc2 += prototype_filter[j*D + i] * st->syn[3][j * 2 * D + (j & 1)*D + i];
		//out_sp[i]=(short)(acc2*D*32);
		out_sp[3][i] = (short)(acc2*D * 10);
		//temp1[i]=(short)(acc2*D*32);
	}

	//ANC
	NS_run(&(st->my_str_anc), out_sp[3], out_sp[4]);
	//AGC
	AGC_Process(&(st->my_str_agc), out_sp[4], NULL, FRM_LEN, out_sp[5], NULL, inMicLevel, &outMicLevel, 0, (uint8_t*)(&saturationWarning));

	st->init_farme_cnt = (st->init_farme_cnt < INIT_LEN) ? (st->init_farme_cnt + 1) : INIT_LEN;
}

void main(void)
{
	Word16 in_sp[NUM][FRM_LEN2], mic_sp[NUM][FRM_LEN], spk[FRM_LEN2], spk_sp[FRM_LEN], out_sp[6][FRM_LEN], out_all_around[FRM_LEN], out_vad[FRM_LEN];
	Word32 i, flag, locate_mid_num, valid_mid_num, locate_moment_num, valid_moment_num;
	Float32 error_mid, error_moment;
	AEC_SRP_ST my_str;
	DOWNSAMPLE_STR down_str[NUM + 1];

	clock_t start, finish, start_file, finish_file;
	double total_time;

	FILE *fq[NUM], *fr, *fs[6], *fp, *fn, *fm, *fvad;

	for (i = 0; i < NUM + 1; i++)
	{
		down_sample_init(&down_str[i]);
	}

	aec_srp_gsc_init(&my_str, 0.0f);

	fp = fopen("azimuth_mid2.txt", "w");
	fn = fopen("azimuth_moment.txt", "w");

	fq[0] = fopen("0.pcm", "rb");
	fq[1] = fopen("1.pcm", "rb");
	fq[2] = fopen("2.pcm", "rb");
	fq[3] = fopen("3.pcm", "rb");
	fq[4] = fopen("4.pcm", "rb");
	fq[5] = fopen("5.pcm", "rb");

	for (i = 0; i < NUM; i++)
	{
		if (fq[i] == NULL) { printf("Cann't open the No.%ld microphone file!!!\n", i);  exit(0); }
	}

	fr = fopen("spk.pcm", "rb");
	if (fr == NULL) { printf("Cann't open spk_sp.pcm file!!!\n");  exit(0); }

	fs[0] = fopen("out0.pcm", "wb");
	fs[1] = fopen("out1.pcm", "wb");
	fs[2] = fopen("out2.pcm", "wb");
	fs[3] = fopen("out3.pcm", "wb");
	fs[4] = fopen("out4.pcm", "wb");
	fs[5] = fopen("out5.pcm", "wb");
	for (i = 0; i < NUM; i++)
	{
		if (fs[i] == NULL) { printf("Cann't open the No.%ld out file!!!\n", i);  exit(0); }
	}

	fm = fopen("out_all_around.pcm", "wb");
	if (fm == NULL) { printf("Cann't open the out_all_around.pcm file!!!\n"); exit(0); }

	fvad = fopen("out_vad.pcm", "wb");
	if (fvad == NULL) { printf("Cann't open the out_vad.pcm file!!!\n"); exit(0); }

	start_file = clock();
	for (total_time = 0, error_mid = error_moment = 0, locate_mid_num = valid_mid_num = locate_moment_num = valid_moment_num = 0, frame = 0; ; frame++)
	{
		for (flag = i = 0; i < NUM; i++)
		{
			//if( fread(in_sp[i], sizeof(short), FRM_LEN2, fq[i]) != FRM_LEN2 )  //48k
			if (fread(mic_sp[i], sizeof(short), FRM_LEN, fq[i]) != FRM_LEN)  //16k
			{
				flag = 1; break;
			}

			//down_sample_run(&down_str[i], in_sp[i], mic_sp[i]);  //48k->16k
		}
		if (flag == 1) break;

		if (fread(spk_sp, sizeof(short), FRM_LEN, fr) != FRM_LEN) break;//16k
		//if( fread(spk, sizeof(short), FRM_LEN2, fr) != FRM_LEN2 ) break;//48k
		//down_sample_run(&down_str[NUM], spk, spk_sp);//48k->16k

		if ((frame & 2047L) == 0) printf("frame=%ld\n", frame);
		if (frame * 160L >= 350752)
		{
			i = i;
			//set_position(&my_str, -1, -1, HUMAN_FACE);
			//get_position(&my_str);
		}

		start = clock();
		aec_srp_gsc(&my_str, mic_sp, spk_sp, out_sp, out_all_around, out_vad);
		finish = clock();
		total_time += (double)(finish - start) / CLOCKS_PER_SEC;


		//娴嬭宸?
		if (get_position(&my_str) > 0)
		{
			locate_mid_num++;
			if ((fabs(get_position(&my_str) - (12 * PI / 6)) <= PI / 18) || (fabs(get_position(&my_str)) <= PI / 18))
				//if(fabs(get_position(&my_str)-(1*PI/6+PI/18))<=PI/36)
			{
				valid_mid_num++;
				//error_mid += fabs(get_position(&my_str)-(1*PI/6+PI/18));
				if (fabs(get_position(&my_str)) <= PI / 36)
				{
					error_mid += fabs(get_position(&my_str));
				}
				else
				{
					error_mid += fabs(get_position(&my_str) - 12 * PI / 6 + PI / 36);
				}
			}
		}

		if (my_str.audio_az[MED_NUM - 1] > 0)
		{
			locate_moment_num++;
			if ((fabs(my_str.audio_az[MED_NUM - 1] - (12 * PI / 6)) <= PI / 18) || (fabs(my_str.audio_az[MED_NUM - 1]) <= PI / 18))
				//if(fabs(my_str.audio_az[MED_NUM-1]-(1*PI/6+PI/18))<=PI/36)
			{
				valid_moment_num++;
				error_moment += fabs(my_str.audio_az[MED_NUM - 1] - 1 * PI / 6);
			}
		}

		fprintf(fp, "%f\n", get_position(&my_str) / PI * 180);
		fprintf(fn, "%f\n", my_str.audio_az[MED_NUM - 1] / PI * 180);

		for (i = 0; i < 6; i++)
		{
			if (i != 2)
			{
				fwrite(out_sp[i], sizeof(short), FRM_LEN, fs[i]);
			}
		}
		if (frame > 0)
		{
			fwrite(out_sp[2], sizeof(short), FRM_LEN, fs[2]);
			fwrite(out_all_around, sizeof(short), FRM_LEN, fm);
			fwrite(out_vad, sizeof(short), FRM_LEN, fvad);
		}
	}

	finish_file = clock();

	printf("中值滤波的定位成功率为%f\n", ((Float32)valid_mid_num) / ((Float32)locate_mid_num));
	printf("中值滤波的平均绝对误差为%f\n", error_mid / valid_mid_num);
	printf("中值滤波的平均相对误差为%f\n", error_mid / valid_mid_num / (12 * PI / 6));

	printf("瞬时值的定位成功率为%f\n", ((Float32)valid_moment_num) / ((Float32)locate_moment_num));
	printf("瞬时值的平均绝对误差为%f\n", error_moment / valid_moment_num);
	printf("瞬时值的平均相对误差为%f\n", error_moment / valid_moment_num / (1 * PI / 6));

	printf("噪声帧数为%d\n", my_str.frame_noise);
	printf("语音帧数为%d\n", my_str.frame_voice);
	printf("大于20dB的语音帧数为%d\n", my_str.frame_success);
	printf("定位成功的帧数为%d\n", my_str.frame_locate);
	printf("真实定位成功率为%f\n", (double)(my_str.frame_locate) / my_str.frame_success);

	aec_srp_gsc_free(&my_str);

	for (i = 0; i < 6; i++) fclose(fs[i]);
	fclose(fr);

	for (i = 0; i < NUM; i++) fclose(fq[i]);
	fclose(fp);
	fclose(fn);
	fclose(fm);
	fclose(fvad);

	printf("程序运行时间为%f秒\n", (double)(finish_file - start_file) / CLOCKS_PER_SEC);
	printf("程序运行时间为%f秒\n", total_time);
	printf("AEC运行时间为%f秒\n", my_str.aec_time);
	printf("SRP运行时间为%f秒\n", my_str.srp_time);
	printf("MVDR运行时间为%f秒\n", my_str.mvdr_time);

	system("pause");
}
