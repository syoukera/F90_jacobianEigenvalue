#include <stdio.h>
#include <stdlib.h>
#include "jacob.h"

#define N_SPECIES 10

// LAPACK関数のプロトタイプ宣言
extern void dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda, 
                   double* wr, double* wi, double* vl, int* ldvl, 
                   double* vr, int* ldvr, double* work, int* lwork, int* info);

int main() {
    const double t = 0.0;
    const double pres = 101325.0;
    
    double y[N_SPECIES] = {1000, 0., 0.04030952, 0., 0.95969048, 0., 0., 0., 0., 0.};
    double jac[N_SPECIES * N_SPECIES] = {0.0};
    
    eval_jacob(t, pres, y, jac);  // ヤコビアンを計算

    // 固有値計算用の変数
    double wr[N_SPECIES], wi[N_SPECIES];  // 実部と虚部
    double vl[N_SPECIES * N_SPECIES];  // 左固有ベクトル（不要ならNULLでも可）
    double vr[N_SPECIES * N_SPECIES];  // 右固有ベクトル（不要ならNULLでも可）
    int lwork = 4 * N_SPECIES;  // 作業配列のサイズ
    double work[4 * N_SPECIES];  // 作業配列
    int info;
    
    char jobvl = 'V';  // 左固有ベクトルは計算しない
    char jobvr = 'V';  // 右固有ベクトルは計算しない
    int n = N_SPECIES;
    int lda = N_SPECIES;
    int ldvl = N_SPECIES;
    int ldvr = N_SPECIES;

    double a_exp[N_SPECIES];
    double b_exp[N_SPECIES];

    double EP[N_SPECIES];
    double EI[N_SPECIES];

    // LAPACK の dgeev を呼び出し（ヤコビアンの固有値を計算）
    dgeev_(&jobvl, &jobvr, &n, jac, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
    
    if (info != 0) {
        printf("Error: dgeev failed with INFO = %d\n", info);
        return 1;
    }

    // 固有値・固有ベクトルの表示
    int i = 2;
    printf("Eigenvalue %d: %f\n", i, wr[i]);
    
    double EP_sum = 0.0;

    printf("EP`:\n");
    // ignore complex value
    for (int j = 0; j < N_SPECIES; j++) {
        a_exp[j] = vr[j + i*N_SPECIES]; // right
        b_exp[j] = vl[j + i*N_SPECIES]; // left
        EP[j] = a_exp[j]*b_exp[j];
        EP_sum += fabs(EP[j]);
        printf("%e\n", EP[j]);  
    }

    printf("EI:\n");
    // ignore complex value
    for (int j = 0; j < N_SPECIES; j++) {
        EI[j] = fabs(EP[j])/EP_sum;
        printf("%e\n", EI[j]);  
    }

    return 0;
}
