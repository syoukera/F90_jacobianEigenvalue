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
    
    char jobvl = 'N';  // 左固有ベクトルは計算しない
    char jobvr = 'N';  // 右固有ベクトルは計算しない
    int n = N_SPECIES;
    int lda = N_SPECIES;
    int ldvl = N_SPECIES;
    int ldvr = N_SPECIES;

    // LAPACK の dgeev を呼び出し（ヤコビアンの固有値を計算）
    dgeev_(&jobvl, &jobvr, &n, jac, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
    
    if (info == 0) {
        printf("固有値（実部 虚部）:\n");
        for (int i = 0; i < N_SPECIES; i++) {
            printf("%f + %fi\n", wr[i], wi[i]);
        }
    } else {
        printf("固有値計算に失敗しました (info = %d)\n", info);
    }

    return 0;
}
