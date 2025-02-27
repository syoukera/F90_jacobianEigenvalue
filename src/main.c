#include <stdio.h>
#include <stdlib.h>
#include "jacob.h"

// LAPACK関数のプロトタイプ宣言
extern void dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda, 
                   double* wr, double* wi, double* vl, int* ldvl, 
                   double* vr, int* ldvr, double* work, int* lwork, int* info);

int main() {
    const double t = 0.0;
    const double pres = 101325.0;
    
    double y[NSP] = {1000, 0., 0.04030952, 0., 0.95969048, 0., 0., 0., 0., 0.};
    double jac[NSP * NSP] = {0.0};
    
    eval_jacob(t, pres, y, jac);  // ヤコビアンを計算

    // 固有値計算用の変数
    double wr[NSP], wi[NSP];  // 実部と虚部
    double vl[NSP * NSP];  // 左固有ベクトル（不要ならNULLでも可）
    double vr[NSP * NSP];  // 右固有ベクトル（不要ならNULLでも可）
    int lwork = 4 * NSP;  // 作業配列のサイズ
    double work[4 * NSP];  // 作業配列
    int info;
    
    char jobvl = 'V';  // calculate left eigenvelcor
    char jobvr = 'V';  // calculate right eigenvelcor
    int n = NSP;
    int lda = NSP;
    int ldvl = NSP;
    int ldvr = NSP;

    double a_exp[NSP]; // right eigenvector
    double b_exp[NSP]; // left eigenvector

    double EP[NSP];
    double EI[NSP];

    // calclate eigenvalue of Jacobian using LAPACK
    dgeev_(&jobvl, &jobvr, &n, jac, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
    
    if (info != 0) {
        printf("Error: dgeev failed with INFO = %d\n", info);
        return 1;
    }

    // find maximum eigenvalue
    double wr_max = 0.0;
    int i_wr = 2;
    for (int i = 0; i < NSP; i++) {
        // printf("Eigenvalue %d: %f\n", i, wr[i]);
        if (wr[i] > wr_max) {
            wr_max = wr[i];
            i_wr = i;
        }
    }
    printf("Maximum Eigenvalue %d: %f\n", i_wr, wr[i_wr]);
    
    // calculate EP
    printf("EP`:\n");
    double EP_sum = 0.0;
    
    for (int j = 0; j < NSP; j++) {
        // ignore complex value in eigenvector
        a_exp[j] = vr[j + i_wr*NSP]; // right
        b_exp[j] = vl[j + i_wr*NSP]; // left
        EP[j] = a_exp[j]*b_exp[j];
        EP_sum += fabs(EP[j]);
        printf("%e\n", EP[j]);  
    }

    // calculate EI
    printf("EI:\n");
    for (int j = 0; j < NSP; j++) {
        EI[j] = fabs(EP[j])/EP_sum;
        printf("%e\n", EI[j]);  
    }

    return 0;
}
