#ifndef Ragnarkinematics_H
#define Ragnarkinematics_H

#include <Arduino.h>
#include <math.h>
#include <Matrix.h>
#include <RMatrix.h>
#include <Trigonometrics.h>

bool ragnarIK(
    float (*params)[4][8], double (*pose)[6], double (*theta)[4], 
    double (*passive)[4][2], double sc[4][6], bool offset = false, 
    bool full = false);
bool ragnarIKf(
    float (*params)[4][8], float (*pose)[6], float (*theta)[4], 
    float (*passive)[4][2], float sc[4][6], bool offset, 
    bool full);
bool ragnarFK(
    float (*params)[4][8], double (*theta)[4], double (*pose)[6], 
    double sc[4][6], bool offset = false);
bool ragnarFKf(
    float (*params)[4][8], float (*theta)[4], float (*pose)[6], 
    float sc[4][6], bool offset = false);
bool ragnarFKpassivef( // computes position for end effector and passive joints 
    float (*params)[4][8], float (*theta)[4], float (*pose)[6], 
    float (*passive)[4][2], float sc[4][6], bool offset= false);
void ragnardAdB(
    double q[7], double dq[7], float param[4][8], double sct[8], 
    double sc[4][6], double (*dA)[4][3], double (*dB)[4][4]);
void ragnardAdBf(
    float q[7], float dq[7], float param[4][8], float sct[8], 
    float sc[4][6], float (*dA)[4][3], float (*dB)[4][4]);
void ragnarAB(
    double q[7], float param[4][8], double sct[8], double sc[4][6], 
    double (*A)[4][3], double (*B)[4][4]);
void ragnarABf(
    float q[7], float param[4][8], float sct[8], float sc[4][6], 
    float (*A)[4][3], float (*B)[4][4]);
//void uvw(double joints[4], float (*params)[4][8], int leg_number, double * u, double * v, double * w);
void Atra(
    double A[3][3], double (*At)[3][3]);
void Atraf(
    float A[3][3], float (*At)[3][3]);
bool Ainv(
    double A[3][3], double (*iA)[3][3]);
bool Ainvf(
    float A[3][3], float (*iA)[3][3]);
bool Ainvd(
    float A[3][3], float (*iA)[3][3]);
double detA33(
    double A[3][3]);
float detA33f(
        float A[3][3]);
void ragnarJpassive(
    float params[4][8], double  theta[4], double passive[4][2],
    double sc[4][6], double sct[8], double scez[4][4], 
    double (*J)[3][3], int leg_num);
void ragnarJpassivef(float params[4][8], float  theta[4], float passive[4][2],
                    float sc[4][6], float sct[8], float scez[4][4], 
                    float (*J)[3][3], int leg_num);
void ragnariJpassive(
    float params[4][8], double  theta[4], double passive[4][2],
    double sc[4][6], double sct[8], double scez[4][4], double (*iJ)[3][3], 
    int leg_num);
void ragnariJpassivef(
    float params[4][8], float  theta[4], float passive[4][2], float sc[4][6], 
    float sct[8], float scez[4][4], float (*iJ)[3][3], int leg_num);
void ragnardJpassive(
    float params[4][8], double  theta[4], double passive[4][2],
    double dtheta[4], double detazeta[4][2], double sc[4][6], 
    double sct[8], double scez[4][4], double (*J)[3][3], 
    int leg_num, double (*dJ)[3][3]);
void ragnardJpassivef(
    float params[4][8], float  theta[4], float passive[4][2],
    float dtheta[4], float detazeta[4][2], float sc[4][6], 
    float sct[8], float scez[4][4], int leg_num, float (*dJ)[3][3]);
bool ragnarpassiveJ(
    float A[4][3], float B[4][4], float iJp[3][3], float (*pJ)[2][3], 
    int option);
bool ragnarpassiveJ3(
    float A[4][3], float B[4][4], float iJp[3][3], float (*pJ)[3][3], 
    int option);
void ragnarParams(float x[8], float (*parameter)[4][8]);

void scfixed(float parameter[4][8], float (*sc)[4][6]);
void scth(float thetaf[4], float (*sctf)[8]);
void scpassive(float passivef[4][2], float (*scez)[4][4]);
#endif