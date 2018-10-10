#include <Arduino.h>
#include <math.h>
#include <Matrix.h>
#include <RMatrix.h>
#include <Trigonometrics.h>
#include "Ragnarkinematics.h"

#define rad2deg 57.2957795
#define deg2rad 0.01745329

void scfixed(float parameter[4][8], float (*sc)[4][6])
{
    // Computes the sine and cosine of fixed parameters 
    for (int i = 0; i < 4; i++) {
    (*sc)[i][0] = sinf(parameter[i][2]);
    (*sc)[i][1] = cosf(parameter[i][2]);
    (*sc)[i][2] = sinf(parameter[i][3]);
    (*sc)[i][3] = cosf(parameter[i][3]);
    (*sc)[i][4] = sinf(parameter[i][7]);
    (*sc)[i][5] = cosf(parameter[i][7]);
  }
}
void scth(float thetaf[4], float (*sctf)[8])
{
    for (int i = 0; i < 4; i++) {
        (*sctf)[i * 2] = sinf(thetaf[i]);
        (*sctf)[(i * 2) + 1] = cosf(thetaf[i]);
    }
}
void scpassive(float passivef[4][2], float (*scez)[4][4])
{
    for (int i = 0; i < 4; i++) {
        (*scez)[i][0] = sinf(passivef[i][0]);
        (*scez)[i][1] = cosf(passivef[i][0]);
        (*scez)[i][2] = sinf(passivef[i][1]);
        (*scez)[i][3] = cosf(passivef[i][1]);
    }
}
void ragnarParams(float x[8], float (*parameter)[4][8])
{

    //This is a first attempt and the calibration results are to be hard coded 
    float axoff[4] = {0.002560630410027, 0.007763510735409, -0.006594553620469,
                     -0.000600044264199};

    float ayoff[4] = {0.003728831051178, 0.004399368285463, 0.004200014859531,
                      0.002670164837209}; 
    
    float alphaoff[4] = {0.051742929564577, -0.001991210100304, 
                         0.038533274830325, 0.051764157170956};
    float betaoff[4] = {-0.017970056331311, -0.051811812995656, 
                        -0.015059011075560, -0.020272267321578};
    float ax[4];
    ax[0] = x[0] + axoff[0];
    ax[1] = -(x[0] + axoff[1]);
    ax[2] = -(x[0] + axoff[2]);
    ax[3] = x[0] + axoff[3];
    
    float ay[4];
    ay[0] = x[1] + ayoff[0];
    ay[1] = x[1] + ayoff[1];
    ay[2] = -(x[1] + ayoff[2]);
    ay[3] = -(x[1] + ayoff[3]);
    
    
    float gama = x[7];
    float gamma[4];
    // gama1 = gama; 
    gamma[0] = gama;
    // gamma[1] = pi / 2 + gama1;
    gamma[1] = pi / 2 + gamma[0];
    // gamma[2] = -gama2;
    gamma[2] = -gamma[1];
    // gamma[3] = -gama1;
    gamma[3] = -gamma[0];
    
    float alpha[4];
    alpha[0] = -(x[2] + alphaoff[0]);
    alpha[1] = x[2] + alphaoff[1];
    alpha[2] = -(x[2] + alphaoff[2]);
    alpha[3] = x[2] + alphaoff[3];

    float beta[4]; 
    beta[0] = -(x[3] + betaoff[0]);
    beta[1] = x[3] + betaoff[1];
    beta[2] = x[3] + betaoff[2];
    beta[3] = -(x[3] + betaoff[3]);

    for(int i=0;i<4;i++){
        (*parameter)[i][0] = ax[i];
        (*parameter)[i][1] = ay[i];
        (*parameter)[i][2] = alpha[i];
        (*parameter)[i][3] = beta[i];
        (*parameter)[i][4] = x[4];
        (*parameter)[i][5] = x[5];
        (*parameter)[i][6] = x[6];
        (*parameter)[i][7] = gamma[i];
    }

}

bool ragnarIK(
    float (*params)[4][8], double (*pose)[6], double (*theta)[4], 
    double (*passive)[4][2], double sc[4][6], bool offset, 
    bool full)
{   
    // Computer the inverse kinematics
    // access them as (*a)[0][0] .. (*a)[2][2]
    // alocate solution array
    // sc[][] contains sin and cosine of alpha beta gamma 
    // sc[i][] = {sin alpha, cos alpha, sin beta, cos beta, sin gama, cos gama}
    // offset is for the offset of the motors 
    // full is to acquire the full inverse kinematics including the passive
    // joints, if false or default then it is not computed 
    double SOL[4][2];
    int sol[4] = {1, 0, 1, 0};
    for (int i = 0; i < 4; i++)
    {
        // Active joint angles are computed to double precision. Float precision
        // result in unsolvable points in non singular workspace 
        // Passive joint angles are computed using float precision. 
        // Computation time is 1 ms in total. 
        float a = (*params)[i][0];
        float b = (*params)[i][1];
        double alpha = (*params)[i][2];
        double beta = (*params)[i][3];
        float l = (*params)[i][4];
        float L = (*params)[i][5];
        float r = (*params)[i][6];

        double sa = sc[i][0]; 
        double ca = sc[i][1];
        double sb = sc[i][2];
        double cb = sc[i][3];
        double sg = sc[i][4];
        double cg = sc[i][5];

        double gama = (*params)[i][7];
        double A[3] = {a, b, 0};
        float Af[3] = {a, b, 0}; 
        Matrix<double> Am(3, 1, (double *)A);
        Matrix<float> Amf(3, 1, (float *)Af);
        double P[3] = {(*pose)[0], (*pose)[1], (*pose)[2]};
        float Pf[3] = {(*pose)[0], (*pose)[1], (*pose)[2]};
        Matrix<double> Pm(3, 1, (double *)P);
        Matrix<float> Pmf(3, 1, (float *)Pf);
        // double C[3] = {cosf(gama) * r, sinf(gama) * r, 0};
        // float Cf[3] = {cosf(gama) * r, sinf(gama) * r, 0};
        double C[3] = {cg * r, sg * r, 0};
        float Cf[3] = {cg * r, sg * r, 0};
        Matrix<double> Cm(3, 1, (double *)C);
        Matrix<float> Cmf(3, 1, (float *)Cf);
        Cm = Cm + Pm;
        
        Cmf = Cmf + Pmf; 
        Matrix<double> CA(3, 1, (double *)C);
        CA = Cm - Am;
        double I, J, K;
        // I = 2 * l * ((Cm._entity[0][0] - a) * sinf(alpha) - (Cm._entity[1][0] - b) * cosf(alpha));
        // J = -2 * l * ((Cm._entity[0][0] - a) * cosf(alpha) * cosf(beta) + (Cm._entity[1][0] - b) * sinf(alpha) * cosf(beta) - Cm._entity[2][0] * sinf(beta));
        I = 2 * l * ((Cm._entity[0][0] - a) * sa - (Cm._entity[1][0] - b) * ca);
        J = -2 * l * ((Cm._entity[0][0] - a) * ca * cb + (Cm._entity[1][0] - b) * sa * cb - Cm._entity[2][0] * sb);

        K = (Matrix<double>::transpose(CA) * CA)._entity[0][0] + l * l - L * L;
        SOL[i][0] = 2 * atan2f((-I + sqrtf(abs(I * I + J * J - K * K))), (K - J));
        SOL[i][1] = 2 * atan2f((-I - sqrtf(abs(I * I + J * J - K * K))), (K - J));
        // this takes 1.1 ms 
        // Matrix<float> R1 = Rzyz(alpha, beta, SOL[i][sol[i]]);
        // Matrix<float> S1_1 = s1_1(alpha, beta, SOL[i][sol[i]],l);
        double zeta = 0.0;
        double eta = 0.0; 
        if (full){
            Matrix<float> R1 = Rzyz(sc, SOL[i][sol[i]], i);
            Matrix<float> S1_1 = s1_1(sc, SOL[i][sol[i]],l, i);
            Matrix<float> temp_term = Matrix<float>::inv(R1) * (Cmf - Amf - S1_1);
            // temp_term = Matrix<float>::inv(R1) * temp_term; 
            // temp_term.show(); 
            zeta = atan2f(temp_term._entity[1][0], temp_term._entity[0][0]);
            double y1 = cosf((zeta )) * (temp_term._entity[2][0]);
            eta = atan2(-y1, temp_term._entity[0][0]);
            // double eta = 0.0; 
            // double zeta = 0.0; 
            if (eta > pi12)
                eta = eta - pi;

            if (eta < -pi12)
                eta = eta + pi;
        }
        
        (*passive)[i][0] = zeta;
        (*passive)[i][1] = eta;
        //*/
        // (*passive)[i][0] = 0.0;
        // (*passive)[i][1] = 0.0;
        if (offset)
        {
            if (i < 2)
                (*theta)[i] = SOL[i][sol[i]] - pi12;
            else
                (*theta)[i] = SOL[i][sol[i]] + pi12;
        }
        else
            (*theta)[i] = SOL[i][sol[i]];

        while ((*theta)[i] > pi)
            (*theta)[i] = (*theta)[i] - 2 * pi;

        while ((*theta)[i] < -pi)
            (*theta)[i] = (*theta)[i] + 2 * pi;
    }
    return true;

} // returns value on pose
bool ragnarIKf(
    float (*params)[4][8], float (*pose)[6], float (*theta)[4], 
    float (*passive)[4][2], float sc[4][6], bool offset, 
    bool full)
{   
    // Computer the inverse kinematics
    // access them as (*a)[0][0] .. (*a)[2][2]
    // alocate solution array
    // sc[][] contains sin and cosine of alpha beta gamma 
    // sc[i][] = {sin alpha, cos alpha, sin beta, cos beta, sin gama, cos gama}
    // offset is for the offset of the motors 
    // full is to acquire the full inverse kinematics including the passive
    // joints, if false or default then it is not computed 
    float SOL[4][2];
    int sol[4] = {1, 0, 1, 0};
    for (int i = 0; i < 4; i++)
    {
        // Active joint angles are computed to double precision. Float precision
        // result in unsolvable points in non singular workspace 
        // Passive joint angles are computed using float precision. 
        // Computation time is 1 ms in total. 
        float a = (*params)[i][0];
        float b = (*params)[i][1];
        float alpha = (*params)[i][2];
        float beta = (*params)[i][3];
        float l = (*params)[i][4];
        float L = (*params)[i][5];
        float r = (*params)[i][6];

        float sa = sc[i][0]; 
        float ca = sc[i][1];
        float sb = sc[i][2];
        float cb = sc[i][3];
        float sg = sc[i][4];
        float cg = sc[i][5];

        float gama = (*params)[i][7];
        float A[3] = {a, b, 0};
        float Af[3] = {a, b, 0}; 
        Matrix<float> Am(3, 1, (float *)A);
        Matrix<float> Amf(3, 1, (float *)Af);
        float P[3] = {(*pose)[0], (*pose)[1], (*pose)[2]};
        float Pf[3] = {(*pose)[0], (*pose)[1], (*pose)[2]};
        Matrix<float> Pm(3, 1, (float *)P);
        Matrix<float> Pmf(3, 1, (float *)Pf);
        // float C[3] = {cosf(gama) * r, sinf(gama) * r, 0};
        // float Cf[3] = {cosf(gama) * r, sinf(gama) * r, 0};
        float C[3] = {cg * r, sg * r, 0};
        float Cf[3] = {cg * r, sg * r, 0};
        Matrix<float> Cm(3, 1, (float *)C);
        Matrix<float> Cmf(3, 1, (float *)Cf);
        Cm = Cm + Pm;
        
        Cmf = Cmf + Pmf; 
        Matrix<float> CA(3, 1, (float *)C);
        CA = Cm - Am;
        float I, J, K;
        // I = 2 * l * ((Cm._entity[0][0] - a) * sinf(alpha) - (Cm._entity[1][0] - b) * cosf(alpha));
        // J = -2 * l * ((Cm._entity[0][0] - a) * cosf(alpha) * cosf(beta) + (Cm._entity[1][0] - b) * sinf(alpha) * cosf(beta) - Cm._entity[2][0] * sinf(beta));
        I = 2 * l * ((Cmf._entity[0][0] - a) * sa - (Cmf._entity[1][0] - b) * ca);
        J = -2 * l * ((Cmf._entity[0][0] - a) * ca * cb + (Cmf._entity[1][0] - b) * sa * cb - Cmf._entity[2][0] * sb);

        K = (Matrix<float>::transpose(CA) * CA)._entity[0][0] + l * l - L * L;
        SOL[i][0] = 2 * atan2f((-I + sqrtf(abs(I * I + J * J - K * K))), (K - J));
        SOL[i][1] = 2 * atan2f((-I - sqrtf(abs(I * I + J * J - K * K))), (K - J));
        // this takes 1.1 ms 
        // Matrix<float> R1 = Rzyz(alpha, beta, SOL[i][sol[i]]);
        // Matrix<float> S1_1 = s1_1(alpha, beta, SOL[i][sol[i]],l);
        float zeta = 0.0;
        float eta = 0.0; 
        if (full){
            Matrix<float> R1 = Rzyzf(sc, SOL[i][sol[i]], i);
            Matrix<float> S1_1 = s1_1f(sc, SOL[i][sol[i]],l, i);
            
            Matrix<float> temp_term = Matrix<float>::inv(R1) * (Cmf - Amf - S1_1);
            
            // temp_term = Matrix<float>::inv(R1) * temp_term; 
            // temp_term.show(); 
            zeta = atan2f(temp_term._entity[1][0], temp_term._entity[0][0]);
            float y1 = cosf((zeta )) * (temp_term._entity[2][0]);
            eta = atan2f(-y1, temp_term._entity[0][0]);
            // float eta = 0.0; 
            // float zeta = 0.0; 
            if (eta > pi12)
                eta = eta - pi;

            if (eta < -pi12)
                eta = eta + pi;
        }
        
        (*passive)[i][0] = zeta;
        (*passive)[i][1] = eta;
        //*/
        // (*passive)[i][0] = 0.0;
        // (*passive)[i][1] = 0.0;
        if (offset)
        {
            if (i < 2)
                (*theta)[i] = SOL[i][sol[i]] - pi12;
            else
                (*theta)[i] = SOL[i][sol[i]] + pi12;
        }
        else
            (*theta)[i] = SOL[i][sol[i]];

        while ((*theta)[i] > pi)
            (*theta)[i] = (*theta)[i] - 2 * pi;

        while ((*theta)[i] < -pi)
            (*theta)[i] = (*theta)[i] + 2 * pi;
    }
    return true;

} // returns value on pose

void uvw(double joints[4], float (*params)[4][8], int leg_number, double *u, double *v, double *w)
{
    double theta1 = joints[leg_number];
    float ax1 = (*params)[leg_number][0];
    float ay1 = (*params)[leg_number][1];
    double alpha1 = (*params)[leg_number][2];
    double beta1 = (*params)[leg_number][3];
    float b1 = (*params)[leg_number][4];
    float r1 = (*params)[leg_number][6];
    double gama1 = (*params)[leg_number][7];

    // Serial.println("theta");
    // Serial.println(theta1);
    // coefficients of polynomials
    *u = -ax1 + r1 * cosf(gama1) + b1 * (sinf(alpha1) * sinf(theta1) - cosf(alpha1) * cosf(beta1) * cosf(theta1));
    *v = -(ay1 - r1 * sinf(gama1) + b1 * (cosf(alpha1) * sinf(theta1) + cosf(beta1) * sinf(alpha1) * (cosf(theta1))));
    *w = b1 * sinf(beta1) * cosf(theta1);    
}

void uvwf(float joints[4], float (*params)[4][8], int leg_number, float *u, float *v, float *w)
{
    float theta1 = joints[leg_number];
    float ax1 = (*params)[leg_number][0];
    float ay1 = (*params)[leg_number][1];
    float alpha1 = (*params)[leg_number][2];
    float beta1 = (*params)[leg_number][3];
    float b1 = (*params)[leg_number][4];
    float r1 = (*params)[leg_number][6];
    float gama1 = (*params)[leg_number][7];

    // Serial.println("theta");
    // Serial.println(theta1);
    // coefficients of polynomials
    *u = -ax1 + r1 * cosf(gama1) + b1 * (sinf(alpha1) * sinf(theta1) - cosf(alpha1) * cosf(beta1) * cosf(theta1));
    *v = -(ay1 - r1 * sinf(gama1) + b1 * (cosf(alpha1) * sinf(theta1) + cosf(beta1) * sinf(alpha1) * (cosf(theta1))));
    *w = b1 * sinf(beta1) * cosf(theta1);    
}
bool ragnarFK(
    float (*params)[4][8], double (*theta)[4], double (*pose)[6], 
    double sc[4][6], bool offset)
{
    // Computation time is 108 microseconds. 
    // receive joints in degrees
    double joints[4];
    if(offset){
        joints[0] = (*theta)[0] + pi12;
        joints[1] = (*theta)[1] + pi12;
        joints[2] = (*theta)[2] - pi12;
        joints[3] = (*theta)[3] - pi12;
    }
    else {
        joints[0] = (*theta)[0];
        joints[1] = (*theta)[1];
        joints[2] = (*theta)[2];
        joints[3] = (*theta)[3];
    }
    double u1, v1, w1;
    double u2, v2, w2;
    double u3, v3, w3;
    double u4, v4, w4;
    double u, v, w;
    double s1x, s1y, s1z;
    double s2x, s2y, s2z;
    double s11, s22;
    double D[4];
    double g[2];
    int r1, r2;
    double xy[2];
    double sxx, sxy, sxz;
    double x1, y1, z1;
    double x2, y2, z2;
    double A, B, C, M;
    double normSol1, normSol2;
    float parameter[32];
    double jointsSol[4];
    double cartSol[3];

    float L = (*params)[0][5];
    
    uvw(joints, params, 0, &u1, &v1, &w1);
    uvw(joints, params, 1, &u2, &v2, &w2);
    uvw(joints, params, 2, &u3, &v3, &w3);
    uvw(joints, params, 3, &u4, &v4, &w4);
    int chosen = 0;
    // Chooses which arms to use for solution
    // If the angle between arms 1 and 2 is larger than arms 3 and 4, use arms 4,2,3 (with 4-2 and 4-3)
    if (fabs(joints[0] - joints[1]) > fabs(joints[2] - joints[3]))
    {
        // Serial.println(u4);
        // Serial.println(v4);
        // Serial.println(w4);
        chosen = 1;
        s1x = 2.0 * u2 - 2.0 * u3;
        s1y = 2.0 * v2 - 2.0 * v3;
        s1z = 2.0 * w2 - 2.0 * w3;
        s11 = ((((powf(u2, 2.0) - powf(u3, 2.0)) + powf(v2, 2.0)) - powf(v3, 2.0)) + powf(w2, 2.0)) - powf(w3, 2.0);
        s2x = 2.0 * u2 - 2.0 * u4;
        s2y = 2.0 * v2 - 2.0 * v4;
        s2z = 2.0 * w2 - 2.0 * w4;
        s22 = ((((powf(u2, 2.0) - powf(u4, 2.0)) + powf(v2, 2.0)) - powf(v4, 2.0)) + powf(w2, 2.0)) - powf(w4, 2.0);
        u = u2;
        v = v2;
        w = w2;
    }
    // Else use arms 2,1,3 (with 2-1 and 2-3)
    else
    {
        chosen = 0;
        s1x = 2.0 * u1 - 2.0 * u2;
        s1y = 2.0 * v1 - 2.0 * v2;
        s1z = 2.0 * w1 - 2.0 * w2;
        s11 = ((((powf(u1, 2.0) - powf(u2, 2.0)) + powf(v1, 2.0)) - powf(v2, 2.0)) + powf(w1, 2.0)) - powf(w2, 2.0);
        s2x = 2.0 * u1 - 2.0 * u3;
        s2y = 2.0 * v1 - 2.0 * v3;
        s2z = 2.0 * w1 - 2.0 * w3;
        s22 = ((((powf(u1, 2.0) - powf(u3, 2.0)) + powf(v1, 2.0)) - powf(v3, 2.0)) + powf(w1, 2.0)) - powf(w3, 2.0);
        u = u1;
        v = v1;
        w = w1;
    }

    if (joints[0 + chosen] == pi12 && joints[1 + chosen] == pi12 && joints[2 + chosen] == -pi12)
    {
        // Serial.println("Extreme option");
        float D[2][2];
        D[0][0] = s1x;
        D[0][1] = s1y;
        D[1][0] = s2x;
        D[1][1] = s2y;
        Matrix<double> Dm(2, 2, (double *)D);
        double g[2] = {-s11, -s22};
        Matrix<double> gm(2, 1, (double *)g);
        Matrix<double> xym = Matrix<double>::inv(Dm) * gm;
        double z;
        z = -sqrtf(fabs(L * L - (xym._entity[0][0] + u) * (xym._entity[0][0] + u) - (xym._entity[1][0] + v) * (xym._entity[1][0] + v)));
        x1 = xym._entity[0][0];
        y1 = xym._entity[1][0];
        z1 = z;
        x2 = x1;
        y2 = y1;
        z2 = z;
    }
    else
    {
        // Serial.println("else");
        sxx = powf((s1x * s2y - s2x * s1y), 2.0) / powf((s1y * s2z - s2y * s1z), 2.0) \
            + powf((s1x * s2z - s2x * s1z), 2.0) / powf((s1y * s2z - s2y * s1z), 2.0) + 1;
        sxy = 2 * u1 - (2 * v1 * (s1x * s2z - s2x * s1z)) / (s1y * s2z - s2y * s1z) \
           + (2 * w1 * (s1x * s2y - s2x * s1y)) / (s1y * s2z - s2y * s1z) \
           + (2 * (s11 * s2y - s22 * s1y) * (s1x * s2y - s2x * s1y)) / powf((s1y * s2z - s2y * s1z), 2.0) \
           + (2 * (s11 * s2z - s22 * s1z) * (s1x * s2z - s2x * s1z)) / powf((s1y * s2z - s2y * s1z), 2.0);
        sxz = powf(u1, 2.0) - powf(L, 2.0) + powf(v1, 2.0) + powf(w1, 2.0) \
            + powf((s11 * s2y - s22 * s1y), 2.0) / powf((s1y * s2z - s2y * s1z), 2.0) \
            + powf((s11 * s2z - s22 * s1z), 2.0) / powf((s1y * s2z - s2y * s1z), 2.0) \
            - (2 * v1 * (s11 * s2z - s22 * s1z)) / (s1y * s2z - s2y * s1z)\
            + (2 * w1 * (s11 * s2y - s22 * s1y)) / (s1y * s2z - s2y * s1z);
        // Serial.println(sxx);
        // Serial.println(sxy);
        // Serial.println(sxz);
        x1 = (-sxy + sqrtf(fabs(powf(sxy, 2.0) - 4 * sxx * sxz))) / (2 * sxx);
        y1 = -(s11 * s2z - s22 * s1z + s1x * s2z * x1 - s2x * s1z * x1) / (s1y * s2z - s2y * s1z);
        z1 = (s11 * s2y - s22 * s1y + s1x * s2y * x1 - s2x * s1y * x1) / (s1y * s2z - s2y * s1z);

        x2 = (-sxy - sqrtf(fabs(powf(sxy, 2.0) - 4 * sxx * sxz))) / (2 * sxx);
        y2 = -(s11 * s2z - s22 * s1z + s1x * s2z * x2 - s2x * s1z * x2) / (s1y * s2z - s2y * s1z);
        z2 = (s11 * s2y - s22 * s1y + s1x * s2y * x2 - s2x * s1y * x2) / (s1y * s2z - s2y * s1z);
    }
    cartSol[0] = x2;
    cartSol[1] = y2;
    cartSol[2] = z2;

    if (z1 < 0)
    {
        (*pose)[0] = x1;
        (*pose)[1] = y1;
        (*pose)[2] = z1;
        (*pose)[3] = 0.0;
        (*pose)[4] = 0.0;
        (*pose)[5] = 0.0;
    }
    else
    {
        (*pose)[0] = x2;
        (*pose)[1] = y2;
        (*pose)[2] = z2;
        (*pose)[3] = 0.0;
        (*pose)[4] = 0.0;
        (*pose)[5] = 0.0;
    }

    return true;
}
bool ragnarFKf(
    float (*params)[4][8], float (*theta)[4], float (*pose)[6], 
    float sc[4][6], bool offset)
{
    // Computation time is 108 microseconds. 
    // receive joints in degrees
    float joints[4];
    if(offset){
        joints[0] = (*theta)[0] + pi12;
        joints[1] = (*theta)[1] + pi12;
        joints[2] = (*theta)[2] - pi12;
        joints[3] = (*theta)[3] - pi12;
    }
    else {
        joints[0] = (*theta)[0];
        joints[1] = (*theta)[1];
        joints[2] = (*theta)[2];
        joints[3] = (*theta)[3];
    }
    float u1, v1, w1;
    float u2, v2, w2;
    float u3, v3, w3;
    float u4, v4, w4;
    float u, v, w;
    float s1x, s1y, s1z;
    float s2x, s2y, s2z;
    float s11, s22;
    float D[4];
    float g[2];
    int r1, r2;
    float xy[2];
    float sxx, sxy, sxz;
    float x1, y1, z1;
    float x2, y2, z2;
    float A, B, C, M;
    float normSol1, normSol2;
    float parameter[32];
    float jointsSol[4];
    float cartSol[3];

    float L = (*params)[0][5];
    
    uvwf(joints, params, 0, &u1, &v1, &w1);
    uvwf(joints, params, 1, &u2, &v2, &w2);
    uvwf(joints, params, 2, &u3, &v3, &w3);
    uvwf(joints, params, 3, &u4, &v4, &w4);
    int chosen = 0;
    // Chooses which arms to use for solution
    // If the angle between arms 1 and 2 is larger than arms 3 and 4, use arms 4,2,3 (with 4-2 and 4-3)
    if (fabs(joints[0] - joints[1]) > fabs(joints[2] - joints[3]))
    {
        // Serial.println(u4);
        // Serial.println(v4);
        // Serial.println(w4);
        chosen = 1;
        s1x = 2.0 * u2 - 2.0 * u3;
        s1y = 2.0 * v2 - 2.0 * v3;
        s1z = 2.0 * w2 - 2.0 * w3;
        s11 = ((((powf(u2, 2.0) - powf(u3, 2.0)) + powf(v2, 2.0)) - powf(v3, 2.0)) + powf(w2, 2.0)) - powf(w3, 2.0);
        s2x = 2.0 * u2 - 2.0 * u4;
        s2y = 2.0 * v2 - 2.0 * v4;
        s2z = 2.0 * w2 - 2.0 * w4;
        s22 = ((((powf(u2, 2.0) - powf(u4, 2.0)) + powf(v2, 2.0)) - powf(v4, 2.0)) + powf(w2, 2.0)) - powf(w4, 2.0);
        u = u2;
        v = v2;
        w = w2;
    }
    // Else use arms 2,1,3 (with 2-1 and 2-3)
    else
    {
        chosen = 0;
        s1x = 2.0 * u1 - 2.0 * u2;
        s1y = 2.0 * v1 - 2.0 * v2;
        s1z = 2.0 * w1 - 2.0 * w2;
        s11 = ((((powf(u1, 2.0) - powf(u2, 2.0)) + powf(v1, 2.0)) - powf(v2, 2.0)) + powf(w1, 2.0)) - powf(w2, 2.0);
        s2x = 2.0 * u1 - 2.0 * u3;
        s2y = 2.0 * v1 - 2.0 * v3;
        s2z = 2.0 * w1 - 2.0 * w3;
        s22 = ((((powf(u1, 2.0) - powf(u3, 2.0)) + powf(v1, 2.0)) - powf(v3, 2.0)) + powf(w1, 2.0)) - powf(w3, 2.0);
        u = u1;
        v = v1;
        w = w1;
    }

    if (joints[0 + chosen] == pi12 && joints[1 + chosen] == pi12 && joints[2 + chosen] == -pi12)
    {
        // Serial.println("Extreme option");
        float D[2][2];
        D[0][0] = s1x;
        D[0][1] = s1y;
        D[1][0] = s2x;
        D[1][1] = s2y;
        Matrix<float> Dm(2, 2, (float *)D);
        float g[2] = {-s11, -s22};
        Matrix<float> gm(2, 1, (float *)g);
        Matrix<float> xym = Matrix<float>::inv(Dm) * gm;
        float z;
        z = -sqrtf(fabs(L * L - (xym._entity[0][0] + u) * (xym._entity[0][0] + u) - (xym._entity[1][0] + v) * (xym._entity[1][0] + v)));
        x1 = xym._entity[0][0];
        y1 = xym._entity[1][0];
        z1 = z;
        x2 = x1;
        y2 = y1;
        z2 = z;
    }
    else
    {
        // Serial.println("else");
        sxx = powf((s1x * s2y - s2x * s1y), 2.0) / powf((s1y * s2z - s2y * s1z), 2.0) \
            + powf((s1x * s2z - s2x * s1z), 2.0) / powf((s1y * s2z - s2y * s1z), 2.0) + 1;
        sxy = 2 * u1 - (2 * v1 * (s1x * s2z - s2x * s1z)) / (s1y * s2z - s2y * s1z) \
           + (2 * w1 * (s1x * s2y - s2x * s1y)) / (s1y * s2z - s2y * s1z) \
           + (2 * (s11 * s2y - s22 * s1y) * (s1x * s2y - s2x * s1y)) / powf((s1y * s2z - s2y * s1z), 2.0) \
           + (2 * (s11 * s2z - s22 * s1z) * (s1x * s2z - s2x * s1z)) / powf((s1y * s2z - s2y * s1z), 2.0);
        sxz = powf(u1, 2.0) - powf(L, 2.0) + powf(v1, 2.0) + powf(w1, 2.0) \
            + powf((s11 * s2y - s22 * s1y), 2.0) / powf((s1y * s2z - s2y * s1z), 2.0) \
            + powf((s11 * s2z - s22 * s1z), 2.0) / powf((s1y * s2z - s2y * s1z), 2.0) \
            - (2 * v1 * (s11 * s2z - s22 * s1z)) / (s1y * s2z - s2y * s1z)\
            + (2 * w1 * (s11 * s2y - s22 * s1y)) / (s1y * s2z - s2y * s1z);
        // Serial.println(sxx);
        // Serial.println(sxy);
        // Serial.println(sxz);
        x1 = (-sxy + sqrtf(fabs(powf(sxy, 2.0) - 4 * sxx * sxz))) / (2 * sxx);
        y1 = -(s11 * s2z - s22 * s1z + s1x * s2z * x1 - s2x * s1z * x1) / (s1y * s2z - s2y * s1z);
        z1 = (s11 * s2y - s22 * s1y + s1x * s2y * x1 - s2x * s1y * x1) / (s1y * s2z - s2y * s1z);

        x2 = (-sxy - sqrtf(fabs(powf(sxy, 2.0) - 4 * sxx * sxz))) / (2 * sxx);
        y2 = -(s11 * s2z - s22 * s1z + s1x * s2z * x2 - s2x * s1z * x2) / (s1y * s2z - s2y * s1z);
        z2 = (s11 * s2y - s22 * s1y + s1x * s2y * x2 - s2x * s1y * x2) / (s1y * s2z - s2y * s1z);
    }
    cartSol[0] = x2;
    cartSol[1] = y2;
    cartSol[2] = z2;

    if (z1 < 0)
    {
        (*pose)[0] = x1;
        (*pose)[1] = y1;
        (*pose)[2] = z1;
        (*pose)[3] = 0.0;
        (*pose)[4] = 0.0;
        (*pose)[5] = 0.0;
    }
    else
    {
        (*pose)[0] = x2;
        (*pose)[1] = y2;
        (*pose)[2] = z2;
        (*pose)[3] = 0.0;
        (*pose)[4] = 0.0;
        (*pose)[5] = 0.0;
    }

    return true;
}

bool ragnarFKpassivef(
    float (*params)[4][8], float (*theta)[4], float (*pose)[6], 
    float (*passive)[4][2], float sc[4][6], bool offset)
{
    // Computation time is 108 microseconds. 
    // receive joints in degrees
    float joints[4];
    if(offset){
        joints[0] = (*theta)[0] + pi12;
        joints[1] = (*theta)[1] + pi12;
        joints[2] = (*theta)[2] - pi12;
        joints[3] = (*theta)[3] - pi12;
    }
    else {
        joints[0] = (*theta)[0];
        joints[1] = (*theta)[1];
        joints[2] = (*theta)[2];
        joints[3] = (*theta)[3];
    }
    float u1, v1, w1;
    float u2, v2, w2;
    float u3, v3, w3;
    float u4, v4, w4;
    float u, v, w;
    float s1x, s1y, s1z;
    float s2x, s2y, s2z;
    float s11, s22;
    float D[4];
    float g[2];
    int r1, r2;
    float xy[2];
    float sxx, sxy, sxz;
    float x1, y1, z1;
    float x2, y2, z2;
    float A, B, C, M;
    float normSol1, normSol2;
    float parameter[32];
    float jointsSol[4];
    float cartSol[3];

    float L = (*params)[0][5];
    
    uvwf(joints, params, 0, &u1, &v1, &w1);
    uvwf(joints, params, 1, &u2, &v2, &w2);
    uvwf(joints, params, 2, &u3, &v3, &w3);
    uvwf(joints, params, 3, &u4, &v4, &w4);
    int chosen = 0;
    // Chooses which arms to use for solution
    // If the angle between arms 1 and 2 is larger than arms 3 and 4, use arms 4,2,3 (with 4-2 and 4-3)
    if (fabs(joints[0] - joints[1]) > fabs(joints[2] - joints[3]))
    {
        // Serial.println(u4);
        // Serial.println(v4);
        // Serial.println(w4);
        chosen = 1;
        s1x = 2.0 * u2 - 2.0 * u3;
        s1y = 2.0 * v2 - 2.0 * v3;
        s1z = 2.0 * w2 - 2.0 * w3;
        s11 = ((((powf(u2, 2.0) - powf(u3, 2.0)) + powf(v2, 2.0)) - powf(v3, 2.0)) + powf(w2, 2.0)) - powf(w3, 2.0);
        s2x = 2.0 * u2 - 2.0 * u4;
        s2y = 2.0 * v2 - 2.0 * v4;
        s2z = 2.0 * w2 - 2.0 * w4;
        s22 = ((((powf(u2, 2.0) - powf(u4, 2.0)) + powf(v2, 2.0)) - powf(v4, 2.0)) + powf(w2, 2.0)) - powf(w4, 2.0);
        u = u2;
        v = v2;
        w = w2;
    }
    // Else use arms 2,1,3 (with 2-1 and 2-3)
    else
    {
        chosen = 0;
        s1x = 2.0 * u1 - 2.0 * u2;
        s1y = 2.0 * v1 - 2.0 * v2;
        s1z = 2.0 * w1 - 2.0 * w2;
        s11 = ((((powf(u1, 2.0) - powf(u2, 2.0)) + powf(v1, 2.0)) - powf(v2, 2.0)) + powf(w1, 2.0)) - powf(w2, 2.0);
        s2x = 2.0 * u1 - 2.0 * u3;
        s2y = 2.0 * v1 - 2.0 * v3;
        s2z = 2.0 * w1 - 2.0 * w3;
        s22 = ((((powf(u1, 2.0) - powf(u3, 2.0)) + powf(v1, 2.0)) - powf(v3, 2.0)) + powf(w1, 2.0)) - powf(w3, 2.0);
        u = u1;
        v = v1;
        w = w1;
    }

    if (joints[0 + chosen] == pi12 && joints[1 + chosen] == pi12 && joints[2 + chosen] == -pi12)
    {
        // Serial.println("Extreme option");
        float D[2][2];
        D[0][0] = s1x;
        D[0][1] = s1y;
        D[1][0] = s2x;
        D[1][1] = s2y;
        Matrix<float> Dm(2, 2, (float *)D);
        float g[2] = {-s11, -s22};
        Matrix<float> gm(2, 1, (float *)g);
        Matrix<float> xym = Matrix<float>::inv(Dm) * gm;
        float z;
        z = -sqrtf(fabs(L * L - (xym._entity[0][0] + u) * (xym._entity[0][0] + u) - (xym._entity[1][0] + v) * (xym._entity[1][0] + v)));
        x1 = xym._entity[0][0];
        y1 = xym._entity[1][0];
        z1 = z;
        x2 = x1;
        y2 = y1;
        z2 = z;
    }
    else
    {
        // Serial.println("else");
        sxx = powf((s1x * s2y - s2x * s1y), 2.0) / powf((s1y * s2z - s2y * s1z), 2.0) \
            + powf((s1x * s2z - s2x * s1z), 2.0) / powf((s1y * s2z - s2y * s1z), 2.0) + 1;
        sxy = 2 * u1 - (2 * v1 * (s1x * s2z - s2x * s1z)) / (s1y * s2z - s2y * s1z) \
           + (2 * w1 * (s1x * s2y - s2x * s1y)) / (s1y * s2z - s2y * s1z) \
           + (2 * (s11 * s2y - s22 * s1y) * (s1x * s2y - s2x * s1y)) / powf((s1y * s2z - s2y * s1z), 2.0) \
           + (2 * (s11 * s2z - s22 * s1z) * (s1x * s2z - s2x * s1z)) / powf((s1y * s2z - s2y * s1z), 2.0);
        sxz = powf(u1, 2.0) - powf(L, 2.0) + powf(v1, 2.0) + powf(w1, 2.0) \
            + powf((s11 * s2y - s22 * s1y), 2.0) / powf((s1y * s2z - s2y * s1z), 2.0) \
            + powf((s11 * s2z - s22 * s1z), 2.0) / powf((s1y * s2z - s2y * s1z), 2.0) \
            - (2 * v1 * (s11 * s2z - s22 * s1z)) / (s1y * s2z - s2y * s1z)\
            + (2 * w1 * (s11 * s2y - s22 * s1y)) / (s1y * s2z - s2y * s1z);
        // Serial.println(sxx);
        // Serial.println(sxy);
        // Serial.println(sxz);
        x1 = (-sxy + sqrtf(fabs(powf(sxy, 2.0) - 4 * sxx * sxz))) / (2 * sxx);
        y1 = -(s11 * s2z - s22 * s1z + s1x * s2z * x1 - s2x * s1z * x1) / (s1y * s2z - s2y * s1z);
        z1 = (s11 * s2y - s22 * s1y + s1x * s2y * x1 - s2x * s1y * x1) / (s1y * s2z - s2y * s1z);

        x2 = (-sxy - sqrtf(fabs(powf(sxy, 2.0) - 4 * sxx * sxz))) / (2 * sxx);
        y2 = -(s11 * s2z - s22 * s1z + s1x * s2z * x2 - s2x * s1z * x2) / (s1y * s2z - s2y * s1z);
        z2 = (s11 * s2y - s22 * s1y + s1x * s2y * x2 - s2x * s1y * x2) / (s1y * s2z - s2y * s1z);
    }
    cartSol[0] = x2;
    cartSol[1] = y2;
    cartSol[2] = z2;

    if (z1 < 0)
    {
        (*pose)[0] = x1;
        (*pose)[1] = y1;
        (*pose)[2] = z1;
        (*pose)[3] = 0.0;
        (*pose)[4] = 0.0;
        (*pose)[5] = 0.0;
    }
    else
    {
        (*pose)[0] = x2;
        (*pose)[1] = y2;
        (*pose)[2] = z2;
        (*pose)[3] = 0.0;
        (*pose)[4] = 0.0;
        (*pose)[5] = 0.0;
    }

    for(int i=0;i<4;i++){
        float a = (*params)[i][0];
        float b = (*params)[i][1];
        float alpha = (*params)[i][2];
        float beta = (*params)[i][3];
        float l = (*params)[i][4];
        float L = (*params)[i][5];
        float r = (*params)[i][6];

        float sa = sc[i][0]; 
        float ca = sc[i][1];
        float sb = sc[i][2];
        float cb = sc[i][3];
        float sg = sc[i][4];
        float cg = sc[i][5];

        float gama = (*params)[i][7];
        
        float Af[3] = {a, b, 0}; 
        Matrix<float> Amf(3, 1, (float *)Af);
        float Pf[3] = {(*pose)[0], (*pose)[1], (*pose)[2]};
        Matrix<float> Pmf(3, 1, (float *)Pf);
        // float C[3] = {cosf(gama) * r, sinf(gama) * r, 0};
        // float Cf[3] = {cosf(gama) * r, sinf(gama) * r, 0};
        float Cf[3] = {cg * r, sg * r, 0};
        Matrix<float> Cmf(3, 1, (float *)Cf);
        
        Cmf = Cmf + Pmf; 
        
        float zeta = 0.0;
        float eta = 0.0; 
        
        Matrix<float> R1 = Rzyzf(sc, joints[i], i);
        Matrix<float> S1_1 = s1_1f(sc, joints[i],l, i);
            
        Matrix<float> temp_term = Matrix<float>::inv(R1) * (Cmf - Amf - S1_1);
            
        // temp_term = Matrix<float>::inv(R1) * temp_term; 
        // temp_term.show(); 
        zeta = atan2f(temp_term._entity[1][0], temp_term._entity[0][0]);
        float yy1 = cosf((zeta )) * (temp_term._entity[2][0]);
        eta = atan2f(-yy1, temp_term._entity[0][0]);
        // float eta = 0.0; 
        // float zeta = 0.0; 
        if (eta > pi12)
            eta = eta - pi;
        if (eta < -pi12)
            eta = eta + pi;
        
            
        (*passive)[i][0] = zeta;
        (*passive)[i][1] = eta;
    }        
    return true;
}


void ragnardAdB(
    double q[7], double dq[7], float param[4][8], double sct[8], 
    double sc[4][6], double (*dA)[4][3], double (*dB)[4][4])
{
    //TODO reduce memory usage
    double sa[4], ca[4], sb[4], cb[4], sg[4], cg[4];
    double st[4], ct[4];
    double x,y,z,dx,dy,dz,th1,th2,th3,th4, dt[4]; 
    float b,r,ax[4],ay[4]; 

    // extract position and derivatives 
    th1 = q[0];
    th2 = q[1];
    th3 = q[2];
    th4 = q[3];
    x = q[4];
    y = q[5];
    z = q[6];

    dx = dq[4];
    dy = dq[5];
    dz = dq[6];

    b = param[0][4]; //same as l 
    r = param[0][6];

    // extract precomputed sines and cosines
    for (int i=0; i<4; i++){
        dt[i] = dq[i];
        sa[i] = sc[i][0]; // sin alpha
        ca[i] = sc[i][1]; // cos alpha 
        sb[i] = sc[i][2]; // sin beta
        cb[i] = sc[i][3]; // cos beta
        sg[i] = sc[i][4]; // sin gama 
        cg[i] = sc[i][5]; // cos gama 
        // they are indexed as {s1, c1, s2, c2, ... 
        st[i] = sct[i*2]; // sin theta 
        ct[i] = sct[(i*2)+1]; // cos theta 
        ax[i] = param[i][0];
        ay[i] = param[i][1]; 
        for (int j=0; j<4; j++){
            (*dB)[i][j] = 0.0; 
            if(j!=3)
                (*dA)[i][j] = 0.0;    
        }
        (*dA)[i][0] = dx + b*dt[i]*(sa[i]*ct[i] + ca[i]*cb[i]*st[i]);
        (*dA)[i][1] = dy - b*dt[i]*(ca[i]*ct[i] - cb[i]*sa[i]*st[i]);
        (*dA)[i][2] = dz - b*dt[i]*sb[i]*st[i];

        (*dB)[i][i] = b*dy*ca[i]*ct[i] - b*dx*sa[i]*ct[i] + b*dz*sb[i]*st[i] + ay[i]*b*dt[i]*ca[i]*st[i] \
                     - ax[i]*b*dt[i]*sa[i]*st[i] - b*dt[i]*y*ca[i]*st[i] + b*dt[i]*z*sb[i]*ct[i] \
                     + b*st[i]*x*sa[i]*st[i] - b*dx*ca[i]*cb[i]*st[i] - b*dy*cb[i]*sa[i]*st[i] - b*dt[i]*x*ca[i]*cb[i]*ct[i] \
                     - b*dt[i]*y*cb[i]*sa[i]*ct[i] - b*dt[i]*r*ca[i]*sg[i]*st[i] + b*dt[i]*r*cg[i]*sa[i]*st[i] \
                     + ax[i]*b*dt[i]*ca[i]*cb[i]*ct[i] + ay[i]*b*dt[i]*cb[i]*sa[i]*ct[i] \
                     - b*dt[i]*r*ca[i]*cb[i]*cg[i]*ct[i] - b*dt[i]*r*cb[i]*sa[i]*ct[i]*sg[i];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0;
    }           
}

void ragnardAdBf(
    float q[7], float dq[7], float param[4][8], float sct[8], 
    float sc[4][6], float (*dA)[4][3], float (*dB)[4][4])
{
    //TODO reduce memory usage
    float sa[4], ca[4], sb[4], cb[4], sg[4], cg[4];
    float st[4], ct[4];
    float x,y,z,dx,dy,dz,th1,th2,th3,th4, dt[4]; 
    float b,r,ax[4],ay[4]; 

    // extract position and derivatives 
    th1 = q[0];
    th2 = q[1];
    th3 = q[2];
    th4 = q[3];
    x = q[4];
    y = q[5];
    z = q[6];

    dx = dq[4];
    dy = dq[5];
    dz = dq[6];

    b = param[0][4]; //same as l 
    r = param[0][6];

    // extract precomputed sines and cosines
    for (int i=0; i<4; i++){
        dt[i] = dq[i];
        sa[i] = sc[i][0]; // sin alpha
        ca[i] = sc[i][1]; // cos alpha 
        sb[i] = sc[i][2]; // sin beta
        cb[i] = sc[i][3]; // cos beta
        sg[i] = sc[i][4]; // sin gama 
        cg[i] = sc[i][5]; // cos gama 
        // they are indexed as {s1, c1, s2, c2, ... 
        st[i] = sct[i*2]; // sin theta 
        ct[i] = sct[(i*2)+1]; // cos theta 
        ax[i] = param[i][0];
        ay[i] = param[i][1]; 
        for (int j=0; j<4; j++){
            (*dB)[i][j] = 0.0; 
            if(j!=3)
                (*dA)[i][j] = 0.0;    
        }
        (*dA)[i][0] = dx + b*dt[i]*(sa[i]*ct[i] + ca[i]*cb[i]*st[i]);
        (*dA)[i][1] = dy - b*dt[i]*(ca[i]*ct[i] - cb[i]*sa[i]*st[i]);
        (*dA)[i][2] = dz - b*dt[i]*sb[i]*st[i];

        (*dB)[i][i] = b*dy*ca[i]*ct[i] - b*dx*sa[i]*ct[i] + b*dz*sb[i]*st[i] + ay[i]*b*dt[i]*ca[i]*st[i] \
                     - ax[i]*b*dt[i]*sa[i]*st[i] - b*dt[i]*y*ca[i]*st[i] + b*dt[i]*z*sb[i]*ct[i] \
                     + b*st[i]*x*sa[i]*st[i] - b*dx*ca[i]*cb[i]*st[i] - b*dy*cb[i]*sa[i]*st[i] - b*dt[i]*x*ca[i]*cb[i]*ct[i] \
                     - b*dt[i]*y*cb[i]*sa[i]*ct[i] - b*dt[i]*r*ca[i]*sg[i]*st[i] + b*dt[i]*r*cg[i]*sa[i]*st[i] \
                     + ax[i]*b*dt[i]*ca[i]*cb[i]*ct[i] + ay[i]*b*dt[i]*cb[i]*sa[i]*ct[i] \
                     - b*dt[i]*r*ca[i]*cb[i]*cg[i]*ct[i] - b*dt[i]*r*cb[i]*sa[i]*ct[i]*sg[i];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0;
    }           
}

void ragnarAB(
    double q[7], float param[4][8], double sct[8], double sc[4][6], 
    double (*A)[4][3], double (*B)[4][4])
{
    // This computes the constraint Jacobian 
    // TODO reduce memory used here. 
    // extract data 
    double x,y,z,b,r;
    double ax[4], ay[4], st[4], ct[4], sa[4], ca[4], sb[4], cb[4], sg[4], cg[4]; 

    x = q[4];
    y = q[5];
    z = q[6];
    b = param[0][4]; //same as l 
    r = param[0][6];

    for(int i=0; i<4; i++){
        st[i] = sct[i*2]; // sin theta 
        ct[i] = sct[(i*2)+1]; // cos theta 
        ax[i] = param[i][0];
        ay[i] = param[i][1]; 
        sa[i] = sc[i][0]; // sin alpha
        ca[i] = sc[i][1]; // cos alpha 
        sb[i] = sc[i][2]; // sin beta
        cb[i] = sc[i][3]; // cos beta
        sg[i] = sc[i][4]; // sin gama 
        cg[i] = sc[i][5]; // cos gama 
        for (int j=0; j<4; j++){
            (*B)[i][j] = 0.0; 
            if(j!=3)
                (*A)[i][j] = 0.0;    
        }
        (*A)[i][0] = x - ax[i] + r*cg[i] + b*(sa[i]*st[i] - ca[i]*cb[i]*ct[i]);
        (*A)[i][1] = y - ay[i] + r*sg[i] - b*(ca[i]*st[i] + cb[i]*sa[i]*ct[i]);
        (*A)[i][2] = z + b*sb[i]*ct[i];

        (*B)[i][i] = (z + b*sb[i]*ct[i])*(b*ca[i]*sb[i]*(ca[i]*st[i] + cb[i]*sa[i]*ct[i]) + b*sa[i]*sb[i]*(sa[i]*st[i] \
                    - ca[i]*cb[i]*ct[i])) - (b*sa[i]*ct[i]*powf(sb[i],2.0) + b*cb[i]*(ca[i]*st[i] \
                    + cb[i]*sa[i]*ct[i]))*(x - ax[i] + r*cg[i] + b*(sa[i]*st[i] - ca[i]*cb[i]*ct[i])) \
                    + (b*cb[i]*(sa[i]*st[i] - ca[i]*cb[i]*ct[i]) - b*ca[i]*powf(sb[i],2.0)*ct[i])*(ay[i] - y - r*sg[i]\
                     + b*(ca[i]*st[i] + cb[i]*sa[i]*ct[i]));
    }
}

void ragnarABf(
    float q[7], float param[4][8], float sct[8], float sc[4][6], 
    float (*A)[4][3], float (*B)[4][4])
{
    // This computes the constraint Jacobian 
    // TODO reduce memory used here. 
    // extract data 
    float x,y,z,b,r;
    float ax[4], ay[4], st[4], ct[4], sa[4], ca[4], sb[4], cb[4], sg[4], cg[4]; 

    x = q[4];
    y = q[5];
    z = q[6];
    b = param[0][4]; //same as l 
    r = param[0][6];

    for(int i=0; i<4; i++){
        st[i] = sct[i*2]; // sin theta 
        ct[i] = sct[(i*2)+1]; // cos theta 
        ax[i] = param[i][0];
        ay[i] = param[i][1]; 
        sa[i] = sc[i][0]; // sin alpha
        ca[i] = sc[i][1]; // cos alpha 
        sb[i] = sc[i][2]; // sin beta
        cb[i] = sc[i][3]; // cos beta
        sg[i] = sc[i][4]; // sin gama 
        cg[i] = sc[i][5]; // cos gama 
        for (int j=0; j<4; j++){
            (*B)[i][j] = 0.0; 
            if(j!=3)
                (*A)[i][j] = 0.0;    
        }
        (*A)[i][0] = x - ax[i] + r*cg[i] + b*(sa[i]*st[i] - ca[i]*cb[i]*ct[i]);
        (*A)[i][1] = y - ay[i] + r*sg[i] - b*(ca[i]*st[i] + cb[i]*sa[i]*ct[i]);
        (*A)[i][2] = z + b*sb[i]*ct[i];

        (*B)[i][i] = (z + b*sb[i]*ct[i])*(b*ca[i]*sb[i]*(ca[i]*st[i] + cb[i]*sa[i]*ct[i]) + b*sa[i]*sb[i]*(sa[i]*st[i] \
                    - ca[i]*cb[i]*ct[i])) - (b*sa[i]*ct[i]*powf(sb[i],2.0) + b*cb[i]*(ca[i]*st[i] \
                    + cb[i]*sa[i]*ct[i]))*(x - ax[i] + r*cg[i] + b*(sa[i]*st[i] - ca[i]*cb[i]*ct[i])) \
                    + (b*cb[i]*(sa[i]*st[i] - ca[i]*cb[i]*ct[i]) - b*ca[i]*powf(sb[i],2.0)*ct[i])*(ay[i] - y - r*sg[i]\
                     + b*(ca[i]*st[i] + cb[i]*sa[i]*ct[i]));
    }
}

void ragnarJpassive(float params[4][8], double  theta[4], double passive[4][2],
                    double sc[4][6], double sct[8], double scez[4][4], 
                    double (*J)[3][3], int leg_num)
{   
    // Computer the inverse kinematics
    // access them as (*a)[0][0] .. (*a)[2][2]
    // alocate solution array
    // sc[][] contains sin and cosine of alpha beta gamma 
    // sc[i][] = {sin alpha, cos alpha, sin beta, cos beta, sin gama, cos gama}
    // offset is for the offset of the motors 
    // full is to acquire the full inverse kinematics including the passive
    // joints, if false or default then it is not computed 
    // scez is defines as 
    // scez[i][] = {sin zeta, cos zeta, sin eta, cos eta }

    double sa = sc[leg_num][0];
    double ca = sc[leg_num][1];
    double sb = sc[leg_num][2];
    double cb = sc[leg_num][3];
    double b = double(params[leg_num][4]);
    double l = double(params[leg_num][5]);
    double st = sct[(2*leg_num)];
    double ct = sct[(2*leg_num)+1];
    double sz = scez[leg_num][0];
    double cz = scez[leg_num][1];
    double se = scez[leg_num][2];
    double ce = scez[leg_num][3];

    (*J)[0][0] = ce*l*sz*(sa*st - ca*cb*ct) - (ct*sa + ca*cb*st)*(b + ce*cz*l);    
    (*J)[0][1] = ce*l*(sz*(sa*st - ca*cb*ct) - cz*(ct*sa + ca*cb*st));
    (*J)[0][2] = l*(cz*se*(sa*st - ca*cb*ct) + se*sz*(ct*sa + ca*cb*st) - ca*ce*sb);
    (*J)[1][0] = (ca*ct - cb*sa*st)*(b + ce*cz*l) - ce*l*sz*(ca*st + cb*ct*sa);
    (*J)[1][1] = -ce*l*(sz*(ca*st + cb*ct*sa) - cz*(ca*ct - cb*sa*st));
    (*J)[1][2] = -l*(cz*se*(ca*st + cb*ct*sa) + se*sz*(ca*ct - cb*sa*st) + ce*sa*sb);
    (*J)[2][0] = sb*st*(b + ce*cz*l) + ce*ct*l*sb*sz;
    // (*J)[2][1] = l*sin(theta + zeta)*ce*sb; Replace for trigonometric identity
    (*J)[2][1] = ce*l*sb*(ct*sz + cz*st); 
    (*J)[2][2] = -l*(cb*ce - ct*cz*sb*se + sb*se*st*sz);
}

void ragnarJpassivef(float params[4][8], float  theta[4], float passive[4][2],
                    float sc[4][6], float sct[8], float scez[4][4], 
                    float (*J)[3][3], int leg_num)
{   
    // Computer the inverse kinematics
    // access them as (*a)[0][0] .. (*a)[2][2]
    // alocate solution array
    // sc[][] contains sin and cosine of alpha beta gamma 
    // sc[i][] = {sin alpha, cos alpha, sin beta, cos beta, sin gama, cos gama}
    // offset is for the offset of the motors 
    // full is to acquire the full inverse kinematics including the passive
    // joints, if false or default then it is not computed 
    // scez is defines as 
    // scez[i][] = {sin zeta, cos zeta, sin eta, cos eta }

    float sa = sc[leg_num][0];
    float ca = sc[leg_num][1];
    float sb = sc[leg_num][2];
    float cb = sc[leg_num][3];
    float b = float(params[leg_num][4]);
    float l = float(params[leg_num][5]);
    float st = sct[(2*leg_num)];
    float ct = sct[(2*leg_num)+1];
    float sz = scez[leg_num][0];
    float cz = scez[leg_num][1];
    float se = scez[leg_num][2];
    float ce = scez[leg_num][3];

    (*J)[0][0] = ce*l*sz*(sa*st - ca*cb*ct) - (ct*sa + ca*cb*st)*(b + ce*cz*l);    
    (*J)[0][1] = ce*l*(sz*(sa*st - ca*cb*ct) - cz*(ct*sa + ca*cb*st));
    (*J)[0][2] = l*(cz*se*(sa*st - ca*cb*ct) + se*sz*(ct*sa + ca*cb*st) - ca*ce*sb);
    (*J)[1][0] = (ca*ct - cb*sa*st)*(b + ce*cz*l) - ce*l*sz*(ca*st + cb*ct*sa);
    (*J)[1][1] = -ce*l*(sz*(ca*st + cb*ct*sa) - cz*(ca*ct - cb*sa*st));
    (*J)[1][2] = -l*(cz*se*(ca*st + cb*ct*sa) + se*sz*(ca*ct - cb*sa*st) + ce*sa*sb);
    (*J)[2][0] = sb*st*(b + ce*cz*l) + ce*ct*l*sb*sz;
    // (*J)[2][1] = l*sin(theta + zeta)*ce*sb; Replace for trigonometric identity
    (*J)[2][1] = ce*l*sb*(ct*sz + cz*st); 
    (*J)[2][2] = -l*(cb*ce - ct*cz*sb*se + sb*se*st*sz);
}

void ragnariJpassive(float params[4][8], double  theta[4], double passive[4][2],
                    double sc[4][6], double sct[8], double scez[4][4], double (*iJ)[3][3], 
                    int leg_num)
{   
    // Computer the inverse kinematics
    // access them as (*a)[0][0] .. (*a)[2][2]
    // alocate solution array
    // sc[][] contains sin and cosine of alpha beta gamma 
    // sc[i][] = {sin alpha, cos alpha, sin beta, cos beta, sin gama, cos gama}
    // offset is for the offset of the motors 
    // full is to acquire the full inverse kinematics including the passive
    // joints, if false or default then it is not computed 
    // scez is defines as 
    // scez[i][] = {sin eta, cos eta, sin zeta, cos zeta }

    double sa = sc[leg_num][0];
    double ca = sc[leg_num][1];
    double sb = sc[leg_num][2];
    double cb = sc[leg_num][3];
    double b = double(params[leg_num][4]);
    double l = double(params[leg_num][5]);
    double st = sct[(2*leg_num)];
    double ct = sct[(2*leg_num)+1];
    double sz = scez[leg_num][0];
    double cz = scez[leg_num][1];
    double se = scez[leg_num][2];
    double ce = scez[leg_num][3];

    //(*iJ)[0][0] = - (sa*ct + ca*cb*st)/b - (ca*sb*tan(eta) + sa*cz*st - ca*cb*ct*cz)/(b*sz);
    (*iJ)[0][0] = - (sa*ct + ca*cb*st)/b - (ca*sb*(se/ce) + sa*cz*st - ca*cb*ct*cz)/(b*sz);
    (*iJ)[0][1] = (ca*ct - cb*sa*st)/b + (ca*cz*st - sa*sb*(se/ce) + cb*sa*ct*cz)/(b*sz);
    (*iJ)[0][2] = (sb*st)/b - (cb*(se/ce) + sb*ct*cz)/(b*sz);
    (*iJ)[1][0] = (l*sa*ct*sz + l*sa*cz*st - l*ca*cb*ct*cz + l*ca*cb*st*sz)/(b*l*sz) + (ce*(b*sa*st - b*ca*cb*ct + l*ca*sb*se) + b*ca*sb*se*cz)/(b*l*powf(ce,2.0)*sz);
    (*iJ)[1][1] = - (l*ca*ct*sz + l*ca*cz*st + l*cb*sa*ct*cz - l*cb*sa*st*sz)/(b*l*sz) - (ce*(b*ca*st + b*cb*sa*ct - l*sa*sb*se) - b*sa*sb*se*cz)/(b*l*powf(ce,2.0)*sz);
    (*iJ)[1][2] = (l*sb*ct*cz - l*sb*st*sz)/(b*l*sz) + (ce*(l*cb*se + b*sb*ct) + b*cb*se*cz)/(b*l*powf(ce,2.0)*sz);
    (*iJ)[2][0] = -(ca*sb)/(l*ce);
    (*iJ)[2][1] = -(sa*sb)/(l*ce);
    (*iJ)[2][2] = -cb/(l*ce);
}

bool ragnarpassiveJ(float A[4][3], float B[4][4], float iJp[3][3], float (*pJ)[2][3], int option)
{
    // A is the forward jacobian 
    // B is the inverse Jacobian
    // iJp is the inverse Jacobian where it is included the passive joints 

    // theta = theta(0:2) or theta = theta(1:3)
    // degrade A to be square by taking a theta out 
    // degrade B to be square 
    // compute the Jacobian to get dx 
    // dx = inv(A)*B*theta 
    // get the Jacobian for a leg and compute the passive joints velocities 
    // dq = inv(J)*dx 
    // dq = inv(J)*inv(A)*B*theta 
    // dq(1:2) = pJ*theta  

    float A1_1, A1_2, A1_3, A2_1, A2_2, A2_3, A3_1, A3_2, A3_3; 
    float B1_1, B1_2, B1_3, B2_1, B2_2, B2_3, B3_1, B3_2, B3_3; 
    float ij1_1, ij1_2, ij1_3, ij2_1, ij2_2, ij2_3, ij3_1, ij3_2, ij3_3; 
    bool returned; 
    if(option == 0){
        A1_1 = A[0][0];
        A1_2 = A[0][1];
        A1_3 = A[0][2];

        A2_1 = A[1][0];
        A2_2 = A[1][1];
        A2_3 = A[1][2];
        
        A3_1 = A[2][0];
        A3_2 = A[2][1];
        A3_3 = A[2][2];
        
        B1_1 = B[0][0];
        B1_2 = B[0][1];
        B1_3 = B[0][2];
        
        B2_1 = B[1][0];
        B2_2 = B[1][1];
        B2_3 = B[1][2];
        
        B3_1 = B[2][0];
        B3_2 = B[2][1];
        B3_3 = B[2][2];
    }
    else {
        A1_1 = A[1][0];
        A1_2 = A[1][1];
        A1_3 = A[1][2];

        A2_1 = A[2][0];
        A2_2 = A[2][1];
        A2_3 = A[2][2];
        
        A3_1 = A[3][0];
        A3_2 = A[3][1];
        A3_3 = A[3][2];
        
        B1_1 = B[1][1];
        B1_2 = B[1][2];
        B1_3 = B[1][3];
        
        B2_1 = B[2][1];
        B2_2 = B[2][2];
        B2_3 = B[2][3];
        
        B3_1 = B[3][1];
        B3_2 = B[3][2];
        B3_3 = B[3][3];
    }
    ij1_1 = iJp[0][0];
    ij1_2 = iJp[0][1];
    ij1_3 = iJp[0][2];
    
    ij2_1 = iJp[1][0];
    ij2_2 = iJp[1][1];
    ij2_3 = iJp[1][2];
    
    ij3_1 = iJp[2][0];
    ij3_2 = iJp[2][1];
    ij3_3 = iJp[2][2];    
 
    float detA = (A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1);
    if(detA == 0)
        returned = false;
    else {
        returned = true; 
        (*pJ)[0][0] = (A1_1*A2_2*B3_1*ij2_3 - A1_1*A2_3*B3_1*ij2_2 - A1_1*A3_2*B2_1*ij2_3 + A1_1*A3_3*B2_1*ij2_2 - A1_2*A2_1*B3_1*ij2_3 + A1_2*A2_3*B3_1*ij2_1 + 
                A1_2*A3_1*B2_1*ij2_3 - A1_2*A3_3*B2_1*ij2_1 + A1_3*A2_1*B3_1*ij2_2 - A1_3*A2_2*B3_1*ij2_1 - A1_3*A3_1*B2_1*ij2_2 + A1_3*A3_2*B2_1*ij2_1 + 
                A2_1*A3_2*B1_1*ij2_3 - A2_1*A3_3*B1_1*ij2_2 - A2_2*A3_1*B1_1*ij2_3 + A2_2*A3_3*B1_1*ij2_1 + A2_3*A3_1*B1_1*ij2_2 - A2_3*A3_2*B1_1*ij2_1)/detA;
                
        (*pJ)[0][1] = (A1_1*A2_2*B3_2*ij2_3 - A1_1*A2_3*B3_2*ij2_2 - A1_1*A3_2*B2_2*ij2_3 + A1_1*A3_3*B2_2*ij2_2 - A1_2*A2_1*B3_2*ij2_3 + A1_2*A2_3*B3_2*ij2_1 + 
                A1_2*A3_1*B2_2*ij2_3 - A1_2*A3_3*B2_2*ij2_1 + A1_3*A2_1*B3_2*ij2_2 - A1_3*A2_2*B3_2*ij2_1 - A1_3*A3_1*B2_2*ij2_2 + A1_3*A3_2*B2_2*ij2_1 + 
                A2_1*A3_2*B1_2*ij2_3 - A2_1*A3_3*B1_2*ij2_2 - A2_2*A3_1*B1_2*ij2_3 + A2_2*A3_3*B1_2*ij2_1 + A2_3*A3_1*B1_2*ij2_2 - A2_3*A3_2*B1_2*ij2_1)/detA;
                
        (*pJ)[0][2] = (A1_1*A2_2*B3_3*ij2_3 - A1_1*A2_3*B3_3*ij2_2 - A1_1*A3_2*B2_3*ij2_3 + A1_1*A3_3*B2_3*ij2_2 - A1_2*A2_1*B3_3*ij2_3 + A1_2*A2_3*B3_3*ij2_1 + 
                A1_2*A3_1*B2_3*ij2_3 - A1_2*A3_3*B2_3*ij2_1 + A1_3*A2_1*B3_3*ij2_2 - A1_3*A2_2*B3_3*ij2_1 - A1_3*A3_1*B2_3*ij2_2 + A1_3*A3_2*B2_3*ij2_1 + 
                A2_1*A3_2*B1_3*ij2_3 - A2_1*A3_3*B1_3*ij2_2 - A2_2*A3_1*B1_3*ij2_3 + A2_2*A3_3*B1_3*ij2_1 + A2_3*A3_1*B1_3*ij2_2 - A2_3*A3_2*B1_3*ij2_1)/detA;
                
        (*pJ)[1][0] = (A1_1*A2_2*B3_1*ij3_3 - A1_1*A2_3*B3_1*ij3_2 - A1_1*A3_2*B2_1*ij3_3 + A1_1*A3_3*B2_1*ij3_2 - A1_2*A2_1*B3_1*ij3_3 + A1_2*A2_3*B3_1*ij3_1 + 
                A1_2*A3_1*B2_1*ij3_3 - A1_2*A3_3*B2_1*ij3_1 + A1_3*A2_1*B3_1*ij3_2 - A1_3*A2_2*B3_1*ij3_1 - A1_3*A3_1*B2_1*ij3_2 + A1_3*A3_2*B2_1*ij3_1 + 
                A2_1*A3_2*B1_1*ij3_3 - A2_1*A3_3*B1_1*ij3_2 - A2_2*A3_1*B1_1*ij3_3 + A2_2*A3_3*B1_1*ij3_1 + A2_3*A3_1*B1_1*ij3_2 - A2_3*A3_2*B1_1*ij3_1)/detA;
                
        (*pJ)[1][1] = (A1_1*A2_2*B3_2*ij3_3 - A1_1*A2_3*B3_2*ij3_2 - A1_1*A3_2*B2_2*ij3_3 + A1_1*A3_3*B2_2*ij3_2 - A1_2*A2_1*B3_2*ij3_3 + A1_2*A2_3*B3_2*ij3_1 + 
                A1_2*A3_1*B2_2*ij3_3 - A1_2*A3_3*B2_2*ij3_1 + A1_3*A2_1*B3_2*ij3_2 - A1_3*A2_2*B3_2*ij3_1 - A1_3*A3_1*B2_2*ij3_2 + A1_3*A3_2*B2_2*ij3_1 + 
                A2_1*A3_2*B1_2*ij3_3 - A2_1*A3_3*B1_2*ij3_2 - A2_2*A3_1*B1_2*ij3_3 + A2_2*A3_3*B1_2*ij3_1 + A2_3*A3_1*B1_2*ij3_2 - A2_3*A3_2*B1_2*ij3_1)/detA;
                
        (*pJ)[1][2] = (A1_1*A2_2*B3_3*ij3_3 - A1_1*A2_3*B3_3*ij3_2 - A1_1*A3_2*B2_3*ij3_3 + A1_1*A3_3*B2_3*ij3_2 - A1_2*A2_1*B3_3*ij3_3 + A1_2*A2_3*B3_3*ij3_1 + 
                A1_2*A3_1*B2_3*ij3_3 - A1_2*A3_3*B2_3*ij3_1 + A1_3*A2_1*B3_3*ij3_2 - A1_3*A2_2*B3_3*ij3_1 - A1_3*A3_1*B2_3*ij3_2 + A1_3*A3_2*B2_3*ij3_1 + 
                A2_1*A3_2*B1_3*ij3_3 - A2_1*A3_3*B1_3*ij3_2 - A2_2*A3_1*B1_3*ij3_3 + A2_2*A3_3*B1_3*ij3_1 + A2_3*A3_1*B1_3*ij3_2 - A2_3*A3_2*B1_3*ij3_1)/detA;
    }
    return returned; 
}
bool ragnarpassiveJ3(float A[4][3], float B[4][4], float iJp[3][3], float (*pJ)[3][3], int option)
{
    // A is the forward jacobian 
    // B is the inverse Jacobian
    // iJp is the inverse Jacobian where it is included the passive joints 

    // theta = theta(0:2) or theta = theta(1:3)
    // degrade A to be square by taking a theta out 
    // degrade B to be square 
    // compute the Jacobian to get dx 
    // dx = inv(A)*B*theta 
    // get the Jacobian for a leg and compute the passive joints velocities 
    // dq = inv(J)*dx 
    // dq = inv(J)*inv(A)*B*theta 
    // dq(1:2) = pJ*theta  

    float A1_1, A1_2, A1_3, A2_1, A2_2, A2_3, A3_1, A3_2, A3_3; 
    float B1_1, B1_2, B1_3, B2_1, B2_2, B2_3, B3_1, B3_2, B3_3; 
    float ij1_1, ij1_2, ij1_3, ij2_1, ij2_2, ij2_3, ij3_1, ij3_2, ij3_3; 
    bool returned; 
    if(option == 0){
        A1_1 = A[0][0];
        A1_2 = A[0][1];
        A1_3 = A[0][2];

        A2_1 = A[1][0];
        A2_2 = A[1][1];
        A2_3 = A[1][2];
        
        A3_1 = A[2][0];
        A3_2 = A[2][1];
        A3_3 = A[2][2];
        
        B1_1 = B[0][0];
        B1_2 = B[0][1];
        B1_3 = B[0][2];
        
        B2_1 = B[1][0];
        B2_2 = B[1][1];
        B2_3 = B[1][2];
        
        B3_1 = B[2][0];
        B3_2 = B[2][1];
        B3_3 = B[2][2];
    }
    else {
        A1_1 = A[1][0];
        A1_2 = A[1][1];
        A1_3 = A[1][2];

        A2_1 = A[2][0];
        A2_2 = A[2][1];
        A2_3 = A[2][2];
        
        A3_1 = A[3][0];
        A3_2 = A[3][1];
        A3_3 = A[3][2];
        
        B1_1 = B[1][1];
        B1_2 = B[1][2];
        B1_3 = B[1][3];
        
        B2_1 = B[2][1];
        B2_2 = B[2][2];
        B2_3 = B[2][3];
        
        B3_1 = B[3][1];
        B3_2 = B[3][2];
        B3_3 = B[3][3];
    }
    ij1_1 = iJp[0][0];
    ij1_2 = iJp[0][1];
    ij1_3 = iJp[0][2];
    
    ij2_1 = iJp[1][0];
    ij2_2 = iJp[1][1];
    ij2_3 = iJp[1][2];
    
    ij3_1 = iJp[2][0];
    ij3_2 = iJp[2][1];
    ij3_3 = iJp[2][2];    
 
    float detA = (A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1);
    if(detA == 0)
        returned = false;
    else {
        returned = true; 
        (*pJ)[0][0] = (A1_1*A2_2*B3_1*ij1_3 - A1_1*A2_3*B3_1*ij1_2 - A1_1*A3_2*B2_1*ij1_3 + A1_1*A3_3*B2_1*ij1_2 - A1_2*A2_1*B3_1*ij1_3 + A1_2*A2_3*B3_1*ij1_1 + A1_2*A3_1*B2_1*ij1_3 - A1_2*A3_3*B2_1*ij1_1 + A1_3*A2_1*B3_1*ij1_2 - A1_3*A2_2*B3_1*ij1_1 - A1_3*A3_1*B2_1*ij1_2 + A1_3*A3_2*B2_1*ij1_1 + A2_1*A3_2*B1_1*ij1_3 - A2_1*A3_3*B1_1*ij1_2 - A2_2*A3_1*B1_1*ij1_3 + A2_2*A3_3*B1_1*ij1_1 + A2_3*A3_1*B1_1*ij1_2 - A2_3*A3_2*B1_1*ij1_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1);
        (*pJ)[0][1] =  (A1_1*A2_2*B3_2*ij1_3 - A1_1*A2_3*B3_2*ij1_2 - A1_1*A3_2*B2_2*ij1_3 + A1_1*A3_3*B2_2*ij1_2 - A1_2*A2_1*B3_2*ij1_3 + A1_2*A2_3*B3_2*ij1_1 + A1_2*A3_1*B2_2*ij1_3 - A1_2*A3_3*B2_2*ij1_1 + A1_3*A2_1*B3_2*ij1_2 - A1_3*A2_2*B3_2*ij1_1 - A1_3*A3_1*B2_2*ij1_2 + A1_3*A3_2*B2_2*ij1_1 + A2_1*A3_2*B1_2*ij1_3 - A2_1*A3_3*B1_2*ij1_2 - A2_2*A3_1*B1_2*ij1_3 + A2_2*A3_3*B1_2*ij1_1 + A2_3*A3_1*B1_2*ij1_2 - A2_3*A3_2*B1_2*ij1_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1);
        (*pJ)[0][2] =   (A1_1*A2_2*B3_3*ij1_3 - A1_1*A2_3*B3_3*ij1_2 - A1_1*A3_2*B2_3*ij1_3 + A1_1*A3_3*B2_3*ij1_2 - A1_2*A2_1*B3_3*ij1_3 + A1_2*A2_3*B3_3*ij1_1 + A1_2*A3_1*B2_3*ij1_3 - A1_2*A3_3*B2_3*ij1_1 + A1_3*A2_1*B3_3*ij1_2 - A1_3*A2_2*B3_3*ij1_1 - A1_3*A3_1*B2_3*ij1_2 + A1_3*A3_2*B2_3*ij1_1 + A2_1*A3_2*B1_3*ij1_3 - A2_1*A3_3*B1_3*ij1_2 - A2_2*A3_1*B1_3*ij1_3 + A2_2*A3_3*B1_3*ij1_1 + A2_3*A3_1*B1_3*ij1_2 - A2_3*A3_2*B1_3*ij1_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1);
        (*pJ)[1][0] = (A1_1*A2_2*B3_1*ij2_3 - A1_1*A2_3*B3_1*ij2_2 - A1_1*A3_2*B2_1*ij2_3 + A1_1*A3_3*B2_1*ij2_2 - A1_2*A2_1*B3_1*ij2_3 + A1_2*A2_3*B3_1*ij2_1 + A1_2*A3_1*B2_1*ij2_3 - A1_2*A3_3*B2_1*ij2_1 + A1_3*A2_1*B3_1*ij2_2 - A1_3*A2_2*B3_1*ij2_1 - A1_3*A3_1*B2_1*ij2_2 + A1_3*A3_2*B2_1*ij2_1 + A2_1*A3_2*B1_1*ij2_3 - A2_1*A3_3*B1_1*ij2_2 - A2_2*A3_1*B1_1*ij2_3 + A2_2*A3_3*B1_1*ij2_1 + A2_3*A3_1*B1_1*ij2_2 - A2_3*A3_2*B1_1*ij2_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1);
        (*pJ)[1][1] = (A1_1*A2_2*B3_2*ij2_3 - A1_1*A2_3*B3_2*ij2_2 - A1_1*A3_2*B2_2*ij2_3 + A1_1*A3_3*B2_2*ij2_2 - A1_2*A2_1*B3_2*ij2_3 + A1_2*A2_3*B3_2*ij2_1 + A1_2*A3_1*B2_2*ij2_3 - A1_2*A3_3*B2_2*ij2_1 + A1_3*A2_1*B3_2*ij2_2 - A1_3*A2_2*B3_2*ij2_1 - A1_3*A3_1*B2_2*ij2_2 + A1_3*A3_2*B2_2*ij2_1 + A2_1*A3_2*B1_2*ij2_3 - A2_1*A3_3*B1_2*ij2_2 - A2_2*A3_1*B1_2*ij2_3 + A2_2*A3_3*B1_2*ij2_1 + A2_3*A3_1*B1_2*ij2_2 - A2_3*A3_2*B1_2*ij2_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1);
        (*pJ)[1][2] = (A1_1*A2_2*B3_3*ij2_3 - A1_1*A2_3*B3_3*ij2_2 - A1_1*A3_2*B2_3*ij2_3 + A1_1*A3_3*B2_3*ij2_2 - A1_2*A2_1*B3_3*ij2_3 + A1_2*A2_3*B3_3*ij2_1 + A1_2*A3_1*B2_3*ij2_3 - A1_2*A3_3*B2_3*ij2_1 + A1_3*A2_1*B3_3*ij2_2 - A1_3*A2_2*B3_3*ij2_1 - A1_3*A3_1*B2_3*ij2_2 + A1_3*A3_2*B2_3*ij2_1 + A2_1*A3_2*B1_3*ij2_3 - A2_1*A3_3*B1_3*ij2_2 - A2_2*A3_1*B1_3*ij2_3 + A2_2*A3_3*B1_3*ij2_1 + A2_3*A3_1*B1_3*ij2_2 - A2_3*A3_2*B1_3*ij2_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1);
        (*pJ)[2][0] = (A1_1*A2_2*B3_1*ij3_3 - A1_1*A2_3*B3_1*ij3_2 - A1_1*A3_2*B2_1*ij3_3 + A1_1*A3_3*B2_1*ij3_2 - A1_2*A2_1*B3_1*ij3_3 + A1_2*A2_3*B3_1*ij3_1 + A1_2*A3_1*B2_1*ij3_3 - A1_2*A3_3*B2_1*ij3_1 + A1_3*A2_1*B3_1*ij3_2 - A1_3*A2_2*B3_1*ij3_1 - A1_3*A3_1*B2_1*ij3_2 + A1_3*A3_2*B2_1*ij3_1 + A2_1*A3_2*B1_1*ij3_3 - A2_1*A3_3*B1_1*ij3_2 - A2_2*A3_1*B1_1*ij3_3 + A2_2*A3_3*B1_1*ij3_1 + A2_3*A3_1*B1_1*ij3_2 - A2_3*A3_2*B1_1*ij3_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1);
        (*pJ)[2][1] = (A1_1*A2_2*B3_2*ij3_3 - A1_1*A2_3*B3_2*ij3_2 - A1_1*A3_2*B2_2*ij3_3 + A1_1*A3_3*B2_2*ij3_2 - A1_2*A2_1*B3_2*ij3_3 + A1_2*A2_3*B3_2*ij3_1 + A1_2*A3_1*B2_2*ij3_3 - A1_2*A3_3*B2_2*ij3_1 + A1_3*A2_1*B3_2*ij3_2 - A1_3*A2_2*B3_2*ij3_1 - A1_3*A3_1*B2_2*ij3_2 + A1_3*A3_2*B2_2*ij3_1 + A2_1*A3_2*B1_2*ij3_3 - A2_1*A3_3*B1_2*ij3_2 - A2_2*A3_1*B1_2*ij3_3 + A2_2*A3_3*B1_2*ij3_1 + A2_3*A3_1*B1_2*ij3_2 - A2_3*A3_2*B1_2*ij3_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1);
        (*pJ)[2][2] = (A1_1*A2_2*B3_3*ij3_3 - A1_1*A2_3*B3_3*ij3_2 - A1_1*A3_2*B2_3*ij3_3 + A1_1*A3_3*B2_3*ij3_2 - A1_2*A2_1*B3_3*ij3_3 + A1_2*A2_3*B3_3*ij3_1 + A1_2*A3_1*B2_3*ij3_3 - A1_2*A3_3*B2_3*ij3_1 + A1_3*A2_1*B3_3*ij3_2 - A1_3*A2_2*B3_3*ij3_1 - A1_3*A3_1*B2_3*ij3_2 + A1_3*A3_2*B2_3*ij3_1 + A2_1*A3_2*B1_3*ij3_3 - A2_1*A3_3*B1_3*ij3_2 - A2_2*A3_1*B1_3*ij3_3 + A2_2*A3_3*B1_3*ij3_1 + A2_3*A3_1*B1_3*ij3_2 - A2_3*A3_2*B1_3*ij3_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1);
    }
    return returned; 
}

void ragnariJpassivef(
    float params[4][8], float  theta[4], float passive[4][2], float sc[4][6], 
    float sct[8], float scez[4][4], float (*iJ)[3][3], int leg_num)
{   
    // Computer the inverse kinematics
    // access them as (*a)[0][0] .. (*a)[2][2]
    // alocate solution array
    // sc[][] contains sin and cosine of alpha beta gamma 
    // sc[i][] = {sin alpha, cos alpha, sin beta, cos beta, sin gama, cos gama}
    // offset is for the offset of the motors 
    // full is to acquire the full inverse kinematics including the passive
    // joints, if false or default then it is not computed 
    // scez is defines as 
    // scez[i][] = {sin eta, cos eta, sin zeta, cos zeta }

    float sa = sc[leg_num][0];
    float ca = sc[leg_num][1];
    float sb = sc[leg_num][2];
    float cb = sc[leg_num][3];
    float b = float(params[leg_num][4]);
    float l = float(params[leg_num][5]);
    float st = sct[(2*leg_num)];
    float ct = sct[(2*leg_num)+1];
    float sz = scez[leg_num][0];
    float cz = scez[leg_num][1];
    float se = scez[leg_num][2];
    float ce = scez[leg_num][3];

    //(*iJ)[0][0] = - (sa*ct + ca*cb*st)/b - (ca*sb*tan(eta) + sa*cz*st - ca*cb*ct*cz)/(b*sz);
    (*iJ)[0][0] = - (sa*ct + ca*cb*st)/b - (ca*sb*(se/ce) + sa*cz*st - ca*cb*ct*cz)/(b*sz);
    (*iJ)[0][1] = (ca*ct - cb*sa*st)/b + (ca*cz*st - sa*sb*(se/ce) + cb*sa*ct*cz)/(b*sz);
    (*iJ)[0][2] = (sb*st)/b - (cb*(se/ce) + sb*ct*cz)/(b*sz);
    (*iJ)[1][0] = (l*sa*ct*sz + l*sa*cz*st - l*ca*cb*ct*cz + l*ca*cb*st*sz)/(b*l*sz) + (ce*(b*sa*st - b*ca*cb*ct + l*ca*sb*se) + b*ca*sb*se*cz)/(b*l*powf(ce,2.0)*sz);
    (*iJ)[1][1] = - (l*ca*ct*sz + l*ca*cz*st + l*cb*sa*ct*cz - l*cb*sa*st*sz)/(b*l*sz) - (ce*(b*ca*st + b*cb*sa*ct - l*sa*sb*se) - b*sa*sb*se*cz)/(b*l*powf(ce,2.0)*sz);
    (*iJ)[1][2] = (l*sb*ct*cz - l*sb*st*sz)/(b*l*sz) + (ce*(l*cb*se + b*sb*ct) + b*cb*se*cz)/(b*l*powf(ce,2.0)*sz);
    (*iJ)[2][0] = -(ca*sb)/(l*ce);
    (*iJ)[2][1] = -(sa*sb)/(l*ce);
    (*iJ)[2][2] = -cb/(l*ce);
}

double detA33(double A[3][3]){
    //Computes the determinant to see if the matrix is invertible 
    double det = 0.0; 
    det = A[1][1]*A[2][2]*A[3][3] - A[1][1]*A[2][3]*A[3][2] \
         -A[1][2]*A[2][1]*A[3][3] + A[1][2]*A[2][3]*A[3][1] \
         +A[1][3]*A[2][1]*A[3][2] - A[1][3]*A[2][2]*A[3][1]; 
    return det;
}
float detA33f(float A[3][3]){
    //Computes the determinant to see if the matrix is invertible 
    float det = 0.0; 
    det = A[1][1]*A[2][2]*A[3][3] - A[1][1]*A[2][3]*A[3][2] \
         -A[1][2]*A[2][1]*A[3][3] + A[1][2]*A[2][3]*A[3][1] \
         +A[1][3]*A[2][1]*A[3][2] - A[1][3]*A[2][2]*A[3][1]; 
    return det;
}

bool Ainv(double A[3][3], double (*iA)[3][3]){
    // Computes the inverse of a 3 by 3 matrix, if it exists, if not it does not
    // do anything. Reported over a returned bool. 
    double A1_1, A1_2, A1_3, A2_1, A2_2, A2_3, A3_1, A3_2, A3_3;
    A1_1 = A[0][0];
    A1_2 = A[0][1];
    A1_3 = A[0][2];
    A2_1 = A[1][0];
    A2_2 = A[1][1];
    A2_3 = A[1][2];
    A3_1 = A[2][0];
    A3_2 = A[2][1];
    A3_3 = A[2][2];

    bool returned; 

    if (detA33(A) != 0.0){
        double iiA[3][3] = {{  (A2_2*A3_3 - A2_3*A3_2)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1), -(A1_2*A3_3 - A1_3*A3_2)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1),  (A1_2*A2_3 - A1_3*A2_2)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1)}, \
                           { -(A2_1*A3_3 - A2_3*A3_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1),  (A1_1*A3_3 - A1_3*A3_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1), -(A1_1*A2_3 - A1_3*A2_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1)}, \
                           {  (A2_1*A3_2 - A2_2*A3_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1), -(A1_1*A3_2 - A1_2*A3_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1),  (A1_1*A2_2 - A1_2*A2_1)/(A1_1*A2_2*A3_3 - A1_1*A2_3*A3_2 - A1_2*A2_1*A3_3 + A1_2*A2_3*A3_1 + A1_3*A2_1*A3_2 - A1_3*A2_2*A3_1)}};
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                (*iA)[i][j] = iiA[i][j];
        returned = true; 
    }
    else
        returned = false;
    return returned;  
}

bool Ainvf(float A[3][3], float (*iA)[3][3]){
    // Computes the inverse of a 3 by 3 matrix, if it exists, if not it does not
    // do anything. Reported over a returned bool. 
    float A1_1, A1_2, A1_3, A2_1, A2_2, A2_3, A3_1, A3_2, A3_3;
    A1_1 = A[0][0];
    A1_2 = A[0][1];
    A1_3 = A[0][2];
    
    A2_1 = A[1][0];
    A2_2 = A[1][1];
    A2_3 = A[1][2];
    
    A3_1 = A[2][0];
    A3_2 = A[2][1];
    A3_3 = A[2][2];

    bool returned; 

    if (detA33f(A) != 0.0){   
        float deto = ((A1_1*A2_2*A3_3) -( A1_1*A2_3*A3_2) - (A1_2*A2_1*A3_3) + (A1_2*A2_3*A3_1) + (A1_3*A2_1*A3_2) - (A1_3*A2_2*A3_1)); 
        float undeto = 1.0/deto; 
        // Serial.println(deto,8); 
        float iiA[3][3] = {{  (A2_2*A3_3 - A2_3*A3_2)*undeto, -(A1_2*A3_3 - A1_3*A3_2)*undeto,  (A1_2*A2_3 - A1_3*A2_2)*undeto}, \
                           { -(A2_1*A3_3 - A2_3*A3_1)*undeto,  (A1_1*A3_3 - A1_3*A3_1)*undeto, -(A1_1*A2_3 - A1_3*A2_1)*undeto}, \
                           {  (A2_1*A3_2 - A2_2*A3_1)*undeto, -(A1_1*A3_2 - A1_2*A3_1)*undeto,  (A1_1*A2_2 - A1_2*A2_1)*undeto}};
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                (*iA)[i][j] = iiA[i][j];
        returned = true; 
    }
    else
        returned = false;
    return returned;  
}

bool Ainvd(float A[3][3], float (*iA)[3][3]){
    // Computes the inverse of a 3 by 3 matrix, if it exists, if not it does not
    // do anything. Reported over a returned bool. 
    double A1_1, A1_2, A1_3, A2_1, A2_2, A2_3, A3_1, A3_2, A3_3;
    A1_1 = double(A[0][0]);
    A1_2 = double(A[0][1]);
    A1_3 = double(A[0][2]);
    
    A2_1 = double(A[1][0]);
    A2_2 = double(A[1][1]);
    A2_3 = double(A[1][2]);
    
    A3_1 = double(A[2][0]);
    A3_2 = double(A[2][1]);
    A3_3 = double(A[2][2]);

    bool returned; 

    if (detA33f(A) != 0.0){   
        double deto = ((A1_1*A2_2*A3_3) -( A1_1*A2_3*A3_2) - (A1_2*A2_1*A3_3) + (A1_2*A2_3*A3_1) + (A1_3*A2_1*A3_2) - (A1_3*A2_2*A3_1)); 
        //double undeto = 1.0/deto; 
        //Serial.println(deto,8); 
        double iiA[3][3] = {{  (A2_2*A3_3 - A2_3*A3_2)/deto, -(A1_2*A3_3 - A1_3*A3_2)/deto,  (A1_2*A2_3 - A1_3*A2_2)/deto}, \
                           { -(A2_1*A3_3 - A2_3*A3_1)/deto,  (A1_1*A3_3 - A1_3*A3_1)/deto, -(A1_1*A2_3 - A1_3*A2_1)/deto}, \
                           {  (A2_1*A3_2 - A2_2*A3_1)/deto, -(A1_1*A3_2 - A1_2*A3_1)/deto,  (A1_1*A2_2 - A1_2*A2_1)/deto}};
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                (*iA)[i][j] = float(iiA[i][j]);
        returned = true; 
    }
    else
        returned = false;
    return returned;  
}

void Atra(double A[3][3], double (*At)[3][3]){
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            (*At)[i][j] = A[j][i];
} 

void Atraf(float A[3][3], float (*At)[3][3]){
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            (*At)[i][j] = A[j][i];
} 

void ragnardJpassive(
    float params[4][8], double  theta[4], double passive[4][2],
    double dtheta[4], double detazeta[4][2], double sc[4][6], 
    double sct[8], double scez[4][4], double (*J)[3][3], 
    int leg_num, double (*dJ)[3][3])
{
    double sa = sc[leg_num][0];
    double ca = sc[leg_num][1];
    double sb = sc[leg_num][2];
    double cb = sc[leg_num][3];
    double b = double(params[leg_num][4]);
    double l = double(params[leg_num][5]);
    double st = sct[(2*leg_num)];
    double ct = sct[(2*leg_num)+1];
    double sz = scez[leg_num][0];
    double cz = scez[leg_num][1];
    double se = scez[leg_num][2];
    double ce = scez[leg_num][3];
    double dth = dtheta[leg_num];
    double dzeta = detazeta[leg_num][0];
    double deta = detazeta[leg_num][1];

    (*dJ)[0][0] = dzeta*(l*ce*cz*(sa*st - ca*cb*ct) + l*ce*sz*(sa*ct + ca*cb*st)) + deta*(l*se*cz*(sa*ct + ca*cb*st) - l*se*sz*(sa*st - ca*cb*ct)) + dth*((b + l*ce*cz)*(sa*st - ca*cb*ct) + l*ce*sz*(sa*ct + ca*cb*st));
    (*dJ)[0][1] = dth*l*ce*(cz*(sa*st - ca*cb*ct) + sz*(sa*ct + ca*cb*st)) + dzeta*l*ce*(cz*(sa*st - ca*cb*ct) + sz*(sa*ct + ca*cb*st)) + deta*l*se*(cz*(sa*ct + ca*cb*st) - sz*(sa*st - ca*cb*ct));
    (*dJ)[0][2] = dth*l*(se*cz*(sa*ct + ca*cb*st) - se*sz*(sa*st - ca*cb*ct)) + dzeta*l*(se*cz*(sa*ct + ca*cb*st) - se*sz*(sa*st - ca*cb*ct)) + deta*l*(ce*cz*(sa*st - ca*cb*ct) + ce*sz*(sa*ct + ca*cb*st) + ca*sb*se);
    (*dJ)[1][0] = - dzeta*(l*ce*cz*(ca*st + cb*sa*ct) + l*ce*sz*(ca*ct - cb*sa*st)) - deta*(l*se*cz*(ca*ct - cb*sa*st) - l*se*sz*(ca*st + cb*sa*ct)) - dth*((b + l*ce*cz)*(ca*st + cb*sa*ct) + l*ce*sz*(ca*ct - cb*sa*st));
    (*dJ)[1][1] = - dth*l*ce*(cz*(ca*st + cb*sa*ct) + sz*(ca*ct - cb*sa*st)) - dzeta*l*ce*(cz*(ca*st + cb*sa*ct) + sz*(ca*ct - cb*sa*st)) - deta*l*se*(cz*(ca*ct - cb*sa*st) - sz*(ca*st + cb*sa*ct));
    (*dJ)[1][2] = - deta*l*(ce*cz*(ca*st + cb*sa*ct) + ce*sz*(ca*ct - cb*sa*st) - sa*sb*se) - dth*l*(se*cz*(ca*ct - cb*sa*st) - se*sz*(ca*st + cb*sa*ct)) - dzeta*l*(se*cz*(ca*ct - cb*sa*st) - se*sz*(ca*st + cb*sa*ct));
    (*dJ)[2][0] = dth*sb*(b*ct + l*ce*ct*cz - l*ce*st*sz) + dzeta*l*(ct*cz-st*sz)*ce*sb - deta*l*(st*cz+ct*sz)*sb*se;
    (*dJ)[2][1] = l*sb*(dth*(ct*cz-st*sz)*ce + dzeta*(ct*cz-st*sz)*ce - deta*(st*cz+ct*sz)*se);
    (*dJ)[2][2] = deta*l*(cb*se + ce*sb*ct*cz - ce*sb*st*sz) - dth*l*(st*cz+ct*sz)*sb*se - dzeta*l*(st*cz+ct*sz)*sb*se;

}


void ragnardJpassivef(
    float params[4][8], float  theta[4], float passive[4][2],
    float dtheta[4], float detazeta[4][2], float sc[4][6], 
    float sct[8], float scez[4][4], int leg_num, float (*dJ)[3][3])
{
    float sa = sc[leg_num][0];
    float ca = sc[leg_num][1];
    float sb = sc[leg_num][2];
    float cb = sc[leg_num][3];
    float b = params[leg_num][4];
    float l = params[leg_num][5];
    float st = sct[(2*leg_num)];
    float ct = sct[(2*leg_num)+1];
    float sz = scez[leg_num][0];
    float cz = scez[leg_num][1];
    float se = scez[leg_num][2];
    float ce = scez[leg_num][3];
    float dth = dtheta[leg_num];
    float dzeta = detazeta[leg_num][0];
    float deta = detazeta[leg_num][1];
    /*
    Serial.println(sa,8);
    Serial.println(ca,8); 
    Serial.println(sb,8);
    Serial.println(cb,8);
    Serial.println(b,8);
    Serial.println(l,8);
    Serial.println(st,8);
    Serial.println(ct,8);
    Serial.println(sz,8);
    Serial.println(cz,8);
    Serial.println(se,8);
    Serial.println(ce,8);
    Serial.println(dth,8);
    Serial.println(dzeta,8);
    Serial.println(deta),8;
    */
    // dj(1,1) = dzeta*(l*ce*cz*(sa*st - ca*cb*ct) + l*ce*sz*(sa*ct + ca*cb*st)) + deta*(l*se*cz*(sa*ct + ca*cb*st) - l*se*sz*(sa*st - ca*cb*ct)) + dth*((b + l*ce*cz)*(sa*st - ca*cb*ct) + l*ce*sz*(sa*ct + ca*cb*st));
    // dj(1,2) = dth*l*ce*(cz*(sa*st - ca*cb*ct) + sz*(sa*ct + ca*cb*st)) + dzeta*l*ce*(cz*(sa*st - ca*cb*ct) + sz*(sa*ct + ca*cb*st)) + deta*l*se*(cz*(sa*ct + ca*cb*st) - sz*(sa*st - ca*cb*ct));
    // dj(1,3) = dth*l*(se*cz*(sa*ct + ca*cb*st) - se*sz*(sa*st - ca*cb*ct)) + dzeta*l*(se*cz*(sa*ct + ca*cb*st) - se*sz*(sa*st - ca*cb*ct)) + deta*l*(ce*cz*(sa*st - ca*cb*ct) + ce*sz*(sa*ct + ca*cb*st) + ca*sb*se);
    // dj(2,1) = - dzeta*(l*ce*cz*(ca*st + cb*sa*ct) + l*ce*sz*(ca*ct - cb*sa*st)) - deta*(l*se*cz*(ca*ct - cb*sa*st) - l*se*sz*(ca*st + cb*sa*ct)) - dth*((b + l*ce*cz)*(ca*st + cb*sa*ct) + l*ce*sz*(ca*ct - cb*sa*st));
    // dj(2,2) = - dth*l*ce*(cz*(ca*st + cb*sa*ct) + sz*(ca*ct - cb*sa*st)) - dzeta*l*ce*(cz*(ca*st + cb*sa*ct) + sz*(ca*ct - cb*sa*st)) - deta*l*se*(cz*(ca*ct - cb*sa*st) - sz*(ca*st + cb*sa*ct));
    // dj(2,3) = - deta*l*(ce*cz*(ca*st + cb*sa*ct) + ce*sz*(ca*ct - cb*sa*st) - sa*sb*se) - dth*l*(se*cz*(ca*ct - cb*sa*st) - se*sz*(ca*st + cb*sa*ct)) - dzeta*l*(se*cz*(ca*ct - cb*sa*st) - se*sz*(ca*st + cb*sa*ct));
    // dj(3,1) = dth*sb*(b*ct + l*ce*ct*cz - l*ce*st*sz) + dzeta*l*cos(theta + zeta)*ce*sb - deta*l*sin(theta + zeta)*sb*se;
    // dj(3,2) = l*sb*(dth*cos(theta + zeta)*ce + dzeta*cos(theta + zeta)*ce - deta*sin(theta + zeta)*se);
    // dj(3,3) = deta*l*(cb*se + ce*sb*ct*cz - ce*sb*st*sz) - dth*l*sin(theta + zeta)*sb*se - dzeta*l*sin(theta + zeta)*sb*se;
    (*dJ)[0][0] = dzeta*(l*ce*cz*(sa*st - ca*cb*ct) + l*ce*sz*(sa*ct + ca*cb*st)) + deta*(l*se*cz*(sa*ct + ca*cb*st) - l*se*sz*(sa*st - ca*cb*ct)) + dth*((b + l*ce*cz)*(sa*st - ca*cb*ct) + l*ce*sz*(sa*ct + ca*cb*st));
    (*dJ)[0][1] = dth*l*ce*(cz*(sa*st - ca*cb*ct) + sz*(sa*ct + ca*cb*st)) + dzeta*l*ce*(cz*(sa*st - ca*cb*ct) + sz*(sa*ct + ca*cb*st)) + deta*l*se*(cz*(sa*ct + ca*cb*st) - sz*(sa*st - ca*cb*ct));
    (*dJ)[0][2] = dth*l*(se*cz*(sa*ct + ca*cb*st) - se*sz*(sa*st - ca*cb*ct)) + dzeta*l*(se*cz*(sa*ct + ca*cb*st) - se*sz*(sa*st - ca*cb*ct)) + deta*l*(ce*cz*(sa*st - ca*cb*ct) + ce*sz*(sa*ct + ca*cb*st) + ca*sb*se);
    (*dJ)[1][0] = - dzeta*(l*ce*cz*(ca*st + cb*sa*ct) + l*ce*sz*(ca*ct - cb*sa*st)) - deta*(l*se*cz*(ca*ct - cb*sa*st) - l*se*sz*(ca*st + cb*sa*ct)) - dth*((b + l*ce*cz)*(ca*st + cb*sa*ct) + l*ce*sz*(ca*ct - cb*sa*st));
    (*dJ)[1][1] = - dth*l*ce*(cz*(ca*st + cb*sa*ct) + sz*(ca*ct - cb*sa*st)) - dzeta*l*ce*(cz*(ca*st + cb*sa*ct) + sz*(ca*ct - cb*sa*st)) - deta*l*se*(cz*(ca*ct - cb*sa*st) - sz*(ca*st + cb*sa*ct));
    (*dJ)[1][2] = - deta*l*(ce*cz*(ca*st + cb*sa*ct) + ce*sz*(ca*ct - cb*sa*st) - sa*sb*se) - dth*l*(se*cz*(ca*ct - cb*sa*st) - se*sz*(ca*st + cb*sa*ct)) - dzeta*l*(se*cz*(ca*ct - cb*sa*st) - se*sz*(ca*st + cb*sa*ct));
    (*dJ)[2][0] = dth*sb*(b*ct + l*ce*ct*cz - l*ce*st*sz) + dzeta*l*(ct*cz-st*sz)*ce*sb - deta*l*(st*cz+ct*sz)*sb*se;
    (*dJ)[2][1] = l*sb*(dth*(ct*cz-st*sz)*ce + dzeta*(ct*cz-st*sz)*ce - deta*(st*cz+ct*sz)*se);
    (*dJ)[2][2] = deta*l*(cb*se + ce*sb*ct*cz - ce*sb*st*sz) - dth*l*(st*cz+ct*sz)*sb*se - dzeta*l*(st*cz+ct*sz)*sb*se;

}

    

    