#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include "read_binary.hpp"

double poly3(const double a0, const double a1, const double a2, const double a3, const double T)
{
    return (a0 + a1*T + a2*T*T + a3*T*T*T);
}

std::array<double,3> zero_offset(const Constants & c, const Data & d)  // 6.Алгоритм «Вычисление смещение нуля АК»
{
    std::array<double,3> da;
    da[0] = poly3(c.da1[0], c.da1[1], c.da1[2], c.da1[3], d.Ta[0]) + c.da1mkd*d.Tmkd;
    da[1] = poly3(c.da2[0], c.da2[1], c.da2[2], c.da2[3], d.Ta[1]) + c.da2mkd*d.Tmkd;
    da[2] = poly3(c.da3[0], c.da3[1], c.da3[2], c.da3[3], d.Ta[2]) + c.da3mkd*d.Tmkd;
    return da;
}

std::array<std::array<double,2>, 3> ak_misalignment_params(const Constants & c, const Data & d) // 5.Алгоритм «Вычисление параметров несоосности АК»
{
    std::array<std::array<double,2>, 3> dGA;                                                    //dGA = dGA12 dGA13
    dGA[0][0] = poly3(c.dGA12[0], c.dGA12[1], c.dGA12[2], c.dGA12[3], d.Ta[0]);                //       dGA21 dGA23
    dGA[0][1] = poly3(c.dGA13[0], c.dGA13[1], c.dGA13[2], c.dGA13[3], d.Ta[0]);               //        dGA31 dGA32
    dGA[1][0] = poly3(c.dGA21[0], c.dGA21[1], c.dGA21[2], c.dGA21[3], d.Ta[1]);
    dGA[1][1] = poly3(c.dGA23[0], c.dGA23[1], c.dGA23[2], c.dGA23[3], d.Ta[1]);
    dGA[2][0] = poly3(c.dGA31[0], c.dGA31[1], c.dGA31[2], c.dGA31[3], d.Ta[2]);
    dGA[2][1] = poly3(c.dGA32[0], c.dGA32[1], c.dGA32[2], c.dGA32[3], d.Ta[2]);
    return dGA;
}

std::array<double, 3> scale_amendments(const Constants & c, const Data & d)  // 4.Алгоритм «Вычисление поправок масштабного коэффициента АК»
{
    std::array<double,3> dKA;
    dKA[0] = poly3(c.dKA1[0], c.dKA1[1], c.dKA1[2], c.dKA1[3], d.Ta[0]);
    dKA[1] = poly3(c.dKA2[0], c.dKA2[1], c.dKA2[2], c.dKA2[3], d.Ta[1]);
    dKA[2] = poly3(c.dKA3[0], c.dKA3[1], c.dKA3[2], c.dKA3[3], d.Ta[2]);
    return dKA;
}

std::array<double, 3> correct_V(const Constants & c, const Data & d)  // 3. Алгоритм коррекции измерений акселерометров
{
    std::array<double,3> V;
    std::array<std::array<double,2>, 3> dGA = ak_misalignment_params(c, d);
    std::array<double,3> da = zero_offset(c, d);
    std::array<double,3> dKA = scale_amendments(c, d);

    V[0] = (1+dKA[0])*d.V[0] + dGA[0][0]*d.V[1] + dGA[0][1]*d.V[2] - da[0]*d.Tsi/1000000;
    V[1] = dGA[1][0]*d.V[0] + (1+dKA[1])*d.V[1] + dGA[1][1]*d.V[2] - da[1]*d.Tsi/1000000;
    V[2] = dGA[2][0]*d.V[0] + dGA[2][1]*d.V[1] + (1+dKA[2])*d.V[2] - da[2]*d.Tsi/1000000;
    return V;
}

std::array<double, 3> scale_coef_amendments(const Constants & c, const Data & d)   // 1.2. Алгоритм «Вычисление поправок масштабного коэффициента ЛГ»
{
    std::array<double,3> dKG;
    dKG[0] = poly3(c.dKG1[0], c.dKG1[1], c.dKG1[2], c.dKG1[3], d.T_lgX[0]);
    dKG[1] = poly3(c.dKG2[0], c.dKG2[1], c.dKG2[2], c.dKG2[3], d.T_lgY[0]);
    dKG[2] = poly3(c.dKG3[0], c.dKG3[1], c.dKG3[2], c.dKG3[3], d.T_lgZ[0]);
    return dKG;
}

std::array<std::array<double,2>, 3> gyro_misalignment_params(const Constants & c, const Data & d) // 1.1.  Алгоритм «Вычисление параметров несоосности строительных осей гироскопов»
{
    std::array<std::array<double,2>, 3> dG;                                                  // dG = dG12 dG13
    dG[0][0] = poly3(c.dG12[0], c.dG12[1], c.dG12[2], c.dG12[3], d.T_lgX[0]);                //      dG21 dG23
    dG[0][1] = poly3(c.dG13[0], c.dG13[1], c.dG13[2], c.dG13[3], d.T_lgX[0]);                //      dG31 dG32
    dG[1][0] = poly3(c.dG21[0], c.dG21[1], c.dG21[2], c.dG21[3], d.T_lgY[0]);
    dG[1][1] = poly3(c.dG23[0], c.dG23[1], c.dG23[2], c.dG23[3], d.T_lgY[0]);
    dG[2][0] = poly3(c.dG31[0], c.dG31[1], c.dG31[2], c.dG31[3], d.T_lgZ[0]);
    dG[2][1] = poly3(c.dG32[0], c.dG32[1], c.dG32[2], c.dG32[3], d.T_lgZ[0]);
    return dG;
}

std::array<double, 3> gyro_drift(const Constants & c, const Data & d) //в строительных осях   2. Алгоритм «Вычисление модельных значений дрейфа».
{
    std::array<double, 3> dw_m; //в измерительных осях
    double Text{(d.T_lgX[0]+d.T_lgY[0]+d.T_lgZ[0])/3};
    dw_m[0] = poly3(c.dw1[0], c.dw1[1], c.dw1[2], c.dw1[3], d.T_lgX[0]) + c.dw1dT*(d.T_lgX[0] - Text - c.dTnom) + c.dwT1dT*d.T_lgX[0]*(d.T_lgX[0] - Text - c.dTnom);  // (d.T_lgX+d.T_lgY+d.T_lgZ)/3=Text
    dw_m[1] = poly3(c.dw2[0], c.dw2[1], c.dw2[2], c.dw2[3], d.T_lgY[0]) + c.dw2dT*(d.T_lgY[0] - Text - c.dTnom) + c.dwT2dT*d.T_lgY[0]*(d.T_lgY[0] - Text - c.dTnom);
    dw_m[2] = poly3(c.dw3[0], c.dw3[1], c.dw3[2], c.dw3[3], d.T_lgZ[0]) + c.dw3dT*(d.T_lgZ[0] - Text - c.dTnom) + c.dwT3dT*d.T_lgZ[0]*(d.T_lgZ[0] - Text - c.dTnom);

    std::array<double, 3> dw;
    dw[0] = c.A1[0]*dw_m[0] + c.A1[1]*dw_m[1] + c.A1[2]*dw_m[2];
    dw[1] = c.A2[0]*dw_m[0] + c.A2[1]*dw_m[1] + c.A2[2]*dw_m[2];
    dw[2] = c.A3[0]*dw_m[0] + c.A3[1]*dw_m[1] + c.A3[2]*dw_m[2];
    return dw;
}

std::array<double,3> correct_angles(const Constants & c, const Data & d)  // 1. Алгоритм коррекции измерений гироскопов 
{
    std::array<double,3> Theta;
    std::array<std::array<double,2>, 3> dG = gyro_misalignment_params(c, d);
    std::array<double,3> dw = gyro_drift(c, d);
    std::array<double,3> dKG = scale_coef_amendments(c, d);

    Theta[0] = (1 + dKG[0])*d.Theta[0] + dG[0][0]*d.Theta[1] + dG[0][1]*d.Theta[2] - dw[0];
    Theta[1] = dG[1][0]*d.Theta[0] + (1 + dKG[1])*d.Theta[1] + dG[1][1]*d.Theta[2] - dw[1];
    Theta[2] = dG[2][0]*d.Theta[0] + dG[2][1]*d.Theta[1] + (1 + dKG[2])*d.Theta[2] - dw[2];
    return Theta;
}

double get_value(std::fstream & fin)
{
    std::string line;
    std::getline(fin, line);
    std::istringstream iss(line);
    std::string str;
    std::getline(iss, str, '=');
    double value;
    iss >> value;
    return value;
}

void read_pcfd(std::fstream & pcfd, Constants & constants)
{
    std::array<double,4> dG12; // [LG] dG
    std::array<double,4> dG13;
    std::array<double,4> dG21;
    std::array<double,4> dG23;
    std::array<double,4> dG31;
    std::array<double,4> dG32;

    std::array<double,4> dKG1;  //  [LG] dG
    std::array<double,4> dKG2;
    std::array<double,4> dKG3;

    std::array<double,4> dw1; double dw1dT; double dwT1dT; // LG dw12, dw1dT
    std::array<double,4> dw2; double dw2dT; double dwT2dT;
    std::array<double,4> dw3; double dw3dT; double dwT3dT;

    double dTnom;

    std::array<double,3> A1;  // [ortho] M
    std::array<double,3> A2; 
    std::array<double,3> A3;

    std::array<double,4> dKA1; // [AC1] dK
    std::array<double,4> dKA2;
    std::array<double,4> dKA3;

    std::array<double,4> dGA12; //[AC1] dG
    std::array<double,4> dGA13;     
    std::array<double,4> dGA21;
    std::array<double,4> dGA23;
    std::array<double,4> dGA31;
    std::array<double,4> dGA32;

    std::array<double,4> da1; double da1mkd;  //[AC1] da... da1mkd=ma
    std::array<double,4> da2; double da2mkd;
    std::array<double,4> da3; double da3mkd;

    std::string line;
    for (auto i = 0; i < 27; i++)
    {
        std::getline(pcfd, line);
    }
    
    dTnom = get_value(pcfd);
    A1[0] = get_value(pcfd); 
    A1[1] = get_value(pcfd);
    A1[2] = get_value(pcfd); 
    A2[0] = get_value(pcfd);
    A2[1] = get_value(pcfd); 
    A2[2] = get_value(pcfd);
    A3[0] = get_value(pcfd); 
    A3[1] = get_value(pcfd);
    A3[2] = get_value(pcfd);

    std::getline(pcfd, line);  //[AK1] 
    dKA1[0] = get_value(pcfd);
    dKA1[1] = get_value(pcfd); 
    dKA1[2] = get_value(pcfd); 
    dKA1[3] = get_value(pcfd); 
    da1[0] = get_value(pcfd); 
    da1[1] = get_value(pcfd);
    da1[2] = get_value(pcfd); 
    da1[3] = get_value(pcfd); 
    da1mkd = get_value(pcfd); 
    dGA12[0] = get_value(pcfd); 
    dGA12[1] = get_value(pcfd); 
    dGA12[2] = get_value(pcfd); 
    dGA12[3] = get_value(pcfd); 
    dGA13[0] = get_value(pcfd); 
    dGA13[1] = get_value(pcfd);
    dGA13[2] = get_value(pcfd); 
    dGA13[3] = get_value(pcfd); 

    std::getline(pcfd, line); //[AK2]
    dKA2[0] = get_value(pcfd); 
    dKA2[1] = get_value(pcfd);
    dKA2[2] = get_value(pcfd);
    dKA2[3] = get_value(pcfd); 
    da2[0] = get_value(pcfd); 
    da2[1] = get_value(pcfd); 
    da2[2] = get_value(pcfd);
    da2[3] = get_value(pcfd); 
    da2mkd = get_value(pcfd); 
    dGA21[0] = get_value(pcfd);
    dGA21[1] = get_value(pcfd);
    dGA21[2] = get_value(pcfd);
    dGA21[3] = get_value(pcfd);
    dGA23[0] = get_value(pcfd);
    dGA23[1] = get_value(pcfd);
    dGA23[2] = get_value(pcfd);
    dGA23[3] = get_value(pcfd);

    std::getline(pcfd, line);  //[AK3]
    dKA3[0] = get_value(pcfd); 
    dKA3[1] = get_value(pcfd); 
    dKA3[2] = get_value(pcfd); 
    dKA3[3] = get_value(pcfd); 
    da3[0] = get_value(pcfd);
    da3[1] = get_value(pcfd);
    da3[2] = get_value(pcfd);
    da3[3] = get_value(pcfd);
    da3mkd = get_value(pcfd);
    dGA31[0] = get_value(pcfd); 
    dGA31[1] = get_value(pcfd);
    dGA31[2] = get_value(pcfd); 
    dGA31[3] = get_value(pcfd);
    dGA32[0] = get_value(pcfd); 
    dGA32[1] = get_value(pcfd); 
    dGA32[2] = get_value(pcfd);
    dGA32[3] = get_value(pcfd);

    std::getline(pcfd, line); //[T_lgX]
    dKG1[0] = get_value(pcfd); 
    dKG1[1] = get_value(pcfd); 
    dKG1[2] = get_value(pcfd); 
    dKG1[3] = get_value(pcfd); 
    dw1[0] = get_value(pcfd); 
    dw1[1] = get_value(pcfd);  
    dw1[2] = get_value(pcfd);  
    dw1[3] = get_value(pcfd);  
    dw1dT = get_value(pcfd); 
    dwT1dT = get_value(pcfd); 
    dG12[0] = get_value(pcfd); 
    dG12[1] = get_value(pcfd); 
    dG12[2] = get_value(pcfd); 
    dG12[3] = get_value(pcfd); 
    dG13[0] = get_value(pcfd); 
    dG13[1] = get_value(pcfd); 
    dG13[2] = get_value(pcfd); 
    dG13[3] = get_value(pcfd); 

    std::getline(pcfd, line);  //[T_lgY]
    dKG2[0] = get_value(pcfd); 
    dKG2[1] = get_value(pcfd); 
    dKG2[2] = get_value(pcfd);
    dKG2[3] = get_value(pcfd); 
    dw2[0] = get_value(pcfd); 
    dw2[1] = get_value(pcfd);
    dw2[2] = get_value(pcfd); 
    dw2[3] = get_value(pcfd); 
    dw2dT = get_value(pcfd); 
    dwT2dT = get_value(pcfd);
    dG21[0] = get_value(pcfd); 
    dG21[1] = get_value(pcfd); 
    dG21[2] = get_value(pcfd); 
    dG21[3] = get_value(pcfd); 
    dG23[0] = get_value(pcfd); 
    dG23[1] = get_value(pcfd); 
    dG23[2] = get_value(pcfd); 
    dG23[3] = get_value(pcfd); 

    std::getline(pcfd, line);  //[T_lgZ]
    dKG3[0] = get_value(pcfd); 
    dKG3[1] = get_value(pcfd); 
    dKG3[2] = get_value(pcfd); 
    dKG3[3] = get_value(pcfd); 
    dw3[0] = get_value(pcfd); 
    dw3[1] = get_value(pcfd);
    dw3[2] = get_value(pcfd); 
    dw3[3] = get_value(pcfd); 
    dw3dT = get_value(pcfd);
    dwT3dT = get_value(pcfd);
    dG31[0] = get_value(pcfd); 
    dG31[1] = get_value(pcfd); 
    dG31[2] = get_value(pcfd); 
    dG31[3] = get_value(pcfd); 
    dG32[0] = get_value(pcfd); 
    dG32[1] = get_value(pcfd); 
    dG32[2] = get_value(pcfd); 
    dG32[3] = get_value(pcfd); 

    constants = Constants(dG12, dG13, dG21, dG23, dG31, dG32, dKG1, dKG2, dKG3, dw1, dw1dT, dwT1dT, dw2, dw2dT, dwT2dT, dw3, dw3dT, dwT3dT, dTnom, A1, A2, A3,
                        dKA1, dKA2, dKA3, dGA12, dGA13, dGA21, dGA23, dGA31, dGA32,da1, da1mkd, da2, da2mkd, da3, da3mkd, 0);
}



