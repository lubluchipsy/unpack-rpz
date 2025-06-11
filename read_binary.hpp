#pragma once
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <map>
#include <array>

//READ_BINARY.CPP
#pragma pack(1)
struct Pack
{
    unsigned short  Npack; //1-2
    int             A[3];
    int             M[4];
    int             Tsi; //31-34
    unsigned int    Tsist; //35-38
    unsigned short  ski; //39-40
    unsigned short  mi; //41-42
};

class Data
{
public: 

    Data(): time{0}, P{0,0,0}, U{0,0,0}, I{0,0,0,0,0,0}, Ta{0,0,0}, Tmkd{0}, T_lgX{0,0}, T_lgY{0,0}, T_lgZ{0,0}, Fi{0,0,0}{};

    Data(Pack & pack);

    void count_L(Data & data_old);

    void count_dFi_r();

    void count_V(Data & data_old);

    void count_W();

    void process_mi();

    void process_mi(Data & data);

    void count_theta();

    unsigned short Npack; 
    std::array<double,3> A; // данные с акселерометров
    std::array<double,4> M; // кватернионы
    std::array<double,4> L; // кватернион поворота
    double         Tsi; 
    double         time;
    unsigned int   Tsist; 
    unsigned short ski; 
    unsigned short mi;
    std::array<double,3> dFi_r; // углы поворота за 5000 мкс (Tsi)
    std::array<double,3> V;
    std::array<double,3> W;
    std::array<double,3> Theta;

    std::array<unsigned short,3> P;
    std::array<unsigned short,3> U;
    std::array<unsigned short,6> I;

    std::array<double,3>  Ta;
    double  Tmkd;
    std::array<double,2>  T_lgX;
    std::array<double,2>  T_lgY;
    std::array<double,2>  T_lgZ;
    std::array<double,3>  Fi;  // накопленный угол, считается вместе с dFi_r
    
private:

    void copy_mi(Data & data);

    unsigned short sw1; 
    unsigned short sw2;
    unsigned short swd1;
    unsigned short swd2;
};
//READ_BINARY.CPP

//MODEL.CPP
struct Constants
{
    Constants():dG12{0,0,0,0}, dG13{0,0,0,0}, dG21{0,0,0,0}, dG23{0,0,0,0}, dG31{0,0,0,0}, dG32{0,0,0,0}, dKG1{0,0,0,0}, dKG2{0,0,0,0}, dKG3{0,0,0,0},
              dw1{0,0,0,0}, dw1dT(0), dwT1dT(0), dw2{0,0,0,0}, dw2dT(0), dwT2dT(0), dw3{0,0,0,0}, dw3dT(0), dwT3dT(0), dTnom(0), A1{1,0,0}, A2{0,1,0}, A3{0,0,1},
              dKA1{0,0,0,0}, dKA2{0,0,0,0}, dKA3{0,0,0,0}, dGA12{0,0,0,0}, dGA13{0,0,0,0}, dGA21{0,0,0,0}, dGA23{0,0,0,0}, dGA31{0,0,0,0}, dGA32{0,0,0,0},
              da1{0,0,0,0}, da1mkd(0), da2{0,0,0,0}, da2mkd(0), da3{0,0,0,0}, da3mkd(0){}

    Constants(std::array<double,4> idG12, std::array<double,4> idG13, std::array<double,4> idG21,
              std::array<double,4> idG23, std::array<double,4> idG31, std::array<double,4> idG32,
              std::array<double,4> idKG1, std::array<double,4> idKG2, std::array<double,4> idKG3,
              std::array<double,4> idw1, double idw1dT, double idwT1dT, std::array<double,4> idw2, double idw2dT,  double idwT2dT, std::array<double,4> idw3, double idw3dT, double idwT3dT,
              double idTnom, std::array<double,3> iA1, std::array<double,3> iA2, std::array<double,3> iA3,
              std::array<double,4> idKA1, std::array<double,4> idKA2, std::array<double,4> idKA3,
              std::array<double,4> idGA12, std::array<double,4> idGA13, std::array<double,4> idGA21,
              std::array<double,4> idGA23, std::array<double,4> idGA31, std::array<double,4> idGA32,
              std::array<double,4> ida1, double ida1mkd, std::array<double,4> ida2, double ida2mkd, std::array<double,4> ida3, double ida3mkd, bool immA):
              dG12(idG12), dG13(idG13), dG21(idG21), dG23(idG23), dG31(idG31), dG32(idG32), dKG1(idKG1), dKG2(idKG2), dKG3(idKG3),
              dw1(idw1), dw1dT(idw1dT), dwT1dT(idwT1dT), dw2(idw2), dw2dT(idw2dT), dwT2dT(idwT2dT), dw3(idw3), dw3dT(idw3dT), dwT3dT(idwT3dT), dTnom(idTnom), A1(iA1), A2(iA2), A3(iA3),
              dKA1(idKA1), dKA2(idKA2), dKA3(idKA3), dGA12(idGA12), dGA13(idGA13), dGA21(idGA21), dGA23(idGA23), dGA31(idGA31), dGA32(idGA32),
              da1(ida1), da1mkd(ida1mkd), da2(ida2), da2mkd(ida2mkd), da3(ida3), da3mkd(ida3mkd){}

    std::array<double,4> dG12; // [LG] dG
    std::array<double,4> dG13;
    std::array<double,4> dG21;
    std::array<double,4> dG23;
    std::array<double,4> dG31;
    std::array<double,4> dG32;

    std::array<double,4> dKG1;  //  [LG] dK
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

    std::array<double,4> dGA12; // [AC1] dG
    std::array<double,4> dGA13;     
    std::array<double,4> dGA21;
    std::array<double,4> dGA23;
    std::array<double,4> dGA31;
    std::array<double,4> dGA32;

    std::array<double,4> da1; double da1mkd;  // [AC1] da... da1mkd=ma
    std::array<double,4> da2; double da2mkd;
    std::array<double,4> da3; double da3mkd;
};
//MODEL.CPP

//READ_BINARY.CPP
class DataSum
{
public:
    DataSum(double taver): Npacks(0), T(0), Tsi(0), time(0), Taver(taver), A{0,0,0}, M{0,0,0,0}, L{0,0,0,0}, dFi_r{0,0,0},
                           V{0,0,0}, W{0,0,0}, Theta{0,0,0}, P{0,0,0}, U{0,0,0}, I{0,0,0,0,0,0}, Ta{0,0,0}, 
                           Tmkd{0}, T_lgX{0,0}, T_lgY{0,0}, T_lgZ{0,0}{};

    void add_data(const Data & data);

    unsigned int Npacks;
    double       T;
    double       Tsi;
    double       time;
    double       Taver;
    double       Tsist;
    unsigned int Npack;

    std::array<double,3> A; // данные с акселерометров
    std::array<double,4> M; // кватернионы
    std::array<double,4> L; // кватернион поворота
    std::array<double,3> dFi_r; // углы поворота за 5000 мкс (Tsi)
    std::array<double,3> V;
    std::array<double,3> W;
    std::array<double,3> Theta;

    std::array<unsigned short,3> P;
    std::array<unsigned short,3> U;
    std::array<unsigned short,6> I; 

    std::array<double,3> Ta; 
    double  Tmkd; 
    std::array<double,2>  T_lgX;
    std::array<double,2>  T_lgY;
    std::array<double,2>  T_lgZ;
};

class DataSum_m
{
public:
    DataSum_m(double taver): Npacks(0), T(0), Tsi(0), time(0), Taver(taver), A{0,0,0}, M{0,0,0,0}, L{0,0,0,0}, dFi_r{0,0,0},
                             V{0,0,0}, W{0,0,0}, Theta{0,0,0}, P{0,0,0}, U{0,0,0}, I{0,0,0,0,0,0}, Ta{0,0,0}, Tmkd{0}, 
                             T_lgX{0,0}, T_lgY{0,0}, T_lgZ{0,0}, dFi_r_corr{0,0,0}, Theta_corr{0,0,0}, V_corr{0,0,0}, W_corr{0,0,0}{};

    void add_data(const Data & data, const Constants & c);

    unsigned int Npacks;
    double       T;
    double       Tsi;
    double       time;
    double       Taver;
    double       Tsist;
    unsigned int Npack;

    std::array<double,3> A; // данные с акселерометров
    std::array<double,4> M; // кватернионы
    std::array<double,4> L; // кватернион поворота
    std::array<double,3> dFi_r; // углы поворота за 5000 мкс (Tsi)
    std::array<double,3> V;
    std::array<double,3> W;
    std::array<double,3> Theta;

    std::array<unsigned short,3> P;
    std::array<unsigned short,3> U;
    std::array<unsigned short,6> I; 

    std::array<double,3> Ta; 
    double  Tmkd; 
    std::array<double,2>  T_lgX;
    std::array<double,2>  T_lgY;
    std::array<double,2>  T_lgZ;

    std::array<double,3> dFi_r_corr;
    std::array<double,3> Theta_corr;
    std::array<double,3> V_corr;
    std::array<double,3> W_corr;
};

unsigned short Crc16(unsigned char * pcBlock, unsigned short len);

void find_header(std::fstream & fin);

bool is_good(std::fstream & fin);

void check_conditions(unsigned short c);

void print_header(std::fstream & fout, std::map<std::string, bool> & output_flags);

void output_corr_dFi(Data & data, std::fstream & fout, const Constants c);

void output_corr_Theta(Data & data, std::fstream & fout, const Constants c);

void output_corr_Omega(Data & data, std::fstream & fout, const Constants c);

void output_corr_V(Data & data, std::fstream & fout, const Constants c);

void output_corr_W(Data & data, std::fstream & fout, const Constants c);

void output_data(Data & data, std::fstream & fout, std::map<std::string, bool> & output_flags);

void count_data(Data & data, Data & data_old);

void read_data(std::fstream & fin, std::fstream & fout, std::pair<std::map<std::string, bool>, std::map<std::string, double>> & output_params, const Constants & c);

void output_average_data(DataSum & datasum, std::fstream & fout, std::map<std::string, bool> & output_flags);

void read_average_data(std::fstream & fin, std::fstream & fout, std::pair<std::map<std::string, bool>, std::map<std::string, double>> & output_params);

void output_average_data(DataSum_m & datasum, std::fstream & fout, std::map<std::string, bool> & output_flags);

void read_average_data(std::fstream & fin, std::fstream & fout, std::pair<std::map<std::string, bool>, std::map<std::string, double>> & output_params, const Constants & c);

void get_flag(std::fstream & config, std::map<std::string, bool> & flags);

std::pair<std::map<std::string, bool>, std::map<std::string, double>> read_config(std::fstream & config);

std::string str_from_config(std::fstream & config);
//READ_BINARY.CPP

//MODEL.CPP
double poly3(const double a0, const double a1, const double a2, const double a3, const double T);

std::array<double,3> zero_offset(const Constants & c, const Data & d);

std::array<std::array<double,2>, 3> ak_misalignment_params(const Constants & c, const Data & d);

std::array<double, 3> scale_amendments(const Constants & c, const Data & d);

std::array<double, 3> correct_V(const Constants & c, const Data & d);

std::array<double, 3> scale_coef_amendments(const Constants & c, const Data & d);

std::array<std::array<double,2>, 3> gyro_misalignment_params(const Constants & c, const Data & d);

std::array<double, 3> gyro_drift(const Constants & c, const Data & d);

std::array<double,3> correct_angles(const Constants & c, const Data & d);

double get_value(std::fstream & fin);

void read_pcfd(std::fstream & pcfd, Constants & constants);

//MODEL.CPP

