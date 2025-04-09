#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <filesystem>
#include <string>
#include <map>
#include <limits>
#include <numbers>
#include <array>
#include <locale>
#include "read_binary.hpp"


namespace LSB 
{
    const double M   = 1 / pow(2,31);
    const double A   = 0.2 / pow(2,19);
    const double Tsi = 1 / 60.0;
    const double t   = 0.05;
}

unsigned int counter = 0;

#pragma pack(1)


Data::Data(Pack & pack): Npack(pack.Npack), Tsi(pack.Tsi * LSB::Tsi), time{0}, Tsist(pack.Tsist), ski(pack.ski), mi(pack.mi),
                         P{0,0,0}, U{0,0,0}, I{0,0,0,0,0,0}, Ta{0,0,0}, Tmkd{0}, T_lgX{0,0}, T_lgY{0,0}, T_lgZ{0,0}
{
    for (auto i = 0; i < 3; i++)
    {
        A[i] = pack.A[i]/1.0;
    }

    for (auto i = 0; i < 4; i++)
    {
        M[i] = pack.M[i] * LSB::M;
    }
}


void Data::count_L(Data & data_old)
{
    std::array<double,4> M_old = data_old.M;

    L[0]=M_old[0]*M[0]+M_old[1]*M[1]+M_old[2]*M[2]+M_old[3]*M[3];
    L[1]=M_old[0]*M[1]-M_old[1]*M[0]-M_old[2]*M[3]+M_old[3]*M[2];
    L[2]=M_old[0]*M[2]-M_old[2]*M[0]-M_old[3]*M[1]+M_old[1]*M[3];
    L[3]=M_old[0]*M[3]-M_old[3]*M[0]-M_old[1]*M[2]+M_old[2]*M[1];

    float Lnorm = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2] + L[3]*L[3]);

    if (Lnorm == 0)
    {
        Lnorm = 1;
    }
    for (auto l : L)
    {
        l /= Lnorm;
    }
}


void Data::count_dFi_r()
{
    for (auto i = 0; i < 3; i++)
    {
        dFi_r[i] = L[i+1]*2;
        Fi[i] += dFi_r[i];
    }
}

void Data::count_V(Data & data_old)
{
    for (auto i = 0; i < 3; i++)
    {
        if ((A[i]-data_old.A[i]) > std::numeric_limits<int>::max())
            V.at(i) = ((A[i] - data_old.A[i]) - pow(2,32)) * LSB::A;
        else if ((A[i]-data_old.A[i]) < std::numeric_limits<int>::min())
            V[i] = ((A[i] - data_old.A[i]) + pow(2,32)) * LSB::A;
        else
            V[i] = (A[i] - data_old.A[i]) * LSB::A;
    }
}

void Data::count_W()
{
    for (auto i = 0; i < 3; i++)
    {
        W[i] = V[i] * 1000000 / Tsi;
    }
}

void Data::count_theta()
{
    for (auto i = 0; i < 3; i++)
    {
        Theta[i] = dFi_r[i] * 1000000 / Tsi;
    }
}


void Data::process_mi()
{
    switch(Npack % 32)
    {
        case 0:
            Ta[0] = static_cast<short>(mi) * LSB::t;
            break;
        case 1:
            Ta[1] = static_cast<short>(mi) * LSB::t;
            break;
        case 2:
            Ta[2] = static_cast<short>(mi) * LSB::t;
            break;
        case 3:
            Tmkd  = static_cast<short>(mi) * LSB::t;
            break;
        case 4:
            T_lgX[0]= static_cast<short>(mi) * LSB::t;
            break;
        case 5:  
            T_lgX[1]= static_cast<short>(mi) * LSB::t;
            break;
        case 6:  
            T_lgY[0]= static_cast<short>(mi) * LSB::t;
            break;
        case 7: 
            T_lgY[1]= static_cast<short>(mi) * LSB::t;
            break;
        case 8: 
            T_lgZ[0]= static_cast<short>(mi) * LSB::t;
            break;
        case 9:  
            T_lgZ[1]= static_cast<short>(mi) * LSB::t;
            break;
        case 13: 
            sw1   = mi;
            break;
        case 14: 
            sw2   = mi;
            break;
        case 15: 
            swd1  = mi;
            break;
        case 16: 
            swd2  = mi;
            break;
        case 17: 
            P[0]  = mi;
            break;
        case 18: 
            P[1]  = mi;
            break;
        case 19: 
            P[2]  = mi;
            break;
        case 20: 
            U[0]  = mi;
            break;
        case 21: 
            U[1]  = mi;
            break;
        case 22:
            U[2]  = mi;
            break;
        case 23: 
            I[0]  = mi;
            break;
        case 24: 
            I[1]  = mi;
            break;
        case 25: 
            I[2]  = mi;
            break;
        case 26: 
            I[3]  = mi;
            break;
        case 27: 
            I[4]  = mi;
            break;
        case 28: 
            I[5]  = mi;
            break;
    } 
}


void Data::process_mi(Data & data)
{
    switch(Npack % 32)
    {
        case 0:
            copy_mi(data);
            Ta[0] = static_cast<short>(mi) * LSB::t;
            break;
        case 1: 
            copy_mi(data);
            Ta[1] = static_cast<short>(mi) * LSB::t;
            break;
        case 2:  
            copy_mi(data);
            Ta[2] = static_cast<short>(mi) * LSB::t;
            break;
        case 3:  
            copy_mi(data);
            Tmkd  = static_cast<short>(mi) * LSB::t;
            break;
        case 4:  
            copy_mi(data);
            T_lgX[0]= static_cast<short>(mi) * LSB::t;
            break;
        case 5:  
            copy_mi(data);
            T_lgX[1]= static_cast<short>(mi) * LSB::t;
            break;
        case 6:  
            copy_mi(data);
            T_lgY[0]= static_cast<short>(mi) * LSB::t;
            break;
        case 7: 
            copy_mi(data); 
            T_lgY[1] = static_cast<short>(mi) * LSB::t;
            break;
        case 8:  
            copy_mi(data);
            T_lgZ[0]= static_cast<short>(mi) * LSB::t;
            break;
        case 9:  
            copy_mi(data);
            T_lgZ[1]= static_cast<short>(mi) * LSB::t;
            break;
        case 13: 
            copy_mi(data);
            sw1   = mi;
            break;
        case 14:
            copy_mi(data); 
            sw2   = mi;
            break;
        case 15: 
            copy_mi(data);
            swd1  = mi;
            break;
        case 16: 
            copy_mi(data);
            swd2  = mi;
            break;
        case 17: 
            copy_mi(data);
            P[0]  = mi;
            break;
        case 18: 
            copy_mi(data);
            P[1]  = mi;
            break;
        case 19: 
            copy_mi(data);
            P[2]  = mi;
            break;
        case 20: 
            copy_mi(data);
            U[0]  = mi;
            break;
        case 21: 
            copy_mi(data);
            U[1]  = mi;
            break;
        case 22:
            copy_mi(data); 
            U[2]  = mi;
            break;
        case 23: 
            copy_mi(data);
            I[0]  = mi;
            break;
        case 24: 
            copy_mi(data);
            I[1]  = mi;
            break;
        case 25: 
            copy_mi(data);
            I[2]  = mi;
            break;
        case 26: 
            copy_mi(data);
            I[3]  = mi;
            break;
        case 27: 
            copy_mi(data);
            I[4]  = mi;
            break;
        case 28: 
            copy_mi(data);
            I[5]  = mi;
            break;
        default:
            copy_mi(data);
    } 
}


void Data::copy_mi(Data & data)
{
    for (auto i = 0; i < 3; i++)
    {
        P[i]  = data.P[i];
        U[i]  = data.U[i];
        Ta[i] = data.Ta[i];
    }
    for (auto i = 0; i < 6; i++)
    {
        I[i]  = data.I[i];
    }
    Tmkd = data.Tmkd;
    T_lgX[0]  = data.T_lgX[0];
    T_lgX[1]  = data.T_lgX[1];
    T_lgY[0]  = data.T_lgY[0];
    T_lgY[1]  = data.T_lgY[1];
    T_lgZ[0]  = data.T_lgZ[0];
    T_lgZ[1]  = data.T_lgZ[1];
}


void DataSum::add_data(const Data & data)
{
    T      += data.Tsi;
    Npacks += 1;
    Tsist   = data.Tsist;
    Tsi    += data.Tsi;
    time    = data.time;
    Npack   = data.Npack; 
    
    for (auto i = 0; i < 3; i++)
    {
        A[i]     += data.A[i];
        V[i]     += data.V[i];
        W[i]     += data.W[i];
        P[i]     += data.P[i];
        dFi_r[i] += data.dFi_r[i];
        Ta[i]    += data.Ta[i];
        Theta[i] += data.Theta[i];
    }

    for (auto i = 0; i < 4; i++)
    {
        M[i] += data.M[i];
        L[i] += data.L[i];
    }

    for (auto i = 0; i < 6; i++)
    {
        I[i] = I[i] + data.I[i];
    } 

    T_lgX[0]+= data.T_lgX[0];
    T_lgX[1]+= data.T_lgX[1];
    T_lgY[0]+= data.T_lgY[0];
    T_lgY[1]+= data.T_lgY[1];
    T_lgZ[0]+= data.T_lgZ[0];
    T_lgZ[1]+= data.T_lgZ[1];
    Tmkd  += data.Tmkd;
}


void DataSum_m::add_data(const Data & data, const Constants & c)
{
    T      += data.Tsi;
    Npacks += 1;
    Tsist   = data.Tsist;
    Tsi    += data.Tsi;
    time    = data.time;
    Npack   = data.Npack; 
    
    for (auto i = 0; i < 3; i++)
    {
        A[i]     += data.A[i];
        V[i]     += data.V[i];
        W[i]     += data.W[i];
        P[i]     += data.P[i];
        dFi_r[i] += data.dFi_r[i];
        Ta[i]    += data.Ta[i];
        Theta[i] += data.Theta[i];
    }

    for (auto i = 0; i < 4; i++)
    {
        M[i] += data.M[i];
        L[i] += data.L[i];
    }

    for (auto i = 0; i < 6; i++)
    {
        I[i] = I[i] + data.I[i];
    } 

    T_lgX[0]+= data.T_lgX[0];
    T_lgX[1]+= data.T_lgX[1];
    T_lgY[0]+= data.T_lgY[0];
    T_lgY[1]+= data.T_lgY[1];
    T_lgZ[0]+= data.T_lgZ[0];
    T_lgZ[1]+= data.T_lgZ[1];
    Tmkd  += data.Tmkd;

    std::array<double,3> data_Theta_corr = correct_angles(c, data);
    std::array<double,3> data_dFi_r_corr;
    std::array<double,3> data_V_corr = correct_V(c, data);
    std::array<double,3> data_W_corr;

    for (auto i = 0; i < 3; i++)
    {
        data_dFi_r_corr[i] = (data.Tsi/1000000) * data_Theta_corr[i];
        data_W_corr[i] = data_V_corr[i] * 1000000 / data.Tsi;
    }

    for (auto i = 0; i < 3; i++)
    {
        dFi_r_corr[i] += data_dFi_r_corr[i];
        Theta_corr[i] += data_Theta_corr[i];
        V_corr[i] += data_V_corr[i];
        W_corr[i] += data_W_corr[i];
    }
}


const unsigned short Crc16Table[256] = {
    0x0000, 0x1021, 0x2042, 0x3063, 0x4084, 0x50A5, 0x60C6, 0x70E7,
    0x8108, 0x9129, 0xA14A, 0xB16B, 0xC18C, 0xD1AD, 0xE1CE, 0xF1EF,
    0x1231, 0x0210, 0x3273, 0x2252, 0x52B5, 0x4294, 0x72F7, 0x62D6,
    0x9339, 0x8318, 0xB37B, 0xA35A, 0xD3BD, 0xC39C, 0xF3FF, 0xE3DE,
    0x2462, 0x3443, 0x0420, 0x1401, 0x64E6, 0x74C7, 0x44A4, 0x5485,
    0xA56A, 0xB54B, 0x8528, 0x9509, 0xE5EE, 0xF5CF, 0xC5AC, 0xD58D,
    0x3653, 0x2672, 0x1611, 0x0630, 0x76D7, 0x66F6, 0x5695, 0x46B4,
    0xB75B, 0xA77A, 0x9719, 0x8738, 0xF7DF, 0xE7FE, 0xD79D, 0xC7BC,
    0x48C4, 0x58E5, 0x6886, 0x78A7, 0x0840, 0x1861, 0x2802, 0x3823,
    0xC9CC, 0xD9ED, 0xE98E, 0xF9AF, 0x8948, 0x9969, 0xA90A, 0xB92B,
    0x5AF5, 0x4AD4, 0x7AB7, 0x6A96, 0x1A71, 0x0A50, 0x3A33, 0x2A12,
    0xDBFD, 0xCBDC, 0xFBBF, 0xEB9E, 0x9B79, 0x8B58, 0xBB3B, 0xAB1A,
    0x6CA6, 0x7C87, 0x4CE4, 0x5CC5, 0x2C22, 0x3C03, 0x0C60, 0x1C41,
    0xEDAE, 0xFD8F, 0xCDEC, 0xDDCD, 0xAD2A, 0xBD0B, 0x8D68, 0x9D49,
    0x7E97, 0x6EB6, 0x5ED5, 0x4EF4, 0x3E13, 0x2E32, 0x1E51, 0x0E70,
    0xFF9F, 0xEFBE, 0xDFDD, 0xCFFC, 0xBF1B, 0xAF3A, 0x9F59, 0x8F78,
    0x9188, 0x81A9, 0xB1CA, 0xA1EB, 0xD10C, 0xC12D, 0xF14E, 0xE16F,
    0x1080, 0x00A1, 0x30C2, 0x20E3, 0x5004, 0x4025, 0x7046, 0x6067,
    0x83B9, 0x9398, 0xA3FB, 0xB3DA, 0xC33D, 0xD31C, 0xE37F, 0xF35E,
    0x02B1, 0x1290, 0x22F3, 0x32D2, 0x4235, 0x5214, 0x6277, 0x7256,
    0xB5EA, 0xA5CB, 0x95A8, 0x8589, 0xF56E, 0xE54F, 0xD52C, 0xC50D,
    0x34E2, 0x24C3, 0x14A0, 0x0481, 0x7466, 0x6447, 0x5424, 0x4405,
    0xA7DB, 0xB7FA, 0x8799, 0x97B8, 0xE75F, 0xF77E, 0xC71D, 0xD73C,
    0x26D3, 0x36F2, 0x0691, 0x16B0, 0x6657, 0x7676, 0x4615, 0x5634,
    0xD94C, 0xC96D, 0xF90E, 0xE92F, 0x99C8, 0x89E9, 0xB98A, 0xA9AB,
    0x5844, 0x4865, 0x7806, 0x6827, 0x18C0, 0x08E1, 0x3882, 0x28A3,
    0xCB7D, 0xDB5C, 0xEB3F, 0xFB1E, 0x8BF9, 0x9BD8, 0xABBB, 0xBB9A,
    0x4A75, 0x5A54, 0x6A37, 0x7A16, 0x0AF1, 0x1AD0, 0x2AB3, 0x3A92,
    0xFD2E, 0xED0F, 0xDD6C, 0xCD4D, 0xBDAA, 0xAD8B, 0x9DE8, 0x8DC9,
    0x7C26, 0x6C07, 0x5C64, 0x4C45, 0x3CA2, 0x2C83, 0x1CE0, 0x0CC1,
    0xEF1F, 0xFF3E, 0xCF5D, 0xDF7C, 0xAF9B, 0xBFBA, 0x8FD9, 0x9FF8,
    0x6E17, 0x7E36, 0x4E55, 0x5E74, 0x2E93, 0x3EB2, 0x0ED1, 0x1EF0
};


unsigned short Crc16(unsigned char * pcBlock, unsigned short len)
{
    unsigned short crc = 0xFFFF;
    while (len--)
        crc = (crc << 8) ^ Crc16Table[(crc >> 8) ^ *pcBlock++];
    return crc;
}


void find_header(std::fstream & fin)
{
        //проверяем, неполный ли пакет
    unsigned short x;

    for (auto pos = 45; pos != 0; pos--) //позиция отн начала плохого пакета 
    {
        fin.seekp(-1, std::ios::cur);
        fin.read(reinterpret_cast<char *>(&x), 2);
        fin.seekp(-2, std::ios::cur);
        if (x == 50624)
        {
            //std::cout << "After pack " << counter << " incomplete pack" << std::endl;
            return;
        }
    }
        //ищем заголовок после плосле плохого пакета
    fin.seekp(46, std::ios::cur); //пришли в конец плохого пакета +1

    while (fin) //пока файл не кончится
    {
        fin.read(reinterpret_cast<char *>(&x), 2);
        if (x == 50624)
            return;
    }
    //больше нет пакетов!!!!!! файл кончился!!!!!
    return;
}
//void find_header(std::fstream & fin)


bool is_good(std::fstream & fin)
{
    unsigned short header;
    char           data[42];
    unsigned short sum;

    fin.read(reinterpret_cast<char*>(&header), 2);
    if (fin.eof())  // проверка на конец файла
    {
        return 0;
    }
    fin.read(data, 42);
    fin.read(reinterpret_cast<char*>(&sum), 2); 
    
    if (header != 50624)
    {
        //вывод в файл с ошибками!!!!!!!!!!!!!!!
        if ((fin.tellg())!=46)   //если нули в начале файла
            std::cout << "After pack " << counter << " header error" << std::endl;
        find_header(fin);
        return 0;
    }
    
    if (sum != Crc16(reinterpret_cast<unsigned char *>(data), 42))
    {
        unsigned short next;
        fin.read(reinterpret_cast<char *>(&next), 2);
        fin.seekp(-2, std::ios::cur);

        if (next == 50624)
        {
            //целый пакет но битый!!!!!!
            std::cout << "After pack " << counter << " control sum error" << std::endl;
            return 0;
        }
        
        else
        {
            //после пакета не идет сразу следующий
            find_header(fin); //перемещаем курсор к следующему заголовку
            return 0;
        }
    }
    return 1;
}
//bool is_good(std::fstream & fin)


void check_conditions(unsigned short c)
{
    //int nbit = 0; // - номер бита
    if (c & 1) // бит D0
    {
        // равен 1
    }
    c = c >> 1;

    if (c & 1) // бит D1
    {}
    c = c >> 1;

    if (c & 1) // бит D2
    {}
    c = c >> 1;

    if (c & 1) // бит D3
    {}
    c = c >> 1;

    if (c & 1) // бит D4
    {}
    c = c >> 1;

    if (c & 1) // бит D5
    {}
    c = c >> 1;

    if (c & 1) // бит D6
    {}
    c = c >> 1;

    if (c & 1) // бит D7
    {}
    c = c >> 1;

    if (c & 1) // бит D8
    {}
    c = c >> 1;

    if (c & 1) // бит D9
    {}
    c = c >> 1;

    if (c & 1) // бит D10
    {}
    c = c >> 1;

    if (c & 1) // бит D11
    {}
    c = c >> 1;

    if (c & 1) // бит D12
    {}
    c = c >> 1;

    if (c & 1) // бит D13
    {}
    c = c >> 1;

    if (c & 1) // бит D14
    {}

}
//void check_conditions(unsigned short c)


void print_header(std::fstream & fout, std::map<std::string, bool> & output_flags)
{
    if (output_flags["decod"])
    {
        fout << std::setw(10) << "Time[s]";
        fout << std::setw(7)  << "Npack";
        fout << std::setw(10) << "Tsi[mks]";

        fout << std::setw(16) << "dFi_x[rad]"
             << std::setw(16) << "dFi_y[rad]"
             << std::setw(16) << "dFi_z[rad]";
    
        fout << std::setw(14) << "v_x[mps]"
             << std::setw(14) << "v_y[mps]"
             << std::setw(14) << "v_z[mps]";
    
        fout << std::setw(10)  << "Ta1[gC]"
             << std::setw(10)  << "Ta2[gC]"
             << std::setw(10)  << "Ta3[gC]"; 

        fout << std::setw(10)  << "TG11[gC]"
             << std::setw(10)  << "TG12[gC]"
             << std::setw(10)  << "TG13[gC]"
             << std::setw(10)  << "TG21[gC]"
             << std::setw(10)  << "TG22[gC]"
             << std::setw(10)  << "TG23[gC]"
             << std::setw(10)  << "TG31[gC]"
             << std::setw(10)  << "TG32[gC]"
             << std::setw(10)  << "TG33[gC]";
    
        fout << std::setw(10)  << "Tmkd[gC]";
        fout << std::setw(10)  << "Text[gC]";
        fout << std::setw(6)   << "ski";
        fout << std::endl;
    }
    else if (output_flags["decod_corr"])
    {
        fout << std::setw(10) << "Time[s]";
        fout << std::setw(7)  << "Npack";
        fout << std::setw(10) << "Tsi[mks]";

        fout << std::setw(16) << "dFi_x_corr[rad]"
             << std::setw(16) << "dFi_y_corr[rad]"
             << std::setw(16) << "dFi_z_corr[rad]";
    
        fout << std::setw(14) << "v_x_corr[mps]"
             << std::setw(14) << "v_y_corr[mps]"
             << std::setw(14) << "v_z_corr[mps]";
    
        fout << std::setw(10)  << "Ta1[gC]"
             << std::setw(10)  << "Ta2[gC]"
             << std::setw(10)  << "Ta3[gC]"; 

        fout << std::setw(10)  << "TG11[gC]"
             << std::setw(10)  << "TG12[gC]"
             << std::setw(10)  << "TG13[gC]"
             << std::setw(10)  << "TG21[gC]"
             << std::setw(10)  << "TG22[gC]"
             << std::setw(10)  << "TG23[gC]"
             << std::setw(10)  << "TG31[gC]"
             << std::setw(10)  << "TG32[gC]"
             << std::setw(10)  << "TG33[gC]";
    
        fout << std::setw(10)  << "Tmkd[gC]";
        fout << std::setw(10)  << "Text[gC]";
        fout << std::setw(6)   << "ski";
        fout << std::endl;
    }
    else if (output_flags["bins"])
    {
        fout << std::setw(16) << "dFi_x_corr[rad]"
             << std::setw(16) << "dFi_y_corr[rad]"
             << std::setw(16) << "dFi_z_corr[rad]";
    
        fout << std::setw(14) << "v_x_corr[mps]"
             << std::setw(14) << "v_y_corr[mps]"
             << std::setw(14) << "v_z_corr[mps]";
        fout << std::endl;
    }
    else
    {
        if (output_flags["Time"])
            fout << std::setw(10) << "Time[s]";
        if (output_flags["Tsist"])
            fout << std::setw(10) << "Tsist[s]";
        if (output_flags["Tsi"])
            fout << std::setw(10) << "Tsi[mks]";
        if (output_flags["Npack"])
            fout << std::setw(7)  << "Npack";
        if (output_flags["A"])
        {
            fout << std::setw(14) << "AK1" 
                 << std::setw(14) << "AK2"
                 << std::setw(14) << "AK3";
        }
        if (output_flags["M"])
        {
            fout << std::setw(16) << "M0"
                 << std::setw(16) << "M1"
                 << std::setw(16) << "M2"
                 << std::setw(16) << "M3";
        }
        if (output_flags["L"])
        {   
            fout << std::setw(16) << "L0"
                 << std::setw(16) << "L1"
                 << std::setw(16) << "L2"
                 << std::setw(16) << "L3";
        }
        if(output_flags["dFi_r"])
        {
            fout << std::setw(16) << "dFi_x[rad]"
                 << std::setw(16) << "dFi_y[rad]"
                 << std::setw(16) << "dFi_z[rad]";
        }
        if ((output_flags["dFi_r_corr"])&&(output_flags["Model"]))
        {
            fout << std::setw(16) << std::right << "dFi_x_corr"
                 << std::setw(16) << "dFi_y_corr"
                 << std::setw(16) << "dFi_z_corr";
        }
        if(output_flags["Theta"])
        {
            fout << std::setw(16) << "Theta_x[rps]"
                 << std::setw(16) << "Theta_y[rps]"
                 << std::setw(16) << "Theta_z[rps]";
        }
        if ((output_flags["Theta_corr"])&&(output_flags["Model"]))
        {
            fout << std::setw(16) << "Theta_x_corr"
                 << std::setw(16) << "Theta_y_corr"
                 << std::setw(16) << "Theta_z_corr";
        }
        if (output_flags["Omega"])
        {
            fout << std::setw(17) << "omega_x[gph]"
                 << std::setw(17) << "omega_y[gph]"
                 << std::setw(17) << "omega_z[gph]";
        }
        if ((output_flags["Omega_corr"])&&(output_flags["Model"]))
        {
            fout << std::setw(17) << "omega_x_corr"
                 << std::setw(17) << "omega_y_corr"
                 << std::setw(17) << "omega_z_corr";
        }
        if (output_flags["V"])
        {
            fout << std::setw(14) << "v_x_[mps]"
                 << std::setw(14) << "v_y_[mps]"
                 << std::setw(14) << "v_z_[mps]";
        }
        if ((output_flags["V_corr"])&&(output_flags["Model"]))
        {
            fout << std::setw(14) << "v_x_corr"
                 << std::setw(14) << "v_y_corr"
                 << std::setw(14) << "v_z_corr";
        }
        if (output_flags["W"])
        {    
            fout << std::setw(14) << "w_x[mpss]"
                 << std::setw(14) << "w_y[mpss]"
                 << std::setw(14) << "w_z[mpss]";
        }
        if ((output_flags["W_corr"])&&(output_flags["Model"]))
        {    
            fout << std::setw(14) << "w_x_corr"
                 << std::setw(14) << "w_y_corr"
                 << std::setw(14) << "w_z_corr";
        }
        if (output_flags["Ta"])
        {
            fout << std::setw(10) << "Ta1_[gC]"
                 << std::setw(10) << "Ta2_[gC]"
                 << std::setw(10) << "Ta3_[gC]"; 
        }
        if (output_flags["Tmkd"])
            fout << std::setw(10)  << "Tmkd[gC]";
        if (output_flags["T_lg"])
        {
            fout << std::setw(10) << "TlgX1[gC]"
                 << std::setw(10) << "TlgX2[gC]"
                 << std::setw(10) << "TlgY1[gC]"
                 << std::setw(10) << "TlgY2[gC]"
                 << std::setw(10) << "TlgZ1[gC]"
                 << std::setw(10) << "TlgZ2[gC]";
        }
        if (output_flags["T_lg0"])
        {
            fout << std::setw(10) << "TlgX[gC]"
                 << std::setw(10) << "TlgY[gC]"
                 << std::setw(10) << "TlgZ[gC]";
        }
        if (output_flags["Text"])
            fout << std::setw(10) << "Text[gC]";
        if (output_flags["ski"])
            fout << std::setw(6) << "ski";
        if (output_flags["P"])
        {
            fout << std::setw(9) << "P_1"
                 << std::setw(9) << "P_2"
                 << std::setw(9) << "P_3";
        }
        if (output_flags["U"])
        {
            fout << std::setw(9) << "U_1"
                 << std::setw(9) << "U_2"
                 << std::setw(9) << "U_3";
        }
        if (output_flags["I"])
        {
            fout << std::setw(9) << "I_1"
                 << std::setw(9) << "I_2"
                 << std::setw(9) << "I_3"
                 << std::setw(9) << "I_4"
                 << std::setw(9) << "I_5"
                 << std::setw(9) << "I_6";
        }
        fout << std::endl;
    }
}
//void print_header(std::fstream & fout)

void output_corr_dFi(Data & data, std::fstream & fout, const Constants c)
{
    std::array<double,3> Theta = correct_angles(c, data);
    std::array<double,3> dFi_r;
    for (auto i = 0; i < 3; i++)
    {
        dFi_r[i] = Theta[i] * data.Tsi / 1000000;
    }
    fout << std::setw(16) << std::setprecision(12) << std::fixed << dFi_r[0]
         << std::setw(16) << dFi_r[1]
         << std::setw(16) << dFi_r[2];
}


void output_corr_Theta(Data & data, std::fstream & fout, const Constants c)
{
    std::array<double,3> Theta = correct_angles(c, data);
    fout << std::setw(16) << std::setprecision(12) << std::fixed << Theta[0]
         << std::setw(16) << Theta[1]
             << std::setw(16) << Theta[2];
}


void output_corr_Omega(Data & data, std::fstream & fout, const Constants c)
{
    std::array<double,3> Theta = correct_angles(c, data);
    fout << std::setw(17) << std::setprecision(8) << std::fixed << (Theta[0]*180/std::numbers::pi) * 3600
             << std::setw(17) << (Theta[1]*180/std::numbers::pi) * 3600
             << std::setw(17) << (Theta[2]*180/std::numbers::pi) * 3600;
}


void output_corr_V(Data & data, std::fstream & fout, const Constants c)
{
    std::array<double,3> V = correct_V(c, data);
    fout << std::setw(14) << std::setprecision(8) << std::fixed << V[0]
            << std::setw(14) << V[1]
            << std::setw(14) << V[2];
}


void output_corr_W(Data & data, std::fstream & fout, const Constants c)
{
    std::array<double,3> V = correct_V(c, data);
    std::array<double,3> W;
    for (auto i = 0; i < 3; i++)
        W[i] = V[i] * 1000000 / data.Tsi;
    fout << std::setw(14) << std::setprecision(8) << std::fixed << W[0]
            << std::setw(14) << W[1]
            << std::setw(14) << W[2];
}


void output_data(Data & data, std::fstream & fout, std::map<std::string, bool> & output_flags, const Constants & c)
{
    if (output_flags["decod"])
    {
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.time;
        fout << std::setw(7)  << data.Npack;
    
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.Tsi;
    
        fout << std::setw(16) << std::setprecision(12) << std::fixed << data.dFi_r[0]
             << std::setw(16) << data.dFi_r[1]
             << std::setw(16) << data.dFi_r[2];
    
        fout << std::setw(14) << std::setprecision(8) << std::fixed << data.V[0]
             << std::setw(14) << data.V[1]
             << std::setw(14) << data.V[2];
    
        fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.Ta[0]
             << std::setw(10)  << data.Ta[1]
             << std::setw(10)  << data.Ta[2]; 
    
        fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.T_lgX[0]
             << std::setw(10)  << data.T_lgX[0]
             << std::setw(10)  << data.T_lgX[0]
             << std::setw(10)  << data.T_lgY[0]
             << std::setw(10)  << data.T_lgY[0]
             << std::setw(10)  << data.T_lgY[0]
             << std::setw(10)  << data.T_lgZ[0]
             << std::setw(10)  << data.T_lgZ[0]
             << std::setw(10)  << data.T_lgZ[0];

        fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.Tmkd;
        fout << std::setw(10) << (data.T_lgX[0] + data.T_lgY[0] + data.T_lgZ[0])/3;
        fout << std::setw(6) << data.ski;
        fout << std::endl;
    }
    else if (output_flags["decod_corr"])
    {
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.time;
        fout << std::setw(7)  << data.Npack;
    
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.Tsi;
    
        output_corr_dFi(data, fout, c);
    
        output_corr_V(data, fout, c);
    
        fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.Ta[0]
             << std::setw(10)  << data.Ta[1]
             << std::setw(10)  << data.Ta[2]; 
    
        fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.T_lgX[0]
             << std::setw(10)  << data.T_lgX[0]
             << std::setw(10)  << data.T_lgX[0]
             << std::setw(10)  << data.T_lgY[0]
             << std::setw(10)  << data.T_lgY[0]
             << std::setw(10)  << data.T_lgY[0]
             << std::setw(10)  << data.T_lgZ[0]
             << std::setw(10)  << data.T_lgZ[0]
             << std::setw(10)  << data.T_lgZ[0];

        fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.Tmkd;
        fout << std::setw(10) << (data.T_lgX[0] + data.T_lgY[0] + data.T_lgZ[0])/3;
        fout << std::setw(6) << data.ski;
        fout << std::endl;
    }
    else if (output_flags["bins"])
    {
        output_corr_dFi(data, fout, c);
        output_corr_V(data, fout, c);
        fout << std::endl;
    }
    else
    {
        if (output_flags["Time"])
            fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.time;
        if (output_flags["Tsist"])
            fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.Tsist / 1000000.0; 
        if (output_flags["Tsi"])
            fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.Tsi;
        if (output_flags["Npack"])
            fout << std::setw(7)  << data.Npack;
        if (output_flags["A"])
        {
            fout << std::setw(14) << std::setprecision(8) << std::fixed << data.A[0] * LSB::A
                 << std::setw(14) << data.A[1] * LSB::A
                 << std::setw(14) << data.A[2] * LSB::A;
        }
        if (output_flags["M"])
        {
            fout << std::setw(16) << std::setprecision(11) << std::fixed << data.M[0]
                 << std::setw(16) << data.M[1]
                 << std::setw(16) << data.M[2]
                 << std::setw(16) << data.M[3];
        }
        if (output_flags["L"])
        {   
            fout << std::setw(16) << std::setprecision(11) << std::fixed << data.L[0]
                 << std::setw(16) << data.L[1]
                 << std::setw(16) << data.L[2]
                 << std::setw(16) << data.L[3];
        }
        if (output_flags["dFi_r"])
        {
            fout << std::setw(16) << std::setprecision(12) << std::fixed << data.dFi_r[0]
                 << std::setw(16) << data.dFi_r[1]
                 << std::setw(16) << data.dFi_r[2];
        }
        if ((output_flags["dFi_r_corr"])&&(output_flags["Model"]))
            output_corr_dFi(data, fout, c);
        if (output_flags["Theta"])
        {
            fout << std::setw(16) << std::setprecision(12) << std::fixed << data.Theta[0]
                 << std::setw(16) << data.Theta[1]
                 << std::setw(16) << data.Theta[2];
        }
        if ((output_flags["Theta_corr"])&&(output_flags["Model"]))
            output_corr_Theta(data, fout, c);
        if (output_flags["Omega"])
        {
            fout << std::setw(17) << std::setprecision(8) << std::fixed << (data.dFi_r[0]*180/std::numbers::pi) * 3600 / (data.Tsi/1000000)
                 << std::setw(17) << (data.dFi_r[1]*180/std::numbers::pi) * 3600 / (data.Tsi/1000000)
                 << std::setw(17) << (data.dFi_r[2]*180/std::numbers::pi) * 3600 / (data.Tsi/1000000);
        }
        if ((output_flags["Omega_corr"])&&(output_flags["Model"]))
            output_corr_Omega(data, fout, c);
        if (output_flags["V"])
        {
            fout << std::setw(14) << std::setprecision(8) << std::fixed << data.V[0]
                 << std::setw(14) << data.V[1]
                 << std::setw(14) << data.V[2];
        }
        if ((output_flags["V_corr"])&&(output_flags["Model"]))
            output_corr_V(data, fout, c);
        if (output_flags["W"])
        {    
            fout << std::setw(14) << std::setprecision(8) << std::fixed << data.W[0]
                 << std::setw(14) << data.W[1]
                 << std::setw(14) << data.W[2];
        }
        if ((output_flags["W_corr"])&&(output_flags["Model"]))
            output_corr_W(data, fout, c);
        if (output_flags["Ta"])
        {
            fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.Ta[0]
                 << std::setw(10)  << data.Ta[1]
                 << std::setw(10)  << data.Ta[2]; 
        }
        if (output_flags["Tmkd"])
            fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.Tmkd;
        if (output_flags["T_lg"])
        {
            fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.T_lgX[0]
                 << std::setw(10)  << data.T_lgX[1]
                 << std::setw(10)  << data.T_lgY[0]
                 << std::setw(10)  << data.T_lgY[1]
                 << std::setw(10)  << data.T_lgZ[0]
                 << std::setw(10)  << data.T_lgZ[1];
        }
        if (output_flags["T_lg0"])
        {
            fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.T_lgX[0]
                 << std::setw(10)  << data.T_lgY[0]
                 << std::setw(10)  << data.T_lgZ[0];
        }
        if (output_flags["Text"])
            fout << std::setw(10) << (data.T_lgX[0] + data.T_lgY[0] + data.T_lgZ[0])/3;
        if (output_flags["ski"])
            fout << std::setw(6) << data.ski;
        if (output_flags["P"])
        {
            fout << std::setw(9)  << std::setprecision(4) << std::fixed << data.P[0]
                 << std::setw(9)  << data.P[1]
                 << std::setw(9)  << data.P[2];
        }
        if (output_flags["U"])
        {
            fout << std::setw(9)  << std::setprecision(4) << std::fixed << data.U[0]
                 << std::setw(9)  << data.U[1]
                 << std::setw(9)  << data.U[2];
        }
        if (output_flags["I"])
        {
            fout << std::setw(9)  << std::setprecision(4) << std::fixed << data.I[0]
                 << std::setw(9)  << data.I[1]
                 << std::setw(9)  << data.I[2]
                 << std::setw(9)  << data.I[3]
                 << std::setw(9)  << data.I[4]
                 << std::setw(9)  << data.I[5];
        }
        fout << std::endl;
    }
}
//void output_data(Data & data, std::fstream & fout)


void count_data(Data & data, Data & data_old)
{
    data.time = data_old.time + (data.Tsi/1000000);
    data.count_L(data_old);
    data.count_dFi_r(); 
    data.count_V(data_old);
    data.count_W();
    data.count_theta();
    data.process_mi(data_old);
}


void read_data(std::fstream & fin, std::fstream & fout, std::pair<std::map<std::string, bool>, std::map<std::string, double>> & output_params, const Constants & c)
{
    std::map<std::string, bool> output_flags{output_params.first};
    std::map<std::string, double> time_params{output_params.second};
    bool begin{(!time_params["T_beg"])};
    bool end{false};

    print_header(fout, output_flags);

    Pack pack;
    Data data_old;
    Data data;

    if (is_good(fin)) // читаем первый номер пакета, чтобы потом все нормально шло
    {
        fin.seekp(-44, std::ios::cur); // к номеру пакета
        fin.read(reinterpret_cast<char *>(&pack), 42);
        fin.seekp(2, std::ios::cur); // конец
        counter = pack.Npack;
    }

    data_old = Data(pack);
    data_old.process_mi();

    while (fin)
    {
        if (is_good(fin))
        {
            fin.seekp(-44, std::ios::cur);
            fin.read(reinterpret_cast<char *>(&pack), 42);

            if (pack.Npack != (counter+1))
            { 
                if (counter+1 == pack.Npack-1)
                {
                    std::cout << "Missing pack " << counter + 1 << std::endl;
                }
                else if (!((pack.Npack - counter) > (std::numeric_limits<unsigned short>::max()/2)))
                    std::cout << "Missing packs " <<  counter + 1 << " - " << pack.Npack - 1 << std::endl;
            }

            data = Data(pack);
            count_data(data, data_old);
            counter = pack.Npack;

            if (!begin)
                begin = (data.time >= time_params["T_beg"]);
            if (begin)
                output_data(data, fout, output_flags, c);

            fin.seekp(2, std::ios::cur);
            data_old = data;

            if (time_params["T_end"])
                end = (data.time >= time_params["T_end"]);
            if (end)
                break;
        }
    }
    std::cout << "File is processed" << std::endl;
}
//void read_data(std::fstream & fin, std::fstream & fout)


void output_average_data(DataSum & datasum, std::fstream & fout, std::map<std::string, bool> & output_flags)  // без модели
{
    if (output_flags["Time"])
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << datasum.time;
    if (output_flags["Tsist"])
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << datasum.Tsist / 1000000;
    if (output_flags["Tsi"])
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << datasum.Tsi / datasum.Npacks;
    if (output_flags["Npack"])
        fout << std::setw(7)  << datasum.Npack;
    if (output_flags["A"])
    {
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.A[0] * LSB::A / datasum.Npacks
             << std::setw(14) << datasum.A[1] * LSB::A / datasum.Npacks
             << std::setw(14) << datasum.A[2] * LSB::A / datasum.Npacks;
    }
    if (output_flags["M"])
    {
        fout << std::setw(16) << std::setprecision(11) << std::fixed << datasum.M[0] / datasum.Npacks
             << std::setw(16) << datasum.M[1] / datasum.Npacks
             << std::setw(16) << datasum.M[2] / datasum.Npacks
             << std::setw(16) << datasum.M[3] / datasum.Npacks;
    }
    if (output_flags["L"])
    {   
        fout << std::setw(16) << std::setprecision(11) << std::fixed << datasum.L[0] / datasum.Npacks
             << std::setw(16) << datasum.L[1] / datasum.Npacks
             << std::setw(16) << datasum.L[2] / datasum.Npacks
             << std::setw(16) << datasum.L[3] / datasum.Npacks;
    }
    if(output_flags["dFi_r"])
    {
        fout << std::setw(16) << std::setprecision(12) << std::fixed << datasum.dFi_r[0] / datasum.Npacks
             << std::setw(16) << datasum.dFi_r[1] / datasum.Npacks
             << std::setw(16) << datasum.dFi_r[2] / datasum.Npacks;
    }
    if(output_flags["Theta"])
    {
        fout << std::setw(16) << std::setprecision(12) << std::fixed << datasum.Theta[0] / datasum.Npacks
             << std::setw(16) << datasum.Theta[1] / datasum.Npacks
             << std::setw(16) << datasum.Theta[2] / datasum.Npacks;
    }
    if (output_flags["Omega"])                                                
    {
        fout << std::setw(17) << std::setprecision(8) << std::fixed << ((datasum.dFi_r[0]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000)
             << std::setw(17) << ((datasum.dFi_r[1]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000)
             << std::setw(17) << ((datasum.dFi_r[2]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000);
    }
    if (output_flags["V"])
    {
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.V[0] / datasum.Npacks
             << std::setw(14) << datasum.V[1] / datasum.Npacks
             << std::setw(14) << datasum.V[2] / datasum.Npacks;
    }
    if (output_flags["W"])
    {    
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.W[0] / datasum.Npacks
             << std::setw(14) << datasum.W[1] / datasum.Npacks
             << std::setw(14) << datasum.W[2] / datasum.Npacks;
    }
    if (output_flags["Ta"])
    {
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.Ta[0] / datasum.Npacks
             << std::setw(10)  << datasum.Ta[1] / datasum.Npacks
             << std::setw(10)  << datasum.Ta[2] / datasum.Npacks; 
    }
    if (output_flags["Tmkd"])
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.Tmkd / datasum.Npacks;
    if (output_flags["T_lg"])
    {
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.T_lgX[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgX[1] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgY[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgY[1] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgZ[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgZ[1] / datasum.Npacks;
    }
    if (output_flags["T_lg0"])
    {
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.T_lgX[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgY[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgZ[0] / datasum.Npacks;
    }
    if (output_flags["Text"])
        fout << std::setw(10) << (datasum.T_lgX[0] + datasum.T_lgY[0] + datasum.T_lgZ[0])/(3*datasum.Npacks);
    if (output_flags["P"])
    {
        fout << std::setw(9)  << std::setprecision(4) << std::fixed << datasum.P[0] / datasum.Npacks
             << std::setw(9)  << datasum.P[1] / datasum.Npacks
             << std::setw(9)  << datasum.P[2] / datasum.Npacks;
    }
    if (output_flags["U"])
    {
        fout << std::setw(9)  << std::setprecision(4) << std::fixed << datasum.U[0] / datasum.Npacks
             << std::setw(9)  << datasum.U[1] / datasum.Npacks
             << std::setw(9)  << datasum.U[2] / datasum.Npacks;
    }
    if (output_flags["I"])
    {
        fout << std::setw(9)  << std::setprecision(4) << std::fixed << datasum.I[0] / datasum.Npacks
             << std::setw(9)  << datasum.I[1] / datasum.Npacks
             << std::setw(9)  << datasum.I[2] / datasum.Npacks
             << std::setw(9)  << datasum.I[3] / datasum.Npacks
             << std::setw(9)  << datasum.I[4] / datasum.Npacks
             << std::setw(9)  << datasum.I[5] / datasum.Npacks;
    }
    fout << std::endl;
}
//void output_average_data(DataSum & datasum, std::fstream & fout)


void read_average_data(std::fstream & fin, std::fstream & fout, std::pair<std::map<std::string, bool>, std::map<std::string, double>> & output_params)  // без модели
{
    std::map<std::string, bool> output_flags{output_params.first};
    std::map<std::string, double> time_params{output_params.second};
    double taver{time_params["Taverage"]*1000000};
    bool begin{!time_params["T_beg"]};
    bool end{false};

    Pack pack;
    Data data_old;
    Data data;

    if (is_good(fin)) // читаем первый номер пакета, чтобы потом все нормально шло
    {
        fin.seekp(-44, std::ios::cur); // к номеру пакета
        fin.read(reinterpret_cast<char *>(&pack), 42);
        fin.seekp(2, std::ios::cur); // конец
        counter = pack.Npack;
    }

    data_old = Data(pack);

    if (taver < data_old.Tsi)
    {
        fin.seekp(-46, std::ios::cur);
        const Constants c;
        read_data(fin, fout, output_params, c);
        return;
    }

    print_header(fout, output_flags);

    DataSum datasum{taver};

    while (fin)
    {
        if (is_good(fin))
        {
            fin.seekp(-44, std::ios::cur);
            fin.read(reinterpret_cast<char *>(&pack), 42);

            if (pack.Npack != (counter+1))
            { 
                if (counter+1 == pack.Npack-1)
                {
                    std::cout << "Missing pack " << counter + 1 << std::endl;
                }
                else if (!((pack.Npack - counter) > (std::numeric_limits<unsigned short>::max()/2)))
                    std::cout << "Missing packs " <<  counter + 1 << " - " << pack.Npack - 1 << std::endl;
            }

            data = Data(pack);
            count_data(data, data_old);
            counter = pack.Npack;

            if (!begin)
                begin = ((data.time > time_params["T_beg"])&&(abs(data.time-time_params["T_beg"])>(data.Tsi/2000000)));
            if (begin)
                datasum.add_data(data);

            if (datasum.T >= datasum.Taver)
            {
                output_average_data(datasum, fout, output_flags);
                datasum = DataSum(taver);
            }

            fin.seekp(2, std::ios::cur);
            data_old = data;

            if (time_params["T_end"])
                end = ((datasum.time >= (time_params["T_end"]))||(abs(data.time-time_params["T_end"])<(data.Tsi/2000000)));
            if (end)
                break;
        }
    }
    if (datasum.T != 0)
        output_average_data(datasum, fout, output_flags);

    std::cout << "File is processed" << std::endl;
}
//void read_average_data(std::fstream & fin, std::fstream & fout, double taver)


void output_average_data(DataSum_m & datasum, std::fstream & fout, std::map<std::string, bool> & output_flags)  // с моделью
{
    if (output_flags["Time"])
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << datasum.time;
    if (output_flags["Tsist"])
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << datasum.Tsist / 1000000;
    if (output_flags["Tsi"])
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << datasum.Tsi / datasum.Npacks;
    if (output_flags["Npack"])
        fout << std::setw(7)  << datasum.Npack;
    if (output_flags["A"])
    {
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.A[0] * LSB::A / datasum.Npacks
             << std::setw(14) << datasum.A[1] * LSB::A / datasum.Npacks
             << std::setw(14) << datasum.A[2] * LSB::A / datasum.Npacks;
    }
    if (output_flags["M"])
    {
        fout << std::setw(16) << std::setprecision(11) << std::fixed << datasum.M[0] / datasum.Npacks
             << std::setw(16) << datasum.M[1] / datasum.Npacks
             << std::setw(16) << datasum.M[2] / datasum.Npacks
             << std::setw(16) << datasum.M[3] / datasum.Npacks;
    }
    if (output_flags["L"])
    {   
        fout << std::setw(16) << std::setprecision(11) << std::fixed << datasum.L[0] / datasum.Npacks
             << std::setw(16) << datasum.L[1] / datasum.Npacks
             << std::setw(16) << datasum.L[2] / datasum.Npacks
             << std::setw(16) << datasum.L[3] / datasum.Npacks;
    }
    if(output_flags["dFi_r"])
    {
        fout << std::setw(16) << std::setprecision(12) << std::fixed << datasum.dFi_r[0] / datasum.Npacks
             << std::setw(16) << datasum.dFi_r[1] / datasum.Npacks
             << std::setw(16) << datasum.dFi_r[2] / datasum.Npacks;
    }
    if (output_flags["dFi_r_corr"])
    {
        fout << std::setw(16) << std::setprecision(12) << std::fixed << datasum.dFi_r_corr[0] / datasum.Npacks
             << std::setw(16) << datasum.dFi_r_corr[1] / datasum.Npacks
             << std::setw(16) << datasum.dFi_r_corr[2] / datasum.Npacks;
    }
    if(output_flags["Theta"])
    {
        fout << std::setw(16) << std::setprecision(12) << std::fixed << datasum.Theta[0] / datasum.Npacks
             << std::setw(16) << datasum.Theta[1] / datasum.Npacks
             << std::setw(16) << datasum.Theta[2] / datasum.Npacks;
    }
    if(output_flags["Theta_corr"])
    {
        fout << std::setw(16) << std::setprecision(12) << std::fixed << datasum.Theta_corr[0] / datasum.Npacks
             << std::setw(16) << datasum.Theta_corr[1] / datasum.Npacks
             << std::setw(16) << datasum.Theta_corr[2] / datasum.Npacks;
    }
    if (output_flags["Omega"])                                                
    {
        fout << std::setw(17) << std::setprecision(8) << std::fixed << ((datasum.dFi_r[0]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000)
             << std::setw(17) << ((datasum.dFi_r[1]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000)
             << std::setw(17) << ((datasum.dFi_r[2]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000);
    }
    if (output_flags["Omega_corr"])                                                
    {
        fout << std::setw(17) << std::setprecision(8) << std::fixed << ((datasum.dFi_r_corr[0]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000)
             << std::setw(17) << ((datasum.dFi_r_corr[1]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000)
             << std::setw(17) << ((datasum.dFi_r_corr[2]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000);
    }
    if (output_flags["V"])
    {
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.V[0] / datasum.Npacks
             << std::setw(14) << datasum.V[1] / datasum.Npacks
             << std::setw(14) << datasum.V[2] / datasum.Npacks;
    }
    if (output_flags["V_corr"])
    {
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.V_corr[0] / datasum.Npacks
             << std::setw(14) << datasum.V_corr[1] / datasum.Npacks
             << std::setw(14) << datasum.V_corr[2] / datasum.Npacks;
    }
    if (output_flags["W"])
    {    
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.W[0] / datasum.Npacks
             << std::setw(14) << datasum.W[1] / datasum.Npacks
             << std::setw(14) << datasum.W[2] / datasum.Npacks;
    }
    if (output_flags["W_corr"])
    {    
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.W_corr[0] / datasum.Npacks
             << std::setw(14) << datasum.W_corr[1] / datasum.Npacks
             << std::setw(14) << datasum.W_corr[2] / datasum.Npacks;
    }
    if (output_flags["Ta"])
    {
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.Ta[0] / datasum.Npacks
             << std::setw(10)  << datasum.Ta[1] / datasum.Npacks
             << std::setw(10)  << datasum.Ta[2] / datasum.Npacks; 
    }
    if (output_flags["Tmkd"])
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.Tmkd / datasum.Npacks;
    if (output_flags["T_lg"])
    {
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.T_lgX[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgX[1] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgY[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgY[1] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgZ[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgZ[1] / datasum.Npacks;
    }
    if (output_flags["T_lg0"])
    {
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.T_lgX[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgY[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgZ[0] / datasum.Npacks;
    }
    if (output_flags["Text"])
        fout << std::setw(10) << (datasum.T_lgX[0] + datasum.T_lgY[0] + datasum.T_lgZ[0])/(3*datasum.Npacks);
    if (output_flags["P"])
    {
        fout << std::setw(9)  << std::setprecision(4) << std::fixed << datasum.P[0] / datasum.Npacks
             << std::setw(9)  << datasum.P[1] / datasum.Npacks
             << std::setw(9)  << datasum.P[2] / datasum.Npacks;
    }
    if (output_flags["U"])
    {
        fout << std::setw(9)  << std::setprecision(4) << std::fixed << datasum.U[0] / datasum.Npacks
             << std::setw(9)  << datasum.U[1] / datasum.Npacks
             << std::setw(9)  << datasum.U[2] / datasum.Npacks;
    }
    if (output_flags["I"])
    {
        fout << std::setw(9)  << std::setprecision(4) << std::fixed << datasum.I[0] / datasum.Npacks
             << std::setw(9)  << datasum.I[1] / datasum.Npacks
             << std::setw(9)  << datasum.I[2] / datasum.Npacks
             << std::setw(9)  << datasum.I[3] / datasum.Npacks
             << std::setw(9)  << datasum.I[4] / datasum.Npacks
             << std::setw(9)  << datasum.I[5] / datasum.Npacks;
    }
    fout << std::endl;
}
//void output_average_data(DataSum & datasum, std::fstream & fout)


void read_average_data(std::fstream & fin, std::fstream & fout, std::pair<std::map<std::string, bool>, std::map<std::string, double>> & output_params, const Constants & c)  // с моделью
{
    std::map<std::string, bool> output_flags{output_params.first};
    std::map<std::string, double> time_params{output_params.second};
    double taver{time_params["Taverage"]*1000000};
    bool begin{!time_params["T_beg"]};
    bool end{false};

    Pack pack;
    Data data_old;
    Data data;

    if (is_good(fin)) // читаем первый номер пакета, чтобы потом все нормально шло
    {
        fin.seekp(-44, std::ios::cur); // к номеру пакета
        fin.read(reinterpret_cast<char *>(&pack), 42);
        fin.seekp(2, std::ios::cur); // конец
        counter = pack.Npack;
    }

    data_old = Data(pack);

    if (taver < data_old.Tsi)
    {
        fin.seekp(-46, std::ios::cur);
        read_data(fin, fout, output_params, c);
        return;
    }

    print_header(fout, output_flags);

    DataSum_m datasum{taver};

    while (fin)
    {
        if (is_good(fin))
        {
            fin.seekp(-44, std::ios::cur);
            fin.read(reinterpret_cast<char *>(&pack), 42);

            if (pack.Npack != (counter+1))
            { 
                if (counter+1 == pack.Npack-1)
                {
                    std::cout << "Missing pack " << counter + 1 << std::endl;
                }
                else if (!((pack.Npack - counter) > (std::numeric_limits<unsigned short>::max()/2)))
                    std::cout << "Missing packs " <<  counter + 1 << " - " << pack.Npack - 1 << std::endl;
            }

            data = Data(pack);
            count_data(data, data_old);
            counter = pack.Npack;

            if (!begin)
                begin = ((data.time > time_params["T_beg"])&&(abs(data.time-time_params["T_beg"])>(data.Tsi/2000000)));
            if (begin)
                datasum.add_data(data, c);

            if (datasum.T >= datasum.Taver)
            {
                output_average_data(datasum, fout, output_flags);
                datasum = DataSum_m(taver);
            }

            fin.seekp(2, std::ios::cur);
            data_old = data;

            if (time_params["T_end"])
                end = ((datasum.time >= (time_params["T_end"]))||(abs(data.time-time_params["T_end"])<(data.Tsi/2000000)));
            if (end)
                break;
        }
    }
    if (datasum.T != 0)
        output_average_data(datasum, fout, output_flags);

    std::cout << "File is processed" << std::endl;
}
//void read_average_data(std::fstream & fin, std::fstream & fout, double taver)

void get_flag(std::string & line, std::map<std::string, bool> & flags)
{
    std::istringstream iss(line);
    bool value;
    std::string name;
    std::getline(iss, name, '=');
    iss >> value;
    flags.emplace(std::make_pair(name, value));
} 

std::pair<std::map<std::string, bool>, std::map<std::string, double>> read_config(std::fstream & config)
{
    double value;
    std::map<std::string, bool> output_flags;
    std::map<std::string, double> time_params;
    std::string line;

    std::getline(config, line); // форматирование
    std::getline(config, line); // ---------  

    std::getline(config, line); // decod
    get_flag(line, output_flags); 

    std::getline(config, line);
    get_flag(line, output_flags); // decod_corr

    std::getline(config, line);
    get_flag(line, output_flags); // bins

    std::getline(config, line); // --------- 
    std::getline(config, line); // параметры модели
    std::getline(config, line); // ---------    

    std::getline(config, line);
    get_flag(line, output_flags); // model

    std::getline(config, line); // --------- 
    std::getline(config, line); // временные параметры
    std::getline(config, line); // ---------    

    std::getline(config, line); // Taverage=
    std::istringstream iss(line);
    std::string name;
    std::getline(iss, name, '=');
    iss >> value;
    time_params.emplace(std::make_pair(name, value));

    std::getline(config, line); // T_beg=
    std::istringstream iss1(line);
    std::string name1;
    std::getline(iss1, name1, '=');
    iss1 >> value;
    time_params.emplace(std::make_pair(name1, value));

    std::getline(config, line); // T_end=
    std::istringstream iss2(line);
    std::string name2;
    std::getline(iss2, name2, '=');
    iss2 >> value;
    time_params.emplace(std::make_pair(name2, value));

    std::getline(config, line); // ----------
    std::getline(config, line); // параметры вывода
    std::getline(config, line); // ----------

    while(std::getline(config, line))
    {
        get_flag(line, output_flags);
    }
    return std::make_pair(output_flags, time_params);
}
//std::pair<std::map<std::string, bool>, double> read_config(std::fstream & config)

std::string str_from_config(std::fstream & config)
{
    std::string line;
    std::getline(config, line);  // pcfd=(file_path=)
    std::istringstream iss(line);
    std::string name;
    std::getline(iss, name, '=');
    std::getline(iss, name, ' ');
    std::getline(config, line); // --------
    return name;
}
