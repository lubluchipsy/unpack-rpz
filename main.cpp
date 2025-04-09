#include <iostream>
#include <fstream>
#include <math.h>
#include <filesystem>
#include <string>
#include <map>
#include <array>
#include <locale>
#include "read_binary.hpp"


int main()
{
    bool has_config(false);

    std::string this_path = ".";
    for (const auto & entry : std::filesystem::directory_iterator(this_path))
    {
        if (entry.path().filename() == "config.inf")
        {
            has_config = true;
        }
    }
    
    if (!has_config)  // создание config.inf по умолчанию
    {
        std::fstream def_config("config.inf", std::ios::out);
        def_config << "file_path=file.rpz [Относительный путь к файлу .rpz для обработки]" << '\n'
                   << "-------------------------------------------------------------------------------------" << '\n'
                   << "pcfd=p.pcfd        [Название файла с паспортными константами]" << '\n'
                   << "-------------------------------------------------------------------------------------" << '\n'
                   << "Форматирование (остальные параметры не учитываются, если выбран формат):" << '\n'
                   << "-------------------------------------------------------------------------------------" << '\n'
                   << "decod=0      [Формат для программы decod_B5.exe]" << '\n'
                   << "decod_corr=0 [Формат для программы decod_B5.exe со скорр. данными]" << '\n'
                   << "bins=0       [Формат для навигации]" << '\n'
                   << "-------------------------------------------------------------------------------------" << '\n'
                   << "Параметры модели:" << '\n'
                   << "-------------------------------------------------------------------------------------" << '\n'
                   << "Model=0      [1-применять мат. модель, 0-не применять]" << '\n'
                   << "-------------------------------------------------------------------------------------" << '\n'
                   << "Временные параметры:" << '\n'
                   << "-------------------------------------------------------------------------------------" << '\n'
                   << "Taverage=0   [Время усреднения,с. 0, если распаковка без усреднения]" << '\n'
                   << "T_beg=0      [Время начала распаковки,с. 0, если распаковка с начала]" << '\n'
                   << "T_end=0      [Время конца распаковки,с. 0, если распаковка до конца]" << '\n'
                   << "-------------------------------------------------------------------------------------" << '\n'
                   << "Параметры вывода (1-выводить данные в файл, 0-не выводить):" <<'\n'
                   << "-------------------------------------------------------------------------------------" << '\n'
                   << "Time=1       [Счетчик времени, с]" << '\n'
                   << "Npack=1      [номер пакета]" << '\n'
                   << "A=0          [интегральная скорость (данные с АК), м/с]" << '\n'
                   << "M=0          [интегральный кватернион]" << '\n'
                   << "L=0          [кватернион поворота]" << '\n'
                   << "Tsi=0        [время СИ, мкс]" << '\n'
                   << "Tsist=0      [системное время, с]" << '\n'
                   << "dFi_r=0      [приращение угла за такт, рад (нескорр.)]" << '\n'
                   << "dFi_r_corr=0 [приращение угла за такт, рад (скорр.)]" << '\n'
                   << "Theta=0      [углы поворота, рад/с (нескорр.)]" << '\n'
                   << "Theta_corr=0 [углы поворота, рад/с (скорр.)]" << '\n'
                   << "Omega=1      [углы поворота, град/ч (нескорр.)]" << '\n'
                   << "Omega_corr=0 [углы поворота, град/ч (скорр.)]" << '\n'
                   << "V=1          [скорость, м/с (нескорр.)]" << '\n'
                   << "V_corr=0     [скорость, м/с (скорр.)]" << '\n'
                   << "W=1          [ускорение, м/с^2 (нескорр.)]" << '\n'
                   << "W_corr=0     [ускорение, м/с^2 (скорр.)]" << '\n'
                   << "Ta=0         [температура от термодатчиков АК, град.C]" << '\n'
                   << "Tmkd=0       [температура от термодатчика МКД-2, град.C]" << '\n'
                   << "T_lg=0       [все температуры (6) от термодатчиков ЛГ, град.C]" << '\n'
                   << "T_lg0=0      [только верные температуры(3) от термодатчиков ЛГ, град.C]" << '\n'
                   << "Text=0       [внешняя температура, град.С]" << '\n'
                   << "ski=0        [состояние каналов измерения]" << '\n'
                   << "P=0          [сигналы мощностных фотоприемников ЛГ, мВ]" << '\n'
                   << "U=0          [напряжения на ПК ЛГ, В]" << '\n'
                   << "I=0          [контроль тока разряда ЛГ1-3, мкА]" << '\n';
        def_config.close();
    }

    std::fstream config("config.inf", std::ios::in);

    std::string file_path;
    file_path = str_from_config(config);
    std::filesystem::path path{file_path};

    std::string pcfd;
    pcfd = str_from_config(config);
    std::fstream passport(pcfd, std::ios::in);

    std::pair<std::map<std::string, bool>, std::map<std::string, double>> config_data = read_config(config);
    std::map<std::string, bool> output_flags = config_data.first;
    std::map<std::string, double> time_params = config_data.second;

    Constants constants;
    
    if (output_flags["Model"] || output_flags["decod_corr"] || output_flags["bins"])
        read_pcfd(passport, constants);

    std::fstream fin(path, std::ios::in | std::ios::binary); 

    if (!fin)
    {
        std::cout << "Cannot open file " << file_path << std::endl;
        return 1;
    }

    std::string filename = path.filename().string();

    if (!time_params["Taverage"])
    {
        std::string new_file(filename.substr(0, std::size(filename)-4));
        if (output_flags["decod"])
            new_file += "_decod";
        if (output_flags["decod_corr"])
            new_file += "_decod_corr";
        if (output_flags["bins"])
             new_file += "_bins";
        if ((output_flags["Model"])&&(!output_flags["decod"])&&(!output_flags["decod_corr"])&&(!output_flags["bins"]))
            new_file += "_model";
        std::fstream fout(new_file + ".txt", std::ios::out);
        read_data(fin, fout, config_data, constants);
        fout.close();
    }
    else
    {
        std::string new_file(filename.substr(0, std::size(filename)-4));
        if (output_flags["Model"])
            new_file += "_model";
        std::fstream fout(new_file + "_aver" + ".txt", std::ios::out);

        if (output_flags["Model"])
            read_average_data(fin, fout, config_data, constants);
        else
            read_average_data(fin, fout, config_data);
        fout.close();
    }
    fin.close();

    while(true){}

    return 0;
}