#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>

// Linux/Unix 系统调用头文件
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

// 定义数据结构
struct ThermalData {
    int id1;
    int id2;
    double x, y, z;
    double value;
};

// 辅助函数：判断文件名是否符合 "Thermal_Data_On_Patch_n.dat" 格式
bool isTargetFile(const std::string& filename) {
    std::string prefix = "Thermal_Data_On_Patch_";
    std::string suffix = ".dat";

    if (filename.length() <= prefix.length() + suffix.length()) return false;
    if (filename.compare(0, prefix.length(), prefix) != 0) return false;
    if (filename.compare(filename.length() - suffix.length(), suffix.length(), suffix) != 0) return false;

    return true;
}

int main() {
    // ------------------- 路径配置 -------------------
    // 输入路径: ../ThermalSolver/build
    std::string inputDir = "../ThermalSolver/build"; 
    
    // 输出路径: ../ThermalSolver/output (根据你的新要求修改)
    std::string outputDir = "../ThermalSolver/output";
    
    // 输出文件: ../ThermalSolver/output/ThermalDataOut.dat
    std::string outputFilePath = outputDir + "/ThermalDataOut.dat";
    // ------------------------------------------------

    std::vector<ThermalData> allData;

    // 1. 打开输入目录
    DIR *dir;
    struct dirent *ent;

    std::cout << "正在扫描目录: " << inputDir << " ..." << std::endl;

    if ((dir = opendir(inputDir.c_str())) != NULL) {
        int fileCount = 0;
        while ((ent = readdir(dir)) != NULL) {
            std::string filename = ent->d_name;

            if (isTargetFile(filename)) {
                // 拼接完整路径用于读取
                std::string fullPath = inputDir + "/" + filename;
                
                std::ifstream inFile(fullPath);
                if (!inFile.is_open()) {
                    std::cerr << " [错误] 无法打开文件: " << filename << std::endl;
                    continue;
                }

                ThermalData temp;
                while (inFile >> temp.id1 >> temp.id2 >> temp.x >> temp.y >> temp.z >> temp.value) {
                    allData.push_back(temp);
                }
                fileCount++;
                inFile.close();
            }
        }
        closedir(dir);
        std::cout << "扫描结束，共读取 " << fileCount << " 个数据文件。" << std::endl;
    } else {
        std::cerr << "[致命错误] 无法打开输入目录: " << inputDir << std::endl;
        std::cerr << "请检查路径是否存在。" << std::endl;
        return 1;
    }

    if (allData.empty()) {
        std::cerr << "警告: 未获取到任何数据，请检查输入文件夹。" << std::endl;
        return 0;
    }

    // 2. 排序 (C++11 Lambda)
    std::cout << "正在排序 " << allData.size() << " 条数据..." << std::endl;
    std::sort(allData.begin(), allData.end(), [](const ThermalData& a, const ThermalData& b) {
        if (a.id1 != b.id1) {
            return a.id1 < b.id1; // ID1 升序
        }
        return a.id2 < b.id2;     // ID2 升序
    });

    // 3. 创建输出目录
    // 尝试创建 output 文件夹 (假设 ../ThermalSolver 文件夹已存在)
    // 如果文件夹已存在，mkdir 会返回 -1，代码会继续执行，不受影响
    mkdir(outputDir.c_str(), 0777);

    // 4. 写入文件
    std::cout << "正在写入文件: " << outputFilePath << " ..." << std::endl;
    std::ofstream outFile(outputFilePath, std::ios::out | std::ios::trunc); 

    if (!outFile.is_open()) {
        std::cerr << "[错误] 无法创建输出文件。" << std::endl;
        std::cerr << "请检查路径 " << outputDir << " 是否存在或有写入权限。" << std::endl;
        return 1;
    }

    outFile << std::fixed << std::setprecision(12);

    for (const auto& d : allData) {
        outFile << d.id1 << "\t" 
                << d.id2 << "\t" 
                << d.x << "\t" 
                << d.y << "\t" 
                << d.z << "\t" 
                << d.value << "\n";
    }

    outFile.close();
    std::cout << "成功！数据已保存至: " << outputFilePath << std::endl;

    return 0;
}
