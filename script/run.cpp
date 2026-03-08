#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <dirent.h>     // 用于遍历目录
#include <sys/stat.h>   // 用于检查文件属性和创建目录
#include <unistd.h>     // 用于文件操作
#include <vector>

// 检查文件或目录是否存在
bool exists(const std::string& path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

// 递归删除目录 (兼容 C API)
int remove_directory(const char *path) {
    DIR *d = opendir(path);
    size_t path_len = std::string(path).length();
    int r = -1;
    if (d) {
        struct dirent *p;
        r = 0;
        while (!r && (p = readdir(d))) {
            int r2 = -1;
            char *buf;
            size_t len;
            if (!std::string(p->d_name).compare(".") || !std::string(p->d_name).compare("..")) continue;
            len = path_len + std::string(p->d_name).length() + 2;
            buf = new char[len];
            struct stat statbuf;
            snprintf(buf, len, "%s/%s", path, p->d_name);
            if (!stat(buf, &statbuf)) {
                if (S_ISDIR(statbuf.st_mode)) r2 = remove_directory(buf);
                else r2 = unlink(buf);
            }
            delete[] buf;
            r = r2;
        }
        closedir(d);
    }
    if (!r) r = rmdir(path);
    return r;
}

// 替换文件中的内容
void replaceInFile(const std::string& filepath, const std::string& mesh_name) {
    if (!exists(filepath)) {
        std::cerr << "错误：找不到输入文件 " << filepath << std::endl;
        return;
    }

    std::vector<std::string> lines;
    std::ifstream in(filepath.c_str());
    std::string line;
    
    // 按行读取并修改
    while (std::getline(in, line)) {
        // 兼容 Windows 格式的文件，剔除末尾的 '\r'
        if (!line.empty() && line[line.length()-1] == '\r') {
            line.erase(line.length()-1);
        }
        
        // 如果该行包含 ".mphtxt"，直接整行替换为新名字
        if (line.find(".mphtxt") != std::string::npos) {
            lines.push_back(mesh_name + ".mphtxt");
        } 
        // 如果该行包含 ".k"，且不是带有 "#" 的注释行
        else if (line.find(".k") != std::string::npos && line.find("#") == std::string::npos) {
            lines.push_back(mesh_name + ".k");
        } 
        // 其他行保持原样
        else {
            lines.push_back(line);
        }
    }
    in.close();

    // 将修改后的内容写回文件
    std::ofstream out(filepath.c_str());
    for (size_t i = 0; i < lines.size(); ++i) {
        out << lines[i] << "\n";
    }
    out.close();
    
    std::cout << "--> 已更新 input 文件中的网格名称为: " << mesh_name << ".mphtxt / " << mesh_name << ".k" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "用法: ./run_simulation <mesh_name> <step> <num>" << std::endl;
        return 1;
    }

    std::string mesh_name = argv[1];
    int step = std::stoi(argv[2]);
    std::string num = argv[3];

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(5) << step;
    std::string padded_step = oss.str();
    
    std::string target_T = "T" + padded_step;             
    std::string target_TimeStep = "TimeStep." + padded_step; 

    // Step 1
    std::cout << "\n[Step 1] 修改 input 文件并运行 MeshGen..." << std::endl;
    std::string input_file = "../mphtxt2k/input";
    replaceInFile(input_file, mesh_name);

    if (std::system("cd ../mphtxt2k && ./mesh.out") != 0) {
        std::cerr << "错误：MeshGen 运行失败！" << std::endl; return 1;
    }

    // Step 2
    std::cout << "\n[Step 2] 运行 MPI 求解器..." << std::endl;
    std::string cmd_step2 = "cd ../LinearElasticSolver/build && mpirun -n " + num + " ./main3d ../input/LinearElasticSolverAdaptive.input";
    if (std::system(cmd_step2.c_str()) != 0) {
        std::cerr << "错误：求解器运行失败！" << std::endl; return 1;
    }

    // Step 3
    std::cout << "\n[Step 3] 整理结果文件..." << std::endl;
    std::string results_source_dir = "../LinearElasticSolver/build/javis_Linear_elastic_solver_adaptive_JAD/results/TimeStep.00001";
    std::string results_dest_base = "../LinearElasticSolver/results";
    std::string final_dest_path = results_dest_base + "/" + target_TimeStep;

    if (!exists(results_source_dir)) {
        std::cerr << "错误：找不到文件夹 " << results_source_dir << std::endl; return 1;
    }

    // 3.1 遍历目录并重命名内部文件
    DIR* dir = opendir(results_source_dir.c_str());
    if (dir != NULL) {
        struct dirent* ent;
        while ((ent = readdir(dir)) != NULL) {
            std::string filename = ent->d_name;
            if (filename == "." || filename == "..") continue;
            
            size_t pos = filename.find("T00001");
            if (pos != std::string::npos) {
                std::string new_filename = filename;
                new_filename.replace(pos, 6, target_T);
                std::string old_filepath = results_source_dir + "/" + filename;
                std::string new_filepath = results_source_dir + "/" + new_filename;
                rename(old_filepath.c_str(), new_filepath.c_str());
            }
        }
        closedir(dir);
    }
    std::cout << "--> 文件内部重命名完成。" << std::endl;

    // 3.2 准备目标文件夹
    if (!exists(results_dest_base)) {
        mkdir(results_dest_base.c_str(), 0777); // 尝试创建父目录
    }

    // 如果目标已经存在，删除它
    if (exists(final_dest_path)) {
        remove_directory(final_dest_path.c_str());
    }

    // 3.3 移动文件夹
    if (rename(results_source_dir.c_str(), final_dest_path.c_str()) == 0) {
        std::cout << "--> 文件夹已重命名并移动至: " << final_dest_path << std::endl;
    } else {
        std::cerr << "错误：文件夹移动失败！" << std::endl;
        return 1;
    }

    std::cout << "\n>>> 全部流程执行完毕！<<<" << std::endl;
    return 0;
}
