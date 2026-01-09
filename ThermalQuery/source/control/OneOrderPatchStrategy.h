//
// 文件名:     OneOrderPatchStrategy.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:45:12 2014 $
// 描述　:     网格片策略类派生类，线性元网格片算法类.
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_OneOrderPatchStrategy
#define included_OneOrderPatchStrategy

#include "PatchStrategy.h"
#include "JaVisDataWriter.h"
#include "DOFInfo.h"

using namespace JAUMIN;

class OneOrderPatchStrategy : public PatchStrategy {
public:
  /*!
   * @brief 构造函数.
   *
   * @param object_name          输入参数, 字符串, 对象名称。
   * @param is_from_restart      输入参数，逻辑型，是否为重启动计算。
   * @param input_db             输入参数，指针，指向输入数据库。
   *
   */
  OneOrderPatchStrategy(const string& object_name, bool is_from_restart,
                        tbox::Pointer<tbox::Database> input_db);

  /*!
   * @brief 析构函数.
   */
  virtual ~OneOrderPatchStrategy();

  /// @name 重载基类 algs::StandardComponentPatchStrategy<NDIM> 的函数:
  // @{

  /*!
   * @brief 初始化指定的积分构件.
   *
   * 注册待填充的数据片或待调度内存空间的数据片到积分构件.
   *
   * @param component 输入参数, 指针, 指向待初始化的积分构件对象.
   */
  void initializeComponent(algs::IntegratorComponent<NDIM>* component) const;

  /**
   * @brief 初始化数据片（支持有限元初值构件）.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示初始化的时刻.
   * @param initial_time   输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
   * @param component_name 输入参数, 字符串, 当前调用该函数的初值构件之名称.
   */
  void initializePatchData(hier::Patch<NDIM>& patch, const double time,
                           const bool initial_time,
                           const string& component_name);

  /*!
   * @brief 支撑指定名称的数值构件, 在单个网格片上完成矩阵组装.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void buildMatrixOnPatch(hier::Patch<NDIM>& patch, const double time,
                          const double dt, const string& component_name);

  /*!
   * @brief 支撑指定名称的数值构件, 在单个网格片上完成有限元积分.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void buildRHSOnPatch(hier::Patch<NDIM>& patch, const double time,
                       const double dt, const string& component_name);

  /*!
   * @brief 支撑指定名称的数值构件, 在单个网格片上完成约束加载.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void applyConstraint(hier::Patch<NDIM>& patch, const double time,
                       const double dt, const string& component_name);

  /*!
   * @brief 支撑指定名称的数值构件, 在单个网格片上完成后处理计算.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void postProcess(hier::Patch<NDIM>& patch, const double time, const double dt,
                   const string& component_name);
  /**
   * @brief 注册模型变量, 该函数调用EntisyStrategy的注册函数,
   *完成应户变量的注册.
   *
   */
  void registerModelVariable();

  /*!
   * @brief 支撑指定名称的有限元数值构件, 在单个网格片上完成后数值计算.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param initial_time   输入参数, BOOL型, 是否为初始时刻.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void computeOnPatch(hier::Patch<NDIM>& patch, const double time,
                      const double dt, const bool initial_time,
                      const string& component_name);

  /**
   * @brief 支撑指定名称的归约构件, 在单个网格片上执行归约计算.
   *
   * @param vector         输入输出参数, 指针, 指向归约向量.
   * @param len            输入参数, 整型, 归约向量的长度.
   * @param patch          输入参数, 网格片类, 待归约的网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示当前时间步长.
   * @param component_name 输入参数, 字符串, 归约构件的名称.
   *
   * @note
   *  vector是输入输出参数: 输入是已遍历网格片的归约结果,
   *  输出是输入值和当前网格片计算结果的归约值.
   */
  void reduceOnPatch(double* vector, int len, hier::Patch<NDIM>& patch,
                     const double time, const double dt,
                     const string& component_name);

  /*!
   * @brief 将数据成员输出到重启动数据库.
   * @param db 输入参数, 指针, 指向重启动数据库.
   */
  void putToDatabase(tbox::Pointer<tbox::Database> db);

  /**
   * @brief 注册可视化数据.
   */
  void registerPlotData(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer);

  /**
   * @brief 从输入数据库中读如参数。
   *
   * @param db 输入参数，指针，指向输入数据库。
   */
  void getFromInput(tbox::Pointer<tbox::Database> db);

  /**
   * @brief  获取矩阵id
   *
   * @return 整型，矩阵id
   */
  int getMatrixID() { return d_matrix_id; }

  /**
   * @brief  获取右端项id
   *
   * @return 整型，右端项id
   */
  int getRHSID() { return d_rhs_id; }

  /**
   * @brief  获取解向量id
   *
   * @return 整型，解向量id
   */
  int getSolutionID() { return d_solution_id; }

  /**
   * @brief 获取自由度信息
   *
   *
   * @return 指针，指向自由度信息
   */
  tbox::Pointer<solv::DOFInfo<NDIM> > getDOFInfo() { return d_dof_info; }

private:
  /*!@brief 对象名.  */
  string d_object_name;

  /// 自由度信息
  tbox::Pointer<solv::DOFInfo<NDIM> > d_dof_info;

  /// 有限元计算的形函数类型, 单元类型, 积分器类型.
  string d_element_type;
  string d_file_name_query;

  /// For TetQuad
  int d_Cell_volume_id;
  int d_Cell_jacobian_id;

  /// 对应于变量的id, 通过id应户可以获取变量数据.
  int d_solution_id;
  int d_plot_id;
  int d_error_id;
  int d_matrix_id;
  int d_rhs_id;
};
#endif
