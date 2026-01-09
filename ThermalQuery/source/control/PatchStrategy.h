//
// 文件名:     PatchStrategy.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:45:12 2014 $
// 描述　:     网格片策略类派生类.
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_PatchStrategy
#define included_PatchStrategy

#include "StandardComponentPatchStrategy.h"
#include "JaVisDataWriter.h"
#include "DOFInfo.h"
#include "Serializable.h"

using namespace JAUMIN;

/**
 * @brief 该类是网格片算法类的基类，其他个性化算法集成该类，实现个性化的片算法。
 */

class PatchStrategy : public algs::StandardComponentPatchStrategy<NDIM>,
                      tbox::Serializable {
public:
  /*!
   * @brief 构造函数.
   */
  PatchStrategy();

  /*!
   * @brief 析构函数.
   */
  virtual ~PatchStrategy();

  /// @name 重载基类 algs::StandardComponentPatchStrategy<NDIM> 的函数:
  // @{

  /*!
   * @brief 初始化指定的积分构件.
   *
   * 注册待填充的数据片或待调度内存空间的数据片到积分构件.
   *
   * @param component 输入参数, 指针, 指向待初始化的积分构件对象.
   */
  virtual void initializeComponent(
      algs::IntegratorComponent<NDIM>* component) const;

  /**
   * @brief 初始化数据片（支持有限元初值构件）.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示初始化的时刻.
   * @param initial_time   输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
   * @param component_name 输入参数, 字符串, 当前调用该函数的初值构件之名称.
   */
  virtual void initializePatchData(hier::Patch<NDIM>& patch, const double time,
                                   const bool initial_time,
                                   const string& component_name);

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
  virtual void computeOnPatch(hier::Patch<NDIM>& patch, const double time,
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
  virtual void reduceOnPatch(double* vector, int len, hier::Patch<NDIM>& patch,
                             const double time, const double dt,
                             const string& component_name);

  /*!
   * @brief 支撑指定名称的有限元步长构件, 在单个网格片上计算时间步长.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param initial_time   输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
   * @param flag_last_dt   输入参数, 整型, 表示前次返回的状态.
   * @param last_dt        输入参数, 双精度浮点型, 表示前次的时间步长.
   * @param component_name 输入参数, 字符串, 表示步长构件的名称.
   *
   * @return 双精度浮点型, 表示时间步长.
   */
  virtual double getPatchDt(hier::Patch<NDIM>& patch, const double time,
                            const bool initial_time, const int flag_last_dt,
                            const double last_dt, const string& component_name);

  /*!
   * @brief 将数据成员输出到重启动数据库.
   * @param db 输入参数, 指针, 指向重启动数据库.
   */
  virtual void putToDatabase(tbox::Pointer<tbox::Database> db);

  /**
   * @brief 注册可视化数据.
   */
  virtual void registerPlotData(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer);

  /**
   * @brief 从输入数据库中读如参数。
   *
   * @param db 输入参数，指针，指向输入数据库。
   */
  virtual void getFromInput(tbox::Pointer<tbox::Database> db);

  /**
   * @brief  获取矩阵id
   *
   * @return 整型，矩阵id
   */
  virtual int getMatrixID() = 0;

  /**
   * @brief  获取右端项id
   *
   * @return 整型，右端项id
   */
  virtual int getRHSID() = 0;

  /**
   * @brief  获取解向量id
   *
   * @return 整型，解向量id
   */
  virtual int getSolutionID() = 0;

  /**
   * @brief 获取自由度信息
   *
   *
   * @return 指针，指向自由度信息
   */
  virtual tbox::Pointer<solv::DOFInfo<NDIM> > getDOFInfo() = 0;

private:
};
#endif
