//
// 文件名:     ElasFlow.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:21:15 2014 $
// 描述　:     求解线弹性方程组的流程控制类.
//
// 根新情况：
// update #6: 4处  -2017-05-03
//分别在：@1：添加头文件  @2：getLevelDt    @3：acceptTimeDependentSolution   @4：添加步长构件声明

#ifndef included_ElasFlow
#define included_ElasFlow

#include "Timer.h"
#include "PatchLevel.h"
#include "NumericalIntegratorComponent.h"
#include "InitializeIntegratorComponent.h"
#include "MemoryIntegratorComponent.h"
#include "BaseLinearSolver.h"
#include "StandardComponentPatchStrategy.h"
#include "ReductionIntegratorComponent.h"
#include "LinearSolverManager.h"
#include "TimeIntegratorLevelStrategy.h"

//update #6
#include "DtIntegratorComponent.h"

using namespace JAUMIN;
/**
 * @brief 该类从网格层时间积分算法策略类 algs::TimeIntegratorLevelStrategy
 * 的派生
 * 类 BaseFlow 派生,实现个性化(各向同性材料的线弹性方程组)的有限元计算的流程.
 */
class ElasFlow : public algs::TimeIntegratorLevelStrategy<NDIM>

{
public:
  /**
   * @brief 构造函数.
   * @param object_name   输入参数, 字符串, 表示对象名称.
   * @param strategy      输入参数, 指针, 网格片算法类对象.
   * @param input_db      输入参数，指针，指向输入数据库。
   */
  ElasFlow(const string& object_name,
           tbox::Pointer<algs::StandardComponentPatchStrategy<NDIM> > strategy,
           tbox::Pointer<tbox::Database> input_db);

  /**
   * @brief 析构函数.
   */
  virtual ~ElasFlow();

  ///@name 重载基类algs::TimeIntegratorLevelStrategy<NDIM>的函数
  //@{

  /**
   * @初始化网格层上的数据.
   *
   * @param level           输入参数, 指针, 指向待初始化的网格层.
   * @param init_data_time  输入参数, 双精度浮点型, 初始化数据的时刻.
   * @param initial_time    输入参数, 是否为计算起始时间.
   */
  void initializeLevelData(
      const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
      const double init_data_time, const bool initial_time = true);

  /**
   * @brief 初始化该积分算法: 创建所有计算需要的积分构件.
   * @param manager 输入参数, 指针, 指向积分构件管理器.
   */
  void initializeLevelIntegrator(
      tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager);

  /** update #6
   * @brief 返回指定网格层的时间步长.
   *
   * @param level        输入参数, 指针,         指向网格层.
   * @param dt_time      输入参数, 双精度浮点型, 表示计算时间步长的当前时刻.
   * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
   * @param flag_last_dt 输入参数, 整型,         表示上个时间步积分返回的状态.
   * @param last_dt      输入参数, 双精度浮点型, 表示上个时间步长.
   *
   * @return 双精度浮点型, 表示网格层的时间步长.
   *
   * @note
   * 该函数调用了1个时间步长构件对象，该对象又进一步自动调用函数
   * LinAdv::getPatchDt(), 逐个网格片地计算稳定时间步长.
   */
  double getLevelDt(const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
                    const double dt_time, const bool initial_time,
                    const int flag_last_dt, const double last_dt);

  /**
   * @brief 网格层向前积分一个时间步.
   *
   * @param level        输入参数, 指针,         指向待积分的网格层.
   * @param current_time 输入参数, 双精度浮点型, 表示时间步的起始时刻.
   * @param predict_dt   输入参数, 双精度浮点型, 表示为该时间步预测的时间步长.
   * @param max_dt       输入参数, 双精度浮点型, 表示时间步允许的最大时间步长.
   * @param min_dt       输入参数, 双精度浮点型, 表示时间步允许的最小时间步长.
   * @param first_step   输入参数, 逻辑型, 真值当前为重构后或时间步序列的第1步.
   * @param step_number  输入参数, 整型,         表示积分步数.
   * @param actual_dt    输出参数, 双精度浮点型, 表示时间步实际采用的时间步长.
   *
   * @return 整型, 表示该时间步积分的状态.
   *
   */
  int advanceLevel(const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
                   const double current_time, const double predict_dt,
                   const double max_dt, const double min_dt,
                   const bool first_step, const int step_number,
                   double& actual_dt);

  /**update #6
   * 这个函数没有具体实现，好像用不到
   * @brief 更新网格层的状态到新的时刻.
   *
   * @param level           输入参数, 指针, 指向网格层.
   * @param new_time        输入参数, 双精度浮点型, 表示新的时刻.
   * @param deallocate_data 输入参数, 逻辑型, 真值表示接收数值解后,
   *释放新值数据片的内存空间.
   *
   * @note
   * 该函数调用复制构件, 将数据从新值复制到当前值上下文的数据片.
   */
  void acceptTimeDependentSolution(
      const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
      const double new_time, const bool deallocate_data);

private:
  /*!@brief 对象名称. */
  string d_object_name;

  /*!@brief 解法器参数数据库,管理解法器名称,收敛阈值,最大迭代步数等控制参数. */
  tbox::Pointer<tbox::Database> d_solver_db;

  /*!@brief 网格片算法类对象. */
  tbox::Pointer<algs::StandardComponentPatchStrategy<NDIM> > d_patch_strategy;

  /*!@brief 初始化构件: 初始化物理量 */
  tbox::Pointer<algs::InitializeIntegratorComponent<NDIM> > d_init_intc;

  /*!@brief 内存构件: 为矩阵向量开辟内存 */
  tbox::Pointer<algs::MemoryIntegratorComponent<NDIM> > d_alloc_data;

  /*!@brief 数值构件: 计算矩阵 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_num_intc_mat;

  /*!@brief 数值构件: 计算右端项 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_num_intc_rhs;

  /*!@brief 数值构件: 计算应力 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_num_intc_stress;

  /*!@brief 数值构件: 更新位移 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> >
  d_num_intc_displacement;

  /*!@brief 数值构件: 更新恢复应力 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_num_intc_recovery;
  /*!@brief 数值构件: 后处理恢复应力 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_num_intc_postprocess;
  /*!@brief 数值构件: 数据采集并作定量分析 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_num_intc_data_explorer;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_num_intc_thermal_post;
  /*!@brief 数值构件: 加载约束 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_num_intc_cons;

  /*!@brief 数值构件: 加载载荷 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_num_intc_load;

  /// 解法器管理器及解法器对象
  tbox::Pointer<solv::LinearSolverManager<NDIM> > d_solver_manager;
  tbox::Pointer<solv::BaseLinearSolver<NDIM> > d_solver_s;
  tbox::Pointer<solv::BaseLinearSolver<NDIM> > d_solver_th;
  tbox::Pointer<solv::BaseLinearSolver<NDIM> > d_solver_E;

  //update #6
  //！@brief 步长构件： 计算时间步长
  tbox::Pointer<algs::DtIntegratorComponent<NDIM> > d_dt_update;

  //update #8 热计算数值构件
  /*!@brief 数值构件: 计算矩阵 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > th_num_intc_mat;
  /*!@brief 数值构件: 计算右端项 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > th_num_intc_rhs;
  /*!@brief 数值构件: 加载约束 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > th_num_intc_cons;
  /*!@brief 数值构件: 加载载荷 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > th_num_intc_load;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > th_num_intc_plot;

  //update #9 热计算电值构件
  /*!@brief 数值构件: 计算矩阵 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > E_num_intc_mat;
  /*!@brief 数值构件: 计算右端项 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > E_num_intc_rhs;
  /*!@brief 数值构件: 加载约束 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > E_num_intc_cons;
  /*!@brief 数值构件: 电后处理计算 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > E_num_intc_plot;
  //归约构件
  tbox::Pointer<algs::ReductionIntegratorComponent<NDIM> > d_Max_T_intc;
  tbox::Pointer<algs::ReductionIntegratorComponent<NDIM> > d_Max_Stress_intc;

  //计时器类，记录执行时间
  tbox::Pointer<tbox::Timer> t_fem_build_matrix;
  tbox::Pointer<tbox::Timer> t_fem_solve;
  tbox::Pointer<tbox::Timer> t_fem_post;
};
#endif
