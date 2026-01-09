//
// 文件名:     PoissonLevelStrategy.C
// 软件包:     JAUMIN fem
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     有限元网格层策略类.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "PoissonLevelStrategy.h"
#include "PatchStrategy.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

/*************************************************************************
 * 构造函数.
 *************************************************************************/
PoissonLevelStrategy::PoissonLevelStrategy(
    const string& object_name,
    tbox::Pointer<algs::StandardComponentPatchStrategy<NDIM> > strategy,
    tbox::Pointer<tbox::Database> input_db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!object_name.empty());
  assert(strategy.getPointer() != NULL);
#endif
  /// 创建网格片积分算法.
  d_patch_strategy = strategy;

  /// 从数据库中获取解法器类型.
  d_solver_db = input_db;
  d_solver_manager = solv::LinearSolverManager<NDIM>::getManager();
  d_solver = d_solver_manager->lookupLinearSolver(
      d_solver_db->getString("solver_name"));
  d_object_name = object_name;
}

/*************************************************************************
 * 析构函数.
 ************************************************************************/
PoissonLevelStrategy::~PoissonLevelStrategy() {}

/*************************************************************************
 * 初始化该积分算法: 创建所有计算需要的积分构件.
 *
 * 该函数创建了1个内存构件, 1个初始化构件, 1个规约构件， 4个数值构件.
 * 这些构件所操作的数据片,
 * 由函数 d_patch_strategy->initializeComponent() 指定.
 *
 *************************************************************************/
void PoissonLevelStrategy::initializeLevelIntegrator(
    tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager) {
  //初值构件: 管理当前值数据片的内存以及初始化
  d_init_set_value = new algs::InitializeIntegratorComponent<NDIM>(
      "INIT", d_patch_strategy, manager);
  /// 内存构件: 管理矩阵片，向量片的内存开辟及释放.
  d_alloc_data = new algs::MemoryIntegratorComponent<NDIM>(
      "ALLOC", d_patch_strategy, manager);
  /// 数值构件：组装矩阵.
  d_mat_intc = new algs::NumericalIntegratorComponent<NDIM>(
      "MAT", d_patch_strategy, manager);
  /// 数值构件：组装右端项.
  d_rhs_intc = new algs::NumericalIntegratorComponent<NDIM>(
      "RHS", d_patch_strategy, manager);
  /// 数值构件：后处理计算.
  d_num_intc = new algs::NumericalIntegratorComponent<NDIM>(
      "POST", d_patch_strategy, manager);
  /// 数值构件：处理约束边界条件.
  d_constraint_intc = new algs::NumericalIntegratorComponent<NDIM>(
      "CONS", d_patch_strategy, manager);
  /// 规约构件：计算L2误差.
  d_reduction_intc = new algs::ReductionIntegratorComponent<NDIM>(
      "RED", MPI_SUM, d_patch_strategy, manager);
}

/*************************************************************************
 * 初始化网格层.
 *
 * 注解: 该函数调用了初始化构件完成初始化。
 ************************************************************************/
void PoissonLevelStrategy::initializeLevelData(
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
    const double init_data_time, const bool initial_time) {
  d_init_set_value->initializeLevelData(level, init_data_time, initial_time);
}

/*************************************************************************
 * 向前积分一个时间步.
 *
 * 注解: 该函数调用了内存构件，数值构件，规约构件及解法器等完成时间步进的计算。
 ************************************************************************/
int PoissonLevelStrategy::advanceLevel(
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
    const double current_time, const double predict_dt, const double max_dt,
    const double min_dt, const bool first_step, const int step_number,
    double& actual_dt) {
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!level.isNull());
#endif
  const tbox::Pointer<hier::PatchLevel<NDIM> > patch_level = level;
  /// 从网格片算法类获取参数
  tbox::Pointer<PatchStrategy> p_strategy = d_patch_strategy;

  actual_dt = predict_dt;

  ///  在这里输出插值点信息
  tbox::pout<<"正在进行插值点计算......"<<endl;
  d_num_intc->computing(patch_level, current_time, actual_dt, false);

  actual_dt = predict_dt;

  return (0);
}
