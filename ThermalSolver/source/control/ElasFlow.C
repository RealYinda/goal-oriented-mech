//
// 文件名:     ElasFlow.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:25:06 2014 $
// 描述　:     静力问题线弹性材料计算流程的实现.
// 类别　:     %Internal File% ( Don't delete this line )
//
// 根新情况：
// update #6: 3处  -2017-05-03
//分别在：@1：initializeLevelIntegrator   @2：getLevelDt    @3：acceptTimeDependentSolution三个函数处各一个


#include "ElasFlow.h"
#include "PatchStrategy.h"
#include<iostream>
#include<fstream>
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

/*************************************************************************
 * 构造函数.
 *************************************************************************/
ElasFlow::ElasFlow(
    const string& object_name,
    tbox::Pointer<algs::StandardComponentPatchStrategy<NDIM> > strategy,
    tbox::Pointer<tbox::Database> input_db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!object_name.empty());
  assert(strategy.getPointer() != NULL);
#endif

  d_patch_strategy = strategy;

  d_solver_db = input_db;
  d_solver_manager = solv::LinearSolverManager<NDIM>::getManager();
  d_solver_s = d_solver_manager->lookupLinearSolver(
      d_solver_db->getDatabase ("SolverT")->getString("solver_name"));
  d_solver_th = d_solver_manager->lookupLinearSolver(
      d_solver_db->getDatabase ("SolverTH")->getString("solver_name"));
  d_solver_E = d_solver_manager->lookupLinearSolver(
      d_solver_db->getDatabase ("SolverE")->getString("solver_name"));
  d_object_name = object_name;
  t_fem_build_matrix =
      tbox::TimerManager::getManager()->getTimer("ELAS::FEM::buildMatrix");
  t_fem_solve =
      tbox::TimerManager::getManager()->getTimer("ELAS::Solver::TOTAL");
  t_fem_post = tbox::TimerManager::getManager()->getTimer("ELAS::FEM::POST");
}

/*************************************************************************
 * 析构函数.
 ************************************************************************/
ElasFlow::~ElasFlow() {}

/*************************************************************************
 * 初始化网格层积分算法: 创建所有计算需要的积分构件.
 *
 * 该函数创建了8个构件. 这些构件所操作的数据片,
 * 由函数 d_patch_strategy->initializeComponent() 指定.
 *
 *************************************************************************/
void ElasFlow::initializeLevelIntegrator(
    tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager) {
  /// 初始化构件.
  d_init_intc = new algs::InitializeIntegratorComponent<NDIM>(
      "INIT", d_patch_strategy, manager);
  // 数值构件: 更新结点坐标.
  d_num_intc_displacement = new algs::NumericalIntegratorComponent<NDIM>(
      "DISPLACEMENT", d_patch_strategy, manager);
  // 数值构件: 计算应力.
  d_num_intc_stress = new algs::NumericalIntegratorComponent<NDIM>(
      "STRESS", d_patch_strategy, manager);
  /// 数值构件：更新恢复应力
  d_num_intc_recovery = new algs::NumericalIntegratorComponent<NDIM>(
      "RECOVERY", d_patch_strategy, manager);
  /// 数值构件：将更新的恢复应力平均后整理输出
  d_num_intc_postprocess = new algs::NumericalIntegratorComponent<NDIM>(
      "POSTPROCESS", d_patch_strategy, manager);
  /// 数值构件：数据采集并做定量分析
  d_num_intc_data_explorer = new algs::NumericalIntegratorComponent<NDIM>(
              "DATAEXPLORE", d_patch_strategy, manager);

  // 数值构件: 计算矩阵.
  d_num_intc_mat = new algs::NumericalIntegratorComponent<NDIM>(
      "MAT", d_patch_strategy, manager);
  // 数值构件: 计算右端项.
  d_num_intc_rhs = new algs::NumericalIntegratorComponent<NDIM>(
      "RHS", d_patch_strategy, manager);
  // 数值构件: 计算载荷.
  d_num_intc_load = new algs::NumericalIntegratorComponent<NDIM>(
      "LOAD", d_patch_strategy, manager);
  // 数值构件: 计算约束.
  d_num_intc_cons = new algs::NumericalIntegratorComponent<NDIM>(
      "CONS", d_patch_strategy, manager);
  // 内存构件： 管理数据片的内存开辟及释放
  d_alloc_data = new algs::MemoryIntegratorComponent<NDIM>(
      "ALLOC", d_patch_strategy, manager);

  //update #6
  // 步长构件： 计算时间步长
  d_dt_update = new algs::DtIntegratorComponent<NDIM>("Dt", d_patch_strategy, manager);

  //update #8 热计算数值构件
  // 数值构件: 计算矩阵.
  th_num_intc_mat = new algs::NumericalIntegratorComponent<NDIM>(
      "TH_MAT", d_patch_strategy, manager);
  // 数值构件: 计算右端项.
  th_num_intc_rhs = new algs::NumericalIntegratorComponent<NDIM>(
      "TH_RHS", d_patch_strategy, manager);
  // 数值构件: 计算载荷.
  th_num_intc_load = new algs::NumericalIntegratorComponent<NDIM>(
      "TH_LOAD", d_patch_strategy, manager);
  // 数值构件: 计算约束.  Thermal_PostProcesing
  th_num_intc_cons = new algs::NumericalIntegratorComponent<NDIM>(
      "TH_CONS", d_patch_strategy, manager);
  th_num_intc_plot = new algs::NumericalIntegratorComponent<NDIM>(
      "TH_PLOT", d_patch_strategy, manager);

  //update #9 电计算数值构件
  // 数值构件: 计算矩阵.
  E_num_intc_mat = new algs::NumericalIntegratorComponent<NDIM>(
      "E_MAT", d_patch_strategy, manager);
  // 数值构件: 计算右端项.
  E_num_intc_rhs = new algs::NumericalIntegratorComponent<NDIM>(
      "E_RHS", d_patch_strategy, manager);
  // 数值构件: 计算约束.  Thermal_PostProcesing
  E_num_intc_cons = new algs::NumericalIntegratorComponent<NDIM>(
      "E_CONS", d_patch_strategy, manager);
  E_num_intc_plot = new algs::NumericalIntegratorComponent<NDIM>(
      "E_PLOT", d_patch_strategy, manager);

 d_Max_T_intc=new algs::ReductionIntegratorComponent<NDIM>(
              "Max_T",MPI_MAX, d_patch_strategy, manager);
d_Max_Stress_intc=new algs::ReductionIntegratorComponent<NDIM>(
              "Stress_T",MPI_MAX, d_patch_strategy, manager);
}

/*************************************************************************
 *  初始化网格层上的数据.
 *
 * 注解: 该函数调用了初值构件（d_init_set_value)，
 * 该构件又进一步自动调用 d_patch_strategy->initializePatchData(),
 * 完成数据片<uval, current>的初始化.
 ************************************************************************/
void ElasFlow::initializeLevelData(
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
    const double init_data_time, const bool initial_time) {
  /// 初始化网格层上的数据.
  d_init_intc->initializeLevelData(level, init_data_time, initial_time);
}

//update #6
/*************************************************************************
 * 计算时间步长.
 *
 * 注解: 该函数调用了步长构件(d_step_size)，
 * 该构件对象又进一步调用 d_patch_strategy->getPatchDt(),
 * 逐个网格片地计算时间步长.
 ************************************************************************/
double ElasFlow::getLevelDt(
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
    const double dt_time, const bool initial_time, const int flag_last_dt,
    const double last_dt) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!level.isNull());
#endif
  return (d_dt_update->getLevelDt(level, dt_time, initial_time,
                                         flag_last_dt, last_dt, false));
}

/*************************************************************************
 * 向前积分一个时间步.
 *
 * 注解: 该函数调用构件，计算矩阵右端项, 设置边界条件, 并求解线性系统.
 * 调用数值构件计算应力,更新结点坐标.
 *
 ************************************************************************/
int ElasFlow::advanceLevel(
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
    const double current_time, const double predict_dt, const double max_dt,
    const double min_dt, const bool first_step, const int step_number,
    double& actual_dt) {
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!level.isNull());
#endif
  const tbox::Pointer<hier::PatchLevel<NDIM> > patch_level = level;

  actual_dt = predict_dt;

  /// 为矩阵向量开辟内存
  d_alloc_data->allocatePatchData(patch_level, current_time + predict_dt);
  /// 获取参数
  tbox::Pointer<PatchStrategy> p_strategy = d_patch_strategy;

  tbox::pout<<"电场方程求解中...... "<<endl;
  t_fem_solve->start();
  double max[3]={0,0,0};
#if 0
  ///////////////////////////////////////////////////////////////////////////////////////////
  //update #9 电计算
  E_num_intc_mat->computing(patch_level, current_time, actual_dt);

  /// 调用数值构件接口函数,计算并组装右端项
  E_num_intc_rhs->computing(patch_level, current_time, actual_dt);

  /// 调用数值构件接口函数,加载约束
  E_num_intc_cons->computing(patch_level, current_time, actual_dt);

  int mat_id_E = p_strategy->getE_MatrixID();
  int vec_id_E = p_strategy->getE_RHSID();
  int sol_id_E = p_strategy->getE_SolutionID();

  d_solver_E->setMatrix(mat_id_E);
  d_solver_E->setRHS(vec_id_E);
  d_solver_E->solve(first_step, sol_id_E, patch_level, d_solver_db->getDatabase ("SolverE"));
  tbox::pout<<"电场方程求解结束，正在进行后处理...... "<<endl;
  E_num_intc_plot->computing(patch_level, current_time, actual_dt, false);
#endif
  /////////////////////////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////////////////////
  tbox::pout<<"热传导方程求解中...... "<<endl;


  //update #8 热计算
  th_num_intc_mat->computing(patch_level, current_time, actual_dt);

  /// 调用数值构件接口函数,计算并组装右端项
  th_num_intc_rhs->computing(patch_level, current_time, actual_dt);

//  /// 调用数值构件接口函数,加载载荷
//  th_num_intc_load->computing(patch_level, current_time, actual_dt);

  /// 调用数值构件接口函数,加载约束
  th_num_intc_cons->computing(patch_level, current_time, actual_dt);

  int mat_id_th = p_strategy->getTh_MatrixID();
  int vec_id_th = p_strategy->getTh_RHSID();
  int sol_id_th = p_strategy->getTh_SolutionID();

  d_solver_th->setMatrix(mat_id_th);
  d_solver_th->setRHS(vec_id_th);
  d_solver_th->solve(first_step, sol_id_th, patch_level, d_solver_db->getDatabase ("SolverTH"));
  tbox::pout<<"结束热传导方程计算，正在进行后处理 "<<endl;
  th_num_intc_plot->computing(patch_level, current_time, actual_dt, false);
  d_Max_T_intc->reduction(&max[0], 1, patch_level, current_time, actual_dt);

  if (tbox::MPI::getRank() == 0){
	ofstream outdata;
	outdata.open("T_max",ios::app);
        outdata<<current_time+actual_dt<<"\t"<<max[0]<<endl;
     }
  /////////////////////////////////////////////////////////////////////////////////////////////
//#if 1
  //计时开始函数
  t_fem_build_matrix->start();
  tbox::pout<<"Compute Cauchy momentum equations...... "<<endl;
  /// 调用数值构件接口函数,计算并组装矩阵
  //该函数会自动调用用户实现的 algs::StandardComponentPatchStrategy::computeOnPatch().
  d_num_intc_mat->computing(patch_level, current_time, actual_dt);
//  cout<<"Matrix is ok "<<endl;
  /// 调用数值构件接口函数,计算并组装右端项
  d_num_intc_rhs->computing(patch_level, current_time, actual_dt);
  //cout<<"RHS is ok "<<endl;
  /// 调用数值构件接口函数,加载载荷
  d_num_intc_load->computing(patch_level, current_time, actual_dt);
//  cout<<"load is ok "<<endl;
  /// 调用数值构件接口函数,加载约束
  d_num_intc_cons->computing(patch_level, current_time, actual_dt);
//  cout<<"constrain is ok "<<endl;
  t_fem_build_matrix->stop();


  int mat_id = p_strategy->getMatrixID();
  int vec_id = p_strategy->getRHSID();
  int sol_id = p_strategy->getSolutionID();

  /// 设置解法器
  d_solver_s->setMatrix(mat_id);
  d_solver_s->setRHS(vec_id);

  tbox::pout<<"solving "<<endl;
  /// 求解
  d_solver_s->solve(first_step, sol_id, patch_level, d_solver_db->getDatabase ("SolverT"));
  t_fem_solve->stop();
  /// 调用数值构件接口函数, 根据位移更新结点坐标.
  d_num_intc_displacement->computing(patch_level, current_time, actual_dt,false);
  /// 调用数值构件接口函数, 计算应力.
  d_num_intc_stress->computing(patch_level, current_time, actual_dt, false);
  tbox::pout<<"recovery "<<endl;
  t_fem_post->start();
  d_num_intc_recovery->computing(patch_level, current_time, actual_dt, false);
  t_fem_post->stop();
  tbox::pout<<"postprocessing "<<endl;
  d_num_intc_postprocess->computing(patch_level, current_time, actual_dt, false);
  tbox::pout<<"dataexplorer "<<endl;
  if(1)
      d_num_intc_data_explorer->computing(patch_level, current_time, actual_dt, false);
  d_Max_Stress_intc->reduction(&max[1], 2, patch_level, current_time, actual_dt);
  if (tbox::MPI::getRank() == 0){
	ofstream outSdata;
	outSdata.open("max_Stress",ios::app);
        outSdata<<current_time+actual_dt<<"\t"<<max[1]<<endl;
	ofstream outDdata;
	outDdata.open("max_disp",ios::app);
	outDdata<<current_time+actual_dt<<"\t"<<max[2]<<endl;
	
     }

//#endif

  actual_dt = predict_dt;
  d_alloc_data->deallocatePatchData(patch_level);

  return (0);
}

//update #6
/*************************************************************************
 * 接收数值解.
 *
 * 注解: 该函数调用复制构件，
 * 将数据片<uval,new>的值复制到数据片<uval,current>中.
 ************************************************************************/
void ElasFlow::acceptTimeDependentSolution(
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
    const double new_time, const bool last_step) {}
