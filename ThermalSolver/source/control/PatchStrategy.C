//
// 文件名:     PatchStrategy.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:26:35 2014 $
// 描述　:     网格片策略类派生类实现.
//

/***************************************************************************************
***************************************更新说明*******************************************
* update #0: 添加自定义文件或其他未分类更新
* update #1: 添加了单元对应的实体编号的数据片：registerModelVariable() initializeComponent()
*            initializePatchData()
* update #2: 向材料库添加材料及测试：buildApplicationLib.h Materail.C and etc
* update #3: 给属于不同材料实体的单元赋不同材料：computeStress() buildMatrixOnPatch()
*            LinearTet.C
* update #4: 添加位移数据片（绘图量），存储节点位移量 registerModelVariable()
* 			 initializeComponent()  initializePatchData()
*            添加函数计算总位移量（暂时不添加，在应力计算函数中计算）
* update #5: 添加第二类边界条件即面积力载荷
* update #6: 该为瞬态求解，需要添加时间项的函数或文件：Newmark-beta方法，
* 			 beta大于等于1/4时，无条件稳定。
* update #7: 温变材料，添加材料温变特性
* update #8: 热力耦合
* test:      测试代码
* by tong
****************************************************************************************
*************************************** 更新情况******************************************
* update #0: 0处
* update #1: 8处  @1：添加头文件  @2：registerModelVariable()中数据片声明
*                @3：registerModelVariable()中变量注册
*                @4：initializePatchData()中变量初始化
*                @5@6：initializePatchData()中变量定义及赋初值
*                @7@8：buildMatrixOnPatch(),其中后面一个为test
* update #2: 0处
* update #3: 3处  @1@2：computeStress()材料相关部分
*                @3：buildMatrixOnPatch()
* update #4: 8处  @1：添加头文件X2  @2：registerModelVariable()中数据片声明
*                @3：registerModelVariable()中变量注册
*                @4：initializePatchData()中变量初始化
*                @5@6：initializePatchData()中变量定义及赋初值
*                @6@7：updateCoordinate()中计算位移量模值
*                @8：registerPlotData()绘图量输出
* update #5: 3处  @1@2@3：applyLoad()中节点载荷和面力载荷部分
* update #6: 处  @1：registerModelVariable()中数据片声明X4
*                @2：registerModelVariable()中变量注册X4
*                @3：initializePatchData()中变量初始化X4
*                @4@5@6@7：buildRHSOnPatch()瞬态时右端项计算
*                @8@9：buildMatrixOnPatch()矩阵计算
*                @10：initializeComponent()步长构件
* test     : 处  @1：updateCoordinate()
****************************************************************************************/

/**************************************记录**********************************************
 *  Time:2017-05-05
 *  Author:Tong
 *  内容：已加入 1）多材料处理， 2）面力载荷， 3）瞬态求解
 *  	 程序测试没问题
 *  问题：1）小形变求解精度
 *       2）瞬态求解准确性没有验证
 *       3）时间离散方法问题：显示or隐式
 *       4）网格节点位移处理问题
 *  计划：1）温变材料
 *       2）热力耦合
 *       3）程序测试及验证
 *  ***************************************************************************************/

#include "Pointer.h"
#include "Array.h"
#include "RestartManager.h"
#include "DoubleVector.h"
#include "CellData.h"
#include "NodeData.h"
#include "EdgeData.h"
#include "NodeVariable.h"
#include "CellVariable.h"//update #1  @1 2017-04-20 by tong
#include "EdgeVariable.h"
#include "FaceVariable.h"
#include "MPI.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "Patch.h"
#include <assert.h>
#include "CellVariable.h"

#include "PatchStrategy.h"
#include "Matrix.h"
#include "VectorVariable.h"
#include "VectorData.h"
#include "CSRMatrixVariable.h"
#include "CSRMatrixData.h"
#include "BaseElement.h"
#include "IntegratorManager.h"
#include "ShapeFunctionManager.h"
#include "MaterialManager.h"
#include "Material.h"
#include "BaseShapeFunction.h"
#include "BaseIntegrator.h"
#include "JAUMIN_Macros.h"
#include "fem/TetQuad.h"
#include "fem/TetGeom.h"
#include "fem/HighOrderNodal.h"

//update #4  @1
#include <cmath>
#include<fstream>
#include <Eigen/Dense>

/*************************************************************************
 * 构造函数.
 *************************************************************************/
PatchStrategy::PatchStrategy(const string& object_name,
                             tbox::Pointer<tbox::Database> input_db,
                             bool register_for_restart) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!object_name.empty());
#endif
  d_object_name = object_name;
  d_registered_for_restart = register_for_restart;
  // 读取从输入文件或重启动文件读入数据.
  bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
  registerModelVariable();
  if (is_from_restart) {
    getFromRestart(input_db);
  } else {
    getFromInput(input_db);
  }
  d_element_manager = ElementManager<NDIM>::getManager();
  /// 设置单元.
  d_element_manager->setElement(d_element_type, d_element_marks);

  // 注册为重启动对象.
  if (d_registered_for_restart) {
    tbox::RestartManager::getManager()->registerRestartItem(d_object_name,
                                                            this);
  }
}

/*************************************************************************
 *
 * 构造函数.
 *
 ************************************************************************/
PatchStrategy::~PatchStrategy() {
  if (d_registered_for_restart) {
    tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
  }
}

/*************************************************************************
 *
 * 注册变量和数据片.
 *
 ************************************************************************/
void PatchStrategy::registerModelVariable() {
  /// 取出有限元变量数据库.
  hier::VariableDatabase<NDIM>* db = hier::VariableDatabase<NDIM>::getDatabase();
  /// 自由度信息，参数介绍：四个bool型表示点，边，面，体上是否有自由度存在
  d_dof_info = new solv::DOFInfo<NDIM>(true, false, false, false);
  //update #8 @1
  d_dof_info_th = new solv::DOFInfo<NDIM>(true, false, false, false);

  /// 取出变量上下文.
  tbox::Pointer<hier::VariableContext> current = db->getContext("CURRENT");

  /// 有限元矩阵型变量.
  tbox::Pointer<pdat::CSRMatrixVariable<NDIM, double> > matrix =
      new pdat::CSRMatrixVariable<NDIM, double>("elas_matrix", d_dof_info);

  /// 有限元向量型变量.
  tbox::Pointer<pdat::VectorVariable<NDIM, double> >
      solution /**< 解向量，位移解 */
      = new pdat::VectorVariable<NDIM, double>("elas_solution", d_dof_info);

  tbox::Pointer<pdat::VectorVariable<NDIM, double> > rhs /**< 右端项 */
      = new pdat::VectorVariable<NDIM, double>("elas_rhs", d_dof_info);

  tbox::Pointer<hier::Variable<NDIM> > stress /**< 应力张量 绘图量*/
      = new pdat::CellVariable<NDIM, double>("elas_stress", 2 * NDIM);

  //update #2  绘图量
   tbox::Pointer<hier::Variable<NDIM> > von_mises_stress /**< von Mises 应力 */
       = new pdat::CellVariable<NDIM, double>("von_mises_stress", 1);

   //绘图量
  tbox::Pointer<hier::Variable<NDIM> > plot /**< 解向量的可视化量 */
      = new pdat::NodeVariable<NDIM, double>("plot", NDIM);

  //update1 @2记录单元所对应的实体编号  不依赖计算结果
  tbox::Pointer<pdat::CellVariable<NDIM, int> > EntityIdOfCell
      = new pdat::CellVariable<NDIM, int>("EntityIdOfCell", 1);

  //update #4: @2 总位移量,绘图量
  tbox::Pointer<pdat::NodeVariable<NDIM, double> > displacement
      = new pdat::NodeVariable<NDIM, double>("displacement", 1);

  //update #6: @1 前两步右端项值、解向量值以及前一步的右端项和解向量值
  //都需要数据通信
  tbox::Pointer<pdat::VectorVariable<NDIM, double> >
      solution_old /**< 前一个时刻的解向量，位移解 */
      = new pdat::VectorVariable<NDIM, double>("elas_solution_old", d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM, double> >
      solution_older /**< 前两个时刻的解向量，位移解 */
      = new pdat::VectorVariable<NDIM, double>("elas_solution_older", d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > rhs_old /**< 右端项前一个时刻的值 */
      = new pdat::VectorVariable<NDIM, double>("elas_rhs_old", d_dof_info);
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > rhs_older /**< 右端项前两个时刻的值 */
      = new pdat::VectorVariable<NDIM, double>("elas_rhs_older", d_dof_info);

  //update #8 @2 定义热求解相关变量  矩阵 右端项 解向量
  /// 矩阵变量， 矩阵
 tbox::Pointer<pdat::CSRMatrixVariable<NDIM, double> > thermal_matrix =
      new pdat::CSRMatrixVariable<NDIM, double>("thermal_matrix", d_dof_info_th);
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > thermal_rhs =
      new pdat::VectorVariable<NDIM, double>("thermal_rhs", d_dof_info_th);
  /// 向量变量， 解向量
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > thermal_solution =
      new pdat::VectorVariable<NDIM, double>("thermal_solution", d_dof_info_th);
  // 结点型变量，可视化数据片 绘图量
  tbox::Pointer<pdat::NodeVariable<NDIM, double> > temperature_plot =
      new pdat::NodeVariable<NDIM, double>("temperature_plot", 1);
  // 前一步温度值  需要数据通信 对于力求解是当前的值
  tbox::Pointer<pdat::NodeVariable<NDIM, double> > temperature_old =
      new pdat::NodeVariable<NDIM, double>("temperature_old", 1);
  // 前两步温度值  需要数据通信  对于力求解时前一步时间的值
  tbox::Pointer<pdat::NodeVariable<NDIM, double> > temperature_older =
      new pdat::NodeVariable<NDIM, double>("temperature_older", 1);

  //update #9  定义电求解相关变量  矩阵 右端项 解向量
  /// 矩阵变量， 矩阵
 tbox::Pointer<pdat::CSRMatrixVariable<NDIM, double> > electric_matrix =
      new pdat::CSRMatrixVariable<NDIM, double>("electric_matrix", d_dof_info_th);
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > electric_rhs =
      new pdat::VectorVariable<NDIM, double>("electric_rhs", d_dof_info_th);
  /// 向量变量， 解向量
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > electric_solution =
      new pdat::VectorVariable<NDIM, double>("electric_solution", d_dof_info_th);
  // 结点型变量，可视化数据片  绘图量
  tbox::Pointer<pdat::NodeVariable<NDIM, double> > voltage_plot =
      new pdat::NodeVariable<NDIM, double>("voltage_plot", 1);
  // 单元型变量，单元上电场的模的平方 需要数据通信
  tbox::Pointer<pdat::CellVariable<NDIM, double> > electric_mag =
      new pdat::CellVariable<NDIM, double>("electric_mag", 1);
  tbox::Pointer<pdat::CellVariable<NDIM, int> > material_num =
      new pdat::CellVariable<NDIM, int>("material_num", 1);

  /// 将变量上下文注册到变量数据库.
  d_solution_id = db->registerVariableAndContext(solution, current, 1);
  d_rhs_id = db->registerVariableAndContext(rhs, current, 1);
  d_matrix_id = db->registerVariableAndContext(matrix, current, 1);

  d_stress_id = db->registerVariableAndContext(stress, current);
  d_plot_id = db->registerVariableAndContext(plot, current);

  //update #2
   d_von_mises_id = db->registerVariableAndContext(von_mises_stress, current);

  //update #1 @3
  d_EntityIdOfCell_id = db->registerVariableAndContext(EntityIdOfCell, current, 1);
  //update #4 @3
  d_displacement_id = db->registerVariableAndContext(displacement, current,1);
  //update #6 @2
  d_solution_old_id = db->registerVariableAndContext(solution_old, current,1);
  d_solution_older_id = db->registerVariableAndContext(solution_older, current,1);
  d_rhs_old_id = db->registerVariableAndContext(rhs_old, current,1);
  d_rhs_older_id = db->registerVariableAndContext(rhs_older, current,1);

  //update #8 @3
  th_solution_id = db->registerVariableAndContext(thermal_solution, current, 1);
  th_rhs_id = db->registerVariableAndContext(thermal_rhs, current, 1);
  th_matrix_id = db->registerVariableAndContext(thermal_matrix, current, 1);
  th_plot_id = db->registerVariableAndContext(temperature_plot, current);
  th_Told_id = db->registerVariableAndContext(temperature_old, current,1);
  th_Tolder_id = db->registerVariableAndContext(temperature_older, current,1);

  //update #9
  E_solution_id = db->registerVariableAndContext(electric_solution, current, 1);
  E_rhs_id = db->registerVariableAndContext(electric_rhs, current, 1);
  E_matrix_id = db->registerVariableAndContext(electric_matrix, current, 1);
  E_plot_id = db->registerVariableAndContext(voltage_plot, current);
  E_mag_id = db->registerVariableAndContext(electric_mag, current,1);
  material_num_id = db->registerVariableAndContext(material_num, current,1);


  tbox::Pointer<pdat::CellVariable<NDIM, double> > improved_coefficient =
      new pdat::CellVariable<NDIM, double>("improved_coefficient", 6*(NDIM+1));
  tbox::Pointer<pdat::CellVariable<NDIM, double> > output_coefficient =
      new pdat::CellVariable<NDIM, double>("output_coefficient", 6*(NDIM+1));
  tbox::Pointer<pdat::CellVariable<NDIM, double> > improved_von_mises =
      new pdat::CellVariable<NDIM, double>("improved_von_mises", 1);
  tbox::Pointer<pdat::CellVariable<NDIM, int> > tool_domain =
      new pdat::CellVariable<NDIM, int>("tool_domain", 5);
  tbox::Pointer<pdat::CellVariable<NDIM, int> > contained_domain =
      new pdat::CellVariable<NDIM, int>("contained_domain", 1);

  d_improved_coefficient_id
          = db->registerVariableAndContext(improved_coefficient, current, 1);
  d_output_coefficient_id
          = db->registerVariableAndContext(output_coefficient, current, 1);
  d_improved_von_mises_id
          = db->registerVariableAndContext(improved_von_mises, current, 1);
  d_tool_domain_id
          = db->registerVariableAndContext(tool_domain, current, 1);
  d_contained_domain_id
          = db->registerVariableAndContext(contained_domain, current, 1);
  tbox::Pointer<pdat::CellVariable<NDIM, double> > Cell_jacobian =
      new pdat::CellVariable<NDIM, double>("Cell_jacobian",(NDIM+1)*(NDIM+1),1);//雅可比矩阵
  tbox::Pointer<pdat::CellVariable<NDIM, double> > Cell_volume =
        new pdat::CellVariable<NDIM, double>("Cell_volume",1,1);//单元体积
  d_Cell_jacobian_id=db->registerVariableAndContext(Cell_jacobian, current,1);
  d_Cell_volume_id=db->registerVariableAndContext(Cell_volume, current,1);


}


/*************************************************************************
 *
 *  初始化指定的积分构件.
 *
 *  注册待填充的数据片或待调度内存空间的数据片到积分构件.
 ************************************************************************/
void PatchStrategy::initializeComponent(
    algs::IntegratorComponent<NDIM>* component) const {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(component);
#endif
  // 读取从输入文件或重启动文件读入数据.
  const string& component_name = component->getName();
  if (component_name == "INIT") {  // 初始化构件.
    component->registerInitPatchData(d_stress_id);
    component->registerInitPatchData(d_plot_id);
    //update #2
    component->registerInitPatchData(d_von_mises_id);
    //update #1 @4
    component->registerInitPatchData(d_EntityIdOfCell_id);
    component->registerInitPatchData(material_num_id);
    //update #4 @4
    component->registerInitPatchData(d_displacement_id);
    //update #8 @4
    component->registerInitPatchData(th_plot_id);
    component->registerInitPatchData(th_Told_id);
    component->registerInitPatchData(th_Tolder_id);
    //update #9
    component->registerInitPatchData(E_plot_id);
    component->registerInitPatchData(E_mag_id);

    /// 注册应力恢复数据片
    component->registerInitPatchData(d_improved_coefficient_id);
    component->registerInitPatchData(d_output_coefficient_id);
    component->registerInitPatchData(d_improved_von_mises_id);
    component->registerInitPatchData(d_tool_domain_id);
    component->registerInitPatchData(d_contained_domain_id);

    component->registerInitPatchData(d_Cell_jacobian_id);
    component->registerInitPatchData(d_Cell_volume_id);
    /// 将dofInfo中的数据片注册到初始化构件。
    d_dof_info->registerToInitComponent(component);
    //update #8 @4
    d_dof_info_th->registerToInitComponent(component);
  } else if (component_name == "ALLOC") {  // 内存构件.
    component->registerPatchData(d_matrix_id);
    component->registerPatchData(d_solution_id);
    component->registerPatchData(d_rhs_id);
    //update #6 @3
    component->registerPatchData(d_solution_old_id);
    component->registerPatchData(d_solution_older_id);
    component->registerPatchData(d_rhs_old_id);
    component->registerPatchData(d_rhs_older_id);
    //update #8 @4
    component->registerPatchData(th_matrix_id);
    component->registerPatchData(th_solution_id);
    component->registerPatchData(th_rhs_id);
    //update #9
    component->registerPatchData(E_matrix_id);
    component->registerPatchData(E_solution_id);
    component->registerPatchData(E_rhs_id);
  } else if (component_name == "LOAD") {  // 数值构件，加载载荷.
  } else if (component_name == "CONS") {  // 数值构件，加载约束.
  } else if (component_name == "MAT") {   // 数值构件，计算矩阵。
      component->registerCommunicationPatchData(th_Told_id, th_Told_id);
  } else if (component_name == "RHS") {   // 数值构件，计算右端项.
	  component->registerCommunicationPatchData(d_solution_id, d_solution_id);
	  component->registerCommunicationPatchData(th_Told_id, th_Told_id);
	  component->registerCommunicationPatchData(th_Tolder_id, th_Tolder_id);
  } else if (component_name == "DISPLACEMENT") {  // 数值构件, 更新位移.
	  component->registerCommunicationPatchData(d_solution_id, d_solution_id);
  } else if (component_name == "STRESS") {        // 数值构件, 计算应力.
	  component->registerCommunicationPatchData(d_solution_id, d_solution_id);
  } else if (component_name == "RECOVERY") {        // 数值构件, 计算应力.
      component->registerCommunicationPatchData(d_solution_id, d_solution_id);
      component->registerCommunicationPatchData(d_contained_domain_id, d_contained_domain_id);
      component->registerCommunicationPatchData(th_Told_id, th_Told_id);
      component->registerCommunicationPatchData(th_Tolder_id, th_Tolder_id);
  } else if (component_name == "POSTPROCESS") {        // 数值构件, 计算应力.
      component->registerCommunicationPatchData(d_improved_coefficient_id, d_improved_coefficient_id);
  }else if (component_name == "DATAEXPLORE") {        // 数值构件, 计算应力.
  }  else if (component_name == "Dt") {        // update #6 @10  步长构件,
	//update #8 @5 添加热求解相关数值构件
  } else if (component_name == "TH_LOAD") {  // 数值构件，加载热源等.
  } else if (component_name == "TH_CONS") {  // 数值构件，加载固定温度边界， 即第一类边界.
  } else if (component_name == "TH_MAT") {   // 数值构件，计算热求解矩阵。
  } else if (component_name == "TH_RHS") {   // 数值构件，计算热求解右端项.
	  component->registerCommunicationPatchData(th_Told_id, th_Told_id);
	  component->registerCommunicationPatchData(E_mag_id, E_mag_id);
  } else if (component_name == "TH_PLOT") {//Thermal_PostProcesing
	  component->registerCommunicationPatchData(th_solution_id, th_solution_id);
	  component->registerCommunicationPatchData(th_Told_id, th_Told_id);
      component->registerCommunicationPatchData(th_Tolder_id, th_Tolder_id);
	//update #9 添加电求解相关数值构件
  } else if (component_name == "E_CONS") {  // 数值构件，加载电压边界， 即第一类边界.
  } else if (component_name == "E_MAT") {   // 数值构件，计算电求解矩阵。
  } else if (component_name == "E_RHS") {
	  component->registerCommunicationPatchData(th_Told_id, th_Told_id);
  } else if (component_name == "E_PLOT") {// 数值构件，计算电求解右端项.Electric_PostProcesing
	  component->registerCommunicationPatchData(E_solution_id, E_solution_id);
  }else if (component_name == "Max_T")        {// 规约构件，计算最大温度
  }else if (component_name == "Stress_T")        {// 规约构件，计算最大yingli
  }	else {
    TBOX_ERROR("\n::initializeComponent() : component "
               << component_name << " is not matched. " << endl);
  }
}

/*************************************************************************
 *  初始化数据片（支持初值构件）.
 ************************************************************************/
void PatchStrategy::initializePatchData(hier::Patch<NDIM>& patch,
                                        const double time,
                                        const bool initial_time,
                                        const string& component_name) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(component_name == "INIT");
#endif
  NULL_USE(time); /**< 初始化中没有用到time */

  if (initial_time) {
    tbox::Pointer<pdat::NodeData<NDIM, double> > plot =
        patch.getPatchData(d_plot_id);
    tbox::Pointer<pdat::CellData<NDIM, double> > str =
        patch.getPatchData(d_stress_id);

    //update #2
    tbox::Pointer<pdat::CellData<NDIM, double> > von_Mises =
        patch.getPatchData(d_von_mises_id);

    //update #8 @6
    tbox::Pointer<pdat::NodeData<NDIM, double> > th_plot =
            patch.getPatchData(th_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> > th_Told =
                patch.getPatchData(th_Told_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> > th_Tolder =
                patch.getPatchData(th_Tolder_id);

    //update #1 @5
    tbox::Pointer<pdat::CellData<NDIM, int> > entityid_data =
        patch.getPatchData(d_EntityIdOfCell_id);
    tbox::Pointer<pdat::CellData<NDIM, int> > materialid_data =
        patch.getPatchData(material_num_id);

    //update #4 @5
    tbox::Pointer<pdat::NodeData<NDIM, double> >disp_node =
        patch.getPatchData(d_displacement_id);

    //update #9
    tbox::Pointer<pdat::NodeData<NDIM, double> >voltage_node =
        patch.getPatchData(E_plot_id);
    tbox::Pointer<pdat::CellData<NDIM, double> > electric_data =
        patch.getPatchData(E_mag_id);
    tbox::Pointer<pdat::CellData<NDIM, double> > tool_domain =
        patch.getPatchData(d_tool_domain_id);
    tbox::Pointer<pdat::CellData<NDIM, double> > improved_coef =
        patch.getPatchData(d_improved_coefficient_id);
    tbox::Pointer<pdat::CellData<NDIM, double> > plot_coef =
        patch.getPatchData(d_output_coefficient_id);

    tbox::Pointer<pdat::CellData<NDIM, int> > is_on_bnd =
        patch.getPatchData(d_contained_domain_id);
    DECLARE_ADJACENCY(patch,cell,face,Cell,Face);
    DECLARE_ADJACENCY(patch,face,cell,Face,Cell);
    DECLARE_ADJACENCY(patch,cell,cell,Cell,Cell);


    tbox::Pointer<pdat::CellData<NDIM, double> > jacobian =
            patch.getPatchData(d_Cell_jacobian_id);
    tbox::Pointer<pdat::CellData<NDIM, double> > volume=
            patch.getPatchData(d_Cell_volume_id);
    JAUMIN::appu::TetGeom tetrahedron(patch);



    // 获取当前网格片的单元，结点数目.
    int num_nodes = patch.getNumberOfNodes(1);
    int num_local_nodes = patch.getNumberOfNodes();
    int num_local_cells = patch.getNumberOfCells();
    int num_cells = patch.getNumberOfCells(1);

    /// 获取自由度分布和映射数组的指针首地址
    int* dis_ptr = d_dof_info->getDOFDistribution(patch);
    //update #8 @6
    int* th_ptr = d_dof_info_th->getDOFDistribution(patch);

    /// 为分布和映射数组赋值
    for (int i = 0; i < num_nodes; ++i) {
      dis_ptr[i] = NDIM;
      th_ptr[i]=1; //update #8 @6
    }

    /// 建立映射
    d_dof_info->buildPatchDOFMapping(patch);

    d_dof_info_th->buildPatchDOFMapping(patch);//update #8 @6

    /// 初始化tool domain信息
    for(int cc = 0; cc < num_local_cells; cc++){
        (*volume)(0, cc) = tetrahedron.volume(cc);
        tetrahedron.jacobian(cc, &((*jacobian)(0, cc)));
        bool is_computed = false;
        /// 遍历需要关注的域
        for(int von_id = 0; von_id < d_improved_stress.size();von_id++){
            if (HAS_ENTITY_SET(patch, d_improved_stress[von_id], CELL, 1)){
                DECLARE_ENTITY_SET(patch, improved_stress_space,
                                   d_improved_stress[von_id], CELL, 1);
                std::sort(improved_stress_space.getPointer(),
                          improved_stress_space.getPointer() + improved_stress_space.size());
                bool is_computed_cell =
                        std::binary_search(improved_stress_space.getPointer(),
                                           improved_stress_space.getPointer()+improved_stress_space.size(),
                                           cc);
                bool is_neighbor_computed = false;
                for(int c2 = 0; c2<cell_cell_ext[cc+1]-cell_cell_ext[cc]; c2++){
                    int cc2 = cell_cell_idx[cell_cell_ext[cc]+c2];
                    is_neighbor_computed
                            = std::binary_search(improved_stress_space.getPointer(),
                                                 improved_stress_space.getPointer()+improved_stress_space.size(),
                                                 cc2);
                    /// 一旦它变成true就直接退出
                    if (is_neighbor_computed == true)
                        break;
                }
                /// 任意一个满足就置true
                if(is_computed_cell||is_neighbor_computed) is_computed = true;
            }

        }
        bool is_tool = true;
        if(cell_cell_ext[cc+1]-cell_cell_ext[cc]<4) is_tool = false;
        /// 可以被当作工具域中心，且需要计算
        if(is_tool && is_computed) (*is_on_bnd)(0,cc) = 1;
        /// 不能被当作工作域中心，但需要计算
        if(!is_tool && is_computed) (*is_on_bnd)(0,cc) = 0;
        /// 不需要计算
        if(!is_computed) (*is_on_bnd)(0,cc) = -1;
    }

    /// 初始化普通数据片
    for (int i = 0; i < num_local_nodes; ++i) {
      for (int j = 0; j < NDIM; ++j) (*plot)(j, i) = 0;
      (*th_plot)(0,i)=293.15;//update #8 @6
      (*voltage_node)(0,i)=0;
    }

    for (int i = 0; i < num_local_cells; ++i) {
      for (int j = 0; j < 6; ++j) (*str)(j, i) = 0.0;
      (*von_Mises)(0,i)=0;
      for(int j = 0; j < 6 * (NDIM+1); j ++){
          (*improved_coef)(j,i) = 0.;
          (*plot_coef)(j,i) = 0.;
      }
    }

    for (int i = 0; i < num_nodes; ++i) {
      (*th_Told)(0,i)=293.15;//update #8 @6
      (*th_Tolder)(0,i)=293.15;//update #8 @6
    }
    //cout<<"initial half2"<<endl;
    for (int iii = 0; iii < num_cells; ++iii) {
       (*electric_data)(0,iii) = 0.0;
    }
   // cout<<"initial half3"<<endl;
    //update #4 @6
    for (int i = 0; i < num_nodes; ++i) {
      (*disp_node)(0, i) = 0.0;
    }

   // cout<<"initial half4"<<endl;
    /// 对EntityIdOfCell数据片初始化   update #1 @6

    //int copper_id_len=sizeof(Copper_id)/sizeof(Copper_id[0]);
    //int silicon_id_len=sizeof(Silicon_id)/sizeof(Silicon_id[0]);
    int SiO2_id_len=sizeof(SiO2_id)/sizeof(SiO2_id[0]);
    //int SiN_id_len=sizeof(BCB_id)/sizeof(BCB_id[0]);
    int Copper_id_len=sizeof(Copper_id)/sizeof(Copper_id[0]);
    int Aluminum_id_len=sizeof(Aluminum_id)/sizeof(Aluminum_id[0]);
    int GaN_id_len=sizeof(GaN_id)/sizeof(GaN_id[0]);
    int Al2O3_id_len=sizeof(Al2O3_id)/sizeof(Al2O3_id[0]);
    int Alloy_id_len=sizeof(Alloy_id)/sizeof(Alloy_id[0]);
    int Silicon_id_len=sizeof(Silicon_id)/sizeof(Silicon_id[0]);
    

    for(int e_id=1;e_id<(ENTITY_NUM+1);e_id++)
    {
  	  if (!patch.getPatchGeometry()->hasEntitySet
  	      (e_id, hier::EntityUtilities::CELL, patch.getNumberOfCells(1)))
  	    continue;
	    tbox::Array<int> cells =  patch.getPatchGeometry()->getEntityIndicesInSet
	    	    (e_id, hier::EntityUtilities::CELL, patch.getNumberOfCells(1));
	    for(int ii=0;ii<cells.size();ii++)
	    {
		(*entityid_data)(0,cells[ii])=e_id;
	    }
    }
    //cout<<"initial half2"<<endl;
    for (int iii = 0; iii < num_cells; ++iii) {
        (*materialid_data)(0,iii)=3;//Au
        // for(int m_id=0;m_id<copper_id_len;m_id++)


        for(int m_id=0;m_id<Aluminum_id_len;m_id++)
            if((*entityid_data)(0,iii)==Aluminum_id[m_id])
            {(*materialid_data)(0,iii)=6;}

        for(int m_id=0;m_id<GaN_id_len;m_id++)
            if((*entityid_data)(0,iii)==GaN_id[m_id])
            {(*materialid_data)(0,iii)=15;}


        for(int m_id=0;m_id<Copper_id_len;m_id++)
            if((*entityid_data)(0,iii)==Copper_id[m_id])
            {(*materialid_data)(0,iii)=2;}



        for(int m_id=0;m_id<Al2O3_id_len;m_id++)
            if((*entityid_data)(0,iii)==Al2O3_id[m_id])
            {(*materialid_data)(0,iii)=16;}

        for(int m_id=0;m_id<Alloy_id_len;m_id++)
            if((*entityid_data)(0,iii)==Alloy_id[m_id])
            {(*materialid_data)(0,iii)=2;}
        for(int m_id=0;m_id<Silicon_id_len;m_id++)
            if((*entityid_data)(0,iii)==Silicon_id[m_id])
            {(*materialid_data)(0,iii)=1;}
        for(int m_id=0;m_id<SiO2_id_len;m_id++)
            if((*entityid_data)(0,iii)==SiO2_id[m_id])
            {(*materialid_data)(0,iii)=4;}

                        

      // for(int m_id=0;m_id<SiN_id_len;m_id++)
	//  if((*entityid_data)(0,iii)==BCB_id[m_id])
          //   {(*materialid_data)(0,iii)=9;break;}
    }
  }
}

/*************************************************************************
 *  输出数据成员到重启动数据库.
 ************************************************************************/
void PatchStrategy::putToDatabase(tbox::Pointer<tbox::Database> db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!db.isNull());
#endif
  d_dof_info->putToDatabase(db);
  d_dof_info_th->putToDatabase(db);//update #8 @7
}

/*************************************************************************
 *
 *  设置载荷条件.
 *
 ************************************************************************/
void PatchStrategy::applyLoad(hier::Patch<NDIM>& patch, const double time,
                              const double dt, const string& component_name) {
  /// 取出本地PatchGeometry.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();
  tbox::Pointer<pdat::VectorData<NDIM, double> > rhs_data =
      patch.getPatchData(d_rhs_id);

  //update #5 @1
  //节点力载荷
#if NODE_LOAD
  int load_size = d_load_types.getSize();
  for (int k = 0; k < load_size; ++k) {
    // 获取指定编号和类型的集合包含的网格实体的索引。
    if (patch.hasEntitySet(d_load_marks[k], hier::EntityUtilities::NODE)) {
      const tbox::Array<int>& entity_idx = patch_geo->getEntityIndicesInSet(
          d_load_marks[k], hier::EntityUtilities::NODE);
      // 获取物理边界上的边或者面的数目.
      int size = entity_idx.getSize();
      /// 载荷点的y方向加载0.001N的力
      for (int i = 0; i < size; ++i) {
        rhs_data->getPointer()[NDIM * entity_idx[i] + 1] -= 0;
        ;
      }
    }
  }
#endif

  /// 面力载荷
  //update #5 @2
#if FACE_LOAD
  //update #5 @3
  //自由度信息中的映射信息
  int* dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
  /// 取出本地Patch的结点坐标数组.
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch_geo->getNodeCoordinates();
  int num_faces =patch.getNumberOfEntities(hier::EntityUtilities::FACE, 0);

  double q[3]={0,1,0};//单位面积表面力
  if (patch_geo->hasEntitySet(3, hier::EntityUtilities::FACE)) {
    // 获取指定编号和类型的集合包含的网格实体的索引。
    const tbox::Array<int>& Face_idx = patch_geo->getEntityIndicesInSet(
        3, hier::EntityUtilities::FACE,num_faces);
	  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
			      patch.getPatchTopology();
	  //单元节点邻接关系数据
	  tbox::Array<int>face_node_ext,face_node_idx;
	  patch.getPatchTopology()->getFaceAdjacencyNodes(face_node_ext,face_node_idx);
	  //单元节点邻接关系数据
	  tbox::Array<int>face_cell_ext,face_cell_idx;
	  patch.getPatchTopology()->getFaceAdjacencyCells(face_cell_ext,face_cell_idx);
    int size = Face_idx.getSize();
    for (int face=0; face<size;face++)
    {
    	int n=3;//面单元节点数
    	int dof_num=n*NDIM;
        tbox::Array<int> node_mapping(dof_num);
    	int node_id[3]={0,0,0};

    	//面上三角形的三个顶点编号
    	node_id[0]=face_node_idx[face_node_ext[Face_idx[face]]];//Face_idx[face]面的编号
    	node_id[1]=face_node_idx[face_node_ext[Face_idx[face]]+1];
    	node_id[2]=face_node_idx[face_node_ext[Face_idx[face]]+2];
    	tbox::Array<hier::DoubleVector<NDIM> > vertex(n);
    	for (int k=0; k<n;k++){
    		for(int j=0; j<n;j++){
        		vertex[k][j]=(*node_coord)(j,node_id[k]);
    		}
    	}

        for (int i1 = 0, j = face_node_ext[Face_idx[face]]; i1 < n; ++i1, ++j) {
          for (int k = 0; k < NDIM; ++k) {
            node_mapping[NDIM * i1 + k] = dof_map[face_node_idx[j]] + k;
          }
        }

        tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
        ele_vec->resize(dof_num);
        for (int i = 0; i < dof_num; ++i) {
          (*ele_vec)[i] = 0.0;
        }

    	double area=sqrt(AREA(vertex[0],vertex[1],vertex[2]))/2.0;
    	for(int ii=0;ii<n;ii++)
    		for(int jj=0;jj<NDIM;jj++)
    			(*ele_vec)[ii*NDIM+jj]=q[jj]*area/3;
    	for(int ii=0;ii<dof_num;ii++)
    		rhs_data->getPointer()[node_mapping[ii]] += (*ele_vec)[ii];
    }
  }
#endif
}

/*************************************************************************
 *
 *  填充物理边界条件.
 *
 *  error1:
 *  此处对角化1法处理约束存在问题：
 *  1、 没有区分对待x，y，z三个方向上的位移约束
 *  2、 添加约束后右端项处理存在问题
 ************************************************************************/
void PatchStrategy::applyConstraint(hier::Patch<NDIM>& patch,
                                    const double time, const double dt,
                                    const string& component_name) {
  /// 取出本地PatchGeometry.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();

  int constraint_size = d_constraint_types.getSize();

  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > mat_data =
      patch.getPatchData(d_matrix_id);
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_data =
      patch.getPatchData(d_rhs_id);

  int* dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  int* row_start = mat_data->getRowStartPointer();
  int* col_idx = mat_data->getColumnIndicesPointer();
  double* mat_val = mat_data->getValuePointer();
  double* vec_val = vec_data->getPointer();
  for (int k = 0; k < constraint_size; ++k) {
    if (patch_geo->hasEntitySet(d_constraint_marks[k],
                                hier::EntityUtilities::NODE)) {
      // 获取指定编号和类型的集合包含的网格实体的索引。
      const tbox::Array<int>& entity_idx = patch_geo->getEntityIndicesInSet(
          d_constraint_marks[k], hier::EntityUtilities::NODE);
      int size = entity_idx.getSize();

      //error1:
      /// 对角化1法处理约束
      for (int i = 0; i < size; ++i) {
        for (int k = 0; k < NDIM; ++k) {
          int index = dof_map[entity_idx[i]] + k;
          vec_val[index] = 0.0;
          for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
            if (col_idx[j] == index) {
              mat_val[j] = 1.0;
            } else {
              mat_val[j] = 0.0;
              (*mat_data)(col_idx[j], index) = 0.0;
            }
          }
        }
      }
    }
  }
}

/********************************************************************************
* 获取时间步长
********************************************************************************/
double PatchStrategy::getPatchDt(hier::Patch<NDIM>& patch, const double time,
                                 const bool initial_time,
                                 const int flag_last_dt, const double last_dt,
                                 const string& component_name) {
  return 1e-6;
}

/*************************************************************************
 * 完成单个网格片上的数值计算（支持数值构件）.
 ************************************************************************/
void PatchStrategy::computeOnPatch(hier::Patch<NDIM>& patch, const double time,
                                   const double dt, const bool initial_time,
                                   const string& component_name) {
  if (component_name == "MAT") {
    buildMatrixOnPatch(patch, time, dt, component_name);
  } else if (component_name == "RHS") {
    buildRHSOnPatch(patch, time, dt, component_name);
  } else if (component_name == "LOAD") {
    applyLoad(patch, time, dt, component_name);
  } else if (component_name == "CONS") {
    applyConstraint(patch, time, dt, component_name);
  } else if (component_name == "DISPLACEMENT") {
    updateCoordinate(patch, time, dt, component_name);
  } else if (component_name == "STRESS") {  // 数值构件, 计算应力.
    computeStress(patch, time, dt, component_name);
  } else if (component_name == "RECOVERY") {  // 数值构件, 计算应力.
    StressRecovery(patch, time, dt, component_name);
  } else if (component_name == "POSTPROCESS") {  // 数值构件, 计算应力.
    PostprocessStress(patch, time, dt, component_name);
  } else if (component_name == "DATAEXPLORE") {  // 数值构件, 计算应力.
    Dataexplorer(patch, time, dt, component_name);
  } else if(component_name == "TH_MAT") {
	buildTh_MatrixOnPatch(patch, time, dt, component_name);
  } else if (component_name == "TH_RHS") {
	buildTh_RHSOnPatch(patch, time, dt, component_name);
  } else if (component_name == "TH_LOAD") {
	applyTh_Load(patch, time, dt, component_name);
  } else if (component_name == "TH_CONS") {
    applyTh_Constraint(patch, time, dt, component_name);
  } else if (component_name == "TH_PLOT") {
	Thermal_PostProcesing(patch, time, dt, component_name);
  } else if(component_name == "E_MAT") {
	buildE_MatrixOnPatch(patch, time, dt, component_name);
  } else if (component_name == "E_RHS") {
	buildE_RHSOnPatch(patch, time, dt, component_name);
  } else if (component_name == "E_CONS") {
    applyE_Constraint(patch, time, dt, component_name);
  } else if (component_name == "E_PLOT") {
	Electric_PostProcesing(patch, time, dt, component_name);
  } else {
    TBOX_ERROR(" PatchStrategy :: component name is error! ");
  }
}

/*************************************************************************
 * 完成单个网格片上的规约计算（支持规约构件）.
 ************************************************************************/
void PatchStrategy::reduceOnPatch(double* vector, int len, hier::Patch<NDIM>& patch,
                            const double time, const double dt,
                            const string& component_name) {

  if(component_name == "Max_T")
  {
      return Thermal_max(vector,len,patch, time, dt, component_name);
  }
  else if(component_name == "Stress_T")
  {
      return Stress_max(vector,len,patch, time, dt, component_name);
  }
  TBOX_ERROR("Component \"" << component_name << "\" not matched.\n");
}

/********************************************************************************
* 更新坐标
********************************************************************************/
void PatchStrategy::updateCoordinate(hier::Patch<NDIM>& patch,
                                     const double time, const double dt,
                                     const string& component_name) {
  /// 取出本地PatchGeometry.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();
  /// 取出本地Patch的结点坐标数组.
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch_geo->getNodeCoordinates();
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_data =
      patch.getPatchData(d_solution_id);

  //update #4 @7
  tbox::Pointer<pdat::NodeData<NDIM, double> >disp_node =
      patch.getPatchData(d_displacement_id);

  tbox::Pointer<pdat::NodeData<NDIM, double> > plot_data =
      patch.getPatchData(d_plot_id);

  int num_nodes = patch.getNumberOfNodes(1);

  //update #4 @8 计算总位移量
  for(int i=0;i<num_nodes;i++)
  {
	  double disp_tmp=0;
	  for(int j=0;j<3;j++)
		  disp_tmp+=(vec_data->getPointer()[i*3+j])*(vec_data->getPointer()[i*3+j]);
	  (*disp_node)(0,i)=sqrt(disp_tmp);
  }
  /*
  //test @1
  ofstream fout;
  fout.open("disp.txt");
  fout.close();
  */

  //模型坐标被改变
  //这个可以看作是处理位移对模型的影响的一种方式吗
  //这种做法会出现无法求解的问题 经调试问题不是这里导致的，而是变量没有初始化的问题

  //test 最后一步结束再更新坐标
  int step=time/dt;
 if(step>200)
  for (int i = 0; i < NDIM * num_nodes; ++i) {
    node_coord->getPointer()[i] += vec_data->getPointer()[i];
  }

  int num_local_nodes = patch.getNumberOfNodes();
  for (int i = 0; i < NDIM * num_local_nodes; ++i) {
    plot_data->getPointer()[i] = vec_data->getPointer()[i];
  }
}

/*******************************************************************************
* 计算应力
* ******************************************************************************/
void PatchStrategy::computeStress(hier::Patch<NDIM>& patch, const double time,
                                  const double dt,
                                  const string& component_name) {
  /// 取出本地PatchGeometry.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();
  /// 取出本地PatchTopology.
  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
      patch.getPatchTopology();
  /// 取出本地Patch的结点坐标数组.
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch_geo->getNodeCoordinates();
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_data =
      patch.getPatchData(d_solution_id);

  tbox::Pointer<pdat::CellData<NDIM, double> > str_data =
      patch.getPatchData(d_stress_id);

  //update #2
  tbox::Pointer<pdat::CellData<NDIM, double> > von_Mises =
      patch.getPatchData(d_von_mises_id);
  //取前一步的温度分布
  tbox::Pointer<pdat::NodeData<NDIM, double> > T_data =
      patch.getPatchData(th_Told_id);//

  /// 获取单元对应实体编号数组对象 update #3 at 2017-04-21 by tong @1
  tbox::Pointer<pdat::CellData<NDIM, int> > entityid_data =
      patch.getPatchData(d_EntityIdOfCell_id);
  tbox::Pointer<pdat::CellData<NDIM, int> > materialid_data =
      patch.getPatchData(material_num_id);

  int* dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  /// 获取单元周围结点的索引关系.
  tbox::Array<int> can_extent, can_indices;
  patch_top->getCellAdjacencyNodes(can_extent, can_indices);




  int num_cells = patch.getNumberOfCells();

  for (int i = 0; i < num_cells; ++i) {
    int n_vertex = can_extent[i + 1] - can_extent[i];
    int num_dof = NDIM * n_vertex;

    /**< 该单元的结点坐标及自由度映射 */
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<int> node_mapping(num_dof);
    tbox::Array<double> T_val(n_vertex);

    /// 下面的循环做两件事情：1. 建立自由度映射数组；2.取出结点坐标。
    for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
       T_val[i1]=T_data->getPointer()[can_indices[j]];
      for (int k = 0; k < NDIM; ++k) {
        node_mapping[NDIM * i1 + k] = dof_map[can_indices[j]] + k;
        vertex[i1][k] = (*node_coord)(k, can_indices[j]);
      }
    }

    /// 取出积分器对象.
    tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
        IntegratorManager<NDIM>::getManager();
    tbox::Pointer<BaseIntegrator<NDIM> > integrator =
        integrator_manager->getIntegrator("LinearTetrahedron");

    /// 取出形函数对象.
    tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
        ShapeFunctionManager<NDIM>::getManager();
    tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
        shape_manager->getShapeFunction("LinearTetrahedron");

    /// 取出材料.
    tbox::Pointer<MaterialManager<NDIM> > material_manager =
        MaterialManager<NDIM>::getManager();

    //update #3: @2 给不同实体中的单元赋不同的材料
    //实体1（或者实体集1）：Copper
    //实体2（或者实体集2）：Silicon
    //默认：Air
    //at 2017-04-21 by tong

    tbox::Pointer<Material> material =
    				material_manager->getMaterial("Gold");

    if((*materialid_data)(0,i)==1)
		 material =	material_manager->getMaterial("Silicon");
    else if((*materialid_data)(0,i)==13)
		 material =	material_manager->getMaterial("MoCu");
    else if((*materialid_data)(0,i)==2)
		 material =	material_manager->getMaterial("Copper");
    else if((*materialid_data)(0,i)==4)
		 material =	material_manager->getMaterial("SiO2");
    else if((*materialid_data)(0,i)==5)
		 material =	material_manager->getMaterial("SiN");     
    else if((*materialid_data)(0,i)==3)
		 material =	material_manager->getMaterial("Gold");
    else if((*materialid_data)(0,i)==15)
		 material =	material_manager->getMaterial("GaN");  
    else if((*materialid_data)(0,i)==16)
		 material =	material_manager->getMaterial("Al2O3"); 
    else if((*materialid_data)(0,i)==17)
		 material =	material_manager->getMaterial("Alloy");
    else if((*materialid_data)(0,i)==6)
		 material =	material_manager->getMaterial("Aluminum");  
    else
		 material =	material_manager->getMaterial("Gold");

   double temp_T=0;
    for(int t=0;t<4;t++)
    {temp_T+=T_val[t]/4;}
    //double Sigma=material->getSigma(temp_T);	

    /// 取出自由度数目.
    int n_dof = shape_func->getNumberOfDof();

    /// 取出积分点数目.
    int num_quad_pnts = integrator->getNumberOfQuadraturePoints();

    /// 取出积分点.
    tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
        integrator->getQuadraturePoints(vertex);

    /// 取出积分点的积分权重.
    tbox::Array<double> weight = integrator->getQuadratureWeights();

    /// 取出基函数在积分点的值和梯度值.
    tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
        shape_func->gradient(vertex, quad_pnt);

    /// 取出模量矩阵
    tbox::Array<tbox::Array<double> > moduli = material->getModuli(temp_T);
    /// \lambda
    double a = moduli[0][0];
    /// G
    double b = moduli[0][1];
    /// \lambda+2G
    double c = moduli[3][3];

    tbox::Array<double> stress(6);

    /// 计算单元应力.
    for (int i1 = 0; i1 < 6; ++i1) stress[i1] = 0.0;
    for (int l = 0; l < num_quad_pnts; ++l) {
      /// 该点的积分权重.
      double w = weight[l];

      /// 初始化位移在积分点的梯度值.
      double u_grad[NDIM][NDIM];
      for (int i1 = 0; i1 < NDIM; ++i1) {
        for (int j1 = 0; j1 < NDIM; ++j1) {
          u_grad[i1][j1] = 0.0;
        }
      }

      /// 计算位移在积分点的梯度值.
      for (int i1 = 0; i1 < n_dof; ++i1) {
        double* tmp_val = &(vec_data->getPointer()[node_mapping[NDIM * i1]]);
        for (int j1 = 0; j1 < NDIM; ++j1) {
          for (int k1 = 0; k1 < NDIM; ++k1) {
            u_grad[j1][k1] += tmp_val[j1] * bas_grad[l][i1][k1];
          }
        }
      }

      /// 计算单元应力.
      stress[0] += w * (a * u_grad[0][0] + b * (u_grad[1][1] + u_grad[2][2]));
      stress[1] += w * (a * u_grad[1][1] + b * (u_grad[0][0] + u_grad[2][2]));
      stress[2] += w * (a * u_grad[2][2] + b * (u_grad[1][1] + u_grad[0][0]));
      stress[3] += w * c * (u_grad[0][1] + u_grad[1][0]);
      stress[4] += w * c * (u_grad[1][2] + u_grad[2][1]);
      stress[5] += w * c * (u_grad[0][2] + u_grad[2][0]);
    }

    /// 将单元应力回填到应力数据片.
    for (int j = 0; j < 6; ++j) {
      (*str_data)(j, i) = stress[j];
    }

    //计算von Mises应力
    double temp_v=0;
    temp_v=pow((stress[0]-stress[1]),2)+pow((stress[1]-stress[2]),2)+
    		pow((stress[2]-stress[0]),2)+6*(stress[3]*stress[3]+stress[4]*stress[4]+stress[5]*stress[5]);

    (*von_Mises)(0,i)=sqrt(temp_v/2);
  }
}

/// 恢复指定体上的应力
void PatchStrategy::StressRecovery(hier::Patch<NDIM>& patch, const double time,
                                   const double dt, const string& component_name){
    /// 取出本地PatchGeometry.
    tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
            patch.getPatchGeometry();
    /// 取出本地PatchTopology.
    tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
            patch.getPatchTopology();
    /// 应力信息
    tbox::Pointer<pdat::CellData<NDIM, double> > str_data =
            patch.getPatchData(d_stress_id);
    tbox::Pointer<pdat::CellData<NDIM, int> > is_on_bnd =
            patch.getPatchData(d_contained_domain_id);
    tbox::Pointer<pdat::CellData<NDIM, int> > materialid_data =
            patch.getPatchData(material_num_id);


    tbox::Pointer<pdat::NodeData<NDIM, double> > T_data =
            patch.getPatchData(th_Told_id);//
    tbox::Pointer<pdat::NodeData<NDIM, double> > Tolder_data =
            patch.getPatchData(th_Tolder_id);//
    /// 单元数目
    int num_cells = patch.getNumberOfCells(0);
    DECLARE_ADJACENCY(patch,cell,cell,Cell,Cell);
    DECLARE_ADJACENCY(patch,cell,face,Cell,Face);
    DECLARE_ADJACENCY(patch,cell,node,Cell,Node);
    DECLARE_ADJACENCY(patch,face,node,Face,Node);

    /// material 设置
    tbox::Pointer<MaterialManager<NDIM> > material_manager =
            MaterialManager<NDIM>::getManager();
    //高斯积分
    appu::TetQuad quad(patch, patch.getPatchData(d_Cell_volume_id),
                       patch.getPatchData(d_Cell_jacobian_id));
    const appu::TetQuad::Quad *facequad = quad.GetFaceQuadTable()[2];
    const appu::TetQuad::Quad *cellquad = quad.GetQuadTable()[2];
    tbox::Pointer<pdat::VectorData<NDIM, double> > vec_data =
            patch.getPatchData(d_solution_id);

    /// 取出本地Patch的结点坐标数组.
    tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
            patch_geo->getNodeCoordinates();
    tbox::Pointer<BaseElement<NDIM> > ele =
            d_element_manager->getElement(d_element_type[0]);
    tbox::Pointer<pdat::CellData<NDIM, double> > improved_coef =
            patch.getPatchData(d_improved_coefficient_id);
    int patchnumber = patch.getIndex();
    double nnn = patch.getNumberOfCells(1);
    for(int cell = 0; cell < num_cells; cell++){
        int entity_id = (*materialid_data)(0,cell);

        //        tbox::pout<<"进入bnd判断前"<<cell<<"\t"<<patchnumber<<endl;
        tbox::Pointer<Material> material =
                material_manager->getMaterial("Gold");
        if(entity_id==1)
            material =	material_manager->getMaterial("Silicon");
        else if(entity_id==13)
            material =	material_manager->getMaterial("MoCu");
        else if(entity_id==2)
            material =	material_manager->getMaterial("Copper");
        else if(entity_id==4)
            material =	material_manager->getMaterial("SiO2");
        else if(entity_id==5)
            material =	material_manager->getMaterial("SiN");
        else if(entity_id==3)
            material =	material_manager->getMaterial("Gold");
        else if(entity_id==15)
            material =	material_manager->getMaterial("GaN");
        else if(entity_id==16)
            material =	material_manager->getMaterial("Al2O3");
        else if(entity_id==17)
            material =	material_manager->getMaterial("Alloy");
        else if(entity_id==6)
            material =	material_manager->getMaterial("Aluminum");
        else
            material =	material_manager->getMaterial("Gold");


        ///  如果此单元是必须要计算的，且邻接单元都齐
        if((*is_on_bnd)(0,cell) == 1){
            /// 计算超定方程
            /// Yin-Da Wang: 2025-01-15
            /// 被拒稿了 很难过
            /// 共24个未知数，120个方程
            /// 60个方程是面牵引方程
            /// 60个方程是体守恒方程
            /// Syms_coef_overdetermined
            /// nx*(d1 + a1*x + b1*y + c1*z) + ny*(d4 + a4*x + b4*y + c4*z) + nz*(d6 + a6*x + b6*y + c6*z)
            /// nx*(d4 + a4*x + b4*y + c4*z) + ny*(d2 + a2*x + b2*y + c2*z) + nz*(d5 + a5*x + b5*y + c5*z)
            /// nx*(d6 + a6*x + b6*y + c6*z) + ny*(d5 + a5*x + b5*y + c5*z) + nz*(d3 + a3*x + b3*y + c3*z)
            double ToolMatrix[120][24];
            for(int rr = 0; rr < 120; rr ++)
                for(int cc = 0; cc < 24; cc ++)
                    ToolMatrix[rr][cc] = 0.;
            Eigen::VectorXd ToolRHSVec(120);
            int elem_list[5]={cell,
                              cell_cell_idx[cell_cell_ext[cell]+0],
                              cell_cell_idx[cell_cell_ext[cell]+1],
                              cell_cell_idx[cell_cell_ext[cell]+2],
                              cell_cell_idx[cell_cell_ext[cell]+3]};
            for(int elem = 0; elem < 5; elem ++){
                int global_cell = elem_list[elem];
                int n_dof = NDIM*(NDIM+1);
                tbox::Pointer<tbox::Matrix<double> > ele_mat = new tbox::Matrix<double>();
                ele_mat->resize(n_dof, n_dof);
                tbox::Vector<double> ele_disp(n_dof);
                /// 右端项机械应力项
                tbox::Vector<double> ele_sol(n_dof);
                /// 右端项热膨胀项
               tbox::Pointer<tbox::Vector<double> > ele_expansion = new tbox::Vector<double>();
               ele_expansion->resize(n_dof);
               for (int i1 = 0; i1 < n_dof; ++i1) {
                   (*ele_expansion)[i1] = 0.;
                   for (int j = 0; j < n_dof; ++j) {
                       (*ele_mat)(i1, j) = 0.0;
                   }
               }
                tbox::Array<hier::DoubleVector<NDIM> > vertex(NDIM+1);
                int n_vertex = 4;
                //存储上一时刻单元节点温度
                tbox::Array<double> T_val(n_vertex); //当前温度
                tbox::Array<double> Tolder_val(n_vertex);
                /// 填充结点上的坐标和位移信息
                for(int nn = cell_node_ext[global_cell]; nn <cell_node_ext[global_cell+1]; nn++){
                    int i1 = nn - cell_node_ext[global_cell];
                    int index = cell_node_idx[nn];
                    T_val[i1]=(*T_data)(0,index);
//                    T_val[i1] = 650.;
                    T_val[i1] = 300.;
                    Tolder_val[i1]=(*Tolder_data)(0,index);
                    vertex[i1][0] = (*node_coord)(0, index);
                    vertex[i1][1] = (*node_coord)(1, index);
                    vertex[i1][2] = (*node_coord)(2, index);
                    /// 右端项的填充
                    ele_disp[i1*NDIM+0] = (vec_data->getPointer()[index*NDIM+0]);
                    ele_disp[i1*NDIM+1] = (vec_data->getPointer()[index*NDIM+1]);
                    ele_disp[i1*NDIM+2] = (vec_data->getPointer()[index*NDIM+2]);
                }

                double elem_T = 293.15;
                /// 局部应力刚度矩阵
                ele->buildStiffElementMatrix(vertex, dt, time,
                                             ele_mat,(*materialid_data)(0,global_cell),elem_T);
                NewmarkData d_newmark[n_vertex];
                /// 局部热应力矩阵
                ele->buildElementRHS(vertex, dt, time, ele_expansion, d_newmark, (*materialid_data)(0,global_cell),T_val,Tolder_val);
                ele_sol = (*ele_mat) * ele_disp;
                /// 前60行
                for(int row = 0; row < n_dof; row++){
                    ToolRHSVec[elem*n_dof+row] = ele_sol[row] - (*ele_expansion)[row];
                }             
                for(int nbf = 0; nbf < 4; nbf ++){
                    int global_face
                            = cell_face_idx[cell_face_ext[global_cell]+nbf];
                    double outer_normal[NDIM], outer_normal_raw[NDIM];
                    // 计算外推面单元面积
                    outerNormal(cell_node_ext[global_cell+1]-cell_node_ext[global_cell],
                            cell_node_idx.getPointer() + cell_node_ext[global_cell],
                            face_node_ext[global_face+1]-face_node_ext[global_face],
                            face_node_idx.getPointer() + face_node_ext[global_face],
                            patch.getNumberOfNodes(1),
                            patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                            NDIM,
                            outer_normal,
                            outer_normal_raw);
                    /// 这个面的外法向正好和函数得到的相反
                    outer_normal[0] = -outer_normal[0];
                    outer_normal[1] = -outer_normal[1];
                    outer_normal[2] = -outer_normal[2];
                    /// 边界面元的面积
                    double this_area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw))/ 2.0;
                    double nx = outer_normal[0];
                    double ny = outer_normal[1];
                    double nz = outer_normal[2];
                    for(int quad_n = 0; quad_n < facequad->npoints; quad_n ++){
                        double bary3d[NDIM+1]={0.,0.,0.,0.};
                        quad.bary2dTo3d
                                (facequad->points + (quad_n*(NDIM+1-1)),global_cell,nbf,bary3d);
                        double bmr[NDIM]={0.,0.,0.};
                        /// 计算积分点的位矢
                        for(int nn = 0; nn < NDIM+1; nn ++){
                            int this_node = cell_node_idx[cell_node_ext[global_cell]+nn];
                            double coordinate_here[NDIM] = {(*node_coord)(0,this_node),
                                                            (*node_coord)(1,this_node),
                                                            (*node_coord)(2,this_node)};
                            for(int dd = 0; dd < NDIM; dd++)
                                bmr[dd] += bary3d[nn]*coordinate_here[dd];
                        }
                        /// ToolMatrix矩阵顺序
                        /// 第一层cell 1; cell 2; cell 3; cell 4; cell 5 (12元素)
                        /// 第二层bas 1; bas 2; bas3; bas 4 (3元素)
                        /// 第三层 x-axis; y-axis; z-axis

                        double q_weight = facequad->weights[quad_n];
                        /// x方向填充
                        /// 146 对于每个元素，每个面每个积分点累加一次
                        for(int nbas = 0; nbas < NDIM+1; nbas ++){
                            ToolMatrix[elem*12+nbas*3+0][0]
                                    += bary3d[nbas]*nx*bmr[0]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+0][1]
                                    += bary3d[nbas]*nx*bmr[1]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+0][2]
                                    += bary3d[nbas]*nx*bmr[2]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+0][3]
                                    += bary3d[nbas]*nx*1.*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+0][12]
                                    += bary3d[nbas]*ny*bmr[0]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+0][13]
                                    += bary3d[nbas]*ny*bmr[1]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+0][14]
                                    += bary3d[nbas]*ny*bmr[2]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+0][15]
                                    += bary3d[nbas]*ny*1.*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+0][20]
                                    += bary3d[nbas]*nz*bmr[0]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+0][21]
                                    += bary3d[nbas]*nz*bmr[1]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+0][22]
                                    += bary3d[nbas]*nz*bmr[2]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+0][23]
                                    += bary3d[nbas]*nz*1.*this_area*q_weight;
                        }
                        /// y方向填充
                        /// 425 对于每个元素，每个面每个积分点累加一次
                        for(int nbas = 0; nbas < NDIM+1; nbas ++){
                            ToolMatrix[elem*12+nbas*3+1][12]
                                    += bary3d[nbas]*nx*bmr[0]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+1][13]
                                    += bary3d[nbas]*nx*bmr[1]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+1][14]
                                    += bary3d[nbas]*nx*bmr[2]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+1][15]
                                    += bary3d[nbas]*nx*1.*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+1][4]
                                    += bary3d[nbas]*ny*bmr[0]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+1][5]
                                    += bary3d[nbas]*ny*bmr[1]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+1][6]
                                    += bary3d[nbas]*ny*bmr[2]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+1][7]
                                    += bary3d[nbas]*ny*1.*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+1][16]
                                    += bary3d[nbas]*nz*bmr[0]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+1][17]
                                    += bary3d[nbas]*nz*bmr[1]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+1][18]
                                    += bary3d[nbas]*nz*bmr[2]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+1][19]
                                    += bary3d[nbas]*nz*1.*this_area*q_weight;
                        }
                        /// z方向填充
                        /// 653 对于每个元素，每个面每个积分点累加一次
                        for(int nbas = 0; nbas < NDIM+1; nbas ++){
                            ToolMatrix[elem*12+nbas*3+2][20]
                                    += bary3d[nbas]*nx*bmr[0]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+2][21]
                                    += bary3d[nbas]*nx*bmr[1]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+2][22]
                                    += bary3d[nbas]*nx*bmr[2]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+2][23]
                                    += bary3d[nbas]*nx*1.*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+2][16]
                                    += bary3d[nbas]*ny*bmr[0]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+2][17]
                                    += bary3d[nbas]*ny*bmr[1]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+2][18]
                                    += bary3d[nbas]*ny*bmr[2]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+2][19]
                                    += bary3d[nbas]*ny*1.*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+2][8]
                                    += bary3d[nbas]*nz*bmr[0]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+2][9]
                                    += bary3d[nbas]*nz*bmr[1]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+2][10]
                                    += bary3d[nbas]*nz*bmr[2]*this_area*q_weight;
                            ToolMatrix[elem*12+nbas*3+2][11]
                                    += bary3d[nbas]*nz*1.*this_area*q_weight;
                        }

                    }
                }

                /// 后60行
                for(int row = 0; row < n_dof; row++){
                    ToolRHSVec[60+elem*n_dof+row] = ele_sol[row];
                }
                tbox::Pointer<tbox::Matrix<double> > vol_mat = new tbox::Matrix<double>();
                vol_mat->resize(n_dof, 24);
                for (int i1 = 0; i1 < n_dof; ++i1) {
                    for (int j = 0; j < 24; ++j) {
                        (*vol_mat)(i1, j) = 0.0;
                    }
                }
                ele->buildRecoveryMatrix(vertex, dt, time,
                                             vol_mat,(*materialid_data)(0,global_cell));
                for (int i1 = 0; i1 < n_dof; ++i1) {
                    for (int j = 0; j < 24; ++j) {
                        /// 后60行，第elem个单元，第i1个自由度
                       ToolMatrix[60+elem*n_dof+i1][j]  = (*vol_mat)(i1, j);
                    }
                }

            }
            Eigen::MatrixXd ToolLHSMat(120, 24);
            /// 前60行
            for (int i = 0; i < 120; ++i) {
                for (int j = 0; j < 24; ++j) {
                    ToolLHSMat(i, j) = ToolMatrix[i][j]; // 填充 Eigen 矩阵
                }
            }
            Eigen::VectorXd ToolSol(24);
            ToolSol = ToolLHSMat.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(ToolRHSVec);

            for(int row = 0; row < 24; row ++){
                (*improved_coef)(row,cell) = ToolSol[row];
            }
        }

        /// 如果这个单元需要计算，但不是中心域，先看它是不是在某个工作域中心的连接单元上
        if((*is_on_bnd)(0,cell) == 0){
            /// 先看它是不是某个工作域中心的连接单元上
            bool should_computed = true;
            for(int local_cell = cell_cell_ext[cell];local_cell < cell_cell_ext[cell+1];local_cell ++){
                int global_cell = cell_cell_idx[local_cell];
                if((*is_on_bnd)(0,global_cell) == 1)
                    should_computed = false;
            }
            /// 如果这个单元不能由任何一个tool region包含
            if(should_computed){
                /// 此时需要找到五个连通区域
                std::vector<int> elem_list;
                elem_list.push_back(cell);
                for(int source_ele = 0; source_ele < elem_list.size(); source_ele ++){
                    int elem_num = elem_list[source_ele];
                    for(int local_adj = cell_cell_ext[elem_num];local_adj < cell_cell_ext[elem_num+1];local_adj ++){
                        /// 候选单元的全局编号
                        int global_adj = cell_cell_idx[local_adj];
                        /// 在elem_list里找这个单元
                        std::vector<int>::iterator
                                it = std::find(elem_list.begin(),elem_list.end(),global_adj);
                        if (it == elem_list.end()) {
                            /// 如果elem_list里没找到这个单元
                            elem_list.push_back(global_adj);
                        }
                        if(elem_list.size() > 4)
                            break;
                    }
                    if(elem_list.size() > 4)
                        break;
                }

                double ToolMatrix[120][24];
                for(int rr = 0; rr < 120; rr ++)
                    for(int cc = 0; cc < 24; cc ++)
                        ToolMatrix[rr][cc] = 0.;
                Eigen::VectorXd ToolRHSVec(120);
                for(int elem = 0; elem < 5; elem ++){
                    int global_cell = elem_list[elem];
                    int n_dof = NDIM*(NDIM+1);
                    tbox::Pointer<tbox::Matrix<double> > ele_mat = new tbox::Matrix<double>();
                    ele_mat->resize(n_dof, n_dof);
                    tbox::Vector<double> ele_disp(n_dof);
                    /// 右端项机械应力项
                    tbox::Vector<double> ele_sol(n_dof);
                    /// 右端项热膨胀项
                    tbox::Pointer<tbox::Vector<double> > ele_expansion = new tbox::Vector<double>();
                    ele_expansion->resize(n_dof);
                    for (int i1 = 0; i1 < n_dof; ++i1) {
                        (*ele_expansion)[i1] = 0.;
                      for (int j = 0; j < n_dof; ++j) {
                        (*ele_mat)(i1, j) = 0.0;
                      }
                    }
                    tbox::Array<hier::DoubleVector<NDIM> > vertex(NDIM+1);
                    int n_vertex = 4;
                    //存储上一时刻单元节点温度
                    tbox::Array<double> T_val(n_vertex); //当前温度
                    tbox::Array<double> Tolder_val(n_vertex);
                    for(int nn = cell_node_ext[global_cell]; nn <cell_node_ext[global_cell+1]; nn++){
                        int i1 = nn - cell_node_ext[global_cell];
                        int index = cell_node_idx[nn];
                        T_val[i1]=(*T_data)(0,index);
                        T_val[i1] = 650.;
                        T_val[i1] = 300.;
                        Tolder_val[i1]=(*Tolder_data)(0,index);
                        vertex[i1][0] = (*node_coord)(0, index);
                        vertex[i1][1] = (*node_coord)(1, index);
                        vertex[i1][2] = (*node_coord)(2, index);
                        /// 右端项的填充
                        ele_disp[i1*NDIM+0] = (vec_data->getPointer()[index*NDIM+0]);
                        ele_disp[i1*NDIM+1] = (vec_data->getPointer()[index*NDIM+1]);
                        ele_disp[i1*NDIM+2] = (vec_data->getPointer()[index*NDIM+2]);
                    }

                    double elem_T = 293.15;
                    /// 局部应力刚度矩阵
                    ele->buildStiffElementMatrix(vertex, dt, time,
                                                 ele_mat,(*materialid_data)(0,global_cell),elem_T);
                    NewmarkData d_newmark[n_vertex];
                    /// 局部热应力矩阵
                    ele->buildElementRHS(vertex, dt, time, ele_expansion, d_newmark, (*materialid_data)(0,global_cell),T_val,Tolder_val);
                    ele_sol = (*ele_mat) * ele_disp;
                    for(int row = 0; row < n_dof; row++){
                        ToolRHSVec[elem*n_dof+row] = ele_sol[row] - (*ele_expansion)[row];
                    }
                    for(int nbf = 0; nbf < 4; nbf ++){
                        int global_face
                                = cell_face_idx[cell_face_ext[global_cell]+nbf];
                        double outer_normal[NDIM], outer_normal_raw[NDIM];
                        // 计算外推面单元面积
                        outerNormal(cell_node_ext[global_cell+1]-cell_node_ext[global_cell],
                                cell_node_idx.getPointer() + cell_node_ext[global_cell],
                                face_node_ext[global_face+1]-face_node_ext[global_face],
                                face_node_idx.getPointer() + face_node_ext[global_face],
                                patch.getNumberOfNodes(1),
                                patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                                NDIM,
                                outer_normal,
                                outer_normal_raw);
                        /// 这个面的外法向正好和函数得到的相反
                        outer_normal[0] = -outer_normal[0];
                        outer_normal[1] = -outer_normal[1];
                        outer_normal[2] = -outer_normal[2];
                        /// 边界面元的面积
                        double this_area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw))/ 2.0;
                        double nx = outer_normal[0];
                        double ny = outer_normal[1];
                        double nz = outer_normal[2];
                        for(int quad_n = 0; quad_n < facequad->npoints; quad_n ++){
                            double bary3d[NDIM+1]={0.,0.,0.,0.};
                            quad.bary2dTo3d
                                    (facequad->points + (quad_n*(NDIM+1-1)),global_cell,nbf,bary3d);
                            double bmr[NDIM]={0.,0.,0.};
                            /// 计算积分点的位矢
                            for(int nn = 0; nn < NDIM+1; nn ++){
                                int this_node = cell_node_idx[cell_node_ext[global_cell]+nn];
                                double coordinate_here[NDIM] = {(*node_coord)(0,this_node),
                                                                (*node_coord)(1,this_node),
                                                                (*node_coord)(2,this_node)};
                                for(int dd = 0; dd < NDIM; dd++)
                                    bmr[dd] += bary3d[nn]*coordinate_here[dd];
                            }
                            /// ToolMatrix矩阵顺序
                            /// 第一层cell 1; cell 2; cell 3; cell 4; cell 5 (12元素)
                            /// 第二层bas 1; bas 2; bas3; bas 4 (3元素)
                            /// 第三层 x-axis; y-axis; z-axis

                            double q_weight = facequad->weights[quad_n];
                            /// x方向填充
                            /// 146 对于每个元素，每个面每个积分点累加一次
                            for(int nbas = 0; nbas < NDIM+1; nbas ++){
                                ToolMatrix[elem*12+nbas*3+0][0]
                                        += bary3d[nbas]*nx*bmr[0]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+0][1]
                                        += bary3d[nbas]*nx*bmr[1]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+0][2]
                                        += bary3d[nbas]*nx*bmr[2]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+0][3]
                                        += bary3d[nbas]*nx*1.*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+0][12]
                                        += bary3d[nbas]*ny*bmr[0]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+0][13]
                                        += bary3d[nbas]*ny*bmr[1]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+0][14]
                                        += bary3d[nbas]*ny*bmr[2]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+0][15]
                                        += bary3d[nbas]*ny*1.*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+0][20]
                                        += bary3d[nbas]*nz*bmr[0]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+0][21]
                                        += bary3d[nbas]*nz*bmr[1]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+0][22]
                                        += bary3d[nbas]*nz*bmr[2]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+0][23]
                                        += bary3d[nbas]*nz*1.*this_area*q_weight;
                            }
                            /// y方向填充
                            /// 425 对于每个元素，每个面每个积分点累加一次
                            for(int nbas = 0; nbas < NDIM+1; nbas ++){
                                ToolMatrix[elem*12+nbas*3+1][12]
                                        += bary3d[nbas]*nx*bmr[0]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+1][13]
                                        += bary3d[nbas]*nx*bmr[1]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+1][14]
                                        += bary3d[nbas]*nx*bmr[2]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+1][15]
                                        += bary3d[nbas]*nx*1.*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+1][4]
                                        += bary3d[nbas]*ny*bmr[0]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+1][5]
                                        += bary3d[nbas]*ny*bmr[1]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+1][6]
                                        += bary3d[nbas]*ny*bmr[2]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+1][7]
                                        += bary3d[nbas]*ny*1.*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+1][16]
                                        += bary3d[nbas]*nz*bmr[0]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+1][17]
                                        += bary3d[nbas]*nz*bmr[1]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+1][18]
                                        += bary3d[nbas]*nz*bmr[2]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+1][19]
                                        += bary3d[nbas]*nz*1.*this_area*q_weight;
                            }
                            /// z方向填充
                            /// 653 对于每个元素，每个面每个积分点累加一次
                            for(int nbas = 0; nbas < NDIM+1; nbas ++){
                                ToolMatrix[elem*12+nbas*3+2][20]
                                        += bary3d[nbas]*nx*bmr[0]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+2][21]
                                        += bary3d[nbas]*nx*bmr[1]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+2][22]
                                        += bary3d[nbas]*nx*bmr[2]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+2][23]
                                        += bary3d[nbas]*nx*1.*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+2][16]
                                        += bary3d[nbas]*ny*bmr[0]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+2][17]
                                        += bary3d[nbas]*ny*bmr[1]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+2][18]
                                        += bary3d[nbas]*ny*bmr[2]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+2][19]
                                        += bary3d[nbas]*ny*1.*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+2][8]
                                        += bary3d[nbas]*nz*bmr[0]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+2][9]
                                        += bary3d[nbas]*nz*bmr[1]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+2][10]
                                        += bary3d[nbas]*nz*bmr[2]*this_area*q_weight;
                                ToolMatrix[elem*12+nbas*3+2][11]
                                        += bary3d[nbas]*nz*1.*this_area*q_weight;
                            }

                        }
                    }
                    for(int row = 0; row < n_dof; row++){
                        ToolRHSVec[60+elem*n_dof+row] = ele_sol[row];
                    }
                    tbox::Pointer<tbox::Matrix<double> > vol_mat = new tbox::Matrix<double>();
                    vol_mat->resize(n_dof, 24);
                    for (int i1 = 0; i1 < n_dof; ++i1) {
                        for (int j = 0; j < 24; ++j) {
                            (*vol_mat)(i1, j) = 0.0;
                        }
                    }
                    ele->buildRecoveryMatrix(vertex, dt, time,
                                             vol_mat,(*materialid_data)(0,global_cell));
                    for (int i1 = 0; i1 < n_dof; ++i1) {
                        for (int j = 0; j < 24; ++j) {
                            /// 后60行，第elem个单元，第i1个自由度
                            ToolMatrix[60+elem*n_dof+i1][j]  = (*vol_mat)(i1, j);
                        }
                    }
                }

                Eigen::MatrixXd ToolLHSMat(120, 24);
                for (int i = 0; i < 120; ++i) {
                    for (int j = 0; j < 24; ++j) {
                        ToolLHSMat(i, j) = ToolMatrix[i][j]; // 填充 Eigen 矩阵
                    }
                }
                Eigen::VectorXd ToolSol(24);
                ToolSol = ToolLHSMat.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(ToolRHSVec);
                for(int row = 0; row < 24; row ++){
                    (*improved_coef)(row,cell) = ToolSol[row];
                }
            }

        }
    }




}
/// 后处理恢复指定体上的应力
void PatchStrategy::PostprocessStress(hier::Patch<NDIM>& patch, const double time,
                                   const double dt, const string& component_name){
    int num_cells = patch.getNumberOfCells(0);
    DECLARE_ADJACENCY(patch,cell,cell,Cell,Cell);
    DECLARE_ADJACENCY(patch,cell,node,Cell,Node);
    /// 最终的输出应力系数
    tbox::Pointer<pdat::CellData<NDIM, double> > final_improved_coef =
        patch.getPatchData(d_output_coefficient_id);
    tbox::Pointer<pdat::CellData<NDIM, double> > improved_coef =
        patch.getPatchData(d_improved_coefficient_id);
    tbox::Pointer<pdat::CellData<NDIM, double> > rec_von =
        patch.getPatchData(d_improved_von_mises_id);
    tbox::Pointer<pdat::CellData<NDIM, int> > is_on_bnd =
            patch.getPatchData(d_contained_domain_id);
    /// 取出本地PatchGeometry.
    tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
            patch.getPatchGeometry();
    tbox::Pointer<pdat::CellData<NDIM, double> > cell_coord =
            patch_geo->getCellCoordinates();
    tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
            patch_geo->getNodeCoordinates();
    tbox::Pointer<pdat::CellData<NDIM, double> > von_Mises =
        patch.getPatchData(d_von_mises_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> > displacement =
        patch.getPatchData(d_displacement_id);
    tbox::Pointer<pdat::CellData<NDIM, double> > volume =
        patch.getPatchData(d_Cell_volume_id);
    for(int cc = 0; cc < num_cells; cc ++){
        bool is_computed_cell = false;
        /// 遍历需要恢复应力的位置，看单元是否在这些位置
        for(int von_id = 0; von_id < d_improved_stress.size();von_id++){
            if (HAS_ENTITY_SET(patch, d_improved_stress[von_id], CELL, 1)){
                DECLARE_ENTITY_SET(patch, improved_stress_space,
                                   d_improved_stress[von_id], CELL, 1);
                std::sort(improved_stress_space.getPointer(),
                          improved_stress_space.getPointer() + improved_stress_space.size());
                if(std::binary_search(improved_stress_space.getPointer(),
                                      improved_stress_space.getPointer()+improved_stress_space.size(),
                                      cc))
                    is_computed_cell = true;
            }
        }
        /// 如果这个单元需要计算的话
        if(is_computed_cell){
            /// 如果这个单元需要计算，但不是中心域，先看它是不是在某个工作域中心的连接单元上
            /// 否则就好办了，直接5单元平均
            if((*is_on_bnd)(0,cc) == 0){
                bool should_computed = true;
                for(int local_cell = cell_cell_ext[cc];local_cell < cell_cell_ext[cc+1];local_cell ++){
                    int global_cell = cell_cell_idx[local_cell];
                    if((*is_on_bnd)(0,global_cell) == 1)
                        should_computed = false;
                }
                /// 如果这个单元不能由任何一个tool region包含
                /// 则在恢复过程中直接使用那个权衡方案得到的系数
                /// 否则一般好办，看有几个单元能帮助他
                if(should_computed){
                    for(int row = 0; row < 24; row ++){
                        (*final_improved_coef)(row,cc) = (*improved_coef)(row,cc);
                    }
                }
                else{
                    /// 注意到此时ast_cell最终值不可能为0
                    double ast_cell = 0.;
                    double coef_list[24];
                    bzero(coef_list,sizeof(coef_list));
                    for(int cac_ext = cell_cell_ext[cc]; cac_ext < cell_cell_ext[cc+1];cac_ext ++){
                        int cac_idx = cell_cell_idx[cac_ext];
                        if((*is_on_bnd)(0,cac_idx) == 1){
                            ast_cell += 1.;
                            for(int row = 0; row < 24; row ++)
                                coef_list[row] += (*improved_coef)(row,cac_idx);

                        }
                    }


                    /// 做个平均
                    for(int row = 0; row < 24; row ++)
                        (*final_improved_coef)(row,cc) = coef_list[row]/ast_cell;

                }
            }
            else{
                /// 注意到此时ast_cell最终值不可能为0
                /// 首先保底它自己的是，所以这里ast_cell要设为1开始
                double ast_cell = 1.;
                double coef_list[24];
                bzero(coef_list,sizeof(coef_list));
                for(int row = 0; row < 24; row ++)
                    coef_list[row] += (*improved_coef)(row,cc);
                for(int cac_ext = cell_cell_ext[cc]; cac_ext < cell_cell_ext[cc+1];cac_ext ++){
                    int cac_idx = cell_cell_idx[cac_ext];
                    if((*is_on_bnd)(0,cac_idx) == 1){
                        ast_cell += 1.;
                        for(int row = 0; row < 24; row ++)
                            coef_list[row] += (*improved_coef)(row,cac_idx);
                    }
                }


                /// 做个平均
                for(int row = 0; row < 24; row ++)
                    (*final_improved_coef)(row,cc) = coef_list[row]/ast_cell;
            }

        }
        if(is_computed_cell){
            double coord_Q[NDIM] = {(*cell_coord)(0,cc),(*cell_coord)(1,cc),(*cell_coord)(2,cc)};
            double rec_stress_0 = coord_Q[0]*(*final_improved_coef)(0,cc)
                    +coord_Q[1]*(*final_improved_coef)(1,cc)
                    +coord_Q[2]*(*final_improved_coef)(2,cc)
                    +1.*(*final_improved_coef)(3,cc);
            double rec_stress_1 = coord_Q[0]*(*final_improved_coef)(4,cc)
                    +coord_Q[1]*(*final_improved_coef)(5,cc)
                    +coord_Q[2]*(*final_improved_coef)(6,cc)
                    +1.*(*final_improved_coef)(7,cc);
            double rec_stress_2 = coord_Q[0]*(*final_improved_coef)(8,cc)
                    +coord_Q[1]*(*final_improved_coef)(9,cc)
                    +coord_Q[2]*(*final_improved_coef)(10,cc)
                    +1.*(*final_improved_coef)(11,cc);
            double rec_stress_3 = coord_Q[0]*(*final_improved_coef)(12,cc)
                    +coord_Q[1]*(*final_improved_coef)(13,cc)
                    +coord_Q[2]*(*final_improved_coef)(14,cc)
                    +1.*(*final_improved_coef)(15,cc);
            double rec_stress_4 = coord_Q[0]*(*final_improved_coef)(16,cc)
                    +coord_Q[1]*(*final_improved_coef)(17,cc)
                    +coord_Q[2]*(*final_improved_coef)(18,cc)
                    +1.*(*final_improved_coef)(19,cc);
            double rec_stress_5 = coord_Q[0]*(*final_improved_coef)(20,cc)
                    +coord_Q[1]*(*final_improved_coef)(21,cc)
                    +coord_Q[2]*(*final_improved_coef)(22,cc)
                    +1.*(*final_improved_coef)(23,cc);
            double temp_v =
                    pow((rec_stress_0-rec_stress_1),2)+pow((rec_stress_1-rec_stress_2),2)
                    +pow((rec_stress_2-rec_stress_0),2)
                    +6*(pow(rec_stress_3,2)+pow(rec_stress_4,2)+pow(rec_stress_5,2));
            double rec_mises = sqrt(temp_v/2);
            (*rec_von)(0,cc) = rec_mises;
        }
        else{
            (*rec_von)(0,cc) = 0.;
        }
        double dis_0 = (*displacement)(0,cell_node_idx[cell_node_ext[cc]+0]);
        double dis_1 = (*displacement)(0,cell_node_idx[cell_node_ext[cc]+1]);
        double dis_2 = (*displacement)(0,cell_node_idx[cell_node_ext[cc]+2]);
        double dis_3 = (*displacement)(0,cell_node_idx[cell_node_ext[cc]+3]);
        double dis = (dis_0+dis_1+dis_2+dis_3)/4;



        double vol;
          double oned6 = 1.0 / (double)6.0;
          double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;

          int *node = &(cell_node_idx[cell_node_ext[cc]]);

          x0 = (*node_coord)(0, node[0]);
          y0 = (*node_coord)(1, node[0]);
          z0 = (*node_coord)(2, node[0]);
          x1 = (*node_coord)(0, node[1]) - x0;
          y1 = (*node_coord)(1, node[1]) - y0;
          z1 = (*node_coord)(2, node[1]) - z0;
          x2 = (*node_coord)(0, node[2]) - x0;
          y2 = (*node_coord)(1, node[2]) - y0;
          z2 = (*node_coord)(2, node[2]) - z0;
          x3 = (*node_coord)(0, node[3]) - x0;
          y3 = (*node_coord)(1, node[3]) - y0;
          z3 = (*node_coord)(2, node[3]) - z0;

          vol = fabs(x1 * y2 * z3 + x2 * y3 * z1 + y1 * z2 * x3 -
                     (z1 * y2 * x3 + y1 * x2 * z3 + z2 * y3 * x1)) *
                oned6;

        stringstream OutputFile;
        int patch_number = patch.getIndex();
        OutputFile<<"TotalMisesOnPatch"<<patch_number;
        string OutFileName = OutputFile.str();
        char ch[40];strcpy(ch,OutFileName.c_str());
        ofstream vonmises_output;
        vonmises_output.open(ch,ios::app);
        vonmises_output<< std::fixed << std::setprecision(12);
        vonmises_output<<(*cell_coord)(0,cc)<<"\t"<<(*cell_coord)(1,cc)<<"\t"<<(*cell_coord)(2,cc)<<"\t"
                      <<vol<<"\t"
                     <<dis<<"\t"
                    <<(*rec_von)(0,cc)<<"\t"
                   <<(*von_Mises)(0,cc)<<endl;
        vonmises_output.close();



        stringstream OutputFile2;
        //            int patch_number = patch.getIndex();
        OutputFile2<<"CoordQuery"<<patch_number;
        string OutFileName2 = OutputFile2.str();
        char ch2[40];strcpy(ch2,OutFileName2.c_str());
        ofstream vonmises_output2;
        vonmises_output2.open(ch2,ios::app);
        vonmises_output2<< std::fixed << std::setprecision(12);
        vonmises_output2<<(*cell_coord)(0,cc)<<"\t"<<(*cell_coord)(1,cc)<<"\t"<<(*cell_coord)(2,cc)<<endl;
        vonmises_output2.close();


    }
}
void PatchStrategy::Dataexplorer(hier::Patch<NDIM> &patch, const double time,
                                 const double dt,
                                 const string &component_name){
    /// 取出本地PatchGeometry.
    tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
            patch.getPatchGeometry();
    tbox::Pointer<pdat::CellData<NDIM, double> > cell_coord =
            patch_geo->getCellCoordinates();
    tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
            patch_geo->getNodeCoordinates();
    /// 最终的输出应力系数
    tbox::Pointer<pdat::CellData<NDIM, double> > final_improved_coef =
            patch.getPatchData(d_output_coefficient_id);
    /// 算出来的Von Mises应力
    tbox::Pointer<pdat::CellData<NDIM, double> > von_Mises =
        patch.getPatchData(d_von_mises_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> > node_temp =
            patch.getPatchData(th_plot_id);

    int num_node = patch.getNumberOfNodes(0);
    int num_face = patch.getNumberOfFaces(0);
    int num_cell = patch.getNumberOfCells(0);

    DECLARE_ADJACENCY(patch, face, node, Face, Node);
    DECLARE_ADJACENCY(patch, face, cell, Face, Cell);
    DECLARE_ADJACENCY(patch, cell, node, Cell, Node);

    bool should_tree = false;
    for(int von_id = 0; von_id < d_improved_stress.size();von_id++){
        if (HAS_ENTITY_SET(patch, d_improved_stress[von_id], CELL, 1)){
            should_tree = true;
        }
    }
    if(should_tree){
        /// 把patch上每一个点的坐标存储到新的J_point向量中
        for(int nn = 0; nn < num_node; nn++){
            double x_nn = (*node_coord)(0,nn);
            double y_nn = (*node_coord)(1,nn);
            double z_nn = (*node_coord)(2,nn);
            J_point xyz(x_nn,y_nn,z_nn);
            Point_on_patch.push_back(xyz);
        }
        std::vector<J_triangle> triangles;
        /// 把patch上的每个三角形存储到新的triangles向量里
        for(int ff = 0; ff < num_face; ff ++){
            int p0 = face_node_idx[face_node_ext[ff]+0];
            int p1 = face_node_idx[face_node_ext[ff]+1];
            int p2 = face_node_idx[face_node_ext[ff]+2];
            triangles.push_back
                    (J_triangle(&Point_on_patch[p0],&Point_on_patch[p1],&Point_on_patch[p2],ff));

        }
        Tree AABB_tree_on_patch(triangles.begin(),triangles.end());
        AABB_tree_on_patch.build();
        /// 建树后进行读取和寻找
        /// 输出文件名
        stringstream OutputFile;
        int patch_number = patch.getIndex();
        OutputFile<<"RecoveredVonMisesOnPatch"<<patch_number;
        string OutFileName = OutputFile.str();
        char ch[40];strcpy(ch,OutFileName.c_str());
        ifstream vonmises_file;
        ofstream vonmises_output;
        vonmises_file.open("../input/MisesQuery.dat",ios::in);
        if (!vonmises_file) {
            TBOX_ERROR("Von Mises传导插值文件未找到！");
        }
        string line;
        /// 逐行读取Von mises插值文件里的坐标点
        while( getline(vonmises_file, line) ){
            stringstream buf(line);
            int index;int Quad;double coord_Q[NDIM];
            buf >> coord_Q[0];buf >> coord_Q[1];buf >> coord_Q[2];
            CGAL_K::Point_3 point_query(coord_Q[0], coord_Q[1], coord_Q[2]);
            /// 找到最近的点
            Tree::Point_and_primitive_id
                    pp = AABB_tree_on_patch.closest_point_and_primitive(point_query);
            /// 最近邻的三角形面
            Tree::Primitive_id cc = pp.second;
            int global_face = cc->id;
            /// 看下它在哪个邻接体中
            for(int local_cell = 0;
                local_cell < face_cell_ext[global_face+1]-face_cell_ext[global_face];
                local_cell ++){
                int global_cell = face_cell_idx[face_cell_ext[global_face]+local_cell];
                tbox::Vector<double> QuadVec(4);
                tbox::Vector<double> QuadSol(4);
                bool incell = false;
                double local_coord[4][3];
                double candidate_weight[NDIM+1];
                double node_temper[NDIM+1]= {0.,0.,0.,0.};
                for(int can = cell_node_ext[global_cell]; can < cell_node_ext[global_cell+1]; can ++){
                    int ln = can - cell_node_ext[global_cell];
                    int gn = cell_node_idx[can];
                    local_coord[ln][0] = (*node_coord)(0,gn);
                    local_coord[ln][1] = (*node_coord)(1,gn);
                    local_coord[ln][2] = (*node_coord)(2,gn);

                }

                incell = PatchPointWeightInCell(coord_Q, &local_coord[0][0],candidate_weight);
                /// 如果在这个体里面的话，输出估计的Von Mises应力和恢复后的Von Mises应力
                if(incell && global_cell < num_cell){
                    for(int nn = 0; nn < NDIM+1; nn ++){
                        int gn = cell_node_idx[cell_node_ext[global_cell]+nn];
                        node_temper[nn] = (*node_temp)(0,gn);
                    }

                    double output_temp =
                            node_temper[0]*candidate_weight[0]
                            +node_temper[1]*candidate_weight[1]
                            +node_temper[2]*candidate_weight[2]
                            +node_temper[3]*candidate_weight[3];
                    double original_mises = (*von_Mises)(0,global_cell);
                    double rec_stress_0 = coord_Q[0]*(*final_improved_coef)(0,global_cell)
                            +coord_Q[1]*(*final_improved_coef)(1,global_cell)
                            +coord_Q[2]*(*final_improved_coef)(2,global_cell)
                            +1.*(*final_improved_coef)(3,global_cell);
                    double rec_stress_1 = coord_Q[0]*(*final_improved_coef)(4,global_cell)
                            +coord_Q[1]*(*final_improved_coef)(5,global_cell)
                            +coord_Q[2]*(*final_improved_coef)(6,global_cell)
                            +1.*(*final_improved_coef)(7,global_cell);
                    double rec_stress_2 = coord_Q[0]*(*final_improved_coef)(8,global_cell)
                            +coord_Q[1]*(*final_improved_coef)(9,global_cell)
                            +coord_Q[2]*(*final_improved_coef)(10,global_cell)
                            +1.*(*final_improved_coef)(11,global_cell);
                    double rec_stress_3 = coord_Q[0]*(*final_improved_coef)(12,global_cell)
                            +coord_Q[1]*(*final_improved_coef)(13,global_cell)
                            +coord_Q[2]*(*final_improved_coef)(14,global_cell)
                            +1.*(*final_improved_coef)(15,global_cell);
                    double rec_stress_4 = coord_Q[0]*(*final_improved_coef)(16,global_cell)
                            +coord_Q[1]*(*final_improved_coef)(17,global_cell)
                            +coord_Q[2]*(*final_improved_coef)(18,global_cell)
                            +1.*(*final_improved_coef)(19,global_cell);
                    double rec_stress_5 = coord_Q[0]*(*final_improved_coef)(20,global_cell)
                            +coord_Q[1]*(*final_improved_coef)(21,global_cell)
                            +coord_Q[2]*(*final_improved_coef)(22,global_cell)
                            +1.*(*final_improved_coef)(23,global_cell);
                    double temp_v =
                            pow((rec_stress_0-rec_stress_1),2)+pow((rec_stress_1-rec_stress_2),2)
                            +pow((rec_stress_2-rec_stress_0),2)
                            +6*(pow(rec_stress_3,2)+pow(rec_stress_4,2)+pow(rec_stress_5,2));
                    double rec_mises = sqrt(temp_v/2);

                    vonmises_output.open(ch,ios::app);
                    vonmises_output<< std::fixed << std::setprecision(12);
                    vonmises_output<<coord_Q[0]<<"\t"<<coord_Q[1]<<"\t"<<coord_Q[2]<<"\t"
                                              <<output_temp<<"\t"
                                             <<original_mises<<"\t"
                                            <<rec_mises<<endl;
                    vonmises_output.close();
                    break;

                }

            }

        }

    }


}

/*************************************************************************
 *  建立网格片上的右端项.
 *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
 *  update #6 @4 更改右端项
 *  基于前两步位移计算右端项
 ************************************************************************/
void PatchStrategy::buildRHSOnPatch(hier::Patch<NDIM>& patch, const double time,
                                    const double dt,
                                    const string& component_name) {

  /// 取出本地PatchGeometry.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();
  /// 取出本地PatchTopology.
  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
      patch.getPatchTopology();
  /// 取出本地Patch的结点坐标数组.
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch_geo->getNodeCoordinates();
  int* dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_data =
      patch.getPatchData(d_rhs_id);

  /// 获取单元对应实体编号数组对象 update #1  @8 -05-05pm
  tbox::Pointer<pdat::CellData<NDIM, int> > entityid_data =
      patch.getPatchData(d_EntityIdOfCell_id);
  tbox::Pointer<pdat::CellData<NDIM, int> > materialid_data =
      patch.getPatchData(material_num_id);

  /////////////////////////////////////////update #6 @5//////////////////////////////////////
  //@param1 vec_DispData:当前时刻单元解向量
  //@param1 vec_DispData_older:前两步单元解向量
  //@param2 vec_DispData_old:前一步单元解向量
  //@param3 vec_RhsData_older:前两步单元右端项
  //@param4 vec_RhsData_old:前一步单元右端项
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_DispData =
      patch.getPatchData(d_solution_id);
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_DispData_older =
      patch.getPatchData(d_solution_older_id);
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_DispData_old =
      patch.getPatchData(d_solution_old_id);
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_RhsData_older =
      patch.getPatchData(d_rhs_older_id);
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_RhsData_old =
      patch.getPatchData(d_rhs_old_id);

  //update #8 取前一步的温度分布
  tbox::Pointer<pdat::NodeData<NDIM, double> > T_data =
      patch.getPatchData(th_Told_id);//
  tbox::Pointer<pdat::NodeData<NDIM, double> > Tolder_data =
      patch.getPatchData(th_Tolder_id);//

 // cout<<" building RHS"<<endl;
  //update #6 @6
  //更新前两步、前一步的解向量和右端项的值
  int num_nodes = patch.getNumberOfNodes(1);//取本地patch节点数目
  for (int i2 = 0; i2 < NDIM * num_nodes; ++i2) {
	  //更新前两步及前一步的解向量值
	  vec_DispData_older->getPointer()[i2] = vec_DispData_old->getPointer()[i2];
	  vec_DispData_older->getPointer()[i2] = vec_DispData->getPointer()[i2];
	  //更新前两步及前一步的右端项值
	  vec_RhsData_older->getPointer()[i2] = vec_RhsData_old->getPointer()[i2];
	  vec_RhsData_old->getPointer()[i2] = vec_data->getPointer()[i2];
  }
  //////////////////////////////////////////////////////////////////////////////////////////////

  /// 获取单元周围结点的索引关系.
  tbox::Array<int> can_extent, can_indices;
  patch_top->getCellAdjacencyNodes(can_extent, can_indices);

  /// 取出本地Patch的单元数目.
  int num_cells = patch.getNumberOfCells(1);
  //  int num_nodes = patch.getNumberOfNodes(0);

  for (int i = 0; i < num_cells; ++i) {
    int n_vertex = can_extent[i + 1] - can_extent[i];
    int n_dof = NDIM * n_vertex;

    //存储上一时刻单元节点温度
    tbox::Array<double> T_val(n_vertex); //当前温度
    tbox::Array<double> Tolder_val(n_vertex);

    /**< 该单元的结点坐标及自由度映射 */
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);

    ////////////////////////////////////////////update #6 @7////////////////////////////////////////
    //定义一个结构体数组，数量为单元的节点数，
    //内容为Newmark-beta方法所需的解向量及右端项值，具体见MacrosManager.h中的相关定义
    //@param d_newmark[i]:单元的第i个节点的NewmarkData
    //其中：
    //@param1 v_solution_older:前两步单元解向量
    //@param2 v_solution_old:前一步单元解向量
    //@param3 v_rhs_older:前两步单元右端项
    //@param4 v_rhs_old:前一步单元右端项
    NewmarkData d_newmark[n_vertex];

    //赋值
    for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
      for (int k = 0; k < NDIM; ++k) {
    	  //VectorData不支持(*v)(i,j)的索引方式
    	  d_newmark[i1].v_solution_older[k] = vec_DispData_older->getPointer()[can_indices[j]*NDIM+k];
    	  d_newmark[i1].v_solution_old[k] = vec_DispData_old->getPointer()[can_indices[j]*NDIM+k];
    	  d_newmark[i1].v_rhs_older[k] = vec_RhsData_older->getPointer()[can_indices[j]*NDIM+k];
    	  d_newmark[i1].v_rhs_old[k] = vec_RhsData_old->getPointer()[can_indices[j]*NDIM+k];
      }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////

    tbox::Array<int> node_mapping(n_dof);

    for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
    	T_val[i1]=T_data->getPointer()[can_indices[j]];
        Tolder_val[i1]=Tolder_data->getPointer()[can_indices[j]];
      for (int k = 0; k < NDIM; ++k) {
        node_mapping[NDIM * i1 + k] = dof_map[can_indices[j]] + k;
        vertex[i1][k] = (*node_coord)(k, can_indices[j]);
      }
    }
    /// 取出单元对象.
    tbox::Pointer<BaseElement<NDIM> > ele =
        d_element_manager->getElement(d_element_type[0]);

    tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
    ele_vec->resize(n_dof);
    for (int j = 0; j < n_dof; ++j) {
      (*ele_vec)[j] = 0.0;
    }
    /// 退火过程有固定温度

    /// 计算单元右端项.
    /// Vinta Yin-Da Wang:静力学
    T_val[0] = 650.;T_val[1] = 650.;T_val[2] = 650.;T_val[3] = 650.;
    T_val[0] = 300.;T_val[1] = 300.;T_val[2] = 300.;T_val[3] = 300.;
    ele->buildElementRHS(vertex, dt, time, ele_vec, d_newmark, (*materialid_data)(0,i),T_val,Tolder_val);
    for (int ii = 0; ii < n_dof; ++ii)
      vec_data->addVectorValue(node_mapping[ii], (*ele_vec)[ii]);
  }
}

/*************************************************************************
 *  建立网格片上的矩阵.
 ************************************************************************/
void PatchStrategy::buildMatrixOnPatch(hier::Patch<NDIM>& patch,
                                       const double time, const double dt,
                                       const string& component_name) {
    /// 取出本地PatchGeometry.
    tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
            patch.getPatchGeometry();
    /// 取出本地PatchTopology.
    tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
            patch.getPatchTopology();
    /// 取出本地Patch的结点坐标数组.
    tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
            patch_geo->getNodeCoordinates();

    /// 获取单元对应实体编号数组对象 update #1  @7
    tbox::Pointer<pdat::CellData<NDIM, int> > entityid_data =
            patch.getPatchData(d_EntityIdOfCell_id);
    tbox::Pointer<pdat::CellData<NDIM, int> > materialid_data =
            patch.getPatchData(material_num_id);
    //update #8 取前一步的温度分布
    tbox::Pointer<pdat::NodeData<NDIM, double> > T_data =
            patch.getPatchData(th_Told_id);//

    /// 获取单元周围结点的索引关系.
    tbox::Array<int> can_extent, can_indices;
    patch_top->getCellAdjacencyNodes(can_extent, can_indices);

    /// 取出Patch的单元数目.
    int num_cells = patch.getNumberOfCells(1);
    int* dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
    tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > mat_data =
            patch.getPatchData(d_matrix_id);

    for (int i = 0; i < num_cells; ++i) {
        int n_vertex = can_extent[i + 1] - can_extent[i];
        int n_dof = NDIM * n_vertex;
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<int> node_mapping(n_dof);
        tbox::Array<double> T_val(n_vertex); //当前温度
        double elem_T = 293.15;
        /// 这个程序是仿线弹性力学的，所以elem_T要换

        for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
            T_val[i1]=T_data->getPointer()[can_indices[j]];
            for (int k = 0; k < NDIM; ++k) {
                node_mapping[NDIM * i1 + k] = dof_map[can_indices[j]] + k;
                vertex[i1][k] = (*node_coord)(k, can_indices[j]);
            }
        }
        elem_T = (T_val[0]+T_val[1]+T_val[2]+T_val[3])/4;
        elem_T = 600.;

        tbox::Pointer<tbox::Matrix<double> > ele_mat = new tbox::Matrix<double>();
        ele_mat->resize(n_dof, n_dof);
        for (int i1 = 0; i1 < n_dof; ++i1) {
            for (int j = 0; j < n_dof; ++j) {
                (*ele_mat)(i1, j) = 0.0;
            }
        }
        /// 取出单元对象.
        tbox::Pointer<BaseElement<NDIM> > ele =
                d_element_manager->getElement(d_element_type[0]);

        //update #6 @8 通过宏定义值选择程序为静态求解还是瞬态求解
        //静态时
#if STATIC_S
        //update #3   @3 函数添加一个形参
        ele->buildStiffElementMatrix(vertex, dt, time, ele_mat,(*materialid_data)(0,i),T_val);
#endif

        //update #6 @9 瞬态时
        /// Vinta Yin-Da Wang: 应力当作静力学问题去做
#if TIME_S
        //    ele->buildElementMatrix(vertex, dt, time, ele_mat,(*materialid_data)(0,i),T_val);
        ele->buildStiffElementMatrix(vertex, dt, time, ele_mat,(*materialid_data)(0,i),elem_T);
#endif

        /// 累加矩阵
        for (int i1 = 0; i1 < n_dof; ++i1) {
            for (int j = 0; j < n_dof; ++j) {
                mat_data->addMatrixValue
                        (node_mapping[i1], node_mapping[j],(*ele_mat)(i1, j));
            }
        }

    }
    mat_data->assemble();
}


void PatchStrategy::Stress_max(double* vector, int len, hier::Patch<NDIM>& patch, const double time,
                      const double dt, const string& component_name)
{
	  /// 取出本地PatchGeometry.
	  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
	      patch.getPatchGeometry();
  tbox::Pointer<pdat::CellData<NDIM, double> > von_Mises =
      patch.getPatchData(d_von_mises_id);
tbox::Pointer<pdat::NodeData<NDIM,double> > disp = patch.getPatchData(d_displacement_id);
	  int num_cells = patch.getNumberOfCells();
	  int num_nodes = patch.getNumberOfNodes();
	  *vector=(*von_Mises)(0,0);
          vector[1] = (*disp)(0,0);
	  for(int i=0;i<num_cells;i++)
	  {
		  if(*vector<(*von_Mises)(0,i)) *vector=(*von_Mises)(0,i);
	          
	  }
	  for(int i =0; i < num_nodes; i++){
		if(vector[1]<(*disp)(0,i)) vector[1] = (*disp)(0,i);
	  }
}


///////////////////////////////////////////////////////update #8 @8///////////////////////////////////////////////////////
//添加热求解相关函数
//applyTh_Load()
//applyTh_Constrain()
//buildTh_MatrixOnPatch()
//buildTh_RHSOnPatch()

/*************************************************************************
 *
 *  设置热源及第二类边界条件.
 *
 *  未改完
 *
 ************************************************************************/
void PatchStrategy::applyTh_Load(hier::Patch<NDIM>& patch, const double time,
                              const double dt, const string& component_name) {
  /// 取出本地PatchGeometry.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();
  tbox::Pointer<pdat::VectorData<NDIM, double> > th_rhs_data =
      patch.getPatchData(th_rhs_id);

  //update #5 @1 -05-08
  //外加体热源
#if 0
//  int load_size = d_load_types.getSize();
//  for (int k = 0; k < load_size; ++k) {
    // 获取指定编号和类型的集合包含的网格实体的索引。
    if (patch.hasEntitySet(1, hier::EntityUtilities::NODE)) {
      const tbox::Array<int>& entity_idx = patch_geo->getEntityIndicesInSet(
          1, hier::EntityUtilities::NODE);
      // 获取物理边界上节点的数目.
      int size = entity_idx.getSize();
      ///

    }
//  }
#endif

  /// 第二类边界条件
  //update #5 @2  -05-08
#if 0
  //update #5 @3
  //自由度信息中的映射信息
  int* dof_map = th_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
  /// 取出本地Patch的结点坐标数组.
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch_geo->getNodeCoordinates();
  int num_faces =patch.getNumberOfEntities(hier::EntityUtilities::FACE, 0);

  double q=0;//单位面积表面力
  if (patch_geo->hasEntitySet(3, hier::EntityUtilities::FACE)) {
    // 获取指定编号和类型的集合包含的网格实体的索引。
    const tbox::Array<int>& Face_idx = patch_geo->getEntityIndicesInSet(
        3, hier::EntityUtilities::FACE,num_faces);
	  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
			      patch.getPatchTopology();
	  //单元节点邻接关系数据
	  tbox::Array<int>face_node_ext,face_node_idx;
	  patch.getPatchTopology()->getFaceAdjacencyNodes(face_node_ext,face_node_idx);
	  //单元节点邻接关系数据
	  tbox::Array<int>face_cell_ext,face_cell_idx;
	  patch.getPatchTopology()->getFaceAdjacencyCells(face_cell_ext,face_cell_idx);
    int size = Face_idx.getSize();
    for (int face=0; face<size;face++)
    {
    	int n=3;//面单元节点数
    	int dof_num=n;
        tbox::Array<int> node_mapping(dof_num);
    	int node_id[3]={0,0,0};

    	//面上三角形的三个顶点编号
    	node_id[0]=face_node_idx[face_node_ext[Face_idx[face]]];//Face_idx[face]面的编号
    	node_id[1]=face_node_idx[face_node_ext[Face_idx[face]]+1];
    	node_id[2]=face_node_idx[face_node_ext[Face_idx[face]]+2];
    	tbox::Array<hier::DoubleVector<NDIM> > vertex(n);
    	for (int k=0; k<n;k++){
    		for(int j=0; j<n;j++){
        		vertex[k][j]=(*node_coord)(j,node_id[k]);
    		}
    	}

        for (int i1 = 0, j = face_node_ext[Face_idx[face]]; i1 < n; ++i1, ++j) {
            node_mapping[i1] = dof_map[face_node_idx[j]];
        }

        tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
        ele_vec->resize(dof_num);
        for (int i = 0; i < dof_num; ++i) {
          (*ele_vec)[i] = 0.0;
        }

    	double area=sqrt(AREA(vertex[0],vertex[1],vertex[2]))/2.0;
    	for(int ii=0;ii<n;ii++)
    			(*ele_vec)[ii]=q*area/3;
    	for(int ii=0;ii<dof_num;ii++)
    		th_rhs_data->getPointer()[node_mapping[ii]] += (*ele_vec)[ii];
    }
  }
#endif
//cout<<"load ok"<<endl;
}

/*************************************************************************
 *
 *  填充物理边界条件.
 *
 *  已改 -05-08
 ************************************************************************/
void PatchStrategy::applyTh_Constraint(hier::Patch<NDIM>& patch,
                                    const double time, const double dt,
                                    const string& component_name) {
	  /// 取出本地PatchGeometry.
	  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
	      patch.getPatchGeometry();
	  /// 获取网格片矩阵对象
	  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > mat_data =
	      patch.getPatchData(th_matrix_id);
	  /// 获取网格片向量对象
	  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_data =
	      patch.getPatchData(th_rhs_id);
	  /// 结点数目
	  int num_nodes = patch.getNumberOfNodes();
      int num_faces = patch.getNumberOfEntities(hier::EntityUtilities::FACE, 1);

	  /// 自由度信息中的映射信息。
	  int* dof_map = d_dof_info_th->getDOFMapping(patch, hier::EntityUtilities::NODE);
	  /// 取出矩阵向量的核心数据结构的指针
	  int* row_start = mat_data->getRowStartPointer();
	  int* col_idx = mat_data->getColumnIndicesPointer();
	  double* mat_val = mat_data->getValuePointer();
	  double* vec_val = vec_data->getPointer();
      tbox::Array<int>face_node_ext,face_node_idx;
      patch.getPatchTopology()->getFaceAdjacencyNodes(face_node_ext,face_node_idx);
      tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
              patch_geo->getNodeCoordinates();

	  /// 边界条件处理.
	  if (patch_geo->hasEntitySet(3, hier::EntityUtilities::NODE)) {
	    // 获取指定编号和类型的集合包含的网格实体的索引。
	    const tbox::Array<int>& entity_idx = patch_geo->getEntityIndicesInSet(
	        3, hier::EntityUtilities::NODE, num_nodes);
	    int size = entity_idx.getSize();
	    /// 对角化 1 法处理约束。
	    for (int i = 0; i < size; ++i) {
	      int index = dof_map[entity_idx[i]];
          vec_val[index] =293.15; /**< 载荷项设为 0 */
	      for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
	        if (col_idx[j] == index) {
	          mat_val[j] = 1.0; /**< 对角线元素 */
	        } else {
	          mat_val[j] = 0.0;                     /**< 非对角线元素 */
	          vec_val[col_idx[j]] = vec_val[col_idx[j]] - (*mat_data)(col_idx[j], index)*vec_val[index];
	          (*mat_data)(col_idx[j], index) = 0.0; /**< 对应的列元素 */
	        }
	      }
	    }
	  }
#if 0
	  /// 边界条件处理.
	  if (patch_geo->hasEntitySet(2, hier::EntityUtilities::NODE)) {
	    // 获取指定编号和类型的集合包含的网格实体的索引。
	    const tbox::Array<int>& entity_idx = patch_geo->getEntityIndicesInSet(
	        2, hier::EntityUtilities::NODE, num_nodes);
	    int size = entity_idx.getSize();
	    /// 对角化 1 法处理约束。
	    for (int i = 0; i < size; ++i) {
	      int index = dof_map[entity_idx[i]];
	      vec_val[index] =300.0; /**< 载荷项设为 0 */
	      for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
	        if (col_idx[j] == index) {
	          mat_val[j] = 1.0; /**< 对角线元素 */
	        } else {
	          mat_val[j] = 0.0;                     /**< 非对角线元素 */
	          vec_val[col_idx[j]] = vec_val[col_idx[j]] - (*mat_data)(col_idx[j], index)*vec_val[index];
	          (*mat_data)(col_idx[j], index) = 0.0; /**< 对应的列元素 */
	        }
	      }
	    }
	  }
#endif
	  //  cout<<"constrain ok"<<endl;
      for (int conv_id = 0; conv_id < d_convection_boundary.size(); conv_id++){
          double conv_coef = 15.;
          double Text = 293.15;
          if (patch_geo->hasEntitySet(d_convection_boundary[conv_id],
                              hier::EntityUtilities::FACE,patch.getNumberOfFaces(1))){
              const tbox::Array<int>& Face_idx =
                                      patch_geo->getEntityIndicesInSet(
                                              d_convection_boundary[conv_id],
                                              hier::EntityUtilities::FACE, num_faces);
              int size = Face_idx.getSize();
              for (int face=0; face<size;face++){
                  int node_id[3]={0,0,0};
                  int n=3;
                  double b[3]={0,0,0};
                  double k[3][3]={{0,0,0},{0,0,0},{0,0,0}};
                  //面上三角形的三个顶点编号
                  node_id[0]=face_node_idx[face_node_ext[Face_idx[face]]];//Face_idx[face]面的编号
                  node_id[1]=face_node_idx[face_node_ext[Face_idx[face]]+1];
                  node_id[2]=face_node_idx[face_node_ext[Face_idx[face]]+2];
                  tbox::Array<hier::DoubleVector<NDIM> > vertex(n);
                  for (int k=0; k<n;k++){
                      for(int j=0; j<n;j++){
                          vertex[k][j]=(*node_coord)(j,node_id[k]);
                      }
                  }
                  double area = sqrt(AREA(vertex[0], vertex[1], vertex[2]))
                                              / 2.0;
                  for (int i = 0; i < 3; i++) {
                      b[i] =  conv_coef * (Text) * area / 3.0;
                      for (int j = 0; j < 3; j++) {
                          if (i == j)
                              k[i][j] =  conv_coef * area / 6.0;
                          else
                              k[i][j] =  conv_coef * area / 12.0;
                      }

                  }
                  for (int i = 0; i < 3; i++) {
                      vec_data->addVectorValue(dof_map[node_id[i]], b[i]*dt);
                      for (int j = 0; j < 3; j++) {
                          mat_data->addMatrixValue(dof_map[node_id[i]], dof_map[node_id[j]],  k[i][j]*dt);
                      }

                  }
              }

          }


      }
}

/*************************************************************************
 *  建立网格片上thermal右端项.
 *
 *  已改
 ************************************************************************/
void PatchStrategy::buildTh_RHSOnPatch(hier::Patch<NDIM>& patch,
                                            const double time, const double dt,
                                            const string& component_name) {
  /// 取出本地PatchGeometry.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();
  /// 取出本地PatchTopology.
  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
      patch.getPatchTopology();
  /// 取出本地Patch的结点坐标数组.
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch_geo->getNodeCoordinates();
  /// 获取本地单元(边, 面)周围结点的索引关系.
  tbox::Array<int> can_extent, can_indices;
  patch_top->getCellAdjacencyNodes(can_extent, can_indices);
  /// 取出向量数据片
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_data =
      patch.getPatchData(th_rhs_id);

  //取前一步的温度分布
  tbox::Pointer<pdat::NodeData<NDIM, double> > T_data =
      patch.getPatchData(th_Told_id);//

  /// 取出单元数目
  int num_cells = patch.getNumberOfCells(1);

  /// 自由度映射信息
  int* dof_map = d_dof_info_th->getDOFMapping(patch, hier::EntityUtilities::NODE);

  //update #1 @8
  tbox::Pointer<pdat::CellData<NDIM, int> > entityid_data =
      patch.getPatchData(d_EntityIdOfCell_id);
  tbox::Pointer<pdat::CellData<NDIM, int> > materialid_data =
      patch.getPatchData(material_num_id);

  //update #9
  tbox::Pointer<pdat::CellData<NDIM, double> > Emag_data =
      patch.getPatchData(E_mag_id);

  /**
   * 遍历单元，将单元右端项加载到右端项向量中。
   */
  for (int i = 0; i < num_cells; ++i) {
	int cell = i;
    int n_vertex = can_extent[i + 1] - can_extent[i];
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<int> mapping(n_vertex);
    tbox::Array<double> T_val(n_vertex);

    /// 下面的循环做两件事情：1. 建立自由度映射数组；2.取出结点坐标。
    for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
      mapping[i1] = dof_map[can_indices[j]];
      T_val[i1]=T_data->getPointer()[can_indices[j]];
      for (int k = 0; k < NDIM; ++k) {
        vertex[i1][k] = (*node_coord)(k, can_indices[j]);
      }
    }

    //update #9
	  tbox::Pointer<MaterialManager<NDIM> > material_manager =
	      MaterialManager<NDIM>::getManager();
    tbox::Pointer<Material> material =
    				material_manager->getMaterial("Gold");

    if((*materialid_data)(0,i)==1)
		 material =	material_manager->getMaterial("Silicon");
    else if((*materialid_data)(0,i)==13)
		 material =	material_manager->getMaterial("MoCu");
    else if((*materialid_data)(0,i)==2)
		 material =	material_manager->getMaterial("Copper");
    else if((*materialid_data)(0,i)==4)
		 material =	material_manager->getMaterial("SiO2");
    else if((*materialid_data)(0,i)==5)
		 material =	material_manager->getMaterial("SiN");     
    else if((*materialid_data)(0,i)==3)
		 material =	material_manager->getMaterial("Gold");
    else if((*materialid_data)(0,i)==15)
		 material =	material_manager->getMaterial("GaN");  
    else if((*materialid_data)(0,i)==16)
		 material =	material_manager->getMaterial("Al2O3"); 
    else if((*materialid_data)(0,i)==17)
		 material =	material_manager->getMaterial("Alloy");
    else if((*materialid_data)(0,i)==6)
		 material =	material_manager->getMaterial("Aluminum");  
    else
		 material =	material_manager->getMaterial("Gold");

   double temp_T=0;
    for(int t=0;t<4;t++)
    {temp_T+=T_val[t]/4;}
    double Sigma=material->getSigma(temp_T);

//暂时没有添加电求解
//#if 0
    //定义每个单元上的电热源
    double u_val=((*Emag_data)(0,cell))*Sigma;
//#endif
    double e_ThermalSource=u_val;//系统单元热源
    double q=0.;//50e8/0.6438;//30e9/1.033;;//10e9/2.5;//q=5e12/2.929;//单位体积热源
//    if((*entityid_data)(0,cell)==70||(*entityid_data)(0,cell)==71)
//        e_ThermalSource=e_ThermalSource+q;


    /// 取出单元对象.
    tbox::Pointer<BaseElement<NDIM> > ele =
        d_element_manager->getElement(d_element_type[0]);

    /// 声明vec，用来存储单元右端项。
    tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
    ele_vec->resize(n_vertex);
    for (int i = 0; i < n_vertex; ++i) {
      (*ele_vec)[i] = 0.0;
    }

    /// 计算单元右端项，并组装。
    ele->buildTh_ElementRHS(vertex, dt, time, ele_vec, (*materialid_data)(0,cell),T_val, e_ThermalSource);
    for (int i2 = 0; i2 < n_vertex; ++i2)
      vec_data->addVectorValue(mapping[i2], (*ele_vec)[i2]);
  }
  //cout<<"rhs ok"<<endl;
}


/*************************************************************************
 *  建立网格片上的thermal矩阵
 *
 *  已改
 ************************************************************************/
void PatchStrategy::buildTh_MatrixOnPatch(hier::Patch<NDIM>& patch,
                                               const double time,
                                               const double dt,
                                               const string& component_name) {
  /// 取出本地PatchGeometry.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();
  /// 取出本地PatchTopology.
  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
      patch.getPatchTopology();
  /// 取出本地Patch的结点坐标数组.

  /// 获取单元对应实体编号数组对象 update #1  @7
  tbox::Pointer<pdat::CellData<NDIM, int> > entityid_data =
      patch.getPatchData(d_EntityIdOfCell_id);
  tbox::Pointer<pdat::CellData<NDIM, int> > materialid_data =
      patch.getPatchData(material_num_id);

  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch_geo->getNodeCoordinates();
  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > mat_data =
      patch.getPatchData(th_matrix_id);
  //取前一步的温度分布
  tbox::Pointer<pdat::NodeData<NDIM, double> > T_data =
      patch.getPatchData(th_Told_id);//

  /// 获取本地单元(边, 面)周围结点的索引关系.
  tbox::Array<int> can_extent, can_indices;
  patch_top->getCellAdjacencyNodes(can_extent, can_indices);
  int num_cells = patch.getNumberOfCells(1);

  /// 取出自由度映射信息
  int* dof_map = d_dof_info_th->getDOFMapping(patch, hier::EntityUtilities::NODE);
  /// 遍历单元
  for (int i = 0; i < num_cells; ++i) {
	  int cell=i;
    int n_vertex = can_extent[i + 1] - can_extent[i];
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<int> mapping(n_vertex);
    tbox::Array<double> T_val(n_vertex);
    /// 取出单元结点坐标，以及填写映射值。//can_indices 节点编号
    for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
      mapping[i1] = dof_map[can_indices[j]];
      T_val[i1]=T_data->getPointer()[can_indices[j]];
      for (int k = 0; k < NDIM; ++k) {
        vertex[i1][k] = (*node_coord)(k, can_indices[j]);
      }
    }
    /// 初始化单元矩阵
    tbox::Pointer<tbox::Matrix<double> > ele_mat = new tbox::Matrix<double>();
    ele_mat->resize(n_vertex, n_vertex);
    for (int i = 0; i < n_vertex; ++i) {
      for (int j = 0; j < n_vertex; ++j) {
        (*ele_mat)(i, j) = 0.0;
      }
    }

    /// 取出单元对象.
    tbox::Pointer<BaseElement<NDIM> > ele =
        d_element_manager->getElement(d_element_type[0]);

    /// 计算单元矩阵
    ele->buildTh_ElementMatrix(vertex, dt, time, ele_mat,(*materialid_data)(0,cell),T_val);

    /// 累加单元矩阵，这里用户根据映射将矩阵挨个添加。
    int row = 0, col = 0;
    for (int i = 0; i < n_vertex; ++i) {
      row = mapping[i];
      for (int j = 0; j < n_vertex; ++j) {
        col = mapping[j];
        mat_data->addMatrixValue(row, col, (*ele_mat)(i, j));
      }
    }
  }
  /// 矩阵组装。
  mat_data->assemble();
}

void PatchStrategy::Thermal_PostProcesing(hier::Patch<NDIM>& patch, const double time,
                    const double dt, const string& component_name)
{
	  /// 取出本地PatchGeometry.
	  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
	      patch.getPatchGeometry();
	  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_data =
	      patch.getPatchData(th_solution_id);
	  tbox::Pointer<pdat::NodeData<NDIM, double> > plot_data =
	      patch.getPatchData(th_plot_id);
	  tbox::Pointer<pdat::NodeData<NDIM, double> > old_Tdata =
	      patch.getPatchData(th_Told_id);
	  tbox::Pointer<pdat::NodeData<NDIM, double> > older_Tdata =
	      patch.getPatchData(th_Tolder_id);
      int* dof_map = d_dof_info_th->getDOFMapping(patch, hier::EntityUtilities::NODE);
      double* vec_pointer = vec_data->getPointer(0);
	  int local_num_nodes = patch.getNumberOfNodes();
	  for(int i=0;i<local_num_nodes;i++)
          (*plot_data)(0,i)=vec_pointer[dof_map[i]];

	  int num_nodes = patch.getNumberOfNodes(1);

      for(int i=0;i<local_num_nodes;i++)
	  {
		  (*older_Tdata)(0,i)=(*old_Tdata)(0,i);
          (*old_Tdata)(0,i)=vec_pointer[dof_map[i]];
		  
	  }

}

void PatchStrategy::Thermal_max(double* vector, int len, hier::Patch<NDIM>& patch, const double time,
                      const double dt, const string& component_name)
{
	  /// 取出本地PatchGeometry.
	  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
	      patch.getPatchGeometry();
	  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_data =
	      patch.getPatchData(th_solution_id);

	  int num_nodes = patch.getNumberOfNodes(1);

	  *vector=vec_data->getPointer()[0];
	  for(int i=0;i<num_nodes;i++)
	  {
		  if(*vector<vec_data->getPointer()[i]) *vector=vec_data->getPointer()[i];
	  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////update #9////////////////////////////////////////////////////
/*!
 * @brief 支撑指定名称的数值构件, 在单个网格片上完成右端项组装.
 *
 * @param patch          输入参数, 网格片类, 表示网格片.
 * @param time           输入参数, 双精度浮点型, 表示当前时刻.
 * @param dt             输入参数, 双精度浮点型, 表示时间步长.
 * @param component_name 输入参数, 字符串, 表示数值构件的名称.
 *
 */
void PatchStrategy::buildE_RHSOnPatch(hier::Patch<NDIM>& patch, const double time,
                     const double dt, const string& component_name)
{
	 /// 取出本地PatchGeometry.
	  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
	      patch.getPatchGeometry();
	  /// 取出本地PatchTopology.
	  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
	      patch.getPatchTopology();
	  /// 取出本地Patch的结点坐标数组.
	  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
	      patch_geo->getNodeCoordinates();
	  /// 获取本地单元(边, 面)周围结点的索引关系.
	  tbox::Array<int> can_extent, can_indices;
	  patch_top->getCellAdjacencyNodes(can_extent, can_indices);
	  /// 取出向量数据片
	  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_data =
	      patch.getPatchData(E_rhs_id);

	  //取前一步的温度分布
	  tbox::Pointer<pdat::NodeData<NDIM, double> > T_data =
	      patch.getPatchData(th_Told_id);//

	  /// 取出单元数目
	  int num_cells = patch.getNumberOfCells(1);

	  /// 自由度映射信息
	  int* dof_map = d_dof_info_th->getDOFMapping(patch, hier::EntityUtilities::NODE);

	  //update #1 @8
	  tbox::Pointer<pdat::CellData<NDIM, int> > entityid_data =
	      patch.getPatchData(d_EntityIdOfCell_id);
          tbox::Pointer<pdat::CellData<NDIM, int> > materialid_data =
              patch.getPatchData(material_num_id);

	  /**
	   * 遍历单元，将单元右端项加载到右端项向量中。
	   */
	  for (int i = 0; i < num_cells; ++i) {
		int cell = i;
	    int n_vertex = can_extent[i + 1] - can_extent[i];
	    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
	    tbox::Array<int> mapping(n_vertex);
	    tbox::Array<double> T_val(n_vertex);

	    /// 下面的循环做两件事情：1. 建立自由度映射数组；2.取出结点坐标。
	    for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
	      mapping[i1] = dof_map[can_indices[j]];
	      T_val[i1]=T_data->getPointer()[can_indices[j]];
	      for (int k = 0; k < NDIM; ++k) {
	        vertex[i1][k] = (*node_coord)(k, can_indices[j]);
	      }
	    }

	    /// 取出单元对象.
	    tbox::Pointer<BaseElement<NDIM> > ele =
	        d_element_manager->getElement(d_element_type[0]);

	    /// 声明vec，用来存储单元右端项。
	    tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
	    ele_vec->resize(n_vertex);
	    for (int i = 0; i < n_vertex; ++i) {
	      (*ele_vec)[i] = 0.0;
	    }

	    /// 计算单元右端项，并组装。
	    ele->buildE_ElementRHS(vertex, dt, time, ele_vec, (*materialid_data)(0,cell),T_val);
	    for (int i2 = 0; i2 < n_vertex; ++i2)
	      vec_data->addVectorValue(mapping[i2], (*ele_vec)[i2]);
	  }
	//  cout<<"rhs ok"<<endl;
}

/*!
 * @brief 在单个网格片上设置电压约束.
 *
 * @param patch          输入参数, 网格片类, 表示网格片.
 * @param time           输入参数, 双精度浮点型, 表示当前时刻.
 * @param dt             输入参数, 双精度浮点型, 表示时间步长.
 * @param component_name 输入参数, 字符串, 表示数值构件的名称.
 *
 */
void PatchStrategy::applyE_Constraint(hier::Patch<NDIM>& patch, const double time,
                     const double dt, const string& component_name)
{
	  /// 取出本地PatchGeometry.
	  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
	      patch.getPatchGeometry();
	  /// 获取网格片矩阵对象
	  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > mat_data =
	      patch.getPatchData(E_matrix_id);
	  /// 获取网格片向量对象
	  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_data =
	      patch.getPatchData(E_rhs_id);
	  /// 结点数目
	  int num_nodes = patch.getNumberOfNodes();

	  /// 自由度信息中的映射信息。
	  int* dof_map = d_dof_info_th->getDOFMapping(patch, hier::EntityUtilities::NODE);
	  /// 取出矩阵向量的核心数据结构的指针
	  int* row_start = mat_data->getRowStartPointer();
	  int* col_idx = mat_data->getColumnIndicesPointer();
	  double* mat_val = mat_data->getValuePointer();
	  double* vec_val = vec_data->getPointer();

	int step=(time+0.5*dt)/dt;
//#if 0
	  /// 边界条件处理.
	  if (patch_geo->hasEntitySet(1, hier::EntityUtilities::NODE)) {
	    // 获取指定编号和类型的集合包含的网格实体的索引。
	    const tbox::Array<int>& entity_idx = patch_geo->getEntityIndicesInSet(
	        1, hier::EntityUtilities::NODE, num_nodes);
	    int size = entity_idx.getSize();
	    /// 对角化 1 法处理约束。
	    for (int i = 0; i < size; ++i) {
	      int index = dof_map[entity_idx[i]];
     double vpp = vpulse2[step+1];
		  //if(vpulse2[step]!=vpulse2[step+1])
			  //vpp = (vpulse2[step]+vpulse2[step+1])/2;
          vec_val[index] =vpp; /**< 载荷项设为 0 */
	      for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
	        if (col_idx[j] == index) {
	          mat_val[j] = 1.0; /**< 对角线元素 */
	        } else {
	          mat_val[j] = 0.0;                     /**< 非对角线元素 */
	          vec_val[col_idx[j]] = vec_val[col_idx[j]] - (*mat_data)(col_idx[j], index)*vec_val[index];
	          (*mat_data)(col_idx[j], index) = 0.0; /**< 对应的列元素 */
	        }
	      }
	    }
	  }
    if (patch_geo->hasEntitySet(4, hier::EntityUtilities::NODE)) {
	    // 获取指定编号和类型的集合包含的网格实体的索引。
	    const tbox::Array<int>& entity_idx = patch_geo->getEntityIndicesInSet(
	        4, hier::EntityUtilities::NODE, num_nodes);
	    int size = entity_idx.getSize();
	    /// 对角化 1 法处理约束。
	    for (int i = 0; i < size; ++i) {
	      int index = dof_map[entity_idx[i]];
          vec_val[index] =0; /**< 载荷项设为 0 */
	      for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
	        if (col_idx[j] == index) {
	          mat_val[j] = 1.0; /**< 对角线元素 */
	        } else {
	          mat_val[j] = 0.0;                     /**< 非对角线元素 */
	          vec_val[col_idx[j]] = vec_val[col_idx[j]] - (*mat_data)(col_idx[j], index)*vec_val[index];
	          (*mat_data)(col_idx[j], index) = 0.0; /**< 对应的列元素 */
	        }
	      }
	    }
	  }
#if 0
	  if (patch_geo->hasEntitySet(3, hier::EntityUtilities::NODE)) {
	    // 获取指定编号和类型的集合包含的网格实体的索引。
	    const tbox::Array<int>& entity_idx = patch_geo->getEntityIndicesInSet(
	        3, hier::EntityUtilities::NODE, num_nodes);
	    int size = entity_idx.getSize();
	    /// 对角化 1 法处理约束。
	    for (int i = 0; i < size; ++i) {
	      int index = dof_map[entity_idx[i]];
	      vec_val[index] =vpulse3[step]; /**< 载荷项设为 0 */
	      for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
	        if (col_idx[j] == index) {
	          mat_val[j] = 1.0; /**< 对角线元素 */
	        } else {
	          mat_val[j] = 0.0;                     /**< 非对角线元素 */
	          vec_val[col_idx[j]] = vec_val[col_idx[j]] - (*mat_data)(col_idx[j], index)*vec_val[index];
	          (*mat_data)(col_idx[j], index) = 0.0; /**< 对应的列元素 */
	        }
	      }
	    }
	  }
#endif
#if 0
	  /// 边界条件处理.
	  if (patch_geo->hasEntitySet(2, hier::EntityUtilities::NODE)) {
	    // 获取指定编号和类型的集合包含的网格实体的索引。
	    const tbox::Array<int>& entity_idx = patch_geo->getEntityIndicesInSet(
	        2, hier::EntityUtilities::NODE, num_nodes);
	    int size = entity_idx.getSize();
	    /// 对角化 1 法处理约束。
	    for (int i = 0; i < size; ++i) {
	      int index = dof_map[entity_idx[i]];
	      vec_val[index] =0; /**< 载荷项设为 0 */
	      for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
	        if (col_idx[j] == index) {
	          mat_val[j] = 1.0; /**< 对角线元素 */
	        } else {
	          mat_val[j] = 0.0;                     /**< 非对角线元素 */
	          vec_val[col_idx[j]] = vec_val[col_idx[j]] - (*mat_data)(col_idx[j], index)*vec_val[index];
	          (*mat_data)(col_idx[j], index) = 0.0; /**< 对应的列元素 */
	        }
	      }
	    }
	  }
#endif
#if 0
  /// 边界条件处理.
	  if (patch_geo->hasEntitySet(4, hier::EntityUtilities::NODE)) {
	    // 获取指定编号和类型的集合包含的网格实体的索引。
	    const tbox::Array<int>& entity_idx = patch_geo->getEntityIndicesInSet(
	        4, hier::EntityUtilities::NODE, num_nodes);
	    int size = entity_idx.getSize();
	    /// 对角化 1 法处理约束。
	    for (int i = 0; i < size; ++i) {
	      int index = dof_map[entity_idx[i]];
	      vec_val[index] =0.1; /**< 载荷项设为 0 */
	      for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
	        if (col_idx[j] == index) {
	          mat_val[j] = 1.0; /**< 对角线元素 */
	        } else {
	          mat_val[j] = 0.0;                     /**< 非对角线元素 */
	          vec_val[col_idx[j]] = vec_val[col_idx[j]] - (*mat_data)(col_idx[j], index)*vec_val[index];
	          (*mat_data)(col_idx[j], index) = 0.0; /**< 对应的列元素 */
	        }
	      }
	    }
	  }

#endif
	  //  cout<<"constrain ok"<<endl;
}

/*!
 * @brief 完成电计算后处理
 * @param patch          输入参数, 网格片类, 表示网格片.
 * @param time           输入参数, 双精度浮点型, 表示当前时刻.
 * @param dt             输入参数, 双精度浮点型, 表示时间步长.
 * @param component_name 输入参数, 字符串, 表示数值构件的名称.
 *
 */
void PatchStrategy::buildE_MatrixOnPatch(hier::Patch<NDIM>& patch, const double time,
                    const double dt, const string& component_name)
{
	/// 取出本地PatchGeometry.
	  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
	      patch.getPatchGeometry();
	  /// 取出本地PatchTopology.
	  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
	      patch.getPatchTopology();
	  /// 取出本地Patch的结点坐标数组.

	  /// 获取单元对应实体编号数组对象 update #1  @7
	  tbox::Pointer<pdat::CellData<NDIM, int> > entityid_data =
	      patch.getPatchData(d_EntityIdOfCell_id);
	  tbox::Pointer<pdat::CellData<NDIM, int> > materialid_data =
              patch.getPatchData(material_num_id);
	  //取前一步的温度分布
	  tbox::Pointer<pdat::NodeData<NDIM, double> > T_data =
	      patch.getPatchData(th_Told_id);//

	  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
	      patch_geo->getNodeCoordinates();
	  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > mat_data =
	      patch.getPatchData(E_matrix_id);

	  /// 获取本地单元(边, 面)周围结点的索引关系.
	  tbox::Array<int> can_extent, can_indices;
	  patch_top->getCellAdjacencyNodes(can_extent, can_indices);
	  int num_cells = patch.getNumberOfCells(1);

	  /// 取出自由度映射信息
	  int* dof_map = d_dof_info_th->getDOFMapping(patch, hier::EntityUtilities::NODE);
	  /// 遍历单元
	  for (int i = 0; i < num_cells; ++i) {
		  int cell=i;
	    int n_vertex = can_extent[i + 1] - can_extent[i];
	    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
	    tbox::Array<int> mapping(n_vertex);
	    tbox::Array<double> T_val(n_vertex);
	    /// 取出单元结点坐标，以及填写映射值。//can_indices 节点编号
	    for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
	      mapping[i1] = dof_map[can_indices[j]];
	      T_val[i1]=T_data->getPointer()[can_indices[j]];
	      for (int k = 0; k < NDIM; ++k) {
	        vertex[i1][k] = (*node_coord)(k, can_indices[j]);
	      }
	    }
	    /// 初始化单元矩阵
	    tbox::Pointer<tbox::Matrix<double> > ele_mat = new tbox::Matrix<double>();
	    ele_mat->resize(n_vertex, n_vertex);
	    for (int i = 0; i < n_vertex; ++i) {
	      for (int j = 0; j < n_vertex; ++j) {
	        (*ele_mat)(i, j) = 0.0;
	      }
	    }

	    /// 取出单元对象.
	    tbox::Pointer<BaseElement<NDIM> > ele =
	        d_element_manager->getElement(d_element_type[0]);

	    /// 计算单元矩阵
	    ele->buildE_ElementMatrix(vertex, dt, time, ele_mat,(*materialid_data)(0,cell),T_val);

	    /// 累加单元矩阵，这里用户根据映射将矩阵挨个添加。
	    int row = 0, col = 0;
	    for (int i = 0; i < n_vertex; ++i) {
	      row = mapping[i];
	      for (int j = 0; j < n_vertex; ++j) {
	        col = mapping[j];
	        mat_data->addMatrixValue(row, col, (*ele_mat)(i, j));
	      }
	    }
	  }
	  /// 矩阵组装。
	  mat_data->assemble();
	 // cout<<"matrix ok"<<endl;
}

void PatchStrategy::Electric_PostProcesing(hier::Patch<NDIM>& patch, const double time,
                    const double dt, const string& component_name)
{
	 /// 取出本地PatchGeometry.
	  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
	      patch.getPatchGeometry();
	  /// 取出本地PatchTopology.
	  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
	      patch.getPatchTopology();
	  /// 取出本地Patch的结点坐标数组.
	  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
	      patch_geo->getNodeCoordinates();
	  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_data =
	      patch.getPatchData(E_solution_id);

	  tbox::Pointer<pdat::CellData<NDIM, double> > electric_Data =
	      patch.getPatchData(E_mag_id);

	  tbox::Pointer<pdat::NodeData<NDIM, double> > plot_data =
	      patch.getPatchData(E_plot_id);

	  int* dof_map = d_dof_info_th->getDOFMapping(patch, hier::EntityUtilities::NODE);
      double* vec_pointer = vec_data->getPointer(0);

	  /// 获取单元周围结点的索引关系.
	  tbox::Array<int> can_extent, can_indices;
	  patch_top->getCellAdjacencyNodes(can_extent, can_indices);

      int num_cells = patch.getNumberOfCells(0);

	  for (int i = 0; i < num_cells; ++i) {
		int cell_n=i;
	    int n_vertex = can_extent[i + 1] - can_extent[i];
	    int num_dof = NDIM * n_vertex;

	    /**< 该单元的结点坐标及自由度映射 */
	    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
	    tbox::Array<int> node_mapping(num_dof);
	    tbox::Array<double> vol_data(n_vertex);

	    for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
	    	node_mapping[i1] = dof_map[can_indices[j]];
            vol_data[i1] = vec_pointer[node_mapping[i1]];
	      for (int k = 0; k < NDIM; ++k) {
	        vertex[i1][k] = (*node_coord)(k, can_indices[j]);
	      }
	    }

	    /// 取出积分器对象.
	    tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
	        IntegratorManager<NDIM>::getManager();
	    tbox::Pointer<BaseIntegrator<NDIM> > integrator =
	        integrator_manager->getIntegrator("LinearTetrahedron");

	    /// 取出形函数对象.
	    tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
	        ShapeFunctionManager<NDIM>::getManager();
	    tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
	        shape_manager->getShapeFunction("LinearTetrahedron");

	    /// 取出材料.
	    tbox::Pointer<MaterialManager<NDIM> > material_manager =
	        MaterialManager<NDIM>::getManager();

	    /// 取出自由度数目.
	    int n_dof = shape_func->getNumberOfDof();

	    /// 取出积分点数目.
	    int num_quad_pnts = integrator->getNumberOfQuadraturePoints();

	    /// 取出积分点.
	    tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
	        integrator->getQuadraturePoints(vertex);

	    /// 取出积分点的积分权重.
	    tbox::Array<double> weight = integrator->getQuadratureWeights();

	    /// 取出基函数在积分点的值和梯度值.
	    tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
	        shape_func->gradient(vertex, quad_pnt);

	    double Ex=0;
	    double Ey=0;
	    double Ez=0;

	    //求解电场的三个分量
        // Vinta Yin-Da Wang
        // 直接算体心数据
	    for(int ii=0;ii<n_dof;ii++)
	    {
	    	for(int ll=0;ll<num_quad_pnts;ll++)
	    	{
	  	      Ex+=-vol_data[ii]* bas_grad[ll][ii][0]*weight[ll];
	  	      Ey+=-vol_data[ii]* bas_grad[ll][ii][1]*weight[ll];
	  	      Ez+=-vol_data[ii]* bas_grad[ll][ii][2]*weight[ll];

	    	}
	    	//cout<<vec_data->getPointer()[ii]<<endl;
	    }
	    //求电场模的平方
	    (*electric_Data)(0,cell_n)=Ex*Ex+Ey*Ey+Ez*Ez;
	  }

	  int step=time/dt;
	  if(step==1)
	  {
		  ofstream file_out;
		  file_out.open("Emag");
		  for(int ii=0;ii<num_cells;ii++)
			  file_out<<(*electric_Data)(0,ii)<<'\n';
		  file_out.close();
	  }


      int local_num_nodes = patch.getNumberOfNodes();
      for(int i=0;i<local_num_nodes;i++){
          int mapping_i = dof_map[i];
          (*plot_data)(0,i)=vec_pointer[mapping_i];

      }


}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PatchStrategy::registerPlotData(
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer) {
  javis_writer->registerPlotQuantity("plot", "VECTOR", d_plot_id);
  javis_writer->registerPlotQuantity("stress", "COMPOSITE", d_stress_id);
  //update #4 2017-04-21
  javis_writer->registerPlotQuantity("displacement","SCALAR",d_displacement_id);
  //update #2
  javis_writer->registerPlotQuantity("vonMises stress","SCALAR",d_von_mises_id);
  javis_writer->registerPlotQuantity("recovered vonMises stress","SCALAR",d_improved_von_mises_id);
  //update #8
  javis_writer->registerPlotQuantity("temperature_plot","SCALAR",th_plot_id);
  //update #9
  javis_writer->registerPlotQuantity("voltage_plot","SCALAR",E_plot_id);
  javis_writer->registerPlotQuantity("Emag","SCALAR",E_mag_id);
  javis_writer->registerPlotQuantity("mat_plot","SCALAR",material_num_id);
}

/*************************************************************************
 *  从输入数据库读入数据.
 ************************************************************************/
void PatchStrategy::getFromInput(tbox::Pointer<tbox::Database> db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!db.isNull());
#endif

  if (db->keyExists("element_type")) {
    d_element_type = db->getStringArray("element_type");
  } else {
    TBOX_ERROR(d_object_name << ": "
                             << " No key `element_type' found in data."
                             << endl);
  }

  if (db->keyExists("element_marks")) {
    d_element_marks = db->getIntegerArray("element_marks");
  } else {
    TBOX_ERROR(d_object_name << ": "
                             << " No key `element_marks' found in data."
                             << endl);
  }

  if (db->keyExists("constraint_types")) {
    d_constraint_types = db->getStringArray("constraint_types");
  } else {
    TBOX_ERROR(d_object_name << ": "
                             << " No key `constraint_types' found in data."
                             << endl);
  }

  if (db->keyExists("constraint_marks")) {
    d_constraint_marks = db->getIntegerArray("constraint_marks");
  } else {
    TBOX_ERROR(d_object_name << ": "
                             << " No key `constraint_marks' found in data."
                             << endl);
  }

  if (db->keyExists("load_types")) {
    d_load_types = db->getStringArray("load_types");
  } else {
    TBOX_ERROR(d_object_name << ": "
                             << " No key `load_types' found in data." << endl);
  }

  if (db->keyExists("load_marks")) {
    d_load_marks = db->getIntegerArray("load_marks");
  } else {
    TBOX_ERROR(d_object_name << ": "
                             << " No key `load_marks' found in data." << endl);
  }
  if (db->keyExists("convection_boundary")) {
    d_convection_boundary = db->getIntegerArray("convection_boundary");
  } else {
    TBOX_ERROR(d_object_name << ": "
                             << " No key `convection_boundary' found in data." << endl);
  }
  if (db->keyExists("improved_stress")) {
    d_improved_stress = db->getIntegerArray("improved_stress");
  } else {
    TBOX_ERROR(d_object_name << ": "
                             << " No key `improved_stress' found in data." << endl);
  }
}

/*************************************************************************
 *  从重启动数据库读取数据.
 ************************************************************************/
void PatchStrategy::getFromRestart(tbox::Pointer<tbox::Database> db) {
  getFromInput(db);
  tbox::Pointer<tbox::Database> root_db =
      tbox::RestartManager::getManager()->getRootDatabase();

  tbox::Pointer<tbox::Database> sub_db = root_db->getDatabase(d_object_name);
  d_dof_info->getFromDatabase(sub_db);
}

/**
***********************************************************************
*  2026-01-04
*  Yin-Da Wang
*  使用新的方法寻找点处于哪个四面体里
************************************************************************/
void PatchStrategy::QueryFieldAtPoints(hier::Patch<NDIM>& patch, const string& input_filename){
    // 1. 获取几何和拓扑信息
    tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo = patch.getPatchGeometry();
    tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord = patch_geo->getNodeCoordinates();
    tbox::Pointer<hier::PatchTopology<NDIM> > patch_top = patch.getPatchTopology();

    // 获取单元-节点邻接关系
    tbox::Array<int> cell_node_ext, cell_node_idx;
    patch_top->getCellAdjacencyNodes(cell_node_ext, cell_node_idx);

    int num_nodes = patch.getNumberOfNodes(0); // 本地节点数
    int num_cells = patch.getNumberOfCells(0); // 本地单元数

    // 2. 获取物理场数据 (例如温度 T_data)
    // 获取温度场数据，用来进行点插值
    tbox::Pointer<pdat::NodeData<NDIM, double> > T_data = patch.getPatchData(th_plot_id);


    // 3. 准备CGAL点数据
    /// 建树流程
    /// 清空并重新填充类成员 Point_on_patch，防止多次调用数据累积
    /// 将patch上所有的结点都保留下来备用
    Point_on_patch.clear();
    Point_on_patch.reserve(num_nodes);
    for(int nn = 0; nn < num_nodes; nn++){
        Point_on_patch.push_back(J_point((*node_coord)(0,nn), (*node_coord)(1,nn), (*node_coord)(2,nn)));
    }

    // 4. 构建四面体集合
    std::vector<J_tetrahedron> tetrahedrons;
    tetrahedrons.reserve(num_cells);

    for(int i = 0; i < num_cells; ++i) {
        // 获取单元的节点索引
        int n0 = cell_node_idx[cell_node_ext[i] + 0];
        int n1 = cell_node_idx[cell_node_ext[i] + 1];
        int n2 = cell_node_idx[cell_node_ext[i] + 2];
        int n3 = cell_node_idx[cell_node_ext[i] + 3];

        tetrahedrons.push_back(J_tetrahedron(
                                   &Point_on_patch[n0], &Point_on_patch[n1], &Point_on_patch[n2], &Point_on_patch[n3], i
                                   ));
        // 5. 构建 AABB Tree
        Tet_Tree tree(tetrahedrons.begin(), tetrahedrons.end());
        tree.build(); // 用于加速查找

        // 6. 读取输入文件并查询
        ifstream infile;
        infile.open(input_filename.c_str(), ios::in);
        if (!infile) {
            TBOX_ERROR("Query input file not found: " << input_filename << endl);
        }



        // 输出文件
        stringstream ss;
        ss << "Query_Result_Patch_" << patch.getIndex() << ".dat";
        ofstream outfile;
        outfile.open(ss.str().c_str(), ios::out);
        outfile << std::fixed << std::setprecision(12);

        string line;
        while(getline(infile, line)){
            stringstream buf(line);
            double q_num,q_x, q_y, q_z;
            // 假设输入文件格式为: x y z
            if(!(buf >> q_num >> q_x >> q_y >> q_z)) continue;
            double query_coord[3] = {q_x, q_y, q_z};
            CGAL_K::Point_3 query_pt(q_x, q_y, q_z);
            // 7. 使用 Tree 查找包含该点的所有四面体（候选）
            // 这里使用 any_intersected_primitive 或者 all_intersected_primitives
            // 注意：Point 与 Tetrahedron 的 intersection 在 CGAL Kernel 中是支持的
            std::vector<Tet_Tree::Primitive_id> candidates;
            tree.all_intersected_primitives(query_pt, std::back_inserter(candidates));
            bool found = false;
            for(size_t k = 0; k < candidates.size(); ++k){
                int cell_id = candidates[k]->id;
            }
        }
    }
}
