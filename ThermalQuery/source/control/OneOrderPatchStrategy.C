//
// 文件名:     OneOrderPatchStrategy.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:45:28 2014 $
// 描述　:     网格片策略类派生类实现.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "VectorVariable.h"
#include "CSRMatrixVariable.h"
#include "VectorData.h"
#include "CSRMatrixData.h"
#include "ElementManager.h"
#include "ShapeFunctionManager.h"
#include "IntegratorManager.h"
#include "BaseElement.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "Patch.h"
#include "OneOrderPatchStrategy.h"
#include "Matrix.h"
#include "Vector.h"
#include "NodeVariable.h"
#include "NodeData.h"
#include "CellVariable.h"
#include "CellData.h"
#include "LinearTriangle.h"
#include "LinearTetrahedron.h"
#include "IntegratorManager.h"
#include "ShapeFunctionManager.h"
#include "ElementManager.h"
#include "ShapeFunctionManager.h"
#include "IntegratorManager.h"
#include "TriangleIntegrator.h"
#include "Triangle_1_ShapeFunction.h"
#include "TetrahedronShapeFunction.h"
#include "TetrahedronIntegrator.h"
#include "LinearTetrahedron.h"

#include "TetQuad.h"
#include "TetGeom.h"
/*************************************************************************
 *
 * 构造函数.
 *
 *************************************************************************/
OneOrderPatchStrategy::OneOrderPatchStrategy(
    const string& object_name, bool is_from_restart,
    tbox::Pointer<tbox::Database> input_db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!object_name.empty());
#endif

  /// (形函数, 积分器, 单元)管理器
  tbox::Pointer<ElementManager<NDIM> > ele_manager =
      ElementManager<NDIM>::getManager();
  tbox::Pointer<ShapeFunctionManager<NDIM> > func_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();

  /// 创建单元.
#if (NDIM == 2)  // 三角形
  tbox::Pointer<LinearTriangle> linear_ele =
      new LinearTriangle("LinearTriangle");
  tbox::Pointer<TriangleIntegrator> linear_int =
      new TriangleIntegrator(1, "Triangle");
  tbox::Pointer<Triangle_1_ShapeFunction> linear_sha =
      new Triangle_1_ShapeFunction("Triangle");
#else  // 四面体
  tbox::Pointer<LinearTetrahedron> linear_ele =
      new LinearTetrahedron("LinearTetrahedron");
  tbox::Pointer<TetrahedronIntegrator> linear_int =
      new TetrahedronIntegrator(1, "Tetrahedron");
  tbox::Pointer<TetrahedronShapeFunction> linear_sha =
      new TetrahedronShapeFunction("Tetrahedron");
#endif
  /// 将单元添加到单元管理器.
  ele_manager->addElement(linear_ele);
  integrator_manager->addIntegrator(linear_int);
  func_manager->addShapeFunction(linear_sha);
  d_object_name = object_name;

  // 读取从输入文件或重启动文件读入数据.
  if (is_from_restart) {
    // getFromRestart();
  } else {
    getFromInput(input_db);
  }

  d_object_name = object_name;
  registerModelVariable();
}

/*************************************************************************
 *
 * 构造函数.
 *
 ************************************************************************/
OneOrderPatchStrategy::~OneOrderPatchStrategy() {}

/*************************************************************************
 *
 * 注册变量和数据片.
 *
 ************************************************************************/
void OneOrderPatchStrategy::registerModelVariable() {
  /// 获取有限元变量的数据库.
  hier::VariableDatabase<NDIM>* variable_db =
      hier::VariableDatabase<NDIM>::getDatabase();

  /// 创建自由度信息，参数信息（四个逻辑型表示点，边，面，体上自由度是否非零。
  d_dof_info = new solv::DOFInfo<NDIM>(true, false, false, false);

  /// 变量上下文.
  tbox::Pointer<hier::VariableContext> current =
      variable_db->getContext("CURRENT");

  /// 向量变量， 解向量
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > solution =
      new pdat::VectorVariable<NDIM, double>("solution", d_dof_info);
  /// 向量变量， 右端项
  tbox::Pointer<pdat::VectorVariable<NDIM, double> > rhs =
      new pdat::VectorVariable<NDIM, double>("rhs", d_dof_info);

  /// 矩阵变量， 矩阵
  tbox::Pointer<pdat::CSRMatrixVariable<NDIM, double> > matrix =
      new pdat::CSRMatrixVariable<NDIM, double>("matrix", d_dof_info);

  /// 结点型变量，可视化数据片
  tbox::Pointer<pdat::NodeVariable<NDIM, double> > plotdata =
      new pdat::NodeVariable<NDIM, double>("plot", 1);
  /// 结点型变量，误差数据片
  tbox::Pointer<pdat::NodeVariable<NDIM, double> > error =
      new pdat::NodeVariable<NDIM, double>("error", 1);

  /// 将变量及上下文注册到有限元变量数据库.
  d_solution_id = variable_db->registerVariableAndContext(solution, current, 1);
  d_rhs_id = variable_db->registerVariableAndContext(rhs, current, 1);
  d_matrix_id = variable_db->registerVariableAndContext(matrix, current, 1);
  d_plot_id = variable_db->registerVariableAndContext(plotdata, current);
  d_error_id = variable_db->registerVariableAndContext(error, current);

  tbox::Pointer<pdat::CellVariable<NDIM, double> > Cell_jacobian =
      new pdat::CellVariable<NDIM, double>("Cell_jacobian",(NDIM+1)*(NDIM+1),1);//雅可比矩阵
  tbox::Pointer<pdat::CellVariable<NDIM, double> > Cell_volume =
      new pdat::CellVariable<NDIM, double>("Cell_volume",1,1);//单元体积
  d_Cell_jacobian_id=variable_db->registerVariableAndContext(Cell_jacobian, current,1);
  d_Cell_volume_id=variable_db->registerVariableAndContext(Cell_volume, current,1);
}

/*************************************************************************
 *
 *  初始化指定的积分构件.
 *
 *  注册待填充的数据片或待调度内存空间的数据片到积分构件.
 ************************************************************************/
void OneOrderPatchStrategy::initializeComponent(
    algs::IntegratorComponent<NDIM>* component) const {
  const string& component_name = component->getName();
  if (component_name == "INIT") {  // 右端项构件.
    component->registerInitPatchData(d_error_id);
    component->registerInitPatchData(d_plot_id);

    /// 将自由度信息中的若干数据片注册到初始化构件。
    d_dof_info->registerToInitComponent(component);
    component->registerInitPatchData(d_Cell_jacobian_id);
    component->registerInitPatchData(d_Cell_volume_id);
  } else if (component_name == "RHS") {    // 数值构件，计算右端项.
  } else if (component_name == "MAT") {    // 数值构件，计算矩阵.
  } else if (component_name == "ALLOC") {  // 内存构件.
    /// 将矩阵向量注册到内存构件
    component->registerPatchData(d_matrix_id);
    component->registerPatchData(d_solution_id);
    component->registerPatchData(d_rhs_id);
  } else if (component_name == "POST") {  // 数值构件，后处理计算
  } else if (component_name == "CONS") {  // 数值构件，加载约束
  } else if (component_name == "RED") {   // 规约构件，规约误差量
  } else {
    TBOX_ERROR("\n::initializeComponent() : component "
               << component_name << " is not matched. " << endl);
  }
}

/*************************************************************************
 *  初始化数据片（支持初值构件）.
 ************************************************************************/
void OneOrderPatchStrategy::initializePatchData(hier::Patch<NDIM>& patch,
                                                const double time,
                                                const bool initial_time,
                                                const string& component_name) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(component_name == "INIT");
#endif
  NULL_USE(time); /**< 初始化中没有用到time */

  if (initial_time) {
    tbox::Pointer<pdat::NodeData<NDIM, double> > plot_data =
        patch.getPatchData(d_plot_id);
    tbox::Pointer<pdat::NodeData<NDIM, double> > error_data =
        patch.getPatchData(d_error_id);

    int num_nodes = patch.getNumberOfNodes(1);
    int local_num_nodes = patch.getNumberOfNodes();
    JAUMIN::appu::TetGeom tetrahedron(patch);
    /// 取出点上的自由度分布和映射数据片的头指针
    int* dis_ptr =
        d_dof_info->getDOFDistribution(patch, hier::EntityUtilities::NODE);

    /// 对自由度信息的分布信息初始化
    for (int i = 0; i < num_nodes; ++i) {
      dis_ptr[i] = 1;
    }
    /// 建立自由度映射信息
    d_dof_info->buildPatchDOFMapping(patch);

    /// 对普通数据片初始化
    for (int i = 0; i < local_num_nodes; ++i) {
      (*plot_data)(0, i) = 0;
      (*error_data)(0, i) = 0;
    }
  }
}

/*************************************************************************
 *  输出数据成员到重启动数据库.
 ************************************************************************/
void OneOrderPatchStrategy::putToDatabase(tbox::Pointer<tbox::Database> db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!db.isNull());
#endif
}

/*************************************************************************
 * 完成单个网格片上的数值计算（支持数值构件）.
 ************************************************************************/
void OneOrderPatchStrategy::computeOnPatch(hier::Patch<NDIM>& patch,
                                           const double time, const double dt,
                                           const bool initial_time,
                                           const string& component_name) {
  if (component_name == "MAT") { /**< 矩阵组装的数值构件 */
    buildMatrixOnPatch(patch, time, dt, component_name);
  } else if (component_name == "RHS") { /**< 右端项组装的数值构件 */
    buildRHSOnPatch(patch, time, dt, component_name);
  } else if (component_name == "POST") { /**< 后处理的数值构件 */
    postProcess(patch, time, dt, component_name);
  } else if (component_name == "CONS") { /**< 添加约束的数值构件 */
    applyConstraint(patch, time, dt, component_name);
  } else {
    TBOX_ERROR(" OneOrderPatchStrategy :: component name is error! ");
  }
}

/*************************************************************************
 *  加载约束到代数系统.
 ************************************************************************/
void OneOrderPatchStrategy::applyConstraint(hier::Patch<NDIM>& patch,
                                            const double time, const double dt,
                                            const string& component_name) {
  /// 取出本地PatchGeometry.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();
  /// 获取网格片矩阵对象
  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > mat_data =
      patch.getPatchData(d_matrix_id);
  /// 获取网格片向量对象
  tbox::Pointer<pdat::VectorData<NDIM, double> > vec_data =
      patch.getPatchData(d_rhs_id);
  /// 结点数目
  int num_nodes = patch.getNumberOfNodes();

  /// 自由度信息中的映射信息。
  int* dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
  /// 取出矩阵向量的核心数据结构的指针
  int* row_start = mat_data->getRowStartPointer();
  int* col_idx = mat_data->getColumnIndicesPointer();
  double* mat_val = mat_data->getValuePointer();
  double* vec_val = vec_data->getPointer();

  /// 边界条件处理.
  if (patch_geo->hasEntitySet(1, hier::EntityUtilities::NODE)) {
    // 获取指定编号和类型的集合包含的网格实体的索引。
    const tbox::Array<int>& entity_idx = patch_geo->getEntityIndicesInSet(
          1, hier::EntityUtilities::NODE, num_nodes);
    int size = entity_idx.getSize();
    /// 对角化 1 法处理约束。
    for (int i = 0; i < size; ++i) {
      int index = dof_map[entity_idx[i]];
      vec_val[index] = 0.0; /**< 载荷项设为 0 */
      for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
        if (col_idx[j] == index) {
          mat_val[j] = 1.0; /**< 对角线元素 */
        } else {
          mat_val[j] = 0.0;                     /**< 非对角线元素 */
          (*mat_data)(col_idx[j], index) = 0.0; /**< 对应的列元素 */
        }
      }
    }
  }
}

/*************************************************************************
 *  建立网格片上的矩阵和右端项.
 ************************************************************************/
void OneOrderPatchStrategy::buildRHSOnPatch(hier::Patch<NDIM>& patch,
                                            const double time, const double dt,
                                            const string& component_name) {
  /// (形函数, 积分器, 单元)管理器
  tbox::Pointer<ElementManager<NDIM> > ele_manager =
      ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele =
      ele_manager->getElement(d_element_type);

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
      patch.getPatchData(d_rhs_id);
  /// 取出单元数目
  int num_cells = patch.getNumberOfCells(1);

  /// 自由度映射信息
  int* dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  /**
   * 遍历单元，将单元右端项加载到右端项向量中。
   */
  for (int i = 0; i < num_cells; ++i) {
    int n_vertex = can_extent[i + 1] - can_extent[i];
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<int> mapping(n_vertex);

    /// 下面的循环做两件事情：1. 建立自由度映射数组；2.取出结点坐标。
    for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
      mapping[i1] = dof_map[can_indices[j]];
      for (int k = 0; k < NDIM; ++k) {
        vertex[i1][k] = (*node_coord)(k, can_indices[j]);
      }
    }

    /// 声明vec，用来存储单元右端项。
    tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
    ele_vec->resize(n_vertex);
    for (int i = 0; i < n_vertex; ++i) {
      (*ele_vec)[i] = 0.0;
    }

    /// 计算单元右端项，并组装。
    ele->buildElementRHS(vertex, time, dt, ele_vec);
    for (int i2 = 0; i2 < n_vertex; ++i2)
      vec_data->addVectorValue(mapping[i2], (*ele_vec)[i2]);
  }
}

/*************************************************************************
 *  建立网格片上的矩阵和右端项.
 ************************************************************************/
void OneOrderPatchStrategy::buildMatrixOnPatch(hier::Patch<NDIM>& patch,
                                               const double time,
                                               const double dt,
                                               const string& component_name) {
  /// (形函数, 积分器, 单元)管理器
  tbox::Pointer<ElementManager<NDIM> > ele_manager =
      ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele =
      ele_manager->getElement(d_element_type);

  /// 取出本地PatchGeometry.
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();
  /// 取出本地PatchTopology.
  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
      patch.getPatchTopology();
  /// 取出本地Patch的结点坐标数组.
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch_geo->getNodeCoordinates();
  tbox::Pointer<pdat::CSRMatrixData<NDIM, double> > mat_data =
      patch.getPatchData(d_matrix_id);
  /// 获取本地单元(边, 面)周围结点的索引关系.
  tbox::Array<int> can_extent, can_indices;
  patch_top->getCellAdjacencyNodes(can_extent, can_indices);
  int num_cells = patch.getNumberOfCells(1);

  /// 取出自由度映射信息
  int* dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
  /// 遍历单元
  for (int i = 0; i < num_cells; ++i) {
    int n_vertex = can_extent[i + 1] - can_extent[i];
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<int> mapping(n_vertex);
    /// 取出单元结点坐标，以及填写映射值。
    for (int i1 = 0, j = can_extent[i]; i1 < n_vertex; ++i1, ++j) {
      mapping[i1] = dof_map[can_indices[j]];
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
    /// 计算单元矩阵
    ele->buildStiffElementMatrix(vertex, time, dt, ele_mat);

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

/*************************************************************************
 *  后处理计算：1. 将vec写入到可视化的结点量数据片；2. 计算误差.
 ************************************************************************/
void OneOrderPatchStrategy::postProcess(hier::Patch<NDIM>& patch,
                                        const double time, const double dt,
                                        const string& component_name) {
  /**
   * Vinta Yin-Da Wang
   * 2023-11-14
   * 输入热应力中的应力网格，输出采样点的序号和网格
   */
  int num_cell = patch.getNumberOfCells(0);
  /// The Quad rule 2023-11-14
  appu::TetQuad quad(patch, patch.getPatchData(d_Cell_volume_id),
                     patch.getPatchData(d_Cell_jacobian_id));
  appu::TetQuad::Quad ** Quad_patch = quad.GetQuadTable();
  appu::TetQuad::Quad *Quad_order = Quad_patch[2];
  std::ofstream intpFile;




}

/*************************************************************************
 *  规约构件对应的网格片函数，计算网格片上的L1，L2误差.
 ************************************************************************/
void OneOrderPatchStrategy::reduceOnPatch(double* vector, int len,
                                          hier::Patch<NDIM>& patch,
                                          const double time, const double dt,
                                          const string& component_name) {
  double local_l2_error = 0;
  double local_l1_error = 0;
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_error_data =
      patch.getPatchData(d_error_id);
  tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
      patch.getPatchGeometry();
  /// 取出本地PatchTopology.
  tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
      patch.getPatchTopology();
  /// 取出本地Patch的结点坐标数组.
  tbox::Array<int> can_extent, can_indices;
  patch_top->getCellAdjacencyNodes(can_extent, can_indices);
  int num_cells = patch.getNumberOfCells();
  /// 取出本地Patch的结点坐标数组.
  tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord =
      patch_geo->getNodeCoordinates();

  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  for (int i = 0; i < num_cells; ++i) {
#if (NDIM == 2)
    tbox::Pointer<BaseIntegrator<NDIM> > integrator =
        integrator_manager->getIntegrator("Triangle");

    tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
        shape_manager->getShapeFunction("Triangle");
#else
    tbox::Pointer<BaseIntegrator<NDIM> > integrator =
        integrator_manager->getIntegrator("Tetrahedron");

    tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
        shape_manager->getShapeFunction("Tetrahedron");
#endif
    /// 取出单元上自由度数目.
    int n_dof = shape_func->getNumberOfDof();
#if (NDIM > 2)
    tbox::Array<hier::DoubleVector<NDIM> > vertex(4);
#else
    tbox::Array<hier::DoubleVector<NDIM> > vertex(3);
#endif
    tbox::Array<double> dof_value(n_dof);
    for (int i1 = 0, j = can_extent[i]; j < can_extent[i + 1]; ++j, ++i1) {
      dof_value[i1] = node_error_data->getPointer()[can_indices[j]];
      for (int k = 0; k < NDIM; ++k) {
        vertex[i1][k] = (*node_coord)(k, can_indices[j]);
      }
    }

    /// 取出积分点数目.
    int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
    /// 取出模板单元的面积.
    double volume = integrator->getElementVolume();
    /// 取出积分点.
    tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
        integrator->getQuadraturePoints(vertex);
    /// 取出jacobian矩阵行列式.
    double jac = integrator->getLocal2GlobalJacobian(vertex);
    /// 取出积分点的积分权重.
    tbox::Array<double> weight = integrator->getQuadratureWeights();

    /// 取出基函数在积分点的值.
    tbox::Array<tbox::Array<double> > bas_val =
        shape_func->value(vertex, quad_pnt);
    for (int l = 0; l < num_quad_pnts; ++l) {
      double u_value = 0.0;
      for (int i1 = 0; i1 < n_dof; ++i1) {
        u_value += dof_value[i1] * bas_val[l][i1];
      }

      double JxW = volume * jac * weight[l];
      local_l2_error += JxW * u_value * u_value;
      local_l1_error += JxW * fabs(u_value);
    }
  }
  vector[0] = local_l2_error;
  vector[1] = local_l1_error;
}

/*************************************************************************
 *  注册可视化数据.
 ************************************************************************/
void OneOrderPatchStrategy::registerPlotData(
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer) {
  javis_data_writer->registerPlotQuantity("plot", "SCALAR", d_plot_id);
}

/*************************************************************************
 *  从输入数据库读入数据.
 ************************************************************************/
void OneOrderPatchStrategy::getFromInput(tbox::Pointer<tbox::Database> db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!db.isNull());
#endif
  /// 从输入数据库中读入有限元计算的若干参数.
  if (db->keyExists("element_type")) {
    d_element_type = db->getString("element_type");
  } else {
    TBOX_ERROR(d_object_name << ": "
               << " No key `element_type' found in data."
               << endl);
  }
  if (db->keyExists("query_file_name")) {
    d_file_name_query = db->getString("query_file_name");
  } else {
    TBOX_ERROR(d_object_name << ": "
               << " No key `query_file_name' found in data." << endl);
  }
}
