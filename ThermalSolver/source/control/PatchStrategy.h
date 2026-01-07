//
// 文件名:     PatchStrategy.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:18:29 2014 $
// 描述　:     网格片策略类派生类.
//
// 根新情况：
// update #0: 1处  @1：添加头文件
// update #1: 1处  @1：添加id声明
// update #4: 1处  @1：添加id声明
// update #5: 1处  @1：添加宏定义，已做备份用（包含3个宏）
// update #6: 4处  @1@2@3@4：添加id声明

#ifndef included_PatchStrategy
#define included_PatchStrategy

#include "StandardComponentPatchStrategy.h"
#include "DOFInfo.h"
#include "GridInfo.h"
#include "ElementManager.h"
#include "JaVisDataWriter.h"

//update #0
//自定义宏管理文件
#include "MacrosManager.h"

/// 使用CGAL里的功能
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>


//update #5
#if 0
//暂时用不到，如果MacrosManager.h完全没问题，这部分可以直接删去，目前做个备份
//update #5 计算三位坐标系中三角形面积
//at 2017-04-25 by tong
//code by zhu
#define AREA(a, b, c)                      \
  ((((b)[1] - (a)[1]) * ((c)[2] - (a)[2]) -((b)[2] - (a)[2]) * ((c)[1] - (a)[1]))*(((b)[1] - (a)[1]) * ((c)[2] - (a)[2]) -((b)[2] - (a)[2]) * ((c)[1] - (a)[1]))+\
  (((b)[2] - (a)[2]) * ((c)[0] - (a)[0]) -((b)[0] - (a)[0]) * ((c)[2] - (a)[2]))*(((b)[2] - (a)[2]) * ((c)[0] - (a)[0]) -((b)[0] - (a)[0]) * ((c)[2] - (a)[2]))+\
  (((b)[0] - (a)[0]) * ((c)[1] - (a)[1]) -((b)[1] - (a)[1]) * ((c)[0] - (a)[0]))*(((b)[0] - (a)[0]) * ((c)[1] - (a)[1]) -((b)[1] - (a)[1]) * ((c)[0] - (a)[0])))

//update #5 2017-04-26
#define NODE_LOAD 0  //1： 存在节点力载荷； 0： 不存在
#define FACE_LOAD 1  //1： 存在表面力载荷； 0： 不存在
#endif

using namespace JAUMIN;
class PatchStrategy : public algs::StandardComponentPatchStrategy<NDIM>,
                      tbox::Serializable {
public:
  /*! @brief 构造函数.
   * @param object_name          输入参数, 字符串, 表示对象名称.
   * @param input_db             输入参数，指针，指向输入数据库。
   * @param register_for_restart 输入参数, 逻辑型, 是否注册到重启动.
   */
  PatchStrategy(const string& object_name,
                tbox::Pointer<tbox::Database> input_db,
                bool register_for_restart);

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

   /*
   * @note
   *  vector是输入输出参数: 输入是已遍历网格片的归约结果,
   *  输出是输入值和当前网格片计算结果的归约值.
   */
  void reduceOnPatch(double* vector, int len, hier::Patch<NDIM>& patch,
                     const double time, const double dt,
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
   * @brief 支撑指定名称的数值构件, 在单个网格片上完成右端项组装.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void buildRHSOnPatch(hier::Patch<NDIM>& patch, const double time,
                       const double dt, const string& component_name);

  /**
   * @brief 注册模型变量, 该函数调用EntisyStrategy的注册函数,
   *完成应户变量的注册.
   *
   */
  void registerModelVariable();

  /*!
   * @brief 在单个网格片上设置载荷.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void applyLoad(hier::Patch<NDIM>& patch, const double time, const double dt,
                 const string& component_name);

  /*!
   * @brief 在单个网格片上设置约束.
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
   * @brief 在单个网格片上更新坐标.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void updateCoordinate(hier::Patch<NDIM>& patch, const double time,
                        const double dt, const string& component_name);

  /*!
   * @brief 在单个网格片上计算应力.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void computeStress(hier::Patch<NDIM>& patch, const double time,
                     const double dt, const string& component_name);

  /*!
   * @brief 支撑指定名称的有限元数值构件, 在单个网格片上完成数值计算.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param initial_time   输入参数, 逻辑型, 是否为初始时刻.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void computeOnPatch(hier::Patch<NDIM>& patch, const double time,
                      const double dt, const bool initial_time,
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
  double getPatchDt(hier::Patch<NDIM>& patch, const double time,
                    const bool initial_time, const int flag_last_dt,
                    const double last_dt, const string& component_name);

  /*!
   * @brief 将数据成员输出到重启动数据库.
   * @param db 输入参数, 指针, 指向重启动数据库.
   */
  void putToDatabase(tbox::Pointer<tbox::Database> db);

  /**
   * @brief 从重启动数据库中读如参数。
   *
   * @param db 输入参数，指针，指向输入数据库。
   */
  void getFromRestart(tbox::Pointer<tbox::Database> db);

  /**
   * @brief 注册可视化数据.
   */
  void registerPlotData(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer);

  /**
   * @brief 获取矩阵id
   *
   * @return 整型，矩阵id
   */
  int getMatrixID() { return d_matrix_id; }

  /**
   * @brief 获取向量id
   *
   * @return 整型，向量id
   */
  int getRHSID() { return d_rhs_id; }

  /**
   * @brief 获取解向量id
   *
   * @return 整型，解向量id
   */
  int getSolutionID() { return d_solution_id; }

  ///////////////////////////////////////////////update #8///////////////////////////////////////////////////

  /*!
   * @brief 支撑指定名称的数值构件, 在单个网格片上完成矩阵组装.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void buildTh_MatrixOnPatch(hier::Patch<NDIM>& patch, const double time,
                          const double dt, const string& component_name);

  /*!
   * @brief 支撑指定名称的数值构件, 在单个网格片上完成右端项组装.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void buildTh_RHSOnPatch(hier::Patch<NDIM>& patch, const double time,
                       const double dt, const string& component_name);

  /*!
   * @brief 在单个网格片上设置温度载荷：热源及第二类边界.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void applyTh_Load(hier::Patch<NDIM>& patch, const double time, const double dt,
                 const string& component_name);

  /*!
   * @brief 在单个网格片上设置温度约束.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void applyTh_Constraint(hier::Patch<NDIM>& patch, const double time,
                       const double dt, const string& component_name);

  /*!
   * @brief 完成热计算后处理
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   *
   */
  void Thermal_PostProcesing(hier::Patch<NDIM>& patch, const double time,
                      const double dt, const string& component_name);

  void Thermal_max(double* max_T, int len, hier::Patch<NDIM>& patch, const double time,
                      const double dt, const string& component_name);
  void Stress_max(double* vector, int len, hier::Patch<NDIM>& patch, const double time,
                      const double dt, const string& component_name);

  /**
   * @brief 获取热求解矩阵id
   *
   * @return 整型，热求解矩阵id
   */
  int getTh_MatrixID() { return th_matrix_id; }

  /**
   * @brief 获取热求解向量id
   *
   * @return 整型，热求解向量id
   */
  int getTh_RHSID() { return th_rhs_id; }

  /**
   * @brief 获取温度解向量id
   *
   * @return 整型，温度解向量id
   */
  int getTh_SolutionID() { return th_solution_id; }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////update #9///////////////////////////////////////////////////
  //电计算相关函数
    /*!
     * @brief 支撑指定名称的数值构件, 在单个网格片上完成矩阵组装.
     *
     * @param patch          输入参数, 网格片类, 表示网格片.
     * @param time           输入参数, 双精度浮点型, 表示当前时刻.
     * @param dt             输入参数, 双精度浮点型, 表示时间步长.
     * @param component_name 输入参数, 字符串, 表示数值构件的名称.
     *
     */
    void buildE_MatrixOnPatch(hier::Patch<NDIM>& patch, const double time,
                            const double dt, const string& component_name);

    /*!
     * @brief 支撑指定名称的数值构件, 在单个网格片上完成右端项组装.
     *
     * @param patch          输入参数, 网格片类, 表示网格片.
     * @param time           输入参数, 双精度浮点型, 表示当前时刻.
     * @param dt             输入参数, 双精度浮点型, 表示时间步长.
     * @param component_name 输入参数, 字符串, 表示数值构件的名称.
     *
     */
    void buildE_RHSOnPatch(hier::Patch<NDIM>& patch, const double time,
                         const double dt, const string& component_name);

    /*!
     * @brief 在单个网格片上设置电压约束.
     *
     * @param patch          输入参数, 网格片类, 表示网格片.
     * @param time           输入参数, 双精度浮点型, 表示当前时刻.
     * @param dt             输入参数, 双精度浮点型, 表示时间步长.
     * @param component_name 输入参数, 字符串, 表示数值构件的名称.
     *
     */
    void applyE_Constraint(hier::Patch<NDIM>& patch, const double time,
                         const double dt, const string& component_name);

    /*!
     * @brief 完成电计算后处理
     * @param patch          输入参数, 网格片类, 表示网格片.
     * @param time           输入参数, 双精度浮点型, 表示当前时刻.
     * @param dt             输入参数, 双精度浮点型, 表示时间步长.
     * @param component_name 输入参数, 字符串, 表示数值构件的名称.
     *
     */
    void Electric_PostProcesing(hier::Patch<NDIM>& patch, const double time,
                        const double dt, const string& component_name);

    /**
     * @brief StressRecovery: 完成应力在体上的恢复
     * @param patch
     * @param time
     * @param dt
     * @param component_name
     */
    void StressRecovery(hier::Patch<NDIM>& patch, const double time,
                        const double dt, const string& component_name);
    /**
     * @brief PostprocessStress: 完成恢复应力的后处理
     * @param patch
     * @param time
     * @param dt
     * @param component_name
     */
    void PostprocessStress(hier::Patch<NDIM>& patch, const double time,
                        const double dt, const string& component_name);

    /**
     * @brief Dataexplorer: 数据处理
     * @param patch
     * @param time
     * @param dt
     * @param component_name
     */
    void Dataexplorer(hier::Patch<NDIM>& patch, const double time,
                        const double dt, const string& component_name);

    /**
     * @brief ThermalPostprocess: 热问题的后处理
     * @param patch
     * @param time
     * @param dt
     * @param component_name
     */
    void ThermalPostprocess(hier::Patch<NDIM>& patch, const double time,
                        const double dt, const string& component_name);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief 获取热求解矩阵id
     *
     * @return 整型，热求解矩阵id
     */
    int getE_MatrixID() { return E_matrix_id; }

    /**
     * @brief 获取热求解向量id
     *
     * @return 整型，热求解向量id
     */
    int getE_RHSID() { return E_rhs_id; }

    /**
     * @brief 获取温度解向量id
     *
     * @return 整型，温度解向量id
     */
    int getE_SolutionID() { return E_solution_id; }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////


  /**
   * @brief 获取自由度信息
   *
   * @return 指针，指向自由度信息
   */
  tbox::Pointer<solv::DOFInfo<NDIM> > getDOFInfo() { return d_dof_info; }


  /// 支持CGAL的数据结构
  typedef CGAL::Simple_cartesian<double> CGAL_K;
  // custom point type
  struct J_point {
      double m_x;
      double m_y;
      double m_z;

      J_point(const double x,
              const double y,
              const double z)
          : m_x(x), m_y(y), m_z(z) {}
  };
  struct J_triangle {
      J_point *m_pa;
      J_point *m_pb;
      J_point *m_pc;
      int id;

      J_triangle(J_point *pa,
                 J_point *pb,
                 J_point *pc,
                 int id_in)
          : m_pa(pa), m_pb(pb), m_pc(pc), id(id_in) {}
  };
  /// 迭代器，存储J_三角形的原型
  typedef std::vector<J_triangle>::const_iterator Iterator;
  struct J_triangle_primitive {
  public:

      // this is the type of data that the queries returns. For this example
      // we imagine that, for some reasons, we do not want to store the iterators
      // of the vector, but raw pointers. This is to show that the Id type
      // does not have to be the same as the one of the input parameter of the
      // constructor.
      typedef const J_triangle* Id;

      // CGAL types returned
      typedef CGAL_K::Point_3    Point; // CGAL 3D point type
      typedef CGAL_K::Triangle_3 Datum; // CGAL 3D triangle type

  private:
      Id m_pt; // this is what the AABB tree stores internally

  public:
      J_triangle_primitive() {} // default constructor needed

      // the following constructor is the one that receives the iterators from the
      // iterator range given as input to the AABB_tree
      J_triangle_primitive(Iterator it)
          : m_pt(&(*it)) {}

      const Id& id() const { return m_pt; }

      // utility function to convert a custom
      // point type to CGAL point type.
      Point convert(const J_point *p) const
      {
          return Point(p->m_x,p->m_y,p->m_z);
      }

      // on the fly conversion from the internal data to the CGAL types
      Datum datum() const
      {
          return Datum(convert(m_pt->m_pa),
                       convert(m_pt->m_pb),
                       convert(m_pt->m_pc));
      }

      // returns a reference point which must be on the primitive
      Point reference_point() const
      { return convert(m_pt->m_pa); }
  };
  typedef CGAL::AABB_traits<CGAL_K, J_triangle_primitive> J_AABB_traits;
  typedef CGAL::AABB_tree<J_AABB_traits> Tree;

  /// 2026-01-04 更新
  ///  Yin-Da Wang
  ///  定义四面体结构，存储四个顶点指针和单元ID
  struct J_tetrahedron {
      J_point *m_p0;
      J_point *m_p1;
      J_point *m_p2;
      J_point *m_p3;
      int id; // 单元在Patch上的编号

      J_tetrahedron(J_point *p0, J_point *p1, J_point *p2, J_point *p3, int id_in)
          : m_p0(p0), m_p1(p1), m_p2(p2), m_p3(p3), id(id_in) {}
  };
  /// 定义四面体图元，用于 AABB Tree
  struct J_tetrahedron_primitive{
  public:
      typedef const J_tetrahedron* Id;
      typedef CGAL_K::Point_3    Point;
      typedef CGAL_K::Tetrahedron_3 Datum;
  private:
      Id m_pt;
  public:
      J_tetrahedron_primitive() {}
      J_tetrahedron_primitive(std::vector<J_tetrahedron>::const_iterator it) : m_pt(&(*it)) {}

      const Id& id() const { return m_pt; }

      Point convert(const J_point *p) const {
          return Point(p->m_x, p->m_y, p->m_z);
      }
      // 返回 CGAL 的四面体对象，用于构建包围盒
      Datum datum() const {
          return Datum(convert(m_pt->m_p0),
                       convert(m_pt->m_p1),
                       convert(m_pt->m_p2),
                       convert(m_pt->m_p3));
      }

      // 参考点
      Point reference_point() const { return convert(m_pt->m_p0); }
  };

  ///  定义针对四面体的 Tree 类型
  typedef CGAL::AABB_traits<CGAL_K, J_tetrahedron_primitive> J_Tet_Traits;
  typedef CGAL::AABB_tree<J_Tet_Traits> Tet_Tree;

  /// 使用Gemini提供的新方法查找点
  void QueryFieldAtPoints(hier::Patch<NDIM>& patch, const string& input_filename);


  /// 网格片上的CGAL类点对应的vector
  std::vector< J_point > Point_on_patch;
  bool PatchPointWeightInCell(const double* pointcoord,
                         double* localnodecoord, double* weight){
      tbox::Matrix<double> CellMat(NDIM+1);
      double candidate_weight[NDIM+1];
      /// coordinfo 应当是4行3列，在矩阵求解时是3行四列转置排列
      double (*coord_info)[NDIM] = (double(*)[NDIM])localnodecoord;
//      /// Test
//      double tmpdata[4][3];
//      for(int row = 0; row < 4; row ++){
//          for(int col = 0; col < 3; col ++){
//              tmpdata[row][col] = (coord_info)[row][col] ;
//          }
//      }
      for(int row =0; row < NDIM; row++){
          for(int col = 0; col < NDIM+1; col ++){
              (CellMat)(row, col) = coord_info[col][row];
          }
      }
      (CellMat)(NDIM, 0) = 1;(CellMat)(NDIM, 1) = 1;
      (CellMat)(NDIM, 2) = 1;(CellMat)(NDIM, 3) = 1;
      tbox::Matrix<double> PLUinvMat(NDIM+1);
      PLUinvMat = CellMat.getInverse();
//      /// Test
//      double matdata[4][4];
//      double invdata[4][4];
//      for(int row = 0; row < 4; row ++){
//          for(int col = 0; col < 4; col ++){
//              matdata[row][col] = (CellMat)(row, col) ;
//              invdata[row][col] = (PLUinvMat)(row, col) ;
//          }
//      }
      tbox::Vector<double> QuadVec(4);
      (QuadVec)[0] = pointcoord[0];(QuadVec)[1] = pointcoord[1];
      (QuadVec)[2] = pointcoord[2];(QuadVec)[3] = 1;
      tbox::Vector<double> QuadSol(4);
      QuadSol = PLUinvMat * QuadVec;
      for(int vdim = 0; vdim < NDIM+1; vdim ++)
          candidate_weight[vdim] = (QuadSol)[vdim];
      if(abs(candidate_weight[0]-0.5)<0.500001
              && abs(candidate_weight[1]-0.5)<0.500001
              && abs(candidate_weight[2]-0.5)<0.500001
              && abs(candidate_weight[3]-0.5)<0.500001){
          for(int vdim = 0; vdim < NDIM+1; vdim ++){
              weight[vdim] = candidate_weight[vdim];
          }
          return true;
      }
      else{
          return false;
      }

  }

private:
  /*!@brief 从输入数据库读入数据.  */
  void getFromInput(tbox::Pointer<tbox::Database> db);

  bool d_registered_for_restart;

  /*!@brief 对象名.  */
  string d_object_name;
  tbox::Pointer<solv::DOFInfo<NDIM> > d_dof_info;

  //update #8 @1
  tbox::Pointer<solv::DOFInfo<NDIM> > d_dof_info_th;

  tbox::Array<string> d_element_type; /**< 单元类型 */
  tbox::Array<int> d_element_marks;   /**< 单元的标示 */

  tbox::Array<string> d_constraint_types; /**< 约束类型数组 */
  tbox::Array<int> d_constraint_marks; /**< 约束作用的点集标识数组 */

  tbox::Array<string> d_load_types; /**< 载荷类型数组 */
  tbox::Array<int> d_load_marks;    /**< 载荷作用的点集标识数组 */
  /// 对流边界条件对应的边界编号
  tbox::Array<int> d_convection_boundary;
  /// 需要应力恢复的体编号
  tbox::Array<int> d_improved_stress;

  /// 有限元计算相关的管理器(单元, 积分, 形函数, 材料, 约束, 载荷).
  tbox::Pointer<ElementManager<NDIM> > d_element_manager;

  /// 需要读取的坐标点信息
  string d_file_name_query;

  /// 输入数据库指针.
  tbox::Pointer<tbox::Database> d_input_db;

  /// 变量的索引.
  int d_solution_id;
  int d_matrix_id;

  int d_rhs_id;
  int d_stress_id;
  int d_plot_id;

  int d_EntityIdOfCell_id;//update #1
  int d_displacement_id;//update #4
  int d_von_mises_id;//update #2

  //update #6
  int d_solution_old_id;
  int d_solution_older_id;
  int d_rhs_old_id;
  int d_rhs_older_id;

  //update #8 @2
  int th_solution_id;
  int th_matrix_id;
  int th_rhs_id;
  int th_plot_id;
  int th_Told_id;
  int th_Tolder_id;

  //update #9 电计算相关变量
  int E_solution_id;
  int E_matrix_id;
  int E_rhs_id;
  int E_plot_id;
  int E_mag_id;//单元电场模值

  int material_num_id;

  /// 存储以它为中心的domain的插值系数
  int d_improved_coefficient_id;
  /// 平均后的插值系数
  int d_output_coefficient_id;
  /// 存储应力准则（恢复后）的画图量
  int d_improved_von_mises_id;
  /// 工具区域的id
  int d_tool_domain_id;
  /// 0:tool region中心,
  /// 1:不在tool region中心,但是邻接单元在
  /// 2:邻接单元也不在tool region中心
  /// 戳啦,2不用设置
  int d_contained_domain_id;


  /// For TetQuad
  int d_Cell_volume_id;
  int d_Cell_jacobian_id;

};
#endif
