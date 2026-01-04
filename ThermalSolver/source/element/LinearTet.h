//
// 文件名:     LinearTet.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:39:28 2014 $
// 描述　:     四面体线弹性单元.
// 类别　:     %Internal File% ( Don't delete this line )
//
////////////////////////////////////////////////////////////////
//update #3:向单元刚度、质量、阻尼举证三个函数添加参量，标识计算单元的材料属性
//2017-04-21 by tong
//

#ifndef include_LinearTet
#define include_LinearTet
#include "BaseElement.h"
#include "BaseIntegrator.h"
#include "BaseMaterial.h"
#include "BaseShapeFunction.h"
#include "BaseMaterial.h"
#include "DoubleVector.h"
#include "Vector.h"
#include "Matrix.h"

#include "IntegratorManager.h"
#include "ShapeFunctionManager.h"
#include "ElementManager.h"
#include "MaterialManager.h"
#include "Material.h"

//update
#include "Pointer.h"
#include "Array.h"
//自定义宏管理文件
#include "MacrosManager.h"

using namespace JAUMIN;

/**
 * @brief 线弹性单元类, 支撑各向同性线弹性材料计算的单元矩阵和右端项的计算.
 *
 */

class LinearTet : public BaseElement<NDIM> {
public:
  /**
   * @brief 构造函数.
   *
   * @param name  单元名字.
   */
  LinearTet(const string& name);

  /**
   * @brief 析构函数.
   */
  ~LinearTet();

  /**
   * @brief 计算单元刚度矩阵.
   *
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_mat        输出参数, 指针, 指向矩阵.
   *
   * update #3
   * @param entity_id      输入参数，整型，单元对应的实体编号
   */
  virtual void buildStiffElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id,double T_val);

  virtual void buildRecoveryMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id);

  /**
   * @brief 计算单元质量矩阵.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_mat        输出参数, 指针, 指向矩阵.
   *
   * update #3
   * @param entity_id      输入参数，整型，单元对应的实体编号
   */
  virtual void buildMassElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id,double T_val);

  /**
   * @brief 计算单元阻尼矩阵.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_mat        输出参数, 指针, 指向矩阵.
   *
   * update #3
   * @param entity_id      输入参数，整型，单元对应的实体编号
   */
  virtual void buildDumpElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id,double T_val);

  //update #6 计算单元矩阵P
  //P=M/(dt^2)+C/(2*dt)+beta*K
  /**
   * @brief 计算单元矩阵.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_mat        输出参数, 指针, 指向矩阵.
   *
   * update #3
   * @param entity_id      输入参数，整型，单元对应的实体编号
   */
  virtual void buildElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id,tbox::Array<double> T_val);

  /**
   * @brief 计算单元右端项.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_vec        输出参数, 指针, 指向单元向量.
   *
   * update #6
   * @param d_newmark      输入参数，结构体NewmarkData，Newmark-beta方法所需数据
   */
  virtual void buildElementRHS(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
	  NewmarkData *d_newmark,int entity_id,tbox::Array<double> T_val,
      tbox::Array<double> Tolder_val);

  /////////////////////////////////////////////////update #8//////////////////////////////////////////////

  /**
   * @brief 计算热求解单元矩阵.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_mat        输出参数, 指针, 指向矩阵.
   *
   * update #3
   * @param entity_id      输入参数，整型，单元对应的实体编号
   */
  virtual void buildTh_ElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id,tbox::Array<double> T_val);

  /**
   * @brief 计算热求解单元右端项.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_vec        输出参数, 指针, 指向单元向量.
   *
   * update #6
   * @param d_newmark      输入参数，结构体NewmarkData，Newmark-beta方法所需数据
   */
  virtual void buildTh_ElementRHS(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
	  int entity_id, tbox::Array<double> T_val, double e_ThermalSource);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////update #9//////////////////////////////////////////////

  /**
   * @brief 计算电求解单元矩阵.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_mat        输出参数, 指针, 指向矩阵.
   *
   * update #3
   * @param entity_id      输入参数，整型，单元对应的实体编号
   */
  virtual void buildE_ElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id,tbox::Array<double> T_val);

  /**
   * @brief 计算电求解单元右端项.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_vec        输出参数, 指针, 指向单元向量.
   *
   * update #6
   * @param d_newmark      输入参数，结构体NewmarkData，Newmark-beta方法所需数据
   */
  virtual void buildE_ElementRHS(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
	  int entity_id, tbox::Array<double> T_val);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief 获取单元求解问题维数.
   *
   * @return 整型, 单元求解问题的维数.
   */
  const int getProblemDim();

  /**
   * @brief 获取单元上每一个网格实体上自由度的数目.
   *
   * @return 整型, 每一个网格实体上自由度的数目.
   */
  tbox::Array<int> getNumberOfDofOnEntity();
};
#endif
