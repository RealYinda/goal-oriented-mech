//
// 文件名:     BaseElement.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     有限元单元类基类
// 类别　:     %Internal File% ( Don't delete this line )
//
////////////////////////////////////////////////////////////////
//update #3:向单元刚度、质量、阻尼举证三个函数添加参量，标识计算单元的材料属性
//2017-04-21 by tong
//

#ifndef included_fem_BaseElement
#define included_fem_BaseElement

#include <vector>
#include "Array.h"
#include "Pointer.h"
#include "Matrix.h"
#include "Vector.h"

//update #6
#include"MacrosManager.h"

using namespace std;
using namespace JAUMIN;
namespace JAUMIN {

namespace hier {
template <int DIM>
class DoubleVector;
}
}
template <int DIM>
class BaseIntegrator;
template <int DIM>
class BaseShapeFunction;

template <int DIM>
class BaseElement {
public:
  /**
   * @brief 构造函数.
   *
   * @param name 输入参数, 字符串, 对象名称.
   */
  BaseElement(const string& name);

  /**
   * @brief 析构函数.
   */
  virtual ~BaseElement();

  /**
   * @brief 计算单元刚度矩阵.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_mat        输出参数, 指针, 指向单元刚度矩阵.
   *
   *update #3
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
   * @param ele_mat        输出参数, 指针, 指向单元质量矩阵.
   *
   * update #3
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
   * @param ele_mat        输出参数, 指针, 指向单元阻尼矩阵.
   *
   * update #3
   */
  virtual void buildDumpElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id,double T_val);

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
   * @param ele_vec        输出参数, 指针, 指向单元右端项.
   */
  virtual void buildElementRHS(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
	NewmarkData *d_newmark, int entity_id, tbox::Array<double> T_val,
      tbox::Array<double> Tolder_val);//


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
   * @brief 判断给定的名字是否为该对象的名字.
   *
   * @param name 输入参数, 字符串, 约束名字.
   *
   * @return 布尔型, 是否一致.
   */
  virtual bool isMatching(const string& name);

  /**
   * @brief 获取对象的名字.
   *
   * @return 字符串, 对象的名字.
   */
  virtual const string& getName();

  /**
   * @brief 获取单元求解问题维数.
   *
   * @return 整型, 单元求解问题的维数.
   */
  virtual const int getProblemDim();

  /**
   * @brief 获取单元上每一个网格实体上自由度的数目.
   *
   * @return 整型, 每一个网格实体上自由度的数目.
   */
  virtual tbox::Array<int> getNumberOfDofOnEntity();

private:
  string d_object_name;
};
#endif
