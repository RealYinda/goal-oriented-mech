//
// 文件名:     BaseElement.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     有限元单元基类
//
////////////////////////////////////////////////////////////////
//update #3:向单元刚度、质量、阻尼举证三个函数添加参量，标识计算单元的材料属性
//2017-04-21 by tong
//

#include "JAUMIN_config.h"
#include "Utilities.h"
#include "BaseIntegrator.h"
#include "DoubleVector.h"
#include "BaseShapeFunction.h"
#include "BaseElement.h"

#include <string>
#include <list>

#ifndef NULL
#define NULL (0)
#endif

/************************************************************************
 * 构造函数
 *************************************************************************/
template <int DIM>
BaseElement<DIM>::BaseElement(const string& name) {
  d_object_name = name;
}

/************************************************************************
 * 析构函数
 *************************************************************************/
template <int DIM>
BaseElement<DIM>::~BaseElement() {}

/************************************************************************
 * 建立单元刚度矩阵
 *
 * update #3
 *************************************************************************/
template <int DIM>
void BaseElement<DIM>::buildStiffElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id,double T_val) {
  TBOX_WARNING(
      "Function buildStiffElementMatrix is not implemented, please make sure "
      "your code is right! ");
}

template <int DIM>
void BaseElement<DIM>::buildRecoveryMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id) {
  TBOX_WARNING(
      "Function buildRecoveryMatrix is not implemented, please make sure "
      "your code is right! ");
}

/************************************************************************
 * 建立单元质量矩阵
 *
 * update #3
 *************************************************************************/
template <int DIM>
void BaseElement<DIM>::buildMassElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id,double T_val) {
  TBOX_WARNING(
      "Function buildMassElementMatrix is not implemented, please make sure "
      "your code is right! ");
}

/************************************************************************
 * 建立单元阻尼矩阵
 *
 * update #3
 *************************************************************************/
template <int DIM>
void BaseElement<DIM>::buildDumpElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id,double T_val) {
  TBOX_WARNING(
      "Function buildDumpElementMatrix is not implemented, please make sure "
      "your code is right! ");
}


/*****************************************************************************
 * @brief 建立单元矩阵.
 *
 ****************************************************************************/
template <int DIM>
void BaseElement<DIM>::buildElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, tbox::Array<double> T_val){
	  TBOX_WARNING(
	      "Function buidElementMatrix is not implemented, please make sure "
	      "your code is right! ");
	}

/************************************************************************
 * 建立单元右端项
 *************************************************************************/
template <int DIM>
void BaseElement<DIM>::buildElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
	  NewmarkData *d_newmark,int entity_id, tbox::Array<double> T_val,
    tbox::Array<double> Tolder_val) {
  TBOX_WARNING(
      "Function buildElementRHS is not implemented, please make sure your code "
      "is right! ");
}

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
template <int DIM>
void BaseElement<DIM>::buildTh_ElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, tbox::Array<double> T_val)
{
  TBOX_WARNING(
      "Function buildTh_ElementMatrix is not implemented, please make sure your code "
      "is right! ");
}

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
template <int DIM>
void BaseElement<DIM>::buildTh_ElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
	  int entity_id, tbox::Array<double> T_val, double e_ThermalSource)
{
  TBOX_WARNING(
      "Function buildTh_ElementRHS is not implemented, please make sure your code "
      "is right! ");
}

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
template <int DIM>
void BaseElement<DIM>::buildE_ElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, tbox::Array<double> T_val)
{
	TBOX_WARNING(
	"Function buildE_ElementMatrix is not implemented, please make sure your code "
	"is right! ");
}

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
template <int DIM>
void BaseElement<DIM>::buildE_ElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
	  int entity_id, tbox::Array<double> T_val)
{
	TBOX_WARNING(
	"Function buildE_ElementRHS is not implemented, please make sure your code "
	"is right! ");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////


/*****************************************************************************
 * 判断给定名字是否以本对象名字匹配.
 *****************************************************************************/
template <int DIM>
bool BaseElement<DIM>::isMatching(const string& name) {
  return (d_object_name == name);
}

/*****************************************************************************
 * 获取对象名字.
 *****************************************************************************/
template <int DIM>
const string& BaseElement<DIM>::getName() {
  return d_object_name;
}

template <int DIM>
const int BaseElement<DIM>::getProblemDim() {
  TBOX_WARNING(
      "Function getProblemDim is not implemented, please make sure your code "
      "is right! ");
  return 0;
}

template <int DIM>
tbox::Array<int> BaseElement<DIM>::getNumberOfDofOnEntity() {
  TBOX_WARNING(
      "Function getNumberOfDofOnEntity is not implemented, please make sure "
      "your code is right! ");
  return 0;
}
template class BaseElement<NDIM>;
