//
// 文件名:     BaseMaterial.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:35:01 2014 $
// 描述　:     材料基类实现.
//

#include "JAUMIN_config.h"
#include "BaseMaterial.h"
#include "Utilities.h"

#include <string>
#include <list>

#ifndef NULL
#define NULL (0)
#endif

/************************************************************************
 * 构造函数
 ************************************************************************/
template <int DIM>
BaseMaterial<DIM>::BaseMaterial(const string& name) {
  d_object_name = name;
}

/************************************************************************
 * 析构函数
 ************************************************************************/
template <int DIM>
BaseMaterial<DIM>::~BaseMaterial() {}

/************************************************************************
 * 获取杨氏模量.
 ************************************************************************/
template <int DIM>
double BaseMaterial<DIM>::getYoungModulus() {
  TBOX_WARNING(
      "Function getYoungModulus is not implemented, please make sure your code "
      "is right! ");
  return 0.0;
}

/************************************************************************
 * 获取波松比.
 ************************************************************************/
template <int DIM>
double BaseMaterial<DIM>::getPossionRatio() {
  TBOX_WARNING(
      "Function getPossionRatio is not implemented, please make sure your code "
      "is right! ");
  return 0.0;
}

/************************************************************************
 * 获取密度.
 ************************************************************************/
template <int DIM>
double BaseMaterial<DIM>::getDensity() {
  TBOX_WARNING(
      "Function getDensity is not implemented, please make sure your code is "
      "right! ");
  return 0.0;
}

/************************************************************************
 * 获取模量矩阵.
 ************************************************************************/
template <int DIM>
tbox::Array<tbox::Array<double> > BaseMaterial<DIM>::getModuli() {
  tbox::Array<tbox::Array<double> > retval;
  TBOX_WARNING(
      "Function getModuli is not implemented, please make sure your code is "
      "right! ");
  return retval;
}

/************************************************************************
 * 计算模量矩阵.
 ************************************************************************/
template <int DIM>
void BaseMaterial<DIM>::computeMuduli() {
  TBOX_WARNING(
      "Function computeModuli is not implemented, please make sure your code "
      "is right! ");
}

/*****************************************************************************
 * 判断给定名字是否以本对象名字匹配.
 *****************************************************************************/
template <int DIM>
bool BaseMaterial<DIM>::isMatching(const string& name) {
  return (d_object_name == name);
}

/*****************************************************************************
 * 获取对象名字.
 *****************************************************************************/
template <int DIM>
const string& BaseMaterial<DIM>::getName() {
  return d_object_name;
}
template class BaseMaterial<NDIM>;
