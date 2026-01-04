//
// 文件名:     BaseMaterial.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:34:43 2014 $
// 描述　:     材料基类
//

#ifndef included_BaseMaterial
#define included_BaseMaterial

#include "Pointer.h"
#include "Array.h"
using namespace JAUMIN;
using namespace std;

/**
 * @brief 该类是力学计算的材料基类, 提供各种材料的共性接口.
 *
 */
template <int DIM>
class BaseMaterial {
public:
  /**
   * @brief 构造函数.
   *
   * @param name 输入参数, 字符串, 对象名称.
   */
  BaseMaterial(const string& name);

  /**
   * @brief 析构函数.
   *
   */
  virtual ~BaseMaterial();

  /**
   * @brief 获取杨氏模量.
   *
   * @return 双精度浮点型, 杨氏模量.
   */
  virtual double getYoungModulus();

  /**
   * @brief 获取波松比.
   *
   *
   * @return 双精度浮点型, 波松比.
   */
  virtual double getPossionRatio();

  /**
   * @brief 获取密度.
   *
   *
   * @return 双精度浮点型, 密度.
   */
  virtual double getDensity();

  /**
   * @brief 获取模量矩阵.
   *
   * @return 二维双精度浮点型数组, 模量矩阵数组.
   */
  virtual tbox::Array<tbox::Array<double> > getModuli();

  /**
   * @brief 计算弹性刚度矩阵.
   *
   */
  virtual void computeMuduli();

  /**
   * @brief 获取对象的名字.
   *
   * @return 字符串, 对象的名字.
   */
  virtual const string& getName();

  /**
   * @brief 判断给定的名字是否为该对象的名字.
   *
   * @param name 输入参数, 字符串, 约束名字.
   *
   * @return 布尔型, 是否一致.
   */
  virtual bool isMatching(const string& name);

private:
  string d_object_name;
};
#endif
