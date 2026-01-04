//
// 文件名:     Material.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:34:03 2014 $
// 描述　:     线弹性材料
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_Material
#define included_Material

#include "MaterialManager.h"
#include "Pointer.h"
#include "Array.h"
#include "BaseMaterial.h"

#include "MacrosManager.h"
#include <cmath>
using namespace JAUMIN;

/**
 * @brief 该类是弹性力学计算的材料类, 提供材料的各种的属性.
 *
 */
class Material : public BaseMaterial<NDIM> {
public:
  /**
   * @brief 构造函数.
   *
   * @param name 材料名称.
   */
  Material(const string name);

  /**
   * @brief 析构函数.
   *
   */
  ~Material();

  /**
   * @brief 获取杨氏模量.
   *
   * @return 双精度浮点型, 杨氏模量.
   */
  double getYoungModulus(double T);

  /**
   * @brief 获取波松比.
   *
   *
   * @return 双精度浮点型, 波松比.
   */
  double getPossionRatio(double T);

  /**
   * @brief 获取密度.
   *
   *
   * @return 双精度浮点型, 密度.
   */
  double getDensity(double T);

  //update #6
  /**
   * @brief 获取阻尼系数.
   *
   *
   * @return 双精度浮点型, 阻尼系数.
   */
  double getDamping();

  ////////////////////////////////////update #7////////////////////////////////////

  /**
   * @brief 获取导热系数.
   *
   *
   * @return 双精度浮点型, 导热系数.
   */
  double getK(double T);

  /**
   * @brief 获取常压热容.
   *
   *
   * @return 双精度浮点型, 常压热容.
   */
  double getCp(double T);

  /**
   * @brief 获取热膨胀系数.
   *
   *
   * @return 双精度浮点型, 热膨胀系数.
   */
  double getAlpha(double T);

  /**
   * @brief 计算温变参数.
   *
   * 这里采用材料参数随温度变化关系为：四次多项式
   * a0+a1T+a2T^2+a3T^3+a4T^4
   * a0 为基础材料参数
   */
  double computeTimeParam(double T, Material_TParam_Type Param_Type);

  /////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////update #9/////////////////////////////////////
  //电学材料相关函数

  /**
   * @brief 获取相对介电常数.
   *
   *
   * @return 双精度浮点型, 热膨胀系数.
   */
  double getEpsilonr(double T);

  /**
   * @brief 获取相对磁导率.
   *
   *
   * @return 双精度浮点型, 热膨胀系数.
   */
  double getMur(double T);

  /**
   * @brief 获取电导率.
   *
   *
   * @return 双精度浮点型, 热膨胀系数.
   */
  double getSigma(double T);

  /////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief 获取弹性矩阵.
   *
   * @return 二维双精度浮点型数组, 弹性矩阵数组.
   */
  tbox::Array<tbox::Array<double> > getModuli(double T);

  /**
   * @brief 计算弹性刚度矩阵.
   *
   */
  void computeMuduli();

private:
  double d_young_modulus;
  double d_possion_ratio;
  double d_density;
  double d_damping;//update #6 阻尼系数

  //update #7 2017-05-06
  //热学相关材料
  double d_K;//导热系数
  double d_Cp;//常压热容
  double d_alpha;//热膨胀系数

  double Param_T[9][5];//材料温变系数

  //update #9 2017-05-06
  //电学学相关材料
  double d_mur;//相对磁导率
  double d_epsilonr;//相对介电常数
  double d_sigma;//电导率

  tbox::Array<tbox::Array<double> > d_moduli;
};

#endif
