//
// 文件名:     LinearTetrahedron.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     单元计算类
// 类别　:     %Internal File% ( Don't delete this line )
//
#ifndef included_LinearTetrahedron
#define included_LinearTetrahedron
#include "BaseElement.h"
#include "BaseIntegrator.h"
#include "BaseShapeFunction.h"
#include "DoubleVector.h"

#include "IntegratorManager.h"
#include "ShapeFunctionManager.h"
#include "ElementManager.h"

using namespace JAUMIN;
class LinearTetrahedron : public BaseElement<NDIM> {
public:
  /**
   * @brief 构造函数.
   *
   * @param name 输入参数, 单元名字.
   *
   */
  LinearTetrahedron(const string& name);

  /**
   * @brief 析构函数.
   *
   */
  ~LinearTetrahedron();

  /**
   * @brief 计算单元质量矩阵.
   *
   * @param ele_mat        输出参数, 二维数组, 单元质量矩阵.
   * @param integrator     输入参数, 指针, 指向积分器.
   * @param shape_func     输入参数, 指针, 指向形函数.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   */
  virtual void buildStiffElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > coord, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat);

  /**
   * @brief 计算单元右端项.
   *
   * @param ele_rhs        输出参数, 数组, 单元右端项.
   * @param integrator     输入参数, 指针, 指向积分器.
   * @param shape_func     输入参数, 指针, 指向形函数.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   */
  virtual void buildElementRHS(tbox::Array<hier::DoubleVector<NDIM> > coord,
                               const double dt, const double time,
                               tbox::Pointer<tbox::Vector<double> > ele_vec);

private:
  /*!@brief 对象名.  */
  string d_object_name;
};
#endif
