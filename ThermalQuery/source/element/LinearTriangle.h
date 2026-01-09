//
// 文件名:     LinearTriangle.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     单元计算类
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "BaseElement.h"
#include "BaseIntegrator.h"
#include "BaseShapeFunction.h"
#include "DoubleVector.h"
#include "Vector.h"
#include "Matrix.h"

#include "IntegratorManager.h"
#include "ShapeFunctionManager.h"
#include "ElementManager.h"

using namespace JAUMIN;
class LinearTriangle : public BaseElement<NDIM> {
public:
  /**
   * @brief 构造函数.
   *
   * @param name 输入参数, 单元名字.
   *
   */
  LinearTriangle(const string& name);

  /**
   * @brief 析构函数.
   *
   */
  ~LinearTriangle();

  /**
   * @brief 计算单元刚度矩阵.
   *
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_mat        输出参数, 指针, 指向矩阵.
   */
  virtual void buildStiffElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > coord, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat);

  /**
   * @brief 计算单元右端项.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_vec        输出参数, 指针, 指向单元向量.
   */
  virtual void buildElementRHS(tbox::Array<hier::DoubleVector<NDIM> > coord,
                               const double dt, const double time,
                               tbox::Pointer<tbox::Vector<double> > ele_vec);
};
