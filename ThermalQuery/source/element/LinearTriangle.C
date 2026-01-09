//
// 文件名:     LinearTriangle.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     单元计算类实现
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "LinearTriangle.h"
using namespace JAUMIN;
LinearTriangle::LinearTriangle(const string& name) : BaseElement<NDIM>(name) {}
LinearTriangle::~LinearTriangle() {}

void LinearTriangle::buildStiffElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat) {
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Triangle");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Triangle");

  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, quad_pnt);

  /// 计算单元刚度矩阵.
  for (int i = 0; i < n_dof; ++i) {
    for (int j = 0; j < n_dof; ++j) {
      for (int l = 0; l < num_quad_pnts; ++l) {
        double JxW = volume * jac * weight[l];

        (*ele_mat)(i, j) += JxW * (bas_grad[l][i][0] * bas_grad[l][j][0] +
                                   bas_grad[l][i][1] * bas_grad[l][j][1]);
      }
    }
  }
}

void LinearTriangle::buildElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_rhs) {
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Triangle");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Triangle");

  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  /// 计算单元右端项.
  for (int i = 0; i < n_dof; ++i) {
    for (int l = 0; l < num_quad_pnts; ++l) {
      double JxW = volume * jac * weight[l];
      const double PI = 4.0 * atan(1.0);
      double f_val = 104.0 * PI * PI * sin(2.0 * PI * quad_pnt[l][0]) *
                     sin(10.0 * PI * quad_pnt[l][1]);
      (*ele_rhs)[i] += JxW * f_val * bas_val[l][i];
    }
  }
}
