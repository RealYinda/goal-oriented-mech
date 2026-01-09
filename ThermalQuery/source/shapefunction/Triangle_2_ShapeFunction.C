//
// 文件名:     Triangle_2_ShapeFunction.C
// 软件包:     JAUMIN fem
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2012-02-16 12:55:07 +0800 (四, 2012-02-16) $
// 描述　:     有限元基函数类.
// 类别　:     %Internal File% ( Don't delete this line )
//

#define AREA(a, b, c)                      \
  (((b)[0] - (a)[0]) * ((c)[1] - (a)[1]) - \
   ((b)[1] - (a)[1]) * ((c)[0] - (a)[0]))

#define dlambda0_dx ((real_vertex[1][1] - real_vertex[2][1]) / area)
#define dlambda0_dy ((real_vertex[2][0] - real_vertex[1][0]) / area)
#define dlambda1_dx ((real_vertex[2][1] - real_vertex[0][1]) / area)
#define dlambda1_dy ((real_vertex[0][0] - real_vertex[2][0]) / area)
#define dlambda2_dx ((real_vertex[0][1] - real_vertex[1][1]) / area)
#define dlambda2_dy ((real_vertex[1][0] - real_vertex[0][0]) / area)

#include "DoubleVector.h"
#include "Triangle_2_ShapeFunction.h"

/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
Triangle_2_ShapeFunction::Triangle_2_ShapeFunction(const string& name)
    : BaseShapeFunction<NDIM>(name) {
  d_num_dof = 6;
}

/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
Triangle_2_ShapeFunction::~Triangle_2_ShapeFunction() {}

void Triangle_2_ShapeFunction::get_lambda(
    const hier::DoubleVector<NDIM> pnt,
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex, double* lambda,
    double* area) {
  area[0] = AREA(real_vertex[0], real_vertex[1], real_vertex[2]);
  lambda[0] = AREA(pnt, real_vertex[1], real_vertex[2]) / area[0];
  lambda[1] = AREA(pnt, real_vertex[2], real_vertex[0]) / area[0];
  lambda[2] = AREA(pnt, real_vertex[0], real_vertex[1]) / area[0];
}

/*****************************************************************************
 * 获取某一个基函数在某一点的值.
 *****************************************************************************/
double Triangle_2_ShapeFunction::value(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  double tmp_val;
  double lambda[3], area;
  get_lambda(pnt, real_vertex, lambda, &area);
  if (j < 3)
    tmp_val = lambda[j] * (2 * lambda[j] - 1.0);
  else
    tmp_val = 4.0 * lambda[(j + 1) % 3] * lambda[(j + 2) % 3];
  return tmp_val;
}

/*****************************************************************************
 * 获取某一个基函数在某些点的值.
 *****************************************************************************/
tbox::Array<double> Triangle_2_ShapeFunction::value(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const tbox::Array<hier::DoubleVector<NDIM> >& pnts) {
  int pnt_size = pnts.getSize();
  tbox::Array<double> tmp_val(pnt_size);
  for (int i = 0; i < pnt_size; ++i) {
    tmp_val[i] = value(j, real_vertex, pnts[i]);
  }
  return tmp_val;
}

/*****************************************************************************
 * 获取基函数在某一点的值.
 *****************************************************************************/
tbox::Array<double> Triangle_2_ShapeFunction::value(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  tbox::Array<double> tmp_val(6);
  double lambda[3], area;
  get_lambda(pnt, real_vertex, lambda, &area);
  tmp_val[0] = lambda[0] * (2 * lambda[0] - 1.0);
  tmp_val[1] = lambda[1] * (2 * lambda[1] - 1.0);
  tmp_val[2] = lambda[2] * (2 * lambda[2] - 1.0);
  tmp_val[3] = 4.0 * lambda[1] * lambda[2];
  tmp_val[4] = 4.0 * lambda[0] * lambda[2];
  tmp_val[5] = 4.0 * lambda[0] * lambda[1];
  return tmp_val;
}

/*****************************************************************************
 * 获取基函数在某些点的值.
 *****************************************************************************/
tbox::Array<tbox::Array<double> > Triangle_2_ShapeFunction::value(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const tbox::Array<hier::DoubleVector<NDIM> >& pnts) {
  int pnt_size = pnts.getSize();
  tbox::Array<tbox::Array<double> > tmp_val(pnt_size);
  for (int i = 0; i < pnt_size; ++i) {
    tbox::Array<double> val_i = value(real_vertex, pnts[i]);
    tmp_val[i] = val_i;
  }
  return tmp_val;
}

/*****************************************************************************
 * 获取某一个基函数在某一点的梯度值.
 *****************************************************************************/
tbox::Array<double> Triangle_2_ShapeFunction::gradient(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  tbox::Array<double> tmp_grad(2);
  return tmp_grad;  ///?????????????? weibo-bug
}

/*****************************************************************************
 * 获取某一个基函数在某些点的梯度值.
 *****************************************************************************/
tbox::Array<tbox::Array<double> > Triangle_2_ShapeFunction::gradient(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const tbox::Array<hier::DoubleVector<NDIM> >& pnts) {
  int pnt_size = pnts.getSize();
  tbox::Array<tbox::Array<double> > tmp_grad(pnt_size);
  for (int i = 0; i < pnt_size; ++i) {
    tbox::Array<double> grad = gradient(j, real_vertex, pnts[i]);
    tmp_grad[i] = grad;
  }
  return tmp_grad;
}

/*****************************************************************************
 * 获取基函数在某一点的梯度值.
 *****************************************************************************/
tbox::Array<tbox::Array<double> > Triangle_2_ShapeFunction::gradient(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  tbox::Array<tbox::Array<double> > tmp_grad(6);
  double lambda[3], area;
  get_lambda(pnt, real_vertex, lambda, &area);

  tbox::Array<double> val_0(2);
  val_0[0] = (4.0 * lambda[0] - 1.0) * dlambda0_dx;
  val_0[1] = (4.0 * lambda[0] - 1.0) * dlambda0_dy;
  tmp_grad[0] = val_0;

  tbox::Array<double> val_1(2);
  val_1[0] = (4.0 * lambda[1] - 1.0) * dlambda1_dx;
  val_1[1] = (4.0 * lambda[1] - 1.0) * dlambda1_dy;
  tmp_grad[1] = val_1;

  tbox::Array<double> val_2(2);
  val_2[0] = (4.0 * lambda[2] - 1.0) * dlambda2_dx;
  val_2[1] = (4.0 * lambda[2] - 1.0) * dlambda2_dy;
  tmp_grad[2] = val_2;

  tbox::Array<double> val_3(2);
  val_3[0] = 4.0 * (lambda[2] * dlambda1_dx + lambda[1] * dlambda2_dx);
  val_3[1] = 4.0 * (lambda[2] * dlambda1_dy + lambda[1] * dlambda2_dy);
  tmp_grad[3] = val_3;

  tbox::Array<double> val_4(2);
  val_4[0] = 4.0 * (lambda[0] * dlambda2_dx + lambda[2] * dlambda0_dx);
  val_4[1] = 4.0 * (lambda[0] * dlambda2_dy + lambda[2] * dlambda0_dy);
  tmp_grad[4] = val_4;

  tbox::Array<double> val_5(2);
  val_5[0] = 4.0 * (lambda[1] * dlambda0_dx + lambda[0] * dlambda1_dx);
  val_5[1] = 4.0 * (lambda[1] * dlambda0_dy + lambda[0] * dlambda1_dy);
  tmp_grad[5] = val_5;
  return tmp_grad;
}

/*****************************************************************************
 * 获取基函数在某些点的梯度值.
 *****************************************************************************/
tbox::Array<tbox::Array<tbox::Array<double> > >
Triangle_2_ShapeFunction::gradient(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const tbox::Array<hier::DoubleVector<NDIM> >& pnts) {
  int pnt_size = pnts.getSize();
  tbox::Array<tbox::Array<tbox::Array<double> > > tmp_grad(pnt_size);
  for (int i = 0; i < pnt_size; ++i) {
    tbox::Array<tbox::Array<double> > grad = gradient(real_vertex, pnts[i]);
    tmp_grad[i] = grad;
  }
  return tmp_grad;
}

/****************************************************************************
 * 获取自由度数目.
 ****************************************************************************/
int Triangle_2_ShapeFunction::getNumberOfDof() { return d_num_dof; }
