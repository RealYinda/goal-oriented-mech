//
// 文件名:     LinearTet.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:39:40 2014 $
// 描述　:     四面体线弹性单元.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "LinearTet.h"

/*************************************************************************
 * 构造函数.
 ************************************************************************/
LinearTet::LinearTet(const string& name) : BaseElement<NDIM>(name) {}

/*************************************************************************
 * 析构函数.
 ************************************************************************/
LinearTet::~LinearTet() {}


/// 恢复型算法的实现
void LinearTet::buildRecoveryMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id){
    double T_val = 300.;
    /// 取出积分器对象.
    tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
        IntegratorManager<NDIM>::getManager();
    tbox::Pointer<BaseIntegrator<NDIM> > integrator =
        integrator_manager->getIntegrator("LinearTetrahedron");

    /// 取出形函数对象.
    tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
        ShapeFunctionManager<NDIM>::getManager();
    tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
        shape_manager->getShapeFunction("LinearTetrahedron");

    /// 取出材料.
    tbox::Pointer<MaterialManager<NDIM> > material_manager =
        MaterialManager<NDIM>::getManager();

    //update #3:
    tbox::Pointer<Material> material =
                  material_manager->getMaterial("Gold");
        if(entity_id==1)
           material =	material_manager->getMaterial("Silicon");
      else if(entity_id==13)
           material =	material_manager->getMaterial("MoCu");
      else if(entity_id==2)
           material =	material_manager->getMaterial("Copper");
      else if(entity_id==4)
           material =	material_manager->getMaterial("SiO2");
      else if(entity_id==5)
           material =	material_manager->getMaterial("SiN");
      else if(entity_id==3)
           material =	material_manager->getMaterial("Gold");
      else if(entity_id==15)
           material =	material_manager->getMaterial("GaN");
      else if(entity_id==16)
           material =	material_manager->getMaterial("Al2O3");
      else if(entity_id==17)
           material =	material_manager->getMaterial("Alloy");
      else if(entity_id==6)
           material =	material_manager->getMaterial("Aluminum");
      else
           material =	material_manager->getMaterial("Gold");
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
        /// 取出基函数在积分点的值和梯度值.
        tbox::Array<tbox::Array<double> > bas_val =
                shape_func->value(real_vertex, quad_pnt);
        ele_mat->resize(n_dof * NDIM, 24);
        for (int i = 0; i < n_dof * NDIM; ++i) {
          for (int j = 0; j < 24; ++j) {
            (*ele_mat)(i, j) = 0.0;
          }
        }
        for(int quad = 0; quad <num_quad_pnts; quad++){
            double JxW = volume * jac * weight[quad];
            /// 本个积分点的坐标值
            double quad_coord[NDIM] = {0.,0.,0.};
            for(int i = 0; i < n_dof; i++){
                quad_coord[0] += real_vertex[i][0]*bas_val[quad][i];
                quad_coord[1] += real_vertex[i][1]*bas_val[quad][i];
                quad_coord[2] += real_vertex[i][2]*bas_val[quad][i];
            }
            for(int i = 0; i < n_dof; i++){
                /// 自由度i的x分量 146
                (*ele_mat)(NDIM * i + 0, 0 ) += bas_grad[quad][i][0]*quad_coord[0]*JxW;
                (*ele_mat)(NDIM * i + 0, 1 ) += bas_grad[quad][i][0]*quad_coord[1]*JxW;
                (*ele_mat)(NDIM * i + 0, 2 ) += bas_grad[quad][i][0]*quad_coord[2]*JxW;
                (*ele_mat)(NDIM * i + 0, 3 ) += bas_grad[quad][i][0]*1.*JxW;
                (*ele_mat)(NDIM * i + 0, 12 ) += bas_grad[quad][i][1]*quad_coord[0]*JxW;
                (*ele_mat)(NDIM * i + 0, 13 ) += bas_grad[quad][i][1]*quad_coord[1]*JxW;
                (*ele_mat)(NDIM * i + 0, 14 ) += bas_grad[quad][i][1]*quad_coord[2]*JxW;
                (*ele_mat)(NDIM * i + 0, 15 ) += bas_grad[quad][i][1]*1.*JxW;
                (*ele_mat)(NDIM * i + 0, 20 ) += bas_grad[quad][i][2]*quad_coord[0]*JxW;
                (*ele_mat)(NDIM * i + 0, 21 ) += bas_grad[quad][i][2]*quad_coord[1]*JxW;
                (*ele_mat)(NDIM * i + 0, 22 ) += bas_grad[quad][i][2]*quad_coord[2]*JxW;
                (*ele_mat)(NDIM * i + 0, 23 ) += bas_grad[quad][i][2]*1.*JxW;
                /// 自由度i的y分量 425
                (*ele_mat)(NDIM * i + 1, 12 ) += bas_grad[quad][i][0]*quad_coord[0]*JxW;
                (*ele_mat)(NDIM * i + 1, 13 ) += bas_grad[quad][i][0]*quad_coord[1]*JxW;
                (*ele_mat)(NDIM * i + 1, 14 ) += bas_grad[quad][i][0]*quad_coord[2]*JxW;
                (*ele_mat)(NDIM * i + 1, 15 ) += bas_grad[quad][i][0]*1.*JxW;
                (*ele_mat)(NDIM * i + 1, 4 ) += bas_grad[quad][i][1]*quad_coord[0]*JxW;
                (*ele_mat)(NDIM * i + 1, 5 ) += bas_grad[quad][i][1]*quad_coord[1]*JxW;
                (*ele_mat)(NDIM * i + 1, 6 ) += bas_grad[quad][i][1]*quad_coord[2]*JxW;
                (*ele_mat)(NDIM * i + 1, 7 ) += bas_grad[quad][i][1]*1.*JxW;
                (*ele_mat)(NDIM * i + 1, 16 ) += bas_grad[quad][i][2]*quad_coord[0]*JxW;
                (*ele_mat)(NDIM * i + 1, 17 ) += bas_grad[quad][i][2]*quad_coord[1]*JxW;
                (*ele_mat)(NDIM * i + 1, 18 ) += bas_grad[quad][i][2]*quad_coord[2]*JxW;
                (*ele_mat)(NDIM * i + 1, 19 ) += bas_grad[quad][i][2]*1.*JxW;
                /// 自由度i的z分量 653
                (*ele_mat)(NDIM * i + 2, 20 ) += bas_grad[quad][i][0]*quad_coord[0]*JxW;
                (*ele_mat)(NDIM * i + 2, 21 ) += bas_grad[quad][i][0]*quad_coord[1]*JxW;
                (*ele_mat)(NDIM * i + 2, 22 ) += bas_grad[quad][i][0]*quad_coord[2]*JxW;
                (*ele_mat)(NDIM * i + 2, 23 ) += bas_grad[quad][i][0]*1.*JxW;
                (*ele_mat)(NDIM * i + 2, 16 ) += bas_grad[quad][i][1]*quad_coord[0]*JxW;
                (*ele_mat)(NDIM * i + 2, 17 ) += bas_grad[quad][i][1]*quad_coord[1]*JxW;
                (*ele_mat)(NDIM * i + 2, 18 ) += bas_grad[quad][i][1]*quad_coord[2]*JxW;
                (*ele_mat)(NDIM * i + 2, 19 ) += bas_grad[quad][i][1]*1.*JxW;
                (*ele_mat)(NDIM * i + 2, 8 ) += bas_grad[quad][i][2]*quad_coord[0]*JxW;
                (*ele_mat)(NDIM * i + 2, 9 ) += bas_grad[quad][i][2]*quad_coord[1]*JxW;
                (*ele_mat)(NDIM * i + 2, 10 ) += bas_grad[quad][i][2]*quad_coord[2]*JxW;
                (*ele_mat)(NDIM * i + 2, 11 ) += bas_grad[quad][i][2]*1.*JxW;

            }

        }




}


/*************************************************************************
 * 计算单元刚度矩阵.
 *
 * update #3:材料相关部分
 ************************************************************************/
void LinearTet::buildStiffElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, double T_val) {
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  //update #3:
  tbox::Pointer<Material> material =
  				material_manager->getMaterial("Gold");
  	  if(entity_id==1)
		 material =	material_manager->getMaterial("Silicon");
    else if(entity_id==13)
		 material =	material_manager->getMaterial("MoCu");
    else if(entity_id==2)
		 material =	material_manager->getMaterial("Copper");
    else if(entity_id==4)
		 material =	material_manager->getMaterial("SiO2");
    else if(entity_id==5)
		 material =	material_manager->getMaterial("SiN");     
    else if(entity_id==3)
		 material =	material_manager->getMaterial("Gold");
    else if(entity_id==15)
		 material =	material_manager->getMaterial("GaN");  
    else if(entity_id==16)
		 material =	material_manager->getMaterial("Al2O3"); 
    else if(entity_id==17)
		 material =	material_manager->getMaterial("Alloy"); 
	else if(entity_id==6)
		 material =	material_manager->getMaterial("Aluminum");  
    else
		 material =	material_manager->getMaterial("Gold");

				/*
	if(entity_id==1)
		 material =	material_manager->getMaterial("Silicon");
    else if(entity_id==13)
		 material =	material_manager->getMaterial("MoCu");
    else if(entity_id==2)
		 material =	material_manager->getMaterial("Copper");
    else if(entity_id==4)
		 material =	material_manager->getMaterial("SiO2");
    else if(entity_id==5)
		 material =	material_manager->getMaterial("SiN");     
    else if(entity_id==3)
		 material =	material_manager->getMaterial("Gold");
    else if(entity_id==15)
		 material =	material_manager->getMaterial("GaN");  
    else if(entity_id==16)
		 material =	material_manager->getMaterial("Al2O3"); 
    else if(entity_id==17)
		 material =	material_manager->getMaterial("Alloy");  
    else
		 material =	material_manager->getMaterial("Gold");
				*/

  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  ele_mat->resize(n_dof * NDIM, n_dof * NDIM);
  for (int i = 0; i < n_dof * NDIM; ++i) {
    for (int j = 0; j < n_dof * NDIM; ++j) {
      (*ele_mat)(i, j) = 0.0;
    }
  }
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

  tbox::Array<tbox::Array<double> > moduli = material->getModuli(T_val);

  double a = moduli[0][0];
  double b = moduli[0][1];
  double c = moduli[3][3];
  /// 计算单元刚度矩阵.
  for (int i = 0; i < n_dof; ++i) {
    for (int j = 0; j < n_dof; ++j) {
      for (int l = 0; l < num_quad_pnts; ++l) {
        double JxW = volume * jac * weight[l];
        (*ele_mat)(NDIM* i, NDIM* j) +=
            JxW * (a * bas_grad[l][i][0] * bas_grad[l][j][0] +
                   c * bas_grad[l][i][1] * bas_grad[l][j][1] +
                   c * bas_grad[l][i][2] * bas_grad[l][j][2]);

        (*ele_mat)(NDIM* i + 1, NDIM * j + 1) +=
            JxW * (c * bas_grad[l][i][0] * bas_grad[l][j][0] +
                   a * bas_grad[l][i][1] * bas_grad[l][j][1] +
                   c * bas_grad[l][i][2] * bas_grad[l][j][2]);

        (*ele_mat)(NDIM* i + 2, NDIM * j + 2) +=
            JxW * (c * bas_grad[l][i][0] * bas_grad[l][j][0] +
                   c * bas_grad[l][i][1] * bas_grad[l][j][1] +
                   a * bas_grad[l][i][2] * bas_grad[l][j][2]);

        (*ele_mat)(NDIM* i, NDIM* j + 1) +=
            JxW * b * bas_grad[l][i][0] * bas_grad[l][j][1] +
            JxW * c * bas_grad[l][j][0] * bas_grad[l][i][1];

        (*ele_mat)(NDIM* i, NDIM* j + 2) +=
            JxW * b * bas_grad[l][i][0] * bas_grad[l][j][2] +
            JxW * c * bas_grad[l][j][0] * bas_grad[l][i][2];

        (*ele_mat)(NDIM* i + 1, NDIM * j + 2) +=
            JxW * b * bas_grad[l][i][1] * bas_grad[l][j][2] +
            JxW * c * bas_grad[l][j][1] * bas_grad[l][i][2];

        (*ele_mat)(NDIM* i + 1, NDIM * j) +=
            JxW * b * bas_grad[l][i][1] * bas_grad[l][j][0] +
            JxW * c * bas_grad[l][j][1] * bas_grad[l][i][0];

        (*ele_mat)(NDIM* i + 2, NDIM * j) +=
            JxW * b * bas_grad[l][i][2] * bas_grad[l][j][0] +
            JxW * c * bas_grad[l][j][2] * bas_grad[l][i][0];

        (*ele_mat)(NDIM* i + 2, NDIM * j + 1) +=
            JxW * b * bas_grad[l][i][2] * bas_grad[l][j][1] +
            JxW * c * bas_grad[l][j][2] * bas_grad[l][i][1];
      }
    }
  }
}

/*************************************************************************
 * 计算单元质量矩阵.
 *
 * update #3:材料相关部分
 ************************************************************************/
void LinearTet::buildMassElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, double T_val) {
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  //update #3:
  //update #3:
  tbox::Pointer<Material> material =
  				material_manager->getMaterial("Gold");
  	  if(entity_id==1)
		 material =	material_manager->getMaterial("Silicon");
    else if(entity_id==13)
		 material =	material_manager->getMaterial("MoCu");
    else if(entity_id==2)
		 material =	material_manager->getMaterial("Copper");
    else if(entity_id==4)
		 material =	material_manager->getMaterial("SiO2");
    else if(entity_id==5)
		 material =	material_manager->getMaterial("SiN");     
    else if(entity_id==3)
		 material =	material_manager->getMaterial("Gold");
    else if(entity_id==15)
		 material =	material_manager->getMaterial("GaN");  
    else if(entity_id==16)
		 material =	material_manager->getMaterial("Al2O3"); 
    else if(entity_id==17)
		 material =	material_manager->getMaterial("Alloy"); 
	else if(entity_id==6)
		 material =	material_manager->getMaterial("Aluminum");  
    else
		 material =	material_manager->getMaterial("Gold");

  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  ele_mat->resize(n_dof * NDIM, n_dof * NDIM);
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);

  double density =material->getDensity(T_val);

  for (int i = 0; i < n_dof; ++i) {
    for (int j = 0; j < NDIM; ++j) {
      (*ele_mat)(NDIM* i + j, NDIM * i + j) = 0.25 * volume * jac * density;
    }
  }
#if 0 
  ///  计算单元质量矩阵, 非集中质量.
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> >
    quad_pnt = integrator->getQuadraturePoints(real_vertex);
  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();
    
  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> >
    bas_val = shape_func->value(real_vertex,quad_pnt);
  for (int i = 0;i < n_dof; ++ i) {
    for (int j = 0;j < n_dof; ++ j) {
      for (int l = 0;l < num_quad_pnts;++ l) {
        double JxW = volume*jac*weight[l];
        for(int m = 0; m < NDIM; ++m){
          (*ele_mat)(NDIM*i+m,NDIM*j+m) +=
            JxW*density*(bas_val[l][i]*bas_val[l][j]);
        }
      }
    }
  }
#endif
}

/*************************************************************************
 * 计算单元阻尼矩阵.
 *
 * update #3:材料相关部分
 * //update #6 计算阻尼矩阵 at 2017-04-26
 ************************************************************************/
void LinearTet::buildDumpElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, double T_val)
{
	/// 取出积分器对象.
	  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
	      IntegratorManager<NDIM>::getManager();
	  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
	      integrator_manager->getIntegrator("LinearTetrahedron");

	  /// 取出形函数对象.
	  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
	      ShapeFunctionManager<NDIM>::getManager();
	  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
	      shape_manager->getShapeFunction("LinearTetrahedron");

	  /// 取出材料.
	  tbox::Pointer<MaterialManager<NDIM> > material_manager =
	      MaterialManager<NDIM>::getManager();

	  //update #3:
	  //update #3:
   tbox::Pointer<Material> material =
  				material_manager->getMaterial("Gold");
  	  if(entity_id==1)
		 material =	material_manager->getMaterial("Silicon");
    else if(entity_id==13)
		 material =	material_manager->getMaterial("MoCu");
    else if(entity_id==2)
		 material =	material_manager->getMaterial("Copper");
    else if(entity_id==4)
		 material =	material_manager->getMaterial("SiO2");
    else if(entity_id==5)
		 material =	material_manager->getMaterial("SiN");     
    else if(entity_id==3)
		 material =	material_manager->getMaterial("Gold");
    else if(entity_id==15)
		 material =	material_manager->getMaterial("GaN");  
    else if(entity_id==16)
		 material =	material_manager->getMaterial("Al2O3"); 
    else if(entity_id==17)
		 material =	material_manager->getMaterial("Alloy"); 
	else if(entity_id==6)
		 material =	material_manager->getMaterial("Aluminum");  
    else
		 material =	material_manager->getMaterial("Gold");

	  /// 取出单元上自由度数目.
	  int n_dof = shape_func->getNumberOfDof();
	  ele_mat->resize(n_dof * NDIM, n_dof * NDIM);
	  for(int i=0;i<n_dof*NDIM;i++)
		  for(int j=0;j<n_dof*NDIM;j++)
			  (*ele_mat)(i,j)=0;
	  /// 取出模板单元的面积.
	  double volume = integrator->getElementVolume();
	  /// 取出jacobian矩阵行列式.
	  double jac = integrator->getLocal2GlobalJacobian(real_vertex);

	  double damping =material->getDamping();
	  /// 取出积分点数目.
	  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
	  /// 取出积分点.
	  tbox::Array<hier::DoubleVector<NDIM> >
	    quad_pnt = integrator->getQuadraturePoints(real_vertex);
	  /// 取出积分点的积分权重.
	  tbox::Array<double> weight = integrator->getQuadratureWeights();

	  /// 取出基函数在积分点的值和梯度值.
	  tbox::Array<tbox::Array<double> >
	    bas_val = shape_func->value(real_vertex,quad_pnt);
	  for (int i = 0;i < n_dof; ++ i) {
	    for (int j = 0;j < n_dof; ++ j) {
	      for (int l = 0;l < num_quad_pnts;++ l) {
	        double JxW = volume*jac*weight[l];
	        for(int m = 0; m < NDIM; ++m){
	          (*ele_mat)(NDIM*i+m,NDIM*j+m) +=
	            JxW*damping*(bas_val[l][i]*bas_val[l][j]);
	        }
	      }
	    }
	  }
}

/*************************************************************************
 * 计算单元矩阵.
 *
 * 计算单元矩阵P
 * P=M/(dt^2)+C/(2*dt)+beta*K
 *
 * //update #6  at 2017-05-02
 * Tested OK -05-05
 ************************************************************************/
 void LinearTet::buildElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id,  tbox::Array<double> T_val)
 {
	  /// 取出积分器对象. 这里好像用不到
	  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
	      IntegratorManager<NDIM>::getManager();
	  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
	      integrator_manager->getIntegrator("LinearTetrahedron");

	  /// 取出形函数对象.
	  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
	      ShapeFunctionManager<NDIM>::getManager();
	  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
	      shape_manager->getShapeFunction("LinearTetrahedron");

	  /// 取出单元上自由度数目.
	  int n_dof = shape_func->getNumberOfDof();

	  //计算单元刚度矩阵、质量矩阵和阻尼矩阵
	  //刚度矩阵
	  tbox::Pointer<tbox::Matrix<double> > ele_mat_K = new tbox::Matrix<double>();
	  ele_mat_K->resize(n_dof* NDIM, n_dof* NDIM);
	  //质量矩阵
	  tbox::Pointer<tbox::Matrix<double> > ele_mat_M = new tbox::Matrix<double>();
	  ele_mat_M->resize(n_dof* NDIM, n_dof* NDIM);
	  //阻尼矩阵
	  tbox::Pointer<tbox::Matrix<double> > ele_mat_C = new tbox::Matrix<double>();
	  ele_mat_C->resize(n_dof* NDIM, n_dof* NDIM);

	  //初始化
	  for (int i1 = 0; i1 < n_dof* NDIM; ++i1) {
	    for (int j = 0; j < n_dof* NDIM; ++j) {
	      (*ele_mat_K)(i1, j) = 0.0;
	      (*ele_mat_M)(i1, j) = 0.0;
	      (*ele_mat_C)(i1, j) = 0.0;
	    }
	  }
          double e_Temperature=0;
 	  for(int i=0;i<4;i++)
	     e_Temperature+=T_val[i]/4;
	  //计算
	  buildStiffElementMatrix(real_vertex, dt, time, ele_mat_K, entity_id,e_Temperature);
	  buildMassElementMatrix(real_vertex, dt, time, ele_mat_M, entity_id,e_Temperature);
	  buildDumpElementMatrix(real_vertex, dt, time, ele_mat_C, entity_id,e_Temperature);

	  ele_mat->resize(n_dof * NDIM, n_dof * NDIM);
	  for (int i = 0; i < n_dof * NDIM; ++i) {
	    for (int j = 0; j < n_dof * NDIM; ++j) {
	      (*ele_mat)(i, j) = 0.0;
	    }
	  }

	  for (int i = 0; i < n_dof * NDIM; ++i) {
	    for (int j = 0; j < n_dof * NDIM; ++j) {
//	      (*ele_mat)(i, j) = (*ele_mat_M)(i, j)/(dt*dt)+(*ele_mat_C)(i, j)/(2*dt)+(*ele_mat_K)(i, j)*BETA;
            (*ele_mat)(i, j) = (*ele_mat_K)(i, j);
	    }
	  }
 }

/*************************************************************************
 * 计算单元右端项..
 * Tested OK -05-05
 ************************************************************************/
void LinearTet::buildElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
	NewmarkData *d_newmark,int entity_id,tbox::Array<double> T_val,
    tbox::Array<double> Tolder_val) {
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");
  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  ele_vec->resize(n_dof * NDIM);
  for (int j = 0; j < n_dof * NDIM; ++j) {
    (*ele_vec)[j] = 0.0;
  }

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();
  //update #3:
  //update #3:
  tbox::Pointer<Material> material =
  				material_manager->getMaterial("Gold");
  	  if(entity_id==1)
		 material =	material_manager->getMaterial("Silicon");
    else if(entity_id==13)
		 material =	material_manager->getMaterial("MoCu");
    else if(entity_id==2)
		 material =	material_manager->getMaterial("Copper");
    else if(entity_id==4)
		 material =	material_manager->getMaterial("SiO2");
    else if(entity_id==5)
		 material =	material_manager->getMaterial("SiN");     
    else if(entity_id==3)
		 material =	material_manager->getMaterial("Gold");
    else if(entity_id==15)
		 material =	material_manager->getMaterial("GaN");  
    else if(entity_id==16)
		 material =	material_manager->getMaterial("Al2O3"); 
    else if(entity_id==17)
		 material =	material_manager->getMaterial("Alloy"); 
	else if(entity_id==6)
		 material =	material_manager->getMaterial("Aluminum");  
    else
		 material =	material_manager->getMaterial("Gold");

  //update #8 计算单元温度 直接线性插值计算
  double e_Temperature=0;
  for(int i=0;i<n_dof;i++)
	  e_Temperature+=T_val[i]/4;
  double e_Temperature_older=0;
  for(int i=0;i<n_dof;i++)
	  e_Temperature_older+=Tolder_val[i]/4;

  //模量矩阵
  tbox::Array<tbox::Array<double> > moduli = material->getModuli(e_Temperature);
  double a = moduli[0][0];
  double b = moduli[0][1];
  //double c = moduli[3][3];

  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);

  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  //update #6 05-05pm
  /// 取出积分点数目.
    int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
    /// 取出模板单元的面积.
    double volume = integrator->getElementVolume();
    /// 取出jacobian矩阵行列式.
    double jac = integrator->getLocal2GlobalJacobian(real_vertex);

  //计算热应变
  //应变矩阵=th_stress[1 1 1 0 0 0]
  double alpha=material->getAlpha(e_Temperature);
//  double th_stress=alpha*(e_Temperature-293.15);
  double th_stress=alpha*(e_Temperature-650.);
  //取出积分点处基函数梯度值
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, quad_pnt);
  //热应变单元右端项
  double matrix_BTD[12];
  for(int j1=0;j1<12;j1++)
	  matrix_BTD[j1]=0;
  //模量矩阵*热应力=termal_strain*[1 1 1 0 0 0]T
  double thermal_strain=(a+2*b)*th_stress;

  for(int i=0;i<n_dof;i++)
  {
	  for(int j=0;j<NDIM;j++)
	  {
	  	  for(int l=0;l<num_quad_pnts;l++)
	  	  {
	  		matrix_BTD[i*NDIM+j]+=thermal_strain*bas_grad[l][i][j]*weight[l]*jac*volume;
	  		//cout<<i*NDIM+j<<"   :"<<bas_grad[l][i][j]<<"   BTD:"<<matrix_BTD[i*NDIM+j]<<" thermal_strain*:"<<(thermal_strain*bas_grad[l][i][j]*weight[l]*jac)<<"  thermal_strain:"<<thermal_strain<<endl;
	  	  }
	  	//  cout<<matrix_BTD[i*NDIM+j]<<"   ";
	  }
	 // cout<<'\n';
  }

  //静态时
#if 1

  for (int i = 0; i < n_dof; ++i) {
    for (int l = 0; l < num_quad_pnts; ++l) {
      double JxWb = volume * jac * weight[l] * bas_val[l][i];
      for (int j = 0; j < NDIM; ++j) {
          (*ele_vec)[i * NDIM + j] += matrix_BTD[i*NDIM+j];
      }
    }
  }
#endif

  //瞬态时
#if 0
  //计算单元刚度矩阵、质量矩阵和阻尼矩阵
  //刚度矩阵
  tbox::Pointer<tbox::Matrix<double> > ele_mat_K = new tbox::Matrix<double>();
  ele_mat_K->resize(n_dof* NDIM, n_dof* NDIM);
  //质量矩阵
  tbox::Pointer<tbox::Matrix<double> > ele_mat_M = new tbox::Matrix<double>();
  ele_mat_M->resize(n_dof* NDIM, n_dof* NDIM);
  //阻尼矩阵
  tbox::Pointer<tbox::Matrix<double> > ele_mat_C = new tbox::Matrix<double>();
  ele_mat_C->resize(n_dof* NDIM, n_dof* NDIM);

  //初始化
  for (int i1 = 0; i1 < n_dof* NDIM; ++i1) {
    for (int j = 0; j < n_dof* NDIM; ++j) {
      (*ele_mat_K)(i1, j) = 0.0;
      (*ele_mat_M)(i1, j) = 0.0;
      (*ele_mat_C)(i1, j) = 0.0;
    }
  }

  //计算
  buildStiffElementMatrix(real_vertex, dt, time, ele_mat_K, entity_id,e_Temperature);
  buildMassElementMatrix(real_vertex, dt, time, ele_mat_M, entity_id,e_Temperature);
  buildDumpElementMatrix(real_vertex, dt, time, ele_mat_C, entity_id,e_Temperature);


  //计算Q、R矩阵
  //Q=-[-(2*M)/(Dt^2)+(1-2*beta)*K]
  //R=-[M/(Dt^2)-c/(2*Dt)+beta*K]
  //初始化矩阵Q、R
  tbox::Pointer<tbox::Matrix<double> > ele_mat_Q = new tbox::Matrix<double>();
  ele_mat_Q->resize(n_dof* NDIM, n_dof* NDIM);
  tbox::Pointer<tbox::Matrix<double> > ele_mat_R = new tbox::Matrix<double>();
  ele_mat_R->resize(n_dof* NDIM, n_dof* NDIM);

  for (int i1 = 0; i1 < n_dof* NDIM; ++i1) {
    for (int j = 0; j < n_dof* NDIM; ++j) {
      (*ele_mat_Q)(i1, j) = ((2*(*ele_mat_M)(i1, j))/(dt*dt))-((1-2*BETA)*(*ele_mat_K)(i1, j));
      (*ele_mat_R)(i1, j) = -(((*ele_mat_M)(i1, j)/(dt*dt))-((*ele_mat_C)(i1, j)/(2*dt))+(BETA*(*ele_mat_K)(i1, j)));
    }
  }

  //计算Q*solution（t-Dt）和R*solution(t-2*Dt)
  //定义一个向量用于存储计算结果、
  //update #6 -05-05  指针分配内存
  tbox::Pointer<tbox::Vector<double> > solutions_vec= new tbox::Vector<double>();
  solutions_vec->resize(n_dof * NDIM);
  for (int j = 0; j < n_dof * NDIM; ++j) {
    (*solutions_vec)[j] = 0.0;
  }

  for (int i1 = 0; i1 < n_dof* NDIM; ++i1) {
    for (int j = 0; j < n_dof* NDIM; ++j) {
    	(*solutions_vec)[i1]+=(((*ele_mat_Q)(i1, j))*((d_newmark[j/3]).v_solution_old[j%3]))+
    			(((*ele_mat_R)(i1, j))*((d_newmark[j/3]).v_solution_older[j%3]));
    }
  }

  //
  for (int i = 0; i < n_dof; ++i) {
    for (int l = 0; l < num_quad_pnts; ++l) {
      double JxWb = volume * jac * weight[l] * bas_val[l][i];
      for (int j = 0; j < NDIM; ++j) {
   //   if (j % 3 == 2) {
          (*ele_vec)[i * NDIM + j] +=matrix_BTD[i*NDIM+j]-JxWb * 0.0e-3;
   //    }
      }
    }
  }
  for (int i1 = 0; i1 < n_dof; ++i1) {
	  for(int j1=0;j1<NDIM;j1++)
	  (*ele_vec)[i1*NDIM+j1]=BETA*((*ele_vec)[i1*NDIM+j1])+(1-2*BETA)*((d_newmark[i1]).v_rhs_old[j1])+
			  BETA*((d_newmark[i1]).v_rhs_older[j1])+(*solutions_vec)[i1*NDIM+j1];
  }
#endif
}

/////////////////////////////////////////////////update #8//////////////////////////////////////////////

/*************************************************************************
 * 计算热计算单元矩阵.
 *
 ************************************************************************/
 void LinearTet::buildTh_ElementMatrix(
     tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
     const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, tbox::Array<double> T_val)
 {
	 /// 取出积分器对象.
	  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
	      IntegratorManager<NDIM>::getManager();
	  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
	      integrator_manager->getIntegrator("LinearTetrahedron");

	  /// 取出形函数对象.
	  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
	      ShapeFunctionManager<NDIM>::getManager();
	  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
	      shape_manager->getShapeFunction("LinearTetrahedron");

	  /// 取出材料.
	  tbox::Pointer<MaterialManager<NDIM> > material_manager =
	      MaterialManager<NDIM>::getManager();

	  //update #3:
	  //update #3:
 double e_Temperature=0;
  for(int i=0;i<4;i++)
	  e_Temperature+=T_val[i]/4;
  tbox::Pointer<Material> material =
  				material_manager->getMaterial("Gold");
  	  if(entity_id==1)
		 material =	material_manager->getMaterial("Silicon");
    else if(entity_id==13)
		 material =	material_manager->getMaterial("MoCu");
    else if(entity_id==2)
		 material =	material_manager->getMaterial("Copper");
    else if(entity_id==4)
		 material =	material_manager->getMaterial("SiO2");
    else if(entity_id==5)
		 material =	material_manager->getMaterial("SiN");     
    else if(entity_id==3)
		 material =	material_manager->getMaterial("Gold");
    else if(entity_id==15)
		 material =	material_manager->getMaterial("GaN");  
    else if(entity_id==16)
		 material =	material_manager->getMaterial("Al2O3"); 
    else if(entity_id==17)
		 material =	material_manager->getMaterial("Alloy"); 
	else if(entity_id==6)
		 material =	material_manager->getMaterial("Aluminum");  
    else
		 material =	material_manager->getMaterial("Gold");

	  /// 取出单元上自由度数目.
	  int n_dof = shape_func->getNumberOfDof();
	  ele_mat->resize(n_dof, n_dof);
	  for (int i = 0; i < n_dof; ++i) {
	    for (int j = 0; j < n_dof; ++j) {
	      (*ele_mat)(i, j) = 0.0;
	    }
	  }
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
	  tbox::Array<tbox::Array<double> > bas_val =
	  	      shape_func->value(real_vertex, quad_pnt);

	  double density=material->getDensity(e_Temperature);
	  double K=material->getK(e_Temperature);
	  double Cp=material->getCp(e_Temperature);

	  /// 计算单元刚度矩阵.
	  for (int i = 0; i < n_dof; ++i) {
	    for (int j = 0; j < n_dof; ++j) {
	      for (int l = 0; l < num_quad_pnts; ++l) {
	        double JxW = volume * jac * weight[l];
                (*ele_mat)(i,j) += JxW*(density*Cp*bas_val[l][i]*bas_val[l][j]+
	        					   (bas_grad[l][i][0]*bas_grad[l][j][0]+
	        			            bas_grad[l][i][1]*bas_grad[l][j][1]+
								    bas_grad[l][i][2]*bas_grad[l][j][2])*K*dt);
	      }
	    }
	  }
 }

 /*************************************************************************
  * 计算热计算单元右端项.
  *
  ************************************************************************/
 void LinearTet::buildTh_ElementRHS(
     tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
     const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
	 int entity_id, tbox::Array<double> T_val, double e_ThermalSource)
 {
	  /// 取出积分器对象.
	  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
	      IntegratorManager<NDIM>::getManager();
	  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
	      integrator_manager->getIntegrator("LinearTetrahedron");

	  /// 取出形函数对象.
	  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
	      ShapeFunctionManager<NDIM>::getManager();
	  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
	      shape_manager->getShapeFunction("LinearTetrahedron");

	  /// 取出材料.
	  tbox::Pointer<MaterialManager<NDIM> > material_manager =
	      MaterialManager<NDIM>::getManager();
	  //update #3:
	  //update #3:
 double e_Temperature=0;
  for(int i=0;i<4;i++)
	  e_Temperature+=T_val[i]/4;
   tbox::Pointer<Material> material =
  				material_manager->getMaterial("Gold");

  	  if(entity_id==1)
		 material =	material_manager->getMaterial("Silicon");
    else if(entity_id==13)
		 material =	material_manager->getMaterial("MoCu");
    else if(entity_id==2)
		 material =	material_manager->getMaterial("Copper");
    else if(entity_id==4)
		 material =	material_manager->getMaterial("SiO2");
    else if(entity_id==5)
		 material =	material_manager->getMaterial("SiN");     
    else if(entity_id==3)
		 material =	material_manager->getMaterial("Gold");
    else if(entity_id==15)
		 material =	material_manager->getMaterial("GaN");  
    else if(entity_id==16)
		 material =	material_manager->getMaterial("Al2O3"); 
    else if(entity_id==17)
		 material =	material_manager->getMaterial("Alloy"); 
	else if(entity_id==6)
		 material =	material_manager->getMaterial("Aluminum");  
    else
		 material =	material_manager->getMaterial("Gold");

	  /// 取出单元上自由度数目.
	  int n_dof = shape_func->getNumberOfDof();
	  ele_vec->resize(n_dof);
	  for (int j = 0; j < n_dof; ++j) {
	    (*ele_vec)[j] = 0.0;
	  }


	  /// 取出积分点.
	  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
	      integrator->getQuadraturePoints(real_vertex);

	  /// 取出积分点的积分权重.
	  tbox::Array<double> weight = integrator->getQuadratureWeights();

	  /// 取出基函数在积分点的值和梯度值.
	  tbox::Array<tbox::Array<double> > bas_val =
	      shape_func->value(real_vertex, quad_pnt);

	  //静态时
	#if STATIC_S
	//update #6 05-05pm
	/// 取出积分点数目.
	  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
	  /// 取出模板单元的面积.
	  double volume = integrator->getElementVolume();
	  /// 取出jacobian矩阵行列式.
	  double jac = integrator->getLocal2GlobalJacobian(real_vertex);

	  double f1=0.0;//热源
	  for (int i = 0; i < n_dof; ++i) {
	    for (int l = 0; l < num_quad_pnts; ++l) {
	      double JxWb = volume * jac * weight[l] * bas_val[l][i];
	          (*ele_vec)[i] += JxWb * f1 * dt;
	      }
	    }
	  }
	#endif

	  //瞬态时
	#if TIME_S
    //rhs=
	//update #6 05-05pm
	/// 取出积分点数目.
	  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
	  /// 取出模板单元的面积.
	  double volume = integrator->getElementVolume();
	  /// 取出jacobian矩阵行列式.
	  double jac = integrator->getLocal2GlobalJacobian(real_vertex);

	  //先计算单元质量矩阵
	  tbox::Pointer<tbox::Matrix<double> > ele_mat_M = new tbox::Matrix<double>();
	  ele_mat_M->resize(n_dof,n_dof);
	  double density=material->getDensity(e_Temperature);
	  double Cp=material->getCp(e_Temperature);
	  for (int i = 0; i < n_dof; ++i) {
	    for (int j = 0; j < n_dof; ++j) {
	      for (int l = 0; l < num_quad_pnts; ++l) {
	    	  double JxW = volume * jac * weight[l];
              (*ele_mat_M)(i,j) += JxW*density*Cp*(bas_val[l][i]*bas_val[l][j]);
	      }
	    }
	  }

      /// M/dt T^{n-1}
	  for (int i1 = 0; i1 < n_dof; ++i1) {
	    for (int j = 0; j < n_dof; ++j) {
	    	(*ele_vec)[i1]+=((*ele_mat_M)(i1, j))*T_val[j];
	    }
	  }

	  /// 计算单元右端项.
	  for (int i = 0; i < n_dof; ++i) {
	    for (int l = 0; l < num_quad_pnts; ++l) {
	      double JxW = volume * jac * weight[l];
	      (*ele_vec)[i] += JxW * e_ThermalSource * bas_val[l][i] * dt;
	    }
	  }
	#endif
 }

 ////////////////////////////////////////////////////////////////////////////////////////////////////////

 /////////////////////////////////////////////////update #9//////////////////////////////////////////////

 /*************************************************************************
  * 计算电计算单元矩阵.
  *
  ************************************************************************/
  void LinearTet::buildE_ElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, tbox::Array<double> T_val)
  {
 	 /// 取出积分器对象.
 	  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
 	      IntegratorManager<NDIM>::getManager();
 	  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
 	      integrator_manager->getIntegrator("LinearTetrahedron");

 	  /// 取出形函数对象.
 	  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
 	      ShapeFunctionManager<NDIM>::getManager();
 	  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
 	      shape_manager->getShapeFunction("LinearTetrahedron");

 	  /// 取出材料.
 	  tbox::Pointer<MaterialManager<NDIM> > material_manager =
 	      MaterialManager<NDIM>::getManager();

 	  //update #3:
 	  //update #3:
 double e_Temperature=0;
  for(int i=0;i<4;i++)
	  e_Temperature+=T_val[i]/4;
   tbox::Pointer<Material> material =
  				material_manager->getMaterial("Gold");

  	 if(entity_id==1)
		 material =	material_manager->getMaterial("Silicon");
    else if(entity_id==13)
		 material =	material_manager->getMaterial("MoCu");
    else if(entity_id==2)
		 material =	material_manager->getMaterial("Copper");
    else if(entity_id==4)
		 material =	material_manager->getMaterial("SiO2");
    else if(entity_id==5)
		 material =	material_manager->getMaterial("SiN");     
    else if(entity_id==3)
		 material =	material_manager->getMaterial("Gold");
    else if(entity_id==15)
		 material =	material_manager->getMaterial("GaN");  
    else if(entity_id==16)
		 material =	material_manager->getMaterial("Al2O3"); 
    else if(entity_id==17)
		 material =	material_manager->getMaterial("Alloy"); 
	else if(entity_id==6)
		 material =	material_manager->getMaterial("Aluminum");  
    else
		 material =	material_manager->getMaterial("Gold");
 	  /// 取出单元上自由度数目.
 	  int n_dof = shape_func->getNumberOfDof();
 	  ele_mat->resize(n_dof, n_dof);
 	  for (int i = 0; i < n_dof; ++i) {
 	    for (int j = 0; j < n_dof; ++j) {
 	      (*ele_mat)(i, j) = 0.0;
 	    }
 	  }
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

 	  double sigma=material->getSigma(e_Temperature);

 	  /// 计算单元刚度矩阵.
 	  for (int i = 0; i < n_dof; ++i) {
 	    for (int j = 0; j < n_dof; ++j) {
 	      for (int l = 0; l < num_quad_pnts; ++l) {
 	        double JxW = volume * jac * weight[l];
 	        	(*ele_mat)(i,j) += JxW*((bas_grad[l][i][0]*bas_grad[l][j][0]+
 	        			            bas_grad[l][i][1]*bas_grad[l][j][1]+
 								    bas_grad[l][i][2]*bas_grad[l][j][2])*sigma);
 	      }
 	    }
 	  }
  }

  /*************************************************************************
   * 计算电计算单元右端项.
   *
   ************************************************************************/
  void LinearTet::buildE_ElementRHS(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
 	 int entity_id, tbox::Array<double> T_val)
  {
 	  /// 取出积分器对象.
 	  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
 	      IntegratorManager<NDIM>::getManager();
 	  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
 	      integrator_manager->getIntegrator("LinearTetrahedron");

 	  /// 取出形函数对象.
 	  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
 	      ShapeFunctionManager<NDIM>::getManager();
 	  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
 	      shape_manager->getShapeFunction("LinearTetrahedron");

 	  /// 取出材料.
 	  tbox::Pointer<MaterialManager<NDIM> > material_manager =
 	      MaterialManager<NDIM>::getManager();
 	  //update #3:
 	  //update #3:
 double e_Temperature=0;
  for(int i=0;i<4;i++)
	  e_Temperature+=T_val[i]/4;
  tbox::Pointer<Material> material =
  				material_manager->getMaterial("Gold");
  	  if(entity_id==1)
		 material =	material_manager->getMaterial("Silicon");
    else if(entity_id==13)
		 material =	material_manager->getMaterial("MoCu");
    else if(entity_id==2)
		 material =	material_manager->getMaterial("Copper");
    else if(entity_id==4)
		 material =	material_manager->getMaterial("SiO2");
    else if(entity_id==5)
		 material =	material_manager->getMaterial("SiN");     
    else if(entity_id==3)
		 material =	material_manager->getMaterial("Gold");
    else if(entity_id==15)
		 material =	material_manager->getMaterial("GaN");  
    else if(entity_id==16)
		 material =	material_manager->getMaterial("Al2O3"); 
    else if(entity_id==17)
		 material =	material_manager->getMaterial("Alloy"); 
	else if(entity_id==6)
		 material =	material_manager->getMaterial("Aluminum");  
    else
		 material =	material_manager->getMaterial("Gold");

 	  /// 取出单元上自由度数目.
 	  int n_dof = shape_func->getNumberOfDof();
 	  ele_vec->resize(n_dof);
 	  for (int j = 0; j < n_dof; ++j) {
 	    (*ele_vec)[j] = 0.0;
 	  }


 	  /// 取出积分点.
 	  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
 	      integrator->getQuadraturePoints(real_vertex);

 	  /// 取出积分点的积分权重.
 	  tbox::Array<double> weight = integrator->getQuadratureWeights();

 	  /// 取出基函数在积分点的值和梯度值.
 	  tbox::Array<tbox::Array<double> > bas_val =
 	      shape_func->value(real_vertex, quad_pnt);


 	//update #6 05-05pm
 	/// 取出积分点数目.
 	  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
 	  /// 取出模板单元的面积.
 	  double volume = integrator->getElementVolume();
 	  /// 取出jacobian矩阵行列式.
 	  double jac = integrator->getLocal2GlobalJacobian(real_vertex);

 	  double q=0.0;//电荷源
 	  for (int i = 0; i < n_dof; ++i) {
 	    for (int l = 0; l < num_quad_pnts; ++l) {
 	      double JxWb = volume * jac * weight[l] * bas_val[l][i];
 	          (*ele_vec)[i] += JxWb * q;
 	      }
 	    }

  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////


const int LinearTet::getProblemDim() { return NDIM; }

//获取实体自由度
tbox::Array<int> LinearTet::getNumberOfDofOnEntity() {
  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");
  return shape_func->getNumberOfDofOnEntity();
}
