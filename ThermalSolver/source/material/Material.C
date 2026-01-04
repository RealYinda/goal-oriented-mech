//
// 文件名:     Material.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:34:22 2014 $
// 描述　:     线弹性材料
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "Material.h"

/*************************************************************************
 * 构造函数.
 * ***********************************************************************
 * update #2 添加默认材料类型
 * at 2017-04-20 by tong
 * 初始化只是提供室温下的值，默认为300K
 *************************************************************************/
Material::Material(const string name) : BaseMaterial<NDIM>(name) {
	d_damping=0.0;	//update #6
	if(name=="Silicon")
		{
			  d_young_modulus = 170e9;
			  d_possion_ratio = 0.28;
			  d_density = 2329;
			  d_Cp=700;
			  d_K=130;
			  d_alpha=2.6e-6;
			  d_mur=1;
			  d_epsilonr=11.7;
			  d_sigma=1e-12;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  	  Param_T[i][j]=Silicon_Param_T[i][j];
		}
		else if(name=="Copper")
		{
			  d_young_modulus = 110e9;
			  d_possion_ratio = 0.35;
			  d_density = 8700;
			  d_Cp=385;
			  d_K=400;
			  d_alpha=17e-6;
			  d_mur=1;
			  d_epsilonr=1;
			  d_sigma=5.98e7;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=Copper_Param_T[i][j];
		}
		else if(name=="BCB")
		{
			  d_young_modulus = 2.9e9;
			  d_possion_ratio = 0.34;
			  d_density = 1050;
			  d_Cp=2128;
			  d_K=0.29421;
			  d_alpha=42.829e-6;
			  d_mur=1;
			  d_epsilonr=2.65;
			  d_sigma=1200;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=BCB_Param_T[i][j];
		}
		else if(name=="HT_m")
		{
			  d_young_modulus = 2.55e9;
			  d_possion_ratio = 0.4;
			  d_density = 1320;
			  d_Cp=200;
			  d_K=0.21;
			  d_alpha=36e-6;
			  d_mur=1;
			  d_epsilonr=2;
			  d_sigma=1e-12;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=Copper_Param_T[i][j];
		}
		else if(name=="Nickel")
		{
			  d_young_modulus = 210e9;
			  d_possion_ratio = 0.31;
			  d_density = 8900;
			  d_Cp=440;
			  d_K=91;
			  d_alpha=13.4e-6;
			  d_mur=1;
			  d_epsilonr=1;
			  d_sigma=1.43e7;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=Copper_Param_T[i][j];
		}
		else if(name=="Solder_bump")
		{
			  d_young_modulus = 10e9;
			  d_possion_ratio = 0.4;
			  d_density = 9000;
			  d_Cp=150;
			  d_K=50;
			  d_alpha=21e-6;
			  d_mur=1;
			  d_epsilonr=1;
			  d_sigma=6.67e6;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=Copper_Param_T[i][j];
		}
		else if(name=="Underfill")
		{
			  d_young_modulus = 170e9;
			  d_possion_ratio = 0.25;
			  d_density = 2200;
			  d_Cp=7.3;
			  d_K=1.39;
			  d_alpha=2.6e-6;
			  d_mur=1;
			  d_epsilonr=2.1;
			  d_sigma=1e-6;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=Copper_Param_T[i][j];
		}
		else if(name=="Gold")
		{
			  d_young_modulus = 70e9;
			  d_possion_ratio = 0.44;
			  d_density = 19300;
			  d_Cp=129;
			  d_K=317;
			  d_alpha=14.2e-6;
			  d_mur=1;
			  d_epsilonr=1;
			  d_sigma=45.6e6;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=Gold_Param_T[i][j];
		}
		else if(name=="SiO2")
		{
			  d_young_modulus = 70e9;
			  d_possion_ratio = 0.17;
			  d_density = 2200;
			  d_Cp=730;
			  d_K=1.4;
			  d_alpha=0.5e-6;
			  d_mur=1;
			  d_epsilonr=4.2;
			  d_sigma=1e-12;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=SiO2_Param_T[i][j];
		}
		else if(name=="SiN")
		{
			  d_young_modulus = 250e9;
			  d_possion_ratio = 0.23;
			  d_density = 3100;
			  d_Cp=700;
			  d_K=20;
			  d_alpha=2.3e-6;
			  d_mur=1;
			  d_epsilonr=1;
			  d_sigma=1e-12;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=SiN_Param_T[i][j];
     	}
		else if(name=="MoCu") //钼铜
		{
			  d_young_modulus = 200e9;
			  d_possion_ratio = 0.3;//0.32;
			  d_density = 9910; //大于等于9910  MoCu25 9700
			  d_Cp=385;
			  d_K=208; //实际是大于等于150              180
			  d_alpha=6e-6;//5.6+/-1.5           7.9+/-2
			  d_mur=1;
			  d_epsilonr=1;
			  d_sigma=2e7;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=MoCu_Param_T[i][j];
		}
		else if(name=="GaAs")
		{
			  d_young_modulus = 85e9;
			  d_possion_ratio = 0.31;
			  d_density = 5370;
			  d_Cp=550;
			  d_K=35;
			  d_alpha=5.9e-6;//5.8e-6;
			  d_mur=1;
			  d_epsilonr=12.9;
			  d_sigma=1e-6;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=GaAs_Param_T[i][j];
     	}
		else if(name=="GaN")
		{
			  d_young_modulus = 289.50000e9;
			  d_possion_ratio = 0.001;
			  d_density = 6070;
			  d_Cp=200;
			  d_K=0.21;
			  d_alpha=5.9e-6;//5.8e-6;
			  d_mur=1;
			  d_epsilonr=9.0;
			  d_sigma=1e-20;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=GaN_Param_T[i][j];
     	}
		else if(name=="Aluminum")
		{
			  d_young_modulus = 70e9;
			  d_possion_ratio = 0.35;
			  d_density = 2700;
			  d_Cp=904;
			  d_K=237;
			  d_alpha=23.1e-6;
			  d_mur=1;
			  d_epsilonr=1;
			  d_sigma=3.55e7;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=Aluminum_Param_T[i][j];
     	}
		else if(name=="Alloy")//new 未入库
		{
			  d_young_modulus = 143.437e9;
			  d_possion_ratio = 0.3;
			  d_density = 8260;
			  d_Cp=470;
			  d_K=15;
			  d_alpha=5.3e-6;
			  d_mur=1;
			  d_epsilonr=1;
			  d_sigma=1.82e6;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=Alloy_Param_T[i][j];
     	}
		else if(name=="Al2O3")//new 未入库
		{
			  d_young_modulus = 340e9;
			  d_possion_ratio = 0.22;
			  d_density = 3.92e3;
			  d_Cp=840;
			  d_K=26;
			  d_alpha=8.5e-6;
			  d_mur=1;
			  d_epsilonr=10;
			  d_sigma=10e-14;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=Al2O3_Param_T[i][j];
     	}
		else{//默认材料
			  d_young_modulus = 2100;
			  d_possion_ratio = 0.269;
			  d_density = 7850;
			  d_Cp=900;
			  d_K=238;
			  d_alpha=23e-6;
			  d_mur=1;
			  d_epsilonr=1;
			  d_sigma=0;
			  for(int i=0;i<9;i++)
				  for(int j=0;j<5;j++)
					  Param_T[i][j]=Default_Param_T[i][j];
		}

		d_moduli.resizeArray(6);
	    for (int i = 0; i < 6; ++i) {
	    	d_moduli[i].resizeArray(6);
	    	for (int j = 0; j < 6; ++j) {
	    		d_moduli[i][j] = 0.0;
	    	}
	    }
	    computeMuduli();
}

/*************************************************************************
 * 析构函数.
 ************************************************************************/
Material::~Material() {}

//update #7
/*************************************************************************
 * 计算材料温变参数.
 ************************************************************************/
double Material::computeTimeParam(double T, Material_TParam_Type Param_Type)
{
	double a0=Param_T[Param_Type][0];
	double a1=Param_T[Param_Type][1];
	double a2=Param_T[Param_Type][2];
	double a3=Param_T[Param_Type][3];
	double a4=Param_T[Param_Type][4];
	double ans=(a0+a1*T+a2*pow(T,2)+a3*pow(T,3)+a4*pow(T,4));
	return ans;
}

/*************************************************************************
 * 获取杨氏模量.
 ************************************************************************/
double Material::getYoungModulus(double T)
{
	//Material_TParam m_param =e_young_modulus;
	//if(T!=300)
		d_young_modulus= computeTimeParam(T, e_young_modulus);
	return d_young_modulus;
}

/*************************************************************************
 * 获取波松比.
 ************************************************************************/
double Material::getPossionRatio(double T)
{
	//Material_TParam m_param=e_possion_ratio;
	//if(T!=300)
		d_possion_ratio= computeTimeParam(T, e_possion_ratio);
	return d_possion_ratio;
}

/*************************************************************************
 * 获取密度.
 ************************************************************************/
double Material::getDensity(double T)
{
	//Material_TParam m_param=e_density;
	//if(T!=300)
		d_density= computeTimeParam(T, e_density);
	return d_density;
}

//update #6
/*************************************************************************
 * 获取阻尼系数.
 *
 * 暂时，没有考虑温变性质
 ************************************************************************/
double Material::getDamping()
{
	return d_damping;
}

//////////////////////////////update #7///////////////////////////////////

/*************************************************************************
 * 获取导热系数.
 ************************************************************************/
double Material::getK(double T)
{
	//Material_TParam m_param=e_K;
	//if(T!=300)
		d_K= computeTimeParam(T, e_K);
	return d_K;
}

/*************************************************************************
 * 获取常压热容.
 ************************************************************************/
double Material::getCp(double T)
{
	//Material_TParam m_param=e_Cp;
	//if(T!=300)
		d_Cp= computeTimeParam(T, e_Cp);
	return d_Cp;
}

/*************************************************************************
 * 获取热膨胀系数.
 ************************************************************************/
double Material::getAlpha(double T)
{
	//Material_TParam m_param=e_alpha;
	//if(T!=300)
        d_alpha= computeTimeParam(T, e_alpha);
	return d_alpha;
}

//////////////////////////////////////////////////////////////////////////

/////////////////////////////////update #9////////////////////////////////
//电学材料相关函数

/*************************************************************************
 * 获取相对磁导率.
 ************************************************************************/
double Material::getMur(double T)
{
	//if(T!=300)
		d_mur= computeTimeParam(T, e_mur);
	return d_mur;
}

/*************************************************************************
 * 获取相对介电常数.
 ************************************************************************/
double Material::getEpsilonr(double T)
{
	//if(T!=300)
		d_epsilonr= computeTimeParam(T, e_epsilonr);
	return d_epsilonr;
}

/*************************************************************************
 * 获取电导率.
 ************************************************************************/
double Material::getSigma(double T)
{
	//if(T!=300)
		d_sigma= computeTimeParam(T, e_sigma);
	return d_sigma;
}

//////////////////////////////////////////////////////////////////////////

/*************************************************************************
 * 获取模量矩阵.
 ************************************************************************/
tbox::Array<tbox::Array<double> > Material::getModuli(double T)
{
	//7if(T!=300)
{
		//Material_TParam m_param_1=e_young_modulus;
		//Material_TParam m_param=e_possion_ratio;
			d_young_modulus= computeTimeParam(T, e_young_modulus);
			d_possion_ratio= computeTimeParam(T, e_possion_ratio);
		computeMuduli();
	}
	return d_moduli;
}

/*************************************************************************
 * 计算模量矩阵.
 ************************************************************************/
void Material::computeMuduli() {
  double mu = 0.5 * d_young_modulus / (1 + d_possion_ratio);
  double lambda = (d_young_modulus * d_possion_ratio) /
                  ((1 + d_possion_ratio) * (1 - 2 * d_possion_ratio));

  d_moduli[2][2] = d_moduli[1][1] = d_moduli[0][0] = lambda + 2.0 * mu;
  d_moduli[1][2] = d_moduli[0][1] = d_moduli[0][2] = lambda;
  d_moduli[2][1] = d_moduli[1][0] = d_moduli[2][0] = lambda;
  d_moduli[5][5] = d_moduli[4][4] = d_moduli[3][3] = mu;
}
