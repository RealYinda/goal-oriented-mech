//
// 文件名:     PatchStrategy.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:45:28 2014 $
// 描述　:     网格片策略类派生类实现.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "Patch.h"
#include "PatchStrategy.h"

/*************************************************************************
 *
 * 构造函数.
 *
 *************************************************************************/
PatchStrategy::PatchStrategy() {}

/*************************************************************************
 *
 * 构造函数.
 *
 ************************************************************************/
PatchStrategy::~PatchStrategy() {}

/*************************************************************************
 *
 *  初始化指定的积分构件.
 *
 *  注册待填充的数据片或待调度内存空间的数据片到积分构件.
 ************************************************************************/
void PatchStrategy::initializeComponent(
    algs::IntegratorComponent<NDIM>* component) const {
  TBOX_ERROR("PatchStrategy::initializeComponent() need to be implemented！");
}

/*************************************************************************
 *  初始化数据片（支持初值构件）.
 ************************************************************************/
void PatchStrategy::initializePatchData(hier::Patch<NDIM>& patch,
                                        const double time,
                                        const bool initial_time,
                                        const string& component_name) {
  TBOX_ERROR("PatchStrategy::initializePatchData() need to be implemented！");
}

/*************************************************************************
 *  输出数据成员到重启动数据库.
 ************************************************************************/
void PatchStrategy::putToDatabase(tbox::Pointer<tbox::Database> db) {
  TBOX_ERROR("PatchStrategy::putToDatabase() need to be implemented！");
}

double PatchStrategy::getPatchDt(hier::Patch<NDIM>& patch, const double time,
                                 const bool initial_time,
                                 const int flag_last_dt, const double last_dt,
                                 const string& component_name) {
  return 0.0;
}

/*************************************************************************
 * 完成单个网格片上的数值计算（支持数值构件）.
 ************************************************************************/
void PatchStrategy::computeOnPatch(hier::Patch<NDIM>& patch, const double time,
                                   const double dt, const bool initial_time,
                                   const string& component_name) {
  TBOX_ERROR("PatchStrategy::computeOnPatch() need to be implemented！");
}

/*************************************************************************
 * 完成单个网格片上规约值计算（支持规约构件）.
 ************************************************************************/
void PatchStrategy::reduceOnPatch(double* vector, int len,
                                  hier::Patch<NDIM>& patch, const double time,
                                  const double dt,
                                  const string& component_name) {
  TBOX_ERROR("PatchStrategy::reduceOnPatch() need to be implemented！");
}

/*************************************************************************
 * 注册可视化数据片.
 ************************************************************************/
void PatchStrategy::registerPlotData(
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer) {
  TBOX_ERROR("PatchStrategy::registerPlotData() need to be implemented！");
}

/*************************************************************************
 *  从输入数据库读入数据.
 ************************************************************************/
void PatchStrategy::getFromInput(tbox::Pointer<tbox::Database> db) {
  TBOX_ERROR("PatchStrategy::getFromInput() need to be implemented！");
}
