//
// 文件名:     MaterialManager.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:33:51 2014 $
// 描述　:     材料管理器
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "JAUMIN_config.h"
#include "Pointer.h"
#include "BaseMaterial.h"
#include "Utilities.h"
#include <string>
#include <list>
#include "MaterialManager.h"
#ifndef NULL
#define NULL (0)
#endif

template <int DIM>
tbox::Pointer<MaterialManager<DIM> >
    MaterialManager<DIM>::s_material_manager_instance = 0;

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for MaterialManager objects.      *
*                                                                       *
*************************************************************************
*/
template <int DIM>
MaterialManager<DIM>::MaterialManager() {}

template <int DIM>
MaterialManager<DIM>::~MaterialManager() {}

template <int DIM>
tbox::Pointer<MaterialManager<DIM> > MaterialManager<DIM>::getManager() {
  if (s_material_manager_instance.isNull()) {
    s_material_manager_instance = new MaterialManager<DIM>();
  }
  return (s_material_manager_instance);
}

/************************************************************************
*                                                                       *
* Add material to list.                                                     *
*                                                                       *
*************************************************************************/
template <int DIM>
void MaterialManager<DIM>::addMaterial(
    tbox::Pointer<BaseMaterial<DIM> > material) {
  const string &name = material->getName();
  bool exist = false;
  unsigned int size = d_vec_material.size();
  for (unsigned int i = 0; i < size; ++i) {
    exist = d_vec_material[i]->isMatching(name);
    if (exist) break;
  }
  if (!exist)
    d_vec_material.push_back(material);
  else
    TBOX_ERROR("material \"" << name
                             << "\" is existed. please choose another name!!"
                             << endl);
}

/************************************************************************
*                                                                       *
* Search material lists for material matching request. *
*                                                                       *
*************************************************************************/

template <int DIM>
tbox::Pointer<BaseMaterial<DIM> > MaterialManager<DIM>::getMaterial(
    const string &material_name) const {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!material_name.empty());
#endif

  tbox::Pointer<BaseMaterial<DIM> > material = NULL;
  bool found = false;
  int vec_size = d_vec_material.size();
  for (int i = 0; i < vec_size; ++i) {
    tbox::Pointer<BaseMaterial<DIM> > tmp_material = d_vec_material[i];
    if (tmp_material->isMatching(material_name)) {
      material = tmp_material;
      found = true;
    }
  }

  if (!found) {
    TBOX_ERROR("fem material \"" << material_name << "\" is not found. "
                                 << endl);
  }

  if (material.isNull()) {
    TBOX_ERROR("fem material \"" << material_name << "\" is NULL. " << endl);
  }

  return (material);
}

template class MaterialManager<NDIM>;
