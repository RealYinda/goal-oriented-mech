//
// 文件名:     MaterialManager.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:33:18 2014 $
// 描述　:     材料管理器
//

#ifndef included_fem_MaterialManager
#define included_fem_MaterialManager

#include <vector>
#include "Pointer.h"

using namespace std;
using namespace JAUMIN;

template <int DIM>
class BaseMaterial;
/**
 * @brief 有限元材料管理类 adpt::MaterialManager, 为有限元模块管理材料。
 *
 * 该类维护一个数组, 存储用户注册到JAUMIN框架的材料。
 *
 * 该类提供如下函数, 支持用户从数组中获取指定名字的材料。 - getMaterial()
 *
 *
 * @see adpt::BaseMaterial
 */

template <int DIM>
class MaterialManager {
public:
  /**
   * @brief 构造函数。
   */
  MaterialManager();

  /**
   * @brief 虚析构函数。
   */
  virtual ~MaterialManager();

  /**
   * @brief  获取材料管理器.
   *
   * @return 指针, 材料管理器.
   */
  static tbox::Pointer<MaterialManager<DIM> > getManager();

  /**
   * @brief 将材料注册到材料管理器。
   *
   * @param material 输入参数, 指针, 指向添加的材料.
   */
  void addMaterial(tbox::Pointer<BaseMaterial<DIM> > material);

  /**
   * @brief 获取匹配于指定名称的材料.
   * @param material_name 输入参数, 字符串, 材料名字.
   * @return 指针, 指向匹配的材料.
   */
  tbox::Pointer<BaseMaterial<DIM> > getMaterial(
      const string &material_name) const;

private:
  static tbox::Pointer<MaterialManager<DIM> > s_material_manager_instance;

  /**
   * 存储材料的数组。
   */
  std::vector<tbox::Pointer<BaseMaterial<DIM> > > d_vec_material;
};

#endif
