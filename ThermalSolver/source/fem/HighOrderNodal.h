//
// 文件名:      HighOrderNodal.h
// 软件包:
// 版权  :      Zhejiang University
// HighOrder Solver of thermal problem
//
// 版本号:      $Revision$
// 修改  :      $Date$
// 描述  :
//
#ifndef included_appu_HighOrderNodal
#define included_appu_HighOrderNodal
#include "Pointer.h"
#include "Array.h"
#include "Patch.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "EdgeData.h"
#include "CellData.h"
#include "NodeData.h"
namespace JAUMIN {
namespace appu {
class HighOrderNodal{
public:
    HighOrderNodal(const hier::Patch<NDIM>& patch,
                     tbox::Pointer<pdat::EdgeData<NDIM, bool> > edge_order,
                     tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian);
    /**
           * Number of basis functions in a cell.
           */
    int nbas() const;

    /**
           * Dimension of a basis function.
           */
    int dim() const;

    /**
           * mysort function
           */
    static bool mysort(std::pair<int, double*> p1, std::pair<int, double*> p2);

    /**
           * Order of the nodes
           */
    void nodeOrder(const int face, int *order) const;

    /**
           * Basis function values.
           */
    void basis(const int cell,
               const double* lambda,
               double* values_in) const;
    /**
           * Curl of basis functions.
           */
    void gradient(const int cell,
              const double* lambda,
              double* values_in) const;
public:
    void basisnode1(const int cell, const double* lambda, double* values) const;
    void basisedge1(const int cell, const double* lambda, double* values) const;
    void gradientnode1(const int cell, const double* lambda, double* values) const;
    void gradientedge1(const int cell, const double* lambda, double* values) const;
public:
    const hier::Patch<NDIM>& d_patch;
    tbox::Pointer<pdat::EdgeData<NDIM, bool> > d_edge_order;
    tbox::Pointer<pdat::CellData<NDIM, double> > d_cell_jacobian;
    tbox::Pointer<pdat::NodeData<NDIM, double> > d_node_coord;
    tbox::Array<int> d_cell_edge_ext;
    tbox::Array<int> d_cell_edge_idx;
    tbox::Array<int> d_edge_node_ext;
    tbox::Array<int> d_edge_node_idx;
    tbox::Array<int> d_cell_node_ext;
    tbox::Array<int> d_cell_node_idx;
    tbox::Array<int> d_cell_face_ext;
    tbox::Array<int> d_cell_face_idx;
    tbox::Array<int> d_face_node_ext;
    tbox::Array<int> d_face_node_idx;
};
}
}
#endif

