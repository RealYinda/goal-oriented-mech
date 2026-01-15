#include "HighOrderNodal.h"
#include "GridInfo.h"
namespace JAUMIN {
namespace appu{

HighOrderNodal::HighOrderNodal(const hier::Patch<NDIM>& patch,
                                   tbox::Pointer<pdat::EdgeData<NDIM, bool> > edge_order,
                                   tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian) :
    d_patch(patch),
    d_edge_order(edge_order),
    d_cell_jacobian(cell_jacobian)
{
    patch.getPatchTopology()->getCellAdjacencyEdges(d_cell_edge_ext, d_cell_edge_idx);
    patch.getPatchTopology()->getEdgeAdjacencyNodes(d_edge_node_ext, d_edge_node_idx);
    patch.getPatchTopology()->getCellAdjacencyNodes(d_cell_node_ext, d_cell_node_idx);
    patch.getPatchTopology()->getFaceAdjacencyNodes(d_face_node_ext, d_face_node_idx);
    patch.getPatchTopology()->getCellAdjacencyFaces(d_cell_face_ext, d_cell_face_idx);
}
int HighOrderNodal::nbas() const{
    return NDIM == 2? 6 : 10;
}
int HighOrderNodal::dim() const
    { return NDIM; }
void HighOrderNodal::basis(const int cell, const double* lambda, double* values_in) const{
    int nedge = 6; int nnode = 4;
    double *values;
    values = values_in;
    basisnode1(cell,lambda,values);
    values = values + nnode;
    basisedge1(cell,lambda,values);
}
void HighOrderNodal::gradient(const int cell, const double* lambda, double* values_in) const{
    int nedge = 6; int nnode = 4;
    double *values;
    values = values_in;
    gradientnode1(cell,lambda,values);
    values = values + nnode*NDIM;
    gradientedge1(cell,lambda,values);
}
void HighOrderNodal::basisnode1(const int cell, const double* lambda, double* values_in) const{
    double *phi = values_in;
    phi[0] = (2*lambda[0]-1)*lambda[0];
    phi[1] = (2*lambda[1]-1)*lambda[1];
    phi[2] = (2*lambda[2]-1)*lambda[2];
    phi[3] = (2*lambda[3]-1)*lambda[3];
}
void HighOrderNodal::basisedge1(const int cell, const double* lambda, double* values_in) const{
    int nedge = 6;
    double(*nabla)[NDIM + 1] =
            (double(*)[NDIM + 1])(&((*d_cell_jacobian)(0, cell)));
    double *phi = values_in;
    for (int ie = 0; ie < nedge; ie++) {
        int edge = d_cell_edge_idx[d_cell_edge_ext[cell] + ie];//单元的第几条边
        int v0 = d_edge_node_idx[d_edge_node_ext[edge] + 0];//边的一个节点的全局编号
        int v1 = d_edge_node_idx[d_edge_node_ext[edge] + 1];//边的一个节点的全局编号
        bool in_order = (*d_edge_order)(0, edge);
        if (!in_order) {
            int t = v0;
            v0 = v1;
            v1 = t;//交换顺序
        }

        int lv0, lv1;
        lv0 = lv1 = -1;
        //确定边的一个节点的局部编号
        for (int n = 0; n < NDIM + 1; n++) {
            int nidx = d_cell_node_idx[d_cell_node_ext[cell] + n];//单元的节点
            if (nidx == v0) lv0 = n;
            if (nidx == v1) lv1 = n;
        }
        phi[ie] = 4*lambda[lv0]*lambda[lv1];
    }
}
void HighOrderNodal::gradientnode1(const int cell, const double* lambda, double* values_in) const{
    double(*grad)[NDIM] = (double(*)[NDIM]) values_in;
    double(*nabla)[NDIM + 1] =
        (double(*)[NDIM + 1])(&((*d_cell_jacobian)(0, cell)));
    int nnode = 4;
    for (int in = 0; in < nnode; in++) {
        for(int d = 0; d < NDIM; d ++){
            grad[in][d] = nabla[in][d];
        }
    }

}
void HighOrderNodal::gradientedge1(const int cell, const double* lambda, double* values_in) const{
    double(*grad)[NDIM] = (double(*)[NDIM]) values_in;
    double(*nabla)[NDIM + 1] =
        (double(*)[NDIM + 1])(&((*d_cell_jacobian)(0, cell)));
    int nedge =6;
    for (int ie = 0; ie < nedge; ie++) {
        int edge = d_cell_edge_idx[d_cell_edge_ext[cell] + ie];//单元的第几条边
        int v0 = d_edge_node_idx[d_edge_node_ext[edge] + 0];//边的一个节点的全局编号
        int v1 = d_edge_node_idx[d_edge_node_ext[edge] + 1];//边的一个节点的全局编号
        bool in_order = (*d_edge_order)(0, edge);
        if (!in_order) {
            int t = v0;
            v0 = v1;
            v1 = t;//交换顺序
        }

        int lv0, lv1;
        lv0 = lv1 = -1;
        //确定边的一个节点的局部编号
        for (int n = 0; n < NDIM + 1; n++) {
            int nidx = d_cell_node_idx[d_cell_node_ext[cell] + n];//单元的节点
            if (nidx == v0) lv0 = n;
            if (nidx == v1) lv1 = n;
        }
        for(int d =0 ; d < NDIM; d++)
            grad[ie][d] = 4*(lambda[lv0]*nabla[lv1][d]+lambda[lv1]*nabla[lv0][d]);

    }

}
}
}
