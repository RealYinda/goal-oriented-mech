//
// 文件名:      TetQuad.h
// 软件包:
// 版权  :      (c) 2004-2015 北京应用物理与计算数学研究所
//              (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:      $Revision$
// 修改  :      $Date$
// 描述  :
//

#ifndef included_appu_TetQuad
#define included_appu_TetQuad

#if NDIM == 2
#error FOR 3D CASE ONLY
#endif

#include <algorithm>
#include <strings.h>
#include "Pointer.h"
#include "Array.h"
#include "Patch.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "EdgeData.h"
#include "CellData.h"
#include "NodeData.h"
#include "GridInfo.h"
#include "JAUMIN_Macros.h"

namespace JAUMIN {
namespace appu {

class TetQuad
{
public:
    TetQuad(const hier::Patch<NDIM>& patch,
            tbox::Pointer<pdat::CellData<NDIM, double> > cell_volume,
        tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian);

    template<class ShapeFuncSolution, class ShapeFuncField,
             class ShapeFuncCurl, class Function>
    void ComputeGoalOrientedErrorCell
    (const int cell, const ShapeFuncSolution *shapefunc1,
     const ShapeFuncField *shapefunc2, const ShapeFuncCurl *shapefunc3,
     const int order, dcomplex *value, Function func,
     dcomplex *dofs, dcomplex *a_dofs,
     dcomplex *rcurl_dofs, dcomplex *recovery_dofs,
     dcomplex *coef_list){
        /// coef_list: k0^2,$\varepsilor_r$, $\mu_r$
        check(order);
        /// The basis number of
        /// Nedelec: 6, BDM: 12
        const Quad *quad = d_quad_table[order];
        double vol = (*d_cell_volume)(0, cell);

        for (int n = 0; n < quad->npoints; n++){
            dcomplex sum=0.;
            dcomplex e_value[NDIM];
            dcomplex e_curl[NDIM];
            dcomplex v_value[NDIM];
            dcomplex v_curl[NDIM];
            dcomplex vr_value[NDIM];
            dcomplex vr_curl[NDIM];
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               quad->points + n *(NDIM + 1),
                               space_coords);
            dcomplex funcvalues[NDIM];
            func(space_coords[0],space_coords[1],space_coords[2], funcvalues);
            double Nedelec_phi[6][NDIM];
            double Nedelec_curl[6][NDIM];
            shapefunc1->basis(cell, quad->points + n * (NDIM + 1), &(Nedelec_phi[0][0]));
            shapefunc1->curl(cell, quad->points + n * (NDIM + 1), &(Nedelec_curl[0][0]));
            double BDM_phi[12][NDIM];
            shapefunc2->basis(cell, quad->points + n * (NDIM + 1), &(BDM_phi[0][0]));
            double Recovered_curl[6][NDIM];
            shapefunc3->basis(cell, quad->points + n * (NDIM + 1), &(Recovered_curl[0][0]));
            for (int edge = 0; edge < 6; edge ++) {
                e_value[0] += dofs[edge] * Nedelec_phi[edge][0];
                e_value[1] += dofs[edge] * Nedelec_phi[edge][1];
                e_value[2] += dofs[edge] * Nedelec_phi[edge][2];
                v_value[0] += a_dofs[edge] * Nedelec_phi[edge][0];
                v_value[1] += a_dofs[edge] * Nedelec_phi[edge][1];
                v_value[2] += a_dofs[edge] * Nedelec_phi[edge][2];

                e_curl[0] += dofs[edge] * Nedelec_curl[edge][0];
                e_curl[1] += dofs[edge] * Nedelec_curl[edge][1];
                e_curl[2] += dofs[edge] * Nedelec_curl[edge][2];
                v_curl[0] += a_dofs[edge] * Nedelec_curl[edge][0];
                v_curl[1] += a_dofs[edge] * Nedelec_curl[edge][1];
                v_curl[2] += a_dofs[edge] * Nedelec_curl[edge][2];

                vr_curl[0] += rcurl_dofs[edge] * Recovered_curl[edge][0];
                vr_curl[1] += rcurl_dofs[edge] * Recovered_curl[edge][1];
                vr_curl[2] += rcurl_dofs[edge] * Recovered_curl[edge][2];
            }
            for (int face = 0; face < 4; face ++) {
                for (int node = 0; node < 3; node ++){
                    vr_value[0] += recovery_dofs[face*3+node]*BDM_phi[face*3+node][0];
                    vr_value[1] += recovery_dofs[face*3+node]*BDM_phi[face*3+node][1];
                    vr_value[2] += recovery_dofs[face*3+node]*BDM_phi[face*3+node][2];
                }
            }
            /// b(v-v_h)
//            sum += M_PI*M_PI*DOT_PRODUCT(funcvalues,vr_value);
//            sum -= M_PI*M_PI*DOT_PRODUCT(funcvalues,v_value);
            /// curl_curl(v-v_h)
            sum -= DOT_PRODUCT(e_curl,vr_curl);
            sum += DOT_PRODUCT(e_curl,v_curl);
            /// bas_bas(v-v_h)
            sum += coef_list[1]*
                    DOT_PRODUCT(e_value,vr_value);
            sum -= coef_list[1]*
                    DOT_PRODUCT(e_value,v_value);
            sum *= quad->weights[n];
            value[0] += sum;

        }
        value[0] *= vol;


    }
    template<class ShapeFuncSolution, class ShapeFuncField,
             class Function>
    void ComputeHighOrderGoalOrientedErrorCell
    (const int cell, const ShapeFuncSolution *shapefunc1,
     const ShapeFuncField *shapefunc2,
     const int order, dcomplex *value, Function func,
     dcomplex *dofs, dcomplex *adjoint_dofs,
     dcomplex *coef_list){
        /// coef_list: k0^2,$\varepsilor_r$, $\mu_r$
        check(order);
        /// The basis number of
        /// Nedelec: 6, BDM: 12
        const Quad *quad = d_quad_table[order];
        double vol = (*d_cell_volume)(0, cell);

        for (int n = 0; n < quad->npoints; n++){
            dcomplex sum=0.;
            dcomplex e_value[NDIM];
            dcomplex e_curl[NDIM];
            dcomplex v_value[NDIM];
            dcomplex v_curl[NDIM];
            dcomplex vh_value[NDIM];
            dcomplex vh_curl[NDIM];
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               quad->points + n *(NDIM + 1),
                               space_coords);
            dcomplex funcvalues[NDIM];
            func(space_coords[0],space_coords[1],space_coords[2], funcvalues);
            double Nedelec_phi[6][NDIM];
            double Nedelec_curl[6][NDIM];
            shapefunc1->basis(cell, quad->points + n * (NDIM + 1), &(Nedelec_phi[0][0]));
            shapefunc1->curl(cell, quad->points + n * (NDIM + 1), &(Nedelec_curl[0][0]));
            double HighNedelec_phi[20][NDIM];
            double HighNedelec_curl[20][NDIM];
            shapefunc2->basis(cell, quad->points + n * (NDIM + 1), &(HighNedelec_phi[0][0]));
            shapefunc2->curl(cell, quad->points + n * (NDIM + 1), &(HighNedelec_curl[0][0]));
            for (int edge = 0; edge < 6; edge ++) {
                e_value[0] += dofs[edge] * Nedelec_phi[edge][0];
                e_value[1] += dofs[edge] * Nedelec_phi[edge][1];
                e_value[2] += dofs[edge] * Nedelec_phi[edge][2];
                v_value[0] += adjoint_dofs[2*edge] * Nedelec_phi[edge][0];
                v_value[1] += adjoint_dofs[2*edge] * Nedelec_phi[edge][1];
                v_value[2] += adjoint_dofs[2*edge] * Nedelec_phi[edge][2];

                e_curl[0] += dofs[edge] * Nedelec_curl[edge][0];
                e_curl[1] += dofs[edge] * Nedelec_curl[edge][1];
                e_curl[2] += dofs[edge] * Nedelec_curl[edge][2];
                v_curl[0] += adjoint_dofs[2*edge] * Nedelec_curl[edge][0];
                v_curl[1] += adjoint_dofs[2*edge] * Nedelec_curl[edge][1];
                v_curl[2] += adjoint_dofs[2*edge] * Nedelec_curl[edge][2];

                vh_value[0] += adjoint_dofs[2*edge] * HighNedelec_phi[edge][0];
                vh_value[1] += adjoint_dofs[2*edge] * HighNedelec_phi[edge][1];
                vh_value[2] += adjoint_dofs[2*edge] * HighNedelec_phi[edge][2];
                vh_value[0] += adjoint_dofs[2*edge+1] * HighNedelec_phi[6+edge][0];
                vh_value[1] += adjoint_dofs[2*edge+1] * HighNedelec_phi[6+edge][1];
                vh_value[2] += adjoint_dofs[2*edge+1] * HighNedelec_phi[6+edge][2];

                vh_curl[0] += adjoint_dofs[2*edge] * HighNedelec_curl[edge][0];
                vh_curl[1] += adjoint_dofs[2*edge] * HighNedelec_curl[edge][1];
                vh_curl[2] += adjoint_dofs[2*edge] * HighNedelec_curl[edge][2];
                vh_curl[0] += adjoint_dofs[2*edge+1] * HighNedelec_curl[6+edge][0];
                vh_curl[1] += adjoint_dofs[2*edge+1] * HighNedelec_curl[6+edge][1];
                vh_curl[2] += adjoint_dofs[2*edge+1] * HighNedelec_curl[6+edge][2];
            }
            for (int face = 0; face < 4; face ++) {
                vh_value[0] += adjoint_dofs[12+2*face] * HighNedelec_phi[12+face][0];
                vh_value[1] += adjoint_dofs[12+2*face] * HighNedelec_phi[12+face][1];
                vh_value[2] += adjoint_dofs[12+2*face] * HighNedelec_phi[12+face][2];
                vh_value[0] += adjoint_dofs[12+2*face+1] * HighNedelec_phi[16+face][0];
                vh_value[1] += adjoint_dofs[12+2*face+1] * HighNedelec_phi[16+face][1];
                vh_value[2] += adjoint_dofs[12+2*face+1] * HighNedelec_phi[16+face][2];

                vh_curl[0] += adjoint_dofs[12+2*face] * HighNedelec_curl[12+face][0];
                vh_curl[1] += adjoint_dofs[12+2*face] * HighNedelec_curl[12+face][1];
                vh_curl[2] += adjoint_dofs[12+2*face] * HighNedelec_curl[12+face][2];
                vh_curl[0] += adjoint_dofs[12+2*face+1] * HighNedelec_curl[16+face][0];
                vh_curl[1] += adjoint_dofs[12+2*face+1] * HighNedelec_curl[16+face][1];
                vh_curl[2] += adjoint_dofs[12+2*face+1] * HighNedelec_curl[16+face][2];
            }
            /// b(v-v_h)
//            sum += M_PI*M_PI*DOT_PRODUCT(funcvalues,vr_value);
//            sum -= M_PI*M_PI*DOT_PRODUCT(funcvalues,v_value);
            /// curl_curl(v-v_h)
            sum -= DOT_PRODUCT(e_curl,vh_curl);
            sum += DOT_PRODUCT(e_curl,v_curl);
            /// bas_bas(v-v_h)
            sum += coef_list[1]*
                    DOT_PRODUCT(e_value,vh_value);
            sum -= coef_list[1]*
                    DOT_PRODUCT(e_value,v_value);
            sum *= quad->weights[n];
            value[0] += sum;

        }
        value[0] *= vol;


    }

    template <class ShapeFunc>
    void quadBasDotBas(const int cell, const ShapeFunc *shapefunc,
                       const int order, dcomplex *value) const {
      check(order);
      int nbas = shapefunc->nbas();
      int dim = shapefunc->dim();
      dcomplex(*bb)[nbas] = (dcomplex(*)[nbas])value;
      const Quad *quad = d_quad_table[order];

      double vol = (*d_cell_volume)(0, cell);
      bzero(&(bb[0][0]), nbas * nbas * sizeof(*value));

      for (int n = 0; n < quad->npoints; n++) {
        double phi[nbas][dim];
        shapefunc->basis(cell, quad->points + n * (NDIM + 1), &(phi[0][0]));
        for (int i = 0; i < nbas; i++) {
          for (int j = 0; j < nbas; j++) {
            dcomplex v = 0.;
            for (int d = 0; d < dim; d++) v += phi[i][d]*phi[j][d];
            v *= quad->weights[n];
            bb[i][j] += v;
          }
        }
      }
      for (int i = 0; i < nbas; i++)
        for (int j = 0; j < nbas; j++) bb[i][j] *= vol;
    }    

    template<class ShapeFunc>
    void faceQuadBasDotBas(const int cell,
                           const int face,
                   const ShapeFunc* shapefunc,
                   const int order,
               double *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        double (*bb)[nbas] = (double (*)[nbas])value;
        const Quad *quad = d_face_quad_table[order];
        assert(face < NDIM+1 && face >= 0);


        //double vol = (*d_cell_volume)(0, cell);
        bzero(&(bb[0][0]), nbas * nbas * sizeof(*value));

        double outer_normal[NDIM], outer_normal_raw[NDIM];
        // get the id of face in this patch.
        int patch_face = caf_idx[caf_ext[cell] + face];
        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[patch_face+1]-fan_ext[patch_face],
                fan_idx.getPointer() + fan_ext[patch_face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);
        double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw)) / 2.0;


        for(int n = 0; n < quad->npoints; n++) {
            // transform 2d barycentic to 3d
            double bary3d[NDIM+1];
            bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,face,bary3d);
            // get the shape functions' value
            double phi[nbas][dim];
            shapefunc->basis(cell,
                             bary3d,
                             &(phi[0][0]));

            // compute the stiff matrix of face 'face' of element 'cell'
            for(int i = 0; i < nbas; i++) {
                for(int j = 0; j < nbas; j++) {
                    double bas_i[NDIM], bas_j[NDIM];
                    crossProduct(NDIM, outer_normal, phi[i], bas_i);
                    crossProduct(NDIM, outer_normal, phi[j], bas_j);
                    double v = dotProduct(NDIM, bas_i, bas_j);
                    v *= quad->weights[n];
                    bb[i][j] += v;
                }
            }
        }
        for(int i = 0; i < nbas; i++)
            for(int j = 0; j < nbas; j++)
                bb[i][j] *= area;
    }
    template<class ShapeFunc>
    void faceQuadBasDotCurlNormal(const int cell,
                           const int face,
                   const ShapeFunc* shapefunc,
                   const int order,
               double *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        double (*bb)[nbas] = (double (*)[nbas])value;
        const Quad *quad = d_face_quad_table[order];
        assert(face < NDIM+1 && face >= 0);


        //double vol = (*d_cell_volume)(0, cell);
        bzero(&(bb[0][0]), nbas * nbas * sizeof(*value));

        double outer_normal[NDIM], outer_normal_raw[NDIM];
        // get the id of face in this patch.
        int patch_face = caf_idx[caf_ext[cell] + face];
        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[patch_face+1]-fan_ext[patch_face],
                fan_idx.getPointer() + fan_ext[patch_face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);
        double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw)) / 2.0;


        for(int n = 0; n < quad->npoints; n++) {
            // transform 2d barycentic to 3d
            double bary3d[NDIM+1];
            bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,face,bary3d);
            // get the shape functions' value
            double phi[nbas][dim];
            double curlphi[nbas][dim];
            shapefunc->basis(cell,
                             bary3d,
                             &(phi[0][0]));
            shapefunc->curl(cell,
                             bary3d,
                             &(curlphi[0][0]));
            // compute the stiff matrix of face 'face' of element 'cell'
            for(int i = 0; i < nbas; i++) {
                for(int j = 0; j < nbas; j++) {
                    double crossNj[NDIM];
                    crossProduct(NDIM, outer_normal, phi[j], crossNj);
                    double v = dotProduct(NDIM, crossNj, curlphi[i]);
                    v *= quad->weights[n];
                    bb[i][j] += v;
                }
            }
        }
        for(int i = 0; i < nbas; i++)
            for(int j = 0; j < nbas; j++)
                bb[i][j] *= area;
    }

    template<class ShapeFunc, class TYPE, class Function>
    void faceQuadFunctionDotBas(const int cell,
                                const int face,
                                const Function& func,
                                const ShapeFunc* shapefunc,
                                const int order,
                                TYPE *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        TYPE *bb = value;
        const Quad *quad = d_face_quad_table[order];
        assert(face < NDIM+1 && face >= 0);

        //double vol = (*d_cell_volume)(0, cell);
        bzero(&(bb[0]), nbas * sizeof(*value));

        double outer_normal[NDIM], outer_normal_raw[NDIM];
        // get the id of face in this patch.
        int patch_face = caf_idx[caf_ext[cell] + face];
        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[patch_face+1]-fan_ext[patch_face],
                fan_idx.getPointer() + fan_ext[patch_face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);
        double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw))/ 2.0;

        for(int n = 0; n < quad->npoints; n++) {
            // transform this 2d barycentric quadrature point to 3d barycentic on face 'face'
            double bary3d[NDIM+1];
            bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,face,bary3d);
            // get the space coordinate of this quadrature point
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               bary3d,
                               space_coords);
            TYPE funcvalues[NDIM];
            func(space_coords[0],space_coords[1],space_coords[2], funcvalues);

            // get the 'nbas' shape functions' value at this quadrature point
            double phi[nbas][dim];
            shapefunc->basis(cell,
                             bary3d,
                             &(phi[0][0]));

            // compute the load of this element
            for(int i = 0; i < nbas; i++) {
                dcomplex bas_i[dim],bas_value[dim];
                CROSS_PRODUCT(outer_normal, phi[i], bas_value);
                CROSS_PRODUCT(outer_normal, bas_value, bas_i);
                dcomplex v = dotProduct(NDIM, bas_i, funcvalues);
                v *= quad->weights[n];
                bb[i] += v;
            }
        }
        for(int i = 0; i < nbas; i++)
                bb[i] *= area;
    }

    /// When impose the Dirichlet boundary, this should be set strong
    template<class ShapeFunc, class TYPE, class Function>
    void faceQuadCloseFunctionDotCurl(const int cell,
                                const int face,
                                const Function& func,
                                const ShapeFunc* shapefunc,
                                const int order,
                                TYPE *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        TYPE *bb = value;
        const Quad *quad = d_face_quad_table[order];
        assert(face < NDIM+1 && face >= 0);

        //double vol = (*d_cell_volume)(0, cell);
        bzero(&(bb[0]), nbas * sizeof(*value));

        double outer_normal[NDIM], outer_normal_raw[NDIM];
        // get the id of face in this patch.
        int patch_face = caf_idx[caf_ext[cell] + face];
        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[patch_face+1]-fan_ext[patch_face],
                fan_idx.getPointer() + fan_ext[patch_face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);
        double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw))/ 2.0;

        for(int n = 0; n < quad->npoints; n++) {
            // transform this 2d barycentric quadrature point to 3d barycentic on face 'face'
            double bary3d[NDIM+1];
            bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,face,bary3d);
            // get the space coordinate of this quadrature point
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               bary3d,
                               space_coords);
            TYPE funcvalues[NDIM];
            func(space_coords[0],space_coords[1],space_coords[2], funcvalues);

            // get the 'nbas' shape functions' value at this quadrature point
            double phi[nbas][dim];
            shapefunc->basis(cell,
                             bary3d,
                             &(phi[0][0]));

            // compute the load of this element
            for(int i = 0; i < nbas; i++) {
                dcomplex bas_value[dim],curlvalue[dim];
                dcomplex temp[dim];
                CROSS_PRODUCT(outer_normal,phi[i],temp);
                CROSS_PRODUCT(outer_normal,temp,bas_value);
                CROSS_PRODUCT(outer_normal,funcvalues, curlvalue);

                TYPE v = DOT_PRODUCT(bas_value, curlvalue);
                v *= quad->weights[n];
                bb[i] += v;
            }
        }
        for(int i = 0; i < nbas; i++)
                bb[i] *= area;
    }
    template<class ShapeFunc, class TYPE, class Function>
    void faceQuadCloseFunctionDotSolution(const int cell,
                                const int face,
                                const Function& func,
                                const ShapeFunc* shapefunc,
                                const int order,
                                TYPE *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        TYPE *bb = value;
        const Quad *quad = d_face_quad_table[order];
        assert(face < NDIM+1 && face >= 0);

        //double vol = (*d_cell_volume)(0, cell);
        bzero(&(bb[0]), nbas * sizeof(*value));

        double outer_normal[NDIM], outer_normal_raw[NDIM];
        // get the id of face in this patch.
        int patch_face = caf_idx[caf_ext[cell] + face];
        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[patch_face+1]-fan_ext[patch_face],
                fan_idx.getPointer() + fan_ext[patch_face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);
        double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw))/ 2.0;

        for(int n = 0; n < quad->npoints; n++) {
            // transform this 2d barycentric quadrature point to 3d barycentic on face 'face'
            double bary3d[NDIM+1];
            bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,face,bary3d);
            // get the space coordinate of this quadrature point
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               bary3d,
                               space_coords);
            TYPE funcvalues[NDIM];
            func(space_coords[0],space_coords[1],space_coords[2], funcvalues);

            // get the 'nbas' shape functions' value at this quadrature point
            double phi[nbas][dim];
            shapefunc->basis(cell,
                             bary3d,
                             &(phi[0][0]));

            // compute the load of this element
            for(int i = 0; i < nbas; i++) {
                dcomplex bas_value[dim],field_value[dim];
                dcomplex temp[dim];
                CROSS_PRODUCT(outer_normal,phi[i],temp);
                CROSS_PRODUCT(outer_normal,temp,bas_value);
                CROSS_PRODUCT(outer_normal,funcvalues,temp);
                CROSS_PRODUCT(outer_normal,temp,field_value);
                TYPE v = DOT_PRODUCT(bas_value, field_value);
                v *= quad->weights[n];
                bb[i] += v;
            }
        }
        for(int i = 0; i < nbas; i++)
                bb[i] *= area;
    }
    template<class ShapeFunc>
    void testJumpOnFace(const int face,
                        const ShapeFunc* shapefunc,
                        const int order,
                        dcomplex *dofs,
                        dcomplex *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        const Quad *quad = d_face_quad_table[order];
        int cell = fac_idx[fac_ext[face]];
        double outer_normal[NDIM], outer_normal_raw[NDIM];
        // get the id of face in this patch.

        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[face+1]-fan_ext[face],
                fan_idx.getPointer() + fan_ext[face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);

        double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw))/ 2.0;
        int ll_face = 9999;
        for(int ll =0; ll <4; ll ++){
            if(caf_idx[caf_ext[cell]+ll] == face) ll_face = ll;
        }
        assert(ll_face < NDIM+1 && ll_face >= 0);
        for(int n = 0; n < quad->npoints; n++) {
            // transform this 2d barycentric quadrature point to 3d barycentic on face 'face'
            double bary3d[NDIM+1];
            bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,ll_face,bary3d);
            dcomplex gauss_value;
            // get the space coordinate of this quadrature point
            double phi[12][3];
            shapefunc->basis(cell, bary3d, &(phi[0][0]));
            for(int local_face= 0; local_face < 4; local_face ++){
                for(int local_node =0; local_node < 3; local_node ++){
                    dcomplex v = DOT_PRODUCT(phi[local_face*3+local_node],outer_normal);
                    v *= dofs[local_face*3+local_node];
                    gauss_value += v;
                }
            }
            (*value) += gauss_value * quad->weights[n];

        }

    }
    template<class ShapeFunc>
    void computeJumpOnFace(const int face,
                           const ShapeFunc* shapefunc,
                           const int order,
                           dcomplex *dofs,
                           dcomplex *value) const
    {
        /*
         * Vinta Yin-Da Wang
         * 2023-03-23
         * Please notice that this function only used for interior boundaries
         * */
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        const Quad *quad = d_face_quad_table[order];
        int cell = fac_idx[fac_ext[face]];
        double outer_normal[NDIM], outer_normal_raw[NDIM];
        // get the id of face in this patch.

        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[face+1]-fan_ext[face],
                fan_idx.getPointer() + fan_ext[face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);
        // The outer_normal is fixed by the cells

        for(int local_cell = 0; local_cell < 2; local_cell ++){
            int global_cell = fac_idx[fac_ext[face]+local_cell];
            int ll_face = 9999;
            for(int ll =0; ll <4; ll ++){
                if(caf_idx[caf_ext[global_cell]+ll] == face) ll_face = ll;
            }
            assert(ll_face < NDIM+1 && ll_face >= 0);
            for(int n = 0; n < quad->npoints; n++){
                double bary3d[NDIM+1];
                bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,ll_face,bary3d);
                dcomplex gauss_value;
                // get the space coordinate of this quadrature point
                double phi[6][dim];
                shapefunc->basis(cell, bary3d, &(phi[0][0]));
                for (int local_edge = 0;local_edge < 6;local_edge ++) {
                    dcomplex v = DOT_PRODUCT(phi[local_edge],outer_normal);
                    v *= 0.5*dofs[6*local_cell+local_edge];
                    gauss_value += v;
                }
                (*value) += gauss_value * quad->weights[n];
            }
        }
    }
    template<class ShapeFunc,class TYPE>
    void faceQuadDOFcurldotEdge(const int cell,
                            const int face,
                            double* edge_dir,
                            const ShapeFunc* shapefunc,
                            const int order,
                            TYPE *dofs,
                            TYPE *value) const


    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        const Quad *quad = d_face_quad_table[order];
        double outer_normal[NDIM], outer_normal_raw[NDIM];
        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[face+1]-fan_ext[face],
                fan_idx.getPointer() + fan_ext[face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);
        double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw))/ 2.0;
        int local_face=0;
        for(int lf =0; lf < 4; lf ++){
            if(caf_idx[caf_ext[cell]+lf] == face) local_face = lf;
        }
        /**
         * @brief sigmadotte: \sigma * t_e
         * After multiply with area, turning into the integration
         */
        dcomplex sigmadotte(0.,0.);
        for(int n = 0; n < quad->npoints; n++){
            // interior boundary
            double bary3d[NDIM+1];
            bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,local_face,bary3d);
            double curlphi[nbas][NDIM];
            shapefunc->curl(cell,bary3d,&(curlphi[0][0]));
            // basis is arranged in cell->edge local order
            dcomplex dot_gauss(0.,0.);
            for(int bas = 0; bas < nbas; bas ++){
                // dot product in basis
                dcomplex basdot = DOT_PRODUCT(curlphi[bas],edge_dir);
                // multiplex with DOF
                basdot *= dofs[bas];
                // Accumulation on that quad
                dot_gauss += basdot;
            }
            // Weighting on that quad
            dot_gauss *= quad->weights[n];
            sigmadotte += dot_gauss;

        }
        (*value) = sigmadotte * area;

    }
    template<class ShapeFunc,class TYPE>
    void faceQuadDOFbasisdotFace(const int face,
                                 double* outer_normal_fac,
                                 const ShapeFunc* shapefunc,
                                 const int order,
                                 TYPE *dofs,
                                 TYPE *value,
                                 dcomplex *rate) const


    {
        check(order);
        /**
         * Vinta Yin-Da Wang
         * 2023-02-16
         * To check the local index of node
         * lnode should end up the local coeffecient of the control node
         **/
        int nbas = shapefunc->nbas();
        const Quad *quad = d_face_quad_table[order];  
        /**
         * Vinta Yin-Da Wang
         * barycentric: The local face mapping to local volume mapping
         * Because phi[nbas][NDIM] is decided in cells
         * While lambda[3] is decided in faces
         */
        dcomplex jump_value[3] = {0.,0.,0.};
        for (int this_cell = fac_ext[face];this_cell < fac_ext[face+1]; this_cell++) {
            int cell_idx = fac_idx[this_cell];
            int local_cell_num = this_cell - fac_ext[face];
            /// Get the local node number of cell
            int fc_node[3] = {-1, -1, -1};
            for(int inode =0; inode < 4; inode ++){
                if(can_idx[can_ext[cell_idx]+inode] == fan_idx[fan_ext[face]+0])
                    fc_node[0] = inode;
                if(can_idx[can_ext[cell_idx]+inode] == fan_idx[fan_ext[face]+1])
                    fc_node[1] = inode;
                if(can_idx[can_ext[cell_idx]+inode] == fan_idx[fan_ext[face]+2])
                    fc_node[2] = inode;
                // The corresponding element of fc_nodes is the local node of the cell
                // This element can be a index of bary3d
            }
            int local_face=-1;
            for(int lf = 0; lf < 4; lf ++)
                if(caf_idx[caf_ext[cell_idx]+lf] == face) local_face = lf;
            // local_face is the cell local number of the controlled face
            assert(local_face >=0 && local_face < NDIM+1);
            dcomplex cell_related_coef[3] = {dcomplex(0.,0.),dcomplex(0.,0.),dcomplex(0.,0.)};
            for (int n = 0; n < quad->npoints; n++){
                dcomplex vector_jump_on_face[3] = {dcomplex(0.,0.),dcomplex(0.,0.),dcomplex(0.,0.)};
                double bary3d[NDIM+1];
                bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell_idx,local_face,bary3d);
                double phi[nbas][NDIM];
                // Nedelec basis function here
                shapefunc->basis(cell_idx,bary3d,&(phi[0][0]));
                dcomplex q_normal_component(0.,0.);
                /// Compute tau_Fn1 and tau_Fn2
                for(int bas =0; bas < nbas; bas ++){
                    double basdot = DOT_PRODUCT(phi[bas],outer_normal_fac);
                    q_normal_component += basdot * dofs[nbas*local_cell_num+bas];
                }
                for(int f_n =0; f_n < 3; f_n ++){
                    vector_jump_on_face[f_n] = *(rate+local_cell_num)*q_normal_component*bary3d[fc_node[f_n]]*quad->weights[n];
                    cell_related_coef[f_n] += vector_jump_on_face[f_n];
                }

            }
            for(int f_n =0; f_n < 3; f_n ++){
//                cout<<"The Cell "<<local_cell_num<<"jump for node "<<f_n<<
//                      "is equal to"<<cell_related_coef[f_n]<<endl;
                value[f_n] += cell_related_coef[f_n];
            }
        }

        /// Quad first, cell second

/*
        for (int n = 0; n < quad->npoints; n++) {
            dcomplex q_jump_on_face[3]={dcomplex(0.,0.),dcomplex(0.,0.),dcomplex(0.,0.)};

            for (int cc = fac_ext[face];cc < fac_ext[face+1]; cc++) {
                int cell = fac_idx[cc];
                int fc_node[3] = {999,999,999};
                /// The local number of nodes in cell
                for(int inode =0; inode < 4; inode ++){
                    if(can_idx[can_ext[cell]+inode] == fan_idx[fan_ext[face]+0])
                        fc_node[0] = inode;
                    if(can_idx[can_ext[cell]+inode] == fan_idx[fan_ext[face]+1])
                        fc_node[1] = inode;
                    if(can_idx[can_ext[cell]+inode] == fan_idx[fan_ext[face]+2])
                        fc_node[2] = inode;
                    // The corresponding element of fc_nodes is the local node of the cell
                    // This element can be a index of bary3d
                }
                int local_face=9999;
                for(int lf = 0; lf < 4; lf ++)
                    if(caf_idx[caf_ext[cell]+lf] == face) local_face = lf;
                // local_face is the cell local number of the controlled face
                assert(local_face >=0 && local_face < NDIM+1);


                double bary3d[NDIM+1];
                bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,local_face,bary3d);
                double phi[nbas][NDIM];
                // Nedelec basis function here
                shapefunc->basis(cell,bary3d,&(phi[0][0]));
                dcomplex q_normal_component(0.,0.);
                int cc_local = cc-fac_ext[face];
                for(int bas =0; bas < nbas; bas ++){
                    double basdot = DOT_PRODUCT(phi[bas],outer_normal_fac);
                    q_normal_component += basdot * dofs[nbas*cc_local+bas];
//                    cout<<" With basdot "<<basdot<<" locally"<<endl;
                }
//                cout<<q_normal_component<<" on cell "<<cc_local<<" quad point "<<n<<endl;

                for(int f_n =0; f_n < 3; f_n ++){
                    q_jump_on_face[f_n] += *(rate+cc_local)*q_normal_component*bary3d[fc_node[f_n]]*quad->weights[n];
                }
                value[0] += q_jump_on_face[0];
                value[1] += q_jump_on_face[1];
                value[2] += q_jump_on_face[2];
            }


        }
*/



        /*
         * Vinta Yin-Da Wang
         * 2023-03-18
         * local_face is the local number of the control face
         * fc_node is the cell local number in face local order
         * */

    }
    template<class ShapeFunc,class triShapeFunc, class TYPE, class TYPE2>
    void faceQuadE1minusE2DotE2_modal(const int cell,
                             const int face,
                             const TYPE2 *E_edge_tol,//总场边上插值系数
                             const ShapeFunc* shapefunc,//体形函数
                             const TYPE2 *E_edge_modal,//模式场边上插值系数
                             const TYPE2 *E_node_modal,//模式场节点上插值系数
                             const triShapeFunc* trishapefunc,//面形函数
                             tbox::Array<tbox::Array<double> > nodeGradBas,
                             tbox::Array<double> direction,
                                   int faceorder,
                             const int order,
                             TYPE *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        TYPE *bb = value;
        const Quad *quad = d_face_quad_table[order];
        assert(face < NDIM+1 && face >= 0);
        //double vol = (*d_cell_volume)(0, cell);
        bzero(&(bb[0]), 1* sizeof(*value));
        double outer_normal[NDIM], outer_normal_raw[NDIM];
        // get the id of face in this patch.
        int patch_face = caf_idx[caf_ext[cell] + face];
        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[patch_face+1]-fan_ext[patch_face],
                fan_idx.getPointer() + fan_ext[patch_face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);
        double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw))/ 2.0;
        for(int n = 0; n < quad->npoints; n++) {
            // transform this 2d barycentric quadrature point to 3d barycentic on face 'face'
            double bary3d[NDIM+1];
            bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,face,bary3d);
            // get the space coordinate of this quadrature point
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               bary3d,
                               space_coords);
            // get the 'nbas' shape functions' value at this quadrature point
            double phi[nbas][dim];
            double triphi[3][2];
            shapefunc->basis(cell,
                             bary3d,
                             &(phi[0][0]));
            trishapefunc->basis(patch_face,
                             quad->points + (n*(NDIM+1-1)),
                             &(triphi[0][0]),nodeGradBas,faceorder);
            // compute the load of this element
            TYPE2 Ec[3]={0,0,0};
            for(int i = 0; i < nbas; i++) {
                Ec[0]+=E_edge_tol[i]*phi[i][0];
                Ec[1]+=E_edge_tol[i]*phi[i][1];
                Ec[2]+=E_edge_tol[i]*phi[i][2];
            }
            TYPE2 E2_2d[2]={0,0};
            for(int i=0; i<3; i++){
                E2_2d[0]+=E_edge_modal[i]*triphi[i][0];
                E2_2d[1]+=E_edge_modal[i]*triphi[i][1];
            }

            double *nodephi=quad->points + (n*(NDIM+1-1));
            TYPE2 E2z;
            int i_order[3]={0,1,2};
            if(faceorder==0){
                i_order[1]=2;
                i_order[2]=1;
            }
            for(int i=0;i<3;i++){
                E2z+=E_node_modal[i]*nodephi[i_order[i]];
            }

            TYPE2 E2_3d[3]={0,0,0};
            if(direction[0]==1)
            {
                E2_3d[1]=E2_2d[0];//x
                E2_3d[2]=E2_2d[1];//y
                E2_3d[0]=E2z;//z
            }else if(direction[1]==1){
                E2_3d[0]=E2_2d[0];
                E2_3d[2]=E2_2d[1];
                E2_3d[1]=E2z;
            }else{
                E2_3d[0]=E2_2d[0];
                E2_3d[1]=E2_2d[1];
                E2_3d[2]=E2z;
            }


            TYPE2 E2_conj[3]={0,0,0};
            for(int i = 0; i<3; i++){
                E2_conj[i]=dcomplex(E2_3d[i].real(),-E2_3d[i].imag());
            }
//            cout<<"E2 "<<E2_3d[0]<<"||"<<E2_3d[1]<<"||"<<E2_3d[2]<<endl;
//            cout<<"E2_conj "<<E2_conj[0]<<"||"<<E2_conj[1]<<"||"<<E2_conj[2]<<endl;
//            cout<<"Ec "<<Ec[0]<<"||"<<Ec[1]<<"||"<<Ec[2]<<endl;
            TYPE v = 0;
            TYPE2 E1minusE2[3]={0,0,0};
            for(int i=0; i<3; i++){
                E1minusE2[i]=Ec[i]+E2_conj[i];//方向包含在里面
            }
            v=dotProduct(3,E1minusE2,E2_conj);
            v *= quad->weights[n];
           // cout<<quad->weights[n]<<endl;
            (*bb) =(*bb)+ v;
        }
       // cout<<(*bb)<<endl;
                (*bb) =(*bb)* area;
    }
    template<class ShapeFunc,class triShapeFunc, class TYPE, class TYPE2>
    void faceQuadE2DotE2_modal(const int cell,
                             const int face,
                             const TYPE2 *E_edge_tol,//总场边上插值系数
                             const ShapeFunc* shapefunc,//体形函数
                             const TYPE2 *E_edge_modal,//模式场边上插值系数
                             const TYPE2 *E_node_modal,//模式场节点上插值系数
                             const triShapeFunc* trishapefunc,//面形函数
                             tbox::Array<tbox::Array<double> > nodeGradBas,
                             tbox::Array<double> direction,
                                   int faceorder,
                             const int order,
                             TYPE *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        TYPE *bb = value;
        const Quad *quad = d_face_quad_table[order];
        assert(face < NDIM+1 && face >= 0);
        //double vol = (*d_cell_volume)(0, cell);
        bzero(&(bb[0]), 1* sizeof(*value));
        double outer_normal[NDIM], outer_normal_raw[NDIM];
        // get the id of face in this patch.
        int patch_face = caf_idx[caf_ext[cell] + face];
        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[patch_face+1]-fan_ext[patch_face],
                fan_idx.getPointer() + fan_ext[patch_face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);
        double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw))/ 2.0;
        for(int n = 0; n < quad->npoints; n++) {
            // transform this 2d barycentric quadrature point to 3d barycentic on face 'face'
            double bary3d[NDIM+1];
            bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,face,bary3d);
            // get the space coordinate of this quadrature point
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               bary3d,
                               space_coords);
            // get the 'nbas' shape functions' value at this quadrature point
            double phi[nbas][dim];
            double triphi[3][2];
            shapefunc->basis(cell,
                             bary3d,
                             &(phi[0][0]));
            trishapefunc->basis(patch_face,
                             quad->points + (n*(NDIM+1-1)),
                             &(triphi[0][0]),nodeGradBas,faceorder);
            // compute the load of this element
            TYPE2 Ec[3]={0,0,0};
            for(int i = 0; i < nbas; i++) {
                Ec[0]+=E_edge_tol[i]*phi[i][0];
                Ec[1]+=E_edge_tol[i]*phi[i][1];
                Ec[2]+=E_edge_tol[i]*phi[i][2];
            }
            TYPE2 E2_2d[2]={0,0};
            for(int i=0; i<3; i++){
                E2_2d[0]+=E_edge_modal[i]*triphi[i][0];
                E2_2d[1]+=E_edge_modal[i]*triphi[i][1];
            }
            double *nodephi=quad->points + (n*(NDIM+1-1));
            TYPE2 E2z;
            int i_order[3]={0,1,2};
            if(faceorder==0){
                i_order[1]=2;
                i_order[2]=1;
            }
            for(int i=0;i<3;i++){
                E2z+=E_node_modal[i]*nodephi[i_order[i]];
            }

            TYPE2 E2_3d[3]={0,0,0};
            if(direction[0]==1)
            {
                E2_3d[1]=E2_2d[0];//x
                E2_3d[2]=E2_2d[1];//y
                E2_3d[0]=E2z;//z
            }else if(direction[1]==1){
                E2_3d[0]=E2_2d[0];
                E2_3d[2]=E2_2d[1];
                E2_3d[1]=E2z;
            }else{
                E2_3d[0]=E2_2d[0];
                E2_3d[1]=E2_2d[1];
                E2_3d[2]=E2z;
            }
            TYPE2 E2_conj[3]={0,0,0};
            for(int i = 0; i<3; i++){
                E2_conj[i]=dcomplex(E2_3d[i].real(),-E2_3d[i].imag());
            }
            TYPE v = 0;
            v=dotProduct(3,E2_3d,E2_conj);
            v *= quad->weights[n];
           // cout<<quad->weights[n]<<endl;
            (*bb) =(*bb)+ v;
        }
       // cout<<(*bb)<<endl;
                (*bb) =(*bb)* area;
    }
    template<class ShapeFunc,class triShapeFunc, class TYPE, class TYPE2>
    void faceQuadE1DotE2_modal(const int cell,
                             const int face,
                             const TYPE2 *E_edge_tol,//总场边上插值系数
                             const ShapeFunc* shapefunc,//体形函数
                             const TYPE2 *E_edge_modal,//模式场边上插值系数
                             const TYPE2 *E_node_modal,//模式场节点上插值系数
                             const triShapeFunc* trishapefunc,//面形函数
                             tbox::Array<tbox::Array<double> > nodeGradBas,
                             tbox::Array<double> direction,
                                   int faceorder,
                             const int order,
                             TYPE *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        TYPE *bb = value;
        const Quad *quad = d_face_quad_table[order];
        assert(face < NDIM+1 && face >= 0);
        //double vol = (*d_cell_volume)(0, cell);
        bzero(&(bb[0]), 1* sizeof(*value));
        double outer_normal[NDIM], outer_normal_raw[NDIM];
        // get the id of face in this patch.
        int patch_face = caf_idx[caf_ext[cell] + face];
        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[patch_face+1]-fan_ext[patch_face],
                fan_idx.getPointer() + fan_ext[patch_face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);
        double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw))/ 2.0;
        for(int n = 0; n < quad->npoints; n++) {
            // transform this 2d barycentric quadrature point to 3d barycentic on face 'face'
            double bary3d[NDIM+1];
            bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,face,bary3d);
            // get the space coordinate of this quadrature point
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               bary3d,
                               space_coords);
            // get the 'nbas' shape functions' value at this quadrature point
            double phi[nbas][dim];
            double triphi[3][2];
            shapefunc->basis(cell,
                             bary3d,
                             &(phi[0][0]));
            trishapefunc->basis(patch_face,
                             quad->points + (n*(NDIM+1-1)),
                             &(triphi[0][0]),nodeGradBas,faceorder);
            // compute the load of this element
            TYPE2 Ec[3]={0,0,0};
            for(int i = 0; i < nbas; i++) {
                Ec[0]+=E_edge_tol[i]*phi[i][0];
                Ec[1]+=E_edge_tol[i]*phi[i][1];
                Ec[2]+=E_edge_tol[i]*phi[i][2];
            }
            TYPE2 E2_2d[2]={0,0};
            for(int i=0; i<3; i++){
                E2_2d[0]+=E_edge_modal[i]*triphi[i][0];
                E2_2d[1]+=E_edge_modal[i]*triphi[i][1];
            }

            double *nodephi=quad->points + (n*(NDIM+1-1));
            TYPE2 E2z;
            int i_order[3]={0,1,2};
            if(faceorder==0){
                i_order[1]=2;
                i_order[2]=1;
            }
            for(int i=0;i<3;i++){
                E2z+=E_node_modal[i]*nodephi[i_order[i]];
            }

            TYPE2 E2_3d[3]={0,0,0};
            if(direction[0]==1)
            {
                E2_3d[1]=E2_2d[0];//x
                E2_3d[2]=E2_2d[1];//y
                E2_3d[0]=E2z;//z
            }else if(direction[1]==1){
                E2_3d[0]=E2_2d[0];
                E2_3d[2]=E2_2d[1];
                E2_3d[1]=E2z;
            }else{
                E2_3d[0]=E2_2d[0];
                E2_3d[1]=E2_2d[1];
                E2_3d[2]=E2z;
            }


            TYPE2 E2_conj[3]={0,0,0};
            for(int i = 0; i<3; i++){
                E2_conj[i]=dcomplex(E2_3d[i].real(),-E2_3d[i].imag());
            }
//            cout<<"E2 "<<E2_3d[0]<<"||"<<E2_3d[1]<<"||"<<E2_3d[2]<<endl;
//            cout<<"E2_conj "<<E2_conj[0]<<"||"<<E2_conj[1]<<"||"<<E2_conj[2]<<endl;
//            cout<<"Ec "<<Ec[0]<<"||"<<Ec[1]<<"||"<<Ec[2]<<endl;
            TYPE v = 0;
            v=dotProduct(3,Ec,E2_conj);
            v *= quad->weights[n];
           // cout<<quad->weights[n]<<endl;
            (*bb) =(*bb)+ v;
        }
       // cout<<(*bb)<<endl;
                (*bb) =(*bb)* area;
    }
    template<class ShapeFunc, class TYPE, class TYPE2>
    void faceQuadEnorm(const int cell,
                                const int face,
                                const TYPE2 *E_edge,
                                const ShapeFunc* shapefunc,
                                const int order,
                                TYPE *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        TYPE *bb = value;
        const Quad *quad = d_face_quad_table[order];
        assert(face < NDIM+1 && face >= 0);
        //double vol = (*d_cell_volume)(0, cell);
        bzero(&(bb[0]), 1* sizeof(*value));
        double outer_normal[NDIM], outer_normal_raw[NDIM];
        // get the id of face in this patch.
        int patch_face = caf_idx[caf_ext[cell] + face];
        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[patch_face+1]-fan_ext[patch_face],
                fan_idx.getPointer() + fan_ext[patch_face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);
        double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw))/ 2.0;
        for(int n = 0; n < quad->npoints; n++) {
            // transform this 2d barycentric quadrature point to 3d barycentic on face 'face'
            double bary3d[NDIM+1];
            bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,face,bary3d);
            // get the space coordinate of this quadrature point
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               bary3d,
                               space_coords);
            // get the 'nbas' shape functions' value at this quadrature point
            double phi[nbas][dim];
            shapefunc->basis(cell,
                             bary3d,
                             &(phi[0][0]));
            // compute the load of this element
            TYPE2 temp[3]={0,0,0};
            for(int i = 0; i < nbas; i++) {
                temp[0]+=E_edge[i]*phi[i][0];
                temp[1]+=E_edge[i]*phi[i][1];
                temp[2]+=E_edge[i]*phi[i][2];
            }
            TYPE v = 0;
            v=sqrt(norm(temp[0])+norm(temp[1])+norm(temp[2]));
            v *= quad->weights[n];
           // cout<<quad->weights[n]<<endl;
            (*bb) =(*bb)+ v;
        }
       // cout<<(*bb)<<endl;
                (*bb) =(*bb)* area;
    }

    void bary2dTo3d(const double *lambda_2d,int cell, int face,double* lambda_3d)const{
        int patch_face = caf_idx[caf_ext[cell]+face];
        int cnt = 0;
        for (int i=0; i<4; i++){
            lambda_3d[i]=0;
        }
        for (int j = fan_ext[patch_face];j<fan_ext[patch_face+1];j++){
            int face_node = fan_idx[j];
            for (int i=can_ext[cell];i<can_ext[cell+1];i++){
                int cell_node = can_idx[i];
                int local_node = i-can_ext[cell];
                if (cell_node == face_node){
                    lambda_3d[local_node]=lambda_2d[cnt++];
                }
            }
        }
    }
    void bary2dTo2dMap(double *lambda_2d1,int* Mapping1_2,double* lambda_2d2)const{
        for(int cell1_local = 0; cell1_local < NDIM+1; cell1_local ++){
            if(Mapping1_2[cell1_local] == -1) continue;
            int cell2_local = Mapping1_2[cell1_local];
            lambda_2d2[cell2_local] = lambda_2d1[cell1_local];
        }
    }
    template<class ShapeFunc, class TYPE, class TYPE2>
    void lineQuadbas(const int cell,
                                const int edge,
                                const TYPE2 *E_edge,
                                const ShapeFunc* shapefunc,
                                const int order,
                                TYPE *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        TYPE *bb = value;
        const Quad *quad = d_edge_quad_table[order];
        assert(edge < NDIM+3 && edge >= 0);
        //double vol = (*d_cell_volume)(0, cell);
        bzero(&(bb[0]), 1* sizeof(*value));
        double vecedge[NDIM],unit_vecedge[NDIM];
        //int *nedgenode;
        int *edgenode;
        double *nodecoords;
        // get the id of face in this patch.
        int patch_edge = cae_idx[cae_ext[cell] + edge];
        nodecoords=d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer();
        for(int i=0;i<dim;i++)
        {
            edgenode=ean_idx.getPointer() + ean_ext[patch_edge];
            vecedge[i]=nodecoords[edgenode[1]*dim+i]-
            nodecoords[edgenode[0]*dim+i];
        }
        double area = sqrt(dotProduct(NDIM, vecedge, vecedge));//边的长度
        for(int i=0;i<dim;i++){
            unit_vecedge[i]=abs(vecedge[i]/area);
        }
        for(int n = 0; n < quad->npoints; n++) {
            // transform this 1d barycentric quadrature point to 3d barycentic on edge 'edge'
            double bary3d[NDIM+1];
            bary1dTo3d(quad->points + (n*(NDIM-1)),cell,edge,bary3d);

            // get the space coordinate of this quadrature point
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               bary3d,
                               space_coords);
            // get the 'nbas' shape functions' value at this quadrature point
            double phi[nbas][dim];
            shapefunc->basis(cell,
                             bary3d,
                             &(phi[0][0]));
            // compute the load of this element
            TYPE2 temp[3]={0,0,0};
            for(int i = 0; i < nbas; i++) {
                temp[0]=phi[i][0];
                temp[1]=phi[i][1];
                temp[2]=phi[i][2];
                TYPE v = 0;
                v=dotProduct(NDIM, temp, unit_vecedge);
               // v=1;
                v *= quad->weights[n];
               //cout<<quad->weights[n]<<endl;
                bb[i] =bb[i]+ v;
            }
        }
        //cout<<(*bb)<<endl;
        for(int i = 0; i < nbas; i++)
                bb[i] *= area;
    }
    template<class ShapeFunc, class TYPE, class TYPE2>
    void lineQuadEnorm(const int cell,
                                const int edge,
                                const TYPE2 *E_edge,
                                const ShapeFunc* shapefunc,
                                const int order,
                                TYPE *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        TYPE *bb = value;
        const Quad *quad = d_edge_quad_table[order];
        assert(edge < NDIM+3 && edge >= 0);
        //double vol = (*d_cell_volume)(0, cell);
        bzero(&(bb[0]), 1* sizeof(*value));
        double vecedge[NDIM],unit_vecedge[NDIM];
        //int *nedgenode;
        int *edgenode;
        double *nodecoords;
        // get the id of face in this patch.
        int patch_edge = cae_idx[cae_ext[cell] + edge];
        nodecoords=d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer();
        for(int i=0;i<dim;i++)
        {
            edgenode=ean_idx.getPointer() + ean_ext[patch_edge];
            vecedge[i]=nodecoords[edgenode[1]*dim+i]-
            nodecoords[edgenode[0]*dim+i];
        }
        double area = sqrt(dotProduct(NDIM, vecedge, vecedge));
        //cout<<"area  "<<area<<endl;
        for(int i=0;i<dim;i++){
            unit_vecedge[i]=abs(vecedge[i]/area);
        }
        for(int n = 0; n < quad->npoints; n++) {
            // transform this 1d barycentric quadrature point to 3d barycentic on edge 'edge'
            double bary3d[NDIM+1];
            bary1dTo3d(quad->points + (n*(NDIM-1)),cell,edge,bary3d);

            // get the space coordinate of this quadrature point
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               bary3d,
                               space_coords);
            // get the 'nbas' shape functions' value at this quadrature point
            double phi[nbas][dim];
            shapefunc->basis(cell,
                             bary3d,
                             &(phi[0][0]));
            // compute the load of this element
            TYPE2 temp[3]={0,0,0};
            for(int i = 0; i < nbas; i++) {
                temp[0]+=E_edge[i]*phi[i][0];
                temp[1]+=E_edge[i]*phi[i][1];
                temp[2]+=E_edge[i]*phi[i][2];
            }
            TYPE v = 0;

            v=dotProduct(NDIM, temp, unit_vecedge);
           // v=1;
            v *= quad->weights[n];
           //cout<<quad->weights[n]<<endl;
            (*bb) =(*bb)+ v;
        }
        //cout<<(*bb)<<endl;
                (*bb) =(*bb)* area;
    }

void bary1dTo3d(const double *lambda_1d,int cell, int edge,double* lambda_3d)const{
    int patch_edge = cae_idx[cae_ext[cell]+edge];
    int cnt = 0;
    for (int i=0; i<4; i++){
        lambda_3d[i]=0;
    }
    for (int j = ean_ext[patch_edge];j<ean_ext[patch_edge+1];j++){
        int edge_node = ean_idx[j];
        for (int i=can_ext[cell];i<can_ext[cell+1];i++){
            int cell_node = can_idx[i];
            int local_node = i-can_ext[cell];
            if (cell_node == edge_node){
                lambda_3d[local_node]=lambda_1d[cnt++];
            }
        }
    }
}
    template<class ShapeFunc>
    void quadGradBasDotGradBas(const int cell,
                           const ShapeFunc* shapefunc,
                           const int order,
                       double *value) const
    {
    check(order);
    TBOX_ASSERT(shapefunc->dim() == 1);

    int nbas = shapefunc->nbas();
        double (*gg)[nbas] = (double (*)[nbas])value;
    const Quad *quad = d_quad_table[order];

    double vol = (*d_cell_volume)(0, cell);
    bzero(&(gg[0][0]), nbas * nbas * sizeof(*value));

    for(int n = 0; n < quad->npoints; n++) {
        double gradphi[nbas][NDIM];
        shapefunc->gradient(cell,
                        quad->points + n * (NDIM + 1),
                    &(gradphi[0][0]));
        for(int i = 0; i < nbas; i++) {
        for(int j = 0; j < nbas; j++) {
            double v = 0.;
            for(int d = 0; d < NDIM; d++)
            v += gradphi[i][d] * gradphi[j][d];
            v *= quad->weights[n];
            gg[i][j] += v;
        }
        }
    }
    for(int i = 0; i < nbas; i++)
        for(int j = 0; j < nbas; j++)
        gg[i][j] *= vol;
    }

    template <class ShapeFunc>
    void quadCurlBasDotCurlBas(const int cell, const ShapeFunc *shapefunc,
                               const int order, dcomplex *value) const {
      check(order);
      TBOX_ASSERT(shapefunc->dim() == NDIM);

      int nbas = shapefunc->nbas();
      dcomplex(*cc) [nbas] = (dcomplex(*)[nbas])value;
      const Quad *quad = d_quad_table[order];

      double vol = (*d_cell_volume)(0, cell);
      bzero(&(cc[0][0]), nbas * nbas * sizeof(*value));

      for (int n = 0; n < quad->npoints; n++) {
        double curlphi[nbas][NDIM];
        shapefunc->curl(cell, quad->points + n * (NDIM + 1), &(curlphi[0][0]));

        for (int i = 0; i < nbas; i++) {
          for (int j = 0; j < nbas; j++) {
            dcomplex v = 0.;
            for (int d = 0; d < NDIM; d++) v += curlphi[i][d]*curlphi[j][d];
            v *= quad->weights[n];
            cc[i][j] += v;
          }
        }
      }

      for (int i = 0; i < nbas; i++)
        for (int j = 0; j < nbas; j++) cc[i][j] *= vol;
    }

    template<class ShapeFunc, class TYPE>
    void quadDofDotBas(const int cell,
                   const ShapeFunc* shapefunc,
                   const int order,
               const TYPE *dof_value,
               TYPE *value) const
    {
    check(order);
    int nbas = shapefunc->nbas();
    int dim = shapefunc->dim();
    const Quad *quad = d_quad_table[order];

    double vol = (*d_cell_volume)(0, cell);
    bzero(value, nbas * sizeof(*value));

    for(int n = 0; n < quad->npoints; n++) {
        double phi[nbas][dim];
        shapefunc->basis(cell,
                     quad->points + n * (NDIM + 1),
                 &(phi[0][0]));
        for(int i = 0; i < nbas; i++) {
            TYPE v = 0.;
            for(int j = 0; j < nbas; j++) {
                TYPE u = 0.;
                for(int d = 0; d < dim; d++)
                    u += phi[i][d] * phi[j][d];
                v += u * dof_value[j];
            }
            value[i] += v * quad->weights[n];
        }
    }
    for(int i = 0; i < nbas; i++)
      value[i] *= vol;
    }


    template<class ShapeFunc, class TYPE, class Function>
    void quadFunctionDotBas(
            const int cell,
            const Function& func,
            const ShapeFunc* shapefunc,
            const int order,
            TYPE *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        const Quad *quad = d_quad_table[order];

        double vol = (*d_cell_volume)(0, cell);
        bzero(value, nbas * sizeof(*value));

        for(int n = 0; n < quad->npoints; n++) {
            double phi[nbas][dim];
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               quad->points + n *(NDIM + 1),
                               space_coords);
            TYPE funcvalues[NDIM];
            func(space_coords[0],space_coords[1],space_coords[2], funcvalues);
            shapefunc->basis(cell,
                             quad->points + n * (NDIM + 1),
                             &(phi[0][0]));
            for(int i = 0; i < nbas; i++) {
                TYPE prod = dotProduct(NDIM, funcvalues, phi[i]);
                value[i] += prod * quad->weights[n];
            }
        }
        for(int i = 0; i < nbas; i++)
            value[i] *= vol;
    }

    template<class ShapeFunc, class TYPE>
    void quadCurlDofDotCurlDof(const int cell,
                           const ShapeFunc* shapefunc,
                           const int order,
                       const TYPE *dof_u,
                       const TYPE *dof_v,
                       TYPE *value) const
    {
    check(order);
    int nbas = shapefunc->nbas();
    int dim = shapefunc->dim();
    const Quad *quad = d_quad_table[order];

    double vol = (*d_cell_volume)(0, cell);
        *value = 0.;

    for(int n = 0; n < quad->npoints; n++) {
        double curlphi[nbas][NDIM];
        shapefunc->curl(cell,
                    quad->points + n * (NDIM + 1),
                &(curlphi[0][0]));

            TYPE uv[NDIM], vv[NDIM];
            for(int d = 0; d < NDIM; d++)
                uv[d] = vv[d] = 0.;

        for(int i = 0; i < nbas; i++) {
                for(int d = 0; d < dim; d++) {
                    uv[d] += dof_u[i] * curlphi[i][d];
                    vv[d] += dof_v[i] * curlphi[i][d];
                }
            }

            TYPE q = dotProduct(dim, uv, vv);

            //TYPE q = 0.;
            //for(int d = 0; d < dim; d++)
            //    q += uv[d] * vv[d];

            *value += q * quad->weights[n];
    }

        *value *= vol;
    }

    template <class ShapeFunc1, class ShapeFunc2>
    void computeL2norm(const int cell, const ShapeFunc1 *shapefunc,
                       const ShapeFunc2 *shapefunc2, const int order,
                       dcomplex *e_value, dcomplex *r_value,double *L2_errors) const {
        // e_value is the DOF in Nedelec, one edge each
        // r_value is the DOF in Brezzi-Douglas-Thomas one face three
        check(order);
        (*L2_errors) = 0.;
        const Quad *quad = d_quad_table[order];
        for(int n = 0; n < quad->npoints; n++){
            double nedelec_phi[6][NDIM];
            double bdm_phi[12][NDIM];
            shapefunc->basis(cell, quad->points + n * (NDIM + 1), &(nedelec_phi[0][0]));
            shapefunc2->basis(cell, quad->points + n * (NDIM + 1), &(bdm_phi[0][0]));
            dcomplex v_0_gauss[3] = {dcomplex(0.,0.), dcomplex(0.,0.), dcomplex(0.,0.)};
            dcomplex v_s_gauss[3] = {dcomplex(0.,0.), dcomplex(0.,0.), dcomplex(0.,0.)};
            dcomplex error_gauss[3] = {dcomplex(0.,0.), dcomplex(0.,0.), dcomplex(0.,0.)};
            for(int ee = 0; ee < 6; ee++){
                for(int dim =0; dim < NDIM; dim ++) v_0_gauss[dim] += 0.;//e_value[ee]*nedelec_phi[ee][dim];
            }
            for(int ff =0; ff < 4; ff++){
                for(int nn =0; nn < 3; nn++){
                    for(int dim =0; dim < NDIM; dim ++) v_s_gauss[dim] += r_value[ff*3+nn]*bdm_phi[ff*3+nn][dim];
                }
            }
            for(int dim =0; dim < NDIM; dim ++) error_gauss[dim] = v_s_gauss[dim]-v_0_gauss[dim];
            double L2_norm_gauss =
                    error_gauss[0].real()*error_gauss[0].real()+error_gauss[1].real()*error_gauss[1].real()+error_gauss[2].real()*error_gauss[2].real()
                    +error_gauss[0].imag()*error_gauss[0].imag()+error_gauss[1].imag()*error_gauss[1].imag()+error_gauss[2].imag()*error_gauss[2].imag();
            L2_norm_gauss *= quad->weights[n];
            (*L2_errors) += L2_norm_gauss;
        }



    }
    template<class ShapeFunc>
    void testBDMfunction(const int cell,
                            const ShapeFunc* shapefunc,
                            const int order
                         ) const{
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        const Quad *quad = d_face_quad_table[order];
        double bdm_phi[12][NDIM];
        // This function aims at testing every single face to check
        // If the BDM on them with this quad rule
        double outer_normal[NDIM], outer_normal_raw[NDIM];
        //The outer normal should all be in the same direction of one basis function
        //
        double integration_test[4][3][12];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 4; j ++)
                for (int k =0 ;k < 12; k++) integration_test[i][j][k]=0.;
        for(int iface = 0; iface < 4; iface ++){
            int face = caf_idx[caf_ext[cell]+iface];
            int fc_node[3] = {0,0,0};
            for(int inode =0; inode < 4; inode ++){
                if(can_idx[can_ext[cell]+inode] == fan_idx[fan_ext[face]+0])
                    fc_node[0] = inode;
                if(can_idx[can_ext[cell]+inode] == fan_idx[fan_ext[face]+1])
                    fc_node[1] = inode;
                if(can_idx[can_ext[cell]+inode] == fan_idx[fan_ext[face]+2])
                    fc_node[2] = inode;
                // The corresponding element of fc_nodes is the local node of the cell
                // This element can be a index of bary3d

            }
            outerNormal(
                    can_ext[cell+1]-can_ext[cell],
                    can_idx.getPointer() + can_ext[cell],
                    fan_ext[face+1]-fan_ext[face],
                    fan_idx.getPointer() + fan_ext[face],
                    d_patch.getNumberOfNodes(1),
                    d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                    NDIM,
                    outer_normal,
                    outer_normal_raw);
            double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw)) / 2.0;
            if(caf_idx[caf_ext[cell]] != face){
                outer_normal[0] = -outer_normal[0];
                outer_normal[1] = -outer_normal[1];
                outer_normal[2] = -outer_normal[2];
            }
            for(int n = 0; n < quad->npoints; n++){
                // The basis function on face iface and quad point n
                double bary3d[NDIM+1];
                bary2dTo3d(quad->points + (n*(NDIM+1-1)),cell,iface,bary3d);
                shapefunc->basis(cell,bary3d,&(bdm_phi[0][0]));
                for(int bas = 0; bas < nbas; bas ++){
                    for (int nodebas = 0; nodebas < 3; nodebas ++) {
                        double nodal_basdot[3];
                        double basdot = DOT_PRODUCT(bdm_phi[bas*3+nodebas],outer_normal);
                        //basdot is the /psi_i,j where i=bas,j=nodebas
                        for(int lambda_node =0; lambda_node < 3; lambda_node ++){
                            double lambda_node_dot;
                            lambda_node_dot = basdot*bary3d[fc_node[lambda_node]];

                            if(bas == iface){
                                // When the face is corresponding to the basis function
                                cout<<lambda_node_dot<<" : This is the face component computed by BDM"<<endl;
                                // cout<<"face number: "<<iface<<", number node: "<<lambda_node<<endl;
                                double expected_psi = -3*bary3d[fc_node[0]]-3*bary3d[fc_node[1]]-3*bary3d[fc_node[2]];
                                expected_psi += 12*bary3d[fc_node[nodebas]];
                                expected_psi *= bary3d[fc_node[lambda_node]];
                                cout<<expected_psi<<" : This is the standard face component computed by BDM"<<endl;


                            }
                            lambda_node_dot *= quad->weights[n];
                            integration_test[iface][lambda_node][bas*3+nodebas] += lambda_node_dot;
//                            cout<<"This element is added: "<<iface<<","<<lambda_node<<","<<(bas*3+nodebas)<<" by "<<lambda_node_dot<<endl;
//                            cout<<"elemental number"<<"["<<iface<<"]"<<"["<<lambda_node<<"]"<<"["<<bas*3+nodebas<<"]"<<endl;
                            // integration_test[i][j][k]
                            // The k-th element integration on i-th face and j-th node
                            // If success, only k=i*3+j element of integration_test
                        }
                    }

                }
            }

        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j ++) {
                for (int k =0 ;k < 12; k++) {
                    if(k == i*3+j){
                        cout<<"It should be 1: "<<integration_test[i][j][k]<<" "
                           <<i<<","<<j<<","<<k<<endl;
                    }
                    else {
                        cout<<"It should be 0: "<<integration_test[i][j][k]<<" "
                           <<i<<","<<j<<","<<k<<endl;
                    }
                }

            }

        }




    }

public:
    struct Quad {
    const char *name;
    int dim;
    int order;
    int npoints;
    double *points;
    double *weights;
    int id;
    };


private:
    void check(int order) const {
    const int max_order = 3;
    if(order > max_order)
        TBOX_ERROR("Quadrature rules of order "
               << order << " for tetrahedra not implememted.\n");
    }

    const hier::Patch<NDIM>& d_patch;
    tbox::Pointer<pdat::CellData<NDIM, double> > d_cell_volume;
    tbox::Pointer<pdat::CellData<NDIM, double> > d_cell_jacobian;

    tbox::Array<int> can_ext, can_idx;
    tbox::Array<int> fan_ext, fan_idx;
    tbox::Array<int> caf_ext, caf_idx;
    tbox::Array<int> ean_ext, ean_idx;
    tbox::Array<int> cae_ext, cae_idx;
    tbox::Array<int> fac_ext, fac_idx;

    Quad** d_quad_table;
    Quad** d_face_quad_table;
    Quad** d_edge_quad_table;

public:
    Quad** GetQuadTable() {return d_quad_table;}
    Quad** GetFaceQuadTable() {return d_face_quad_table;}
};
} // namespace appu
} // namespace JAUMIN

#endif // included_appu_TetQuad
