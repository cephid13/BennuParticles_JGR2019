///////////////////////////////////////////////////////////////////////////
// GetShapeModelGravity_mex
// Owner: Jay McMahon (The University of Colorado Boulder)
///////////////////////////////////////////////////////////////////////////
//
// Description:
//
//  This function computes the potential, acceleration, and dynamics matrix
//  given the homogeneous density value, the spacecraft position, and shape model
//  for Eros (regular shape model).
//
// Inputs:
//
//     rho                            : [kg/km^3] Asteroid density
//
//     Num_Ptr                        : [n.d.]    Number of vertices, facets, and edges
//
//     SatPos_Ptr                     : [km]      Satellite position
//
//     Vertex_Ptr                     : [n.d.]    Indices and coordinates for each vertex
//
//     FacetVertex_1_Ptr              : [n.d.]    x, y, z coordinate for vertex 1 of each facet
//
//     FacetVertex_2_Ptr              : [n.d.]    x, y, z coordinate for vertex 2 of each facet
//
//     FacetVertex_3_Ptr              : [n.d.]    x, y, z coordinate for vertex 3 of each facet
//
//     EdgesIndex_Ptr                 : [n.d.]    Indices for each edge
//
//     EdgeFacetNumberDirection_1_Ptr : [n.d.]    Facet number and the edge direction for edge 1
//
//     EdgeFacetNumberDirection_2_Ptr : [n.d.]    Facet number and the edge direction for edge 2
//
//     Edges_1_Ptr                    : [n.d.]    x, y, z coordinate for edge 1
//
//     Edges_2_Ptr                    : [n.d.]    x, y, z coordinate for edge 2
//
// Outputs:
//
//     U_Ptr        : [km^2/s^2]  Spacecraft potential
//
//     del_U_Ptr    : [km/s^2]    Spacecraft acceleration
// 
//     deldel_U_Ptr : [1/s^2]     A-matrix for the STM
//
//     del2_U_Ptr   : [str]       signed solid angle
//
// Assumptions/References:
//
//  - Refer to Exterior gravitation of a polyhedron derived and compared
//   with harmonic and mascon gravitation representations of asteroid 4769
//   Castalia
//
// Note:
//
//  - Each facet is triangular
//
// Dependencies:
//
//  - None
//
// Call
//
//  - None
//
// Called by
//
//  - ComputeShapeModelGravity.m
//
// Modification History:
//
//  07May12   Yu Takahashi   original version
//  04Jun12   Yu Takahashi   1st revision
//
///////////////////////////////////////////////////////////////////////////

#include    "mex.h"
#include	<math.h>
#include	<stdarg.h>
#include	<stdio.h>
#include	<stdlib.h>
#include    <string.h>

#define ABS(x) ((x) < 0) ? -(x) : (x)

#define G 6.67384E-20

#define num_vertex_max      30000
#define num_facet_max       50000
#define num_edges_total_max 80000
#define num_edges_per_surface 3

///////////////////
// -- Outputs -- //
///////////////////

double U, del_U [3], deldel_U [3][3], del2_U;
double *U_Ptr, *del_U_Ptr, *deldel_U_Ptr, *del2_U_Ptr, *Omega_Ptr;

/////////////////
// -- Input -- //
/////////////////

double rho;
double *Num_Ptr;
int    num_vertex, num_facet, num_edges_total; 
double *SatPos_Ptr, SatPos [3];
double *Vertex_Ptr, Vertex [num_vertex_max][4];
double *FacetVertex_1_Ptr, FacetVertex_1 [num_facet_max][3];
double *FacetVertex_2_Ptr, FacetVertex_2 [num_facet_max][3];
double *FacetVertex_3_Ptr, FacetVertex_3 [num_facet_max][3];
double *EdgesIndex_Ptr, EdgesIndex [num_edges_total_max][2];
double *EdgeFacetNumberDirection_1_Ptr, EdgeFacetNumberDirection_1[num_edges_total_max][2];
double *EdgeFacetNumberDirection_2_Ptr, EdgeFacetNumberDirection_2[num_edges_total_max][2];
double *Edges_1_Ptr, Edges_1 [num_edges_total_max][3];
double *Edges_2_Ptr, Edges_2 [num_edges_total_max][3];

//////////////////////////////////////////////////////////////////////////////
// -- Number of Elements in the Tetrahedral model & Range of x, y, and z -- //
//////////////////////////////////////////////////////////////////////////////

int xx, yy, zz;

/////////////////
// -- Index -- //
/////////////////

double n, m;

int ee, ff, gg, ii, nn, mm; 
int FacetNumber_1, FacetNumber_2;
int VertexDirection_1, VertexDirection_2;
int IndexVertex_1, IndexVertex_2;

//////////////////////////////
// -- Computation Arrays -- //
//////////////////////////////

double n_hat_f [3][num_facet_max];
double Pi [3], Pj [3], Pk [3];
double rErL, ErL [3],  EL [3][3];
double rFrOmega, FrOmega [3], FOmega [3][3];
double Omega, Omega_temp, L, L_temp;
double E [3][3], F [3][3], F_temp [3][3], E_temp [3][3];
double F_temp_rf [3], E_temp_re [3];
double cos_sin_wf [2], cos_sin_wf_temp [2];
double rj [3], rj_hat [3], ri [3], rk [3], rj_norm;
double rj_ri_cross [3], rj_rk_cross [3];
double sji [3], sji_hat [3], sjk_hat [3], sjk [3], sji_norm, sjk_norm;
double cos_S, sin_S, n_hat_sign[3];
double sign_Omega, Edge_1 [3], Edge_2 [3];
double cross_Edge_12 [3], norm_cross_Edge_12;
double rf [3];
int    SurfaceNumber [2], VertexDirection [2];
double Vertex_1 [3], Vertex_2 [3];
double a_vec [3], b_vec [3], e_vec [3];
double a, b, e;
double n_hat_A [3], EdgeLine_A [3], n_hat_A_12 [3], cross_hat_Edge_A [3], norm_cross_hat_Edge_A;
double n_hat_B [3], EdgeLine_B [3], n_hat_B_21 [3], cross_hat_Edge_B [3], norm_cross_hat_Edge_B;
double re [3];

////////////////////
// -- Function -- //
////////////////////

void ComputePotentialPolygon(double SatPos[3]);

/***************************************************************
 * main program
 ***************************************************************/
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    // Assign Inputs
    
    rho                            = mxGetScalar(prhs[0]); // [kg/km^3] Asteroid density
    Num_Ptr                        = mxGetPr(prhs[1]);     // [n.d.]    Number of vertices, facets, and edges.
    SatPos_Ptr                     = mxGetPr(prhs[2]);     // [km]      Satellite position
    Vertex_Ptr                     = mxGetPr(prhs[3]);     // [n.d.]    Indices for each vertex
    FacetVertex_1_Ptr              = mxGetPr(prhs[4]);     // [n.d.]    x, y, z vertices for surface 1
    FacetVertex_2_Ptr              = mxGetPr(prhs[5]);     // [n.d.]    x, y, z vertices for surface 2
    FacetVertex_3_Ptr              = mxGetPr(prhs[6]);     // [n.d.]    x, y, z vertices for surface 3
    EdgesIndex_Ptr                 = mxGetPr(prhs[7]);     // [n.d.]    Indices for each edge
    EdgeFacetNumberDirection_1_Ptr = mxGetPr(prhs[8]);     // [n.d.]    Facet index and direction of the edges 1
    EdgeFacetNumberDirection_2_Ptr = mxGetPr(prhs[9]);     // [n.d.]    Facet index and direction of the edges 2
    Edges_1_Ptr                    = mxGetPr(prhs[10]);    // [n.d.]    x, y, z vertices for edge 1
    Edges_2_Ptr                    = mxGetPr(prhs[11]);    // [n.d.]    x, y, z vertices for edge 1

    //////////////////////////////
    //// -- Assign Outputs -- ////
    //////////////////////////////
    
    plhs[0]      = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1]      = mxCreateDoubleMatrix(3, 1, mxREAL);
    plhs[2]      = mxCreateDoubleMatrix(3, 3, mxREAL);
    plhs[3]      = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[4]      = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    U_Ptr        = mxGetPr(plhs[0]); // [km^2/s^2]  Spacecraft potential
    del_U_Ptr    = mxGetPr(plhs[1]); // [km/s^2]    Spacecraft acceleration
    deldel_U_Ptr = mxGetPr(plhs[2]); // [1/s^2]     A-matrix for the STM
    del2_U_Ptr   = mxGetPr(plhs[3]); // [1/s^2]     delta^2 U
    Omega_Ptr    = mxGetPr(plhs[4]); // [str]       signed solid angle
    
    ////////////////////////
    // Number of Elements //
    ////////////////////////
    
    num_vertex      = (int) Num_Ptr[0]; // [n.d.] Number of vertices
    num_facet       = (int) Num_Ptr[1]; // [n.d.] Number of facets
    num_edges_total = (int) Num_Ptr[2]; // [n.d.] Number of edges
    
    //////////////////////////
    // Field Point Position //
    //////////////////////////
    
    for (nn = 0; nn < 3; nn++) {
        
        SatPos[nn] = SatPos_Ptr[nn];
        
    } // For nn
    
    ///////////////////////////
    // Allocate the Surfaces //
    ///////////////////////////
    
    for (mm = 0; mm < 3; mm++) {
        
        for (nn = 0; nn < num_facet; nn++) {
            
            // 1st Vertex
            
            FacetVertex_1[nn][mm] = FacetVertex_1_Ptr[num_facet*mm + nn];
            
            // 2nd Vertex
            
            FacetVertex_2[nn][mm] = FacetVertex_2_Ptr[num_facet*mm + nn];
            
            // 3rd Vertex
            
            FacetVertex_3[nn][mm] = FacetVertex_3_Ptr[num_facet*mm + nn];
            
        } // For nn
        
    } // For mm
    
    ////////////////////////
    // Allocate the Edges //
    ////////////////////////
    
    for (mm = 0; mm < 3; mm++) {
        
        for (nn = 0; nn < num_edges_total; nn++) {
            
            // 1st Edge
            
            Edges_1[nn][mm] = Edges_1_Ptr[num_edges_total*mm + nn];
            
            // 2nd Edge
            
            Edges_2[nn][mm] = Edges_2_Ptr[num_edges_total*mm + nn];
            
        } // For nn
        
    } // For mm    
    
    for (mm = 0; mm < 2; mm++) {
        
        for (nn = 0; nn < num_edges_total; nn++) {
            
            // Edge Index
            
            EdgesIndex[nn][mm] = (EdgesIndex_Ptr[num_edges_total*mm + nn] - 1.0);
            
            if (mm == 0) {
            
                EdgeFacetNumberDirection_1[nn][mm] = (EdgeFacetNumberDirection_1_Ptr[num_edges_total*mm + nn] - 1.0);
                EdgeFacetNumberDirection_2[nn][mm] = (EdgeFacetNumberDirection_2_Ptr[num_edges_total*mm + nn] - 1.0);
                    
            } else if (mm == 1) {
                
                EdgeFacetNumberDirection_1[nn][mm] = EdgeFacetNumberDirection_1_Ptr[num_edges_total*mm + nn];
                EdgeFacetNumberDirection_2[nn][mm] = EdgeFacetNumberDirection_2_Ptr[num_edges_total*mm + nn];
                
            } // For if
            
        } // For nn
        
    } // For mm
    
    ///////////////////////////
    // Allocate the Vertices //
    ///////////////////////////
    
    for (mm = 0; mm < 4; mm++) {
        
        for (nn = 0; nn < num_vertex; nn++) {
            
            // Vertex Index
            
            if (mm == 0) {
                
                Vertex[nn][mm] = (Vertex_Ptr[num_vertex*mm + nn] - 1.0);
                
            } else {
                
                Vertex[nn][mm] = Vertex_Ptr[num_vertex*mm + nn];
                
            } // For if
            
        } // For nn
        
    } // For mm
    
    ///////////////////////////////////////////////////////
    // Compute the Potential, Acceleration, and A-matrix //
    ///////////////////////////////////////////////////////
    
    Omega = 0.0; Omega_temp = 0.0; L = 0.0; L_temp = 0.0;
    rFrOmega = 0.0; rErL = 0.0; U = 0.0, del2_U = 0.0;
    
    ComputePotentialPolygon(SatPos);
    
    ////////////////////////
    // Assign the Outputs //
    ////////////////////////
    
    U_Ptr[0]      = U;
    del2_U_Ptr[0] = del2_U;
    Omega_Ptr[0]  = Omega;
    
    for (mm = 0; mm < 3; mm++) {
        
        del_U_Ptr[mm] = del_U[mm];
        
        for (nn = 0; nn < 3; nn++) {
            
            deldel_U_Ptr[3*mm+nn] = deldel_U[nn][mm];
            
        } // For mm
        
        
    } // For nn
    
    return;
    
} // For Main

void ComputePotentialPolygon(double SatPos[3])

{
    
    // Initialization
    
    for (nn = 0; nn < 3; nn++) {
        
        ErL [nn]     = 0.0;
        FrOmega [nn] = 0.0;
        del_U[nn]    = 0.0;
        
        for (mm = 0; mm < 3; mm++) {
            
            EL [nn][mm]      = 0.0;
            FOmega [nn][mm]  = 0.0;
            E [nn][mm]       = 0.0;
            F [nn][mm]       = 0.0;
            
            deldel_U[nn][mm] = 0.0;
            
        } // for mm
        
    } // for nn
    
    for (nn = 0; nn < 3; nn++) {
        
        for (mm = 0; mm < num_facet; mm++) {
            
            n_hat_f[nn][mm] = 0.0;
            
        } // % For mm
        
    } // For nn
    
    /////////////////////////////////////////////
    //// -- Compute the surface integrals -- ////
    /////////////////////////////////////////////
    
    for (ff = 0; ff < num_facet; ++ff) {
        
        cos_sin_wf [0] = 1.0; cos_sin_wf [1] = 0.0;
        
        for (gg = 0; gg < num_edges_per_surface; ++gg) {
                        
            if (gg == 0) { // For the first vertex, take the last and the next vertex
                
                Pj [0] = FacetVertex_1[ff][0]; Pj [1] = FacetVertex_1[ff][1]; Pj [2] = FacetVertex_1[ff][2]; // Middle Vertex
                Pi [0] = FacetVertex_3[ff][0]; Pi [1] = FacetVertex_3[ff][1]; Pi [2] = FacetVertex_3[ff][2]; // Middle Vertex
                Pk [0] = FacetVertex_2[ff][0]; Pk [1] = FacetVertex_2[ff][1]; Pk [2] = FacetVertex_2[ff][2]; // Middle Vertex
                
            } else if (gg == 1) { // For others, take the previous and the next vertex

                Pj [0] = FacetVertex_2[ff][0]; Pj [1] = FacetVertex_2[ff][1]; Pj [2] = FacetVertex_2[ff][2]; // Middle Vertex
                Pi [0] = FacetVertex_1[ff][0]; Pi [1] = FacetVertex_1[ff][1]; Pi [2] = FacetVertex_1[ff][2]; // Middle Vertex
                Pk [0] = FacetVertex_3[ff][0]; Pk [1] = FacetVertex_3[ff][1]; Pk [2] = FacetVertex_3[ff][2]; // Middle Vertex
                
            } else if (gg == 2) { // For the last vertex, take the previous and the last vertex

                Pj [0] = FacetVertex_3[ff][0]; Pj [1] = FacetVertex_3[ff][1]; Pj [2] = FacetVertex_3[ff][2]; // Middle Vertex
                Pi [0] = FacetVertex_2[ff][0]; Pi [1] = FacetVertex_2[ff][1]; Pi [2] = FacetVertex_2[ff][2]; // Middle Vertex
                Pk [0] = FacetVertex_1[ff][0]; Pk [1] = FacetVertex_1[ff][1]; Pk [2] = FacetVertex_1[ff][2]; // Middle Vertex
                
            } // For if
            
            // rj, ri, and rk //
            
            rj[0] = Pj[0] - SatPos[0]; rj[1] = Pj[1] - SatPos[1]; rj[2] = Pj[2] - SatPos[2];
            ri[0] = Pi[0] - SatPos[0]; ri[1] = Pi[1] - SatPos[1]; ri[2] = Pi[2] - SatPos[2];
            rk[0] = Pk[0] - SatPos[0]; rk[1] = Pk[1] - SatPos[1]; rk[2] = Pk[2] - SatPos[2];
            
            // |rj|, rj_hat //
            
            rj_norm   = sqrt(rj[0]*rj[0] + rj[1]*rj[1] + rj[2]*rj[2]);
            rj_hat[0] = rj[0]/rj_norm; rj_hat[1] = rj[1]/rj_norm; rj_hat[2] = rj[2]/rj_norm;
            
            // sji_hat and sjk_hat, and their cross product
            
            rj_ri_cross[0] = rj[1]*ri[2] - rj[2]*ri[1]; rj_ri_cross[1] = rj[2]*ri[0] - rj[0]*ri[2]; rj_ri_cross[2] = rj[0]*ri[1] - rj[1]*ri[0];
            rj_rk_cross[0] = rj[1]*rk[2] - rj[2]*rk[1]; rj_rk_cross[1] = rj[2]*rk[0] - rj[0]*rk[2]; rj_rk_cross[2] = rj[0]*rk[1] - rj[1]*rk[0];
            
            sji[0] = rj_ri_cross[1]*rj[2] - rj_ri_cross[2]*rj[1]; sji[1] = rj_ri_cross[2]*rj[0] - rj_ri_cross[0]*rj[2]; sji[2] = rj_ri_cross[0]*rj[1] - rj_ri_cross[1]*rj[0];
            sjk[0] = rj_rk_cross[1]*rj[2] - rj_rk_cross[2]*rj[1]; sjk[1] = rj_rk_cross[2]*rj[0] - rj_rk_cross[0]*rj[2]; sjk[2] = rj_rk_cross[0]*rj[1] - rj_rk_cross[1]*rj[0];
            
            sji_norm = sqrt(sji[0]*sji[0] + sji[1]*sji[1] + sji[2]*sji[2]);
            sjk_norm = sqrt(sjk[0]*sjk[0] + sjk[1]*sjk[1] + sjk[2]*sjk[2]);
            
            sji_hat[0] = sji[0]/sji_norm; sji_hat[1] = sji[1]/sji_norm; sji_hat[2] = sji[2]/sji_norm;
            sjk_hat[0] = sjk[0]/sjk_norm; sjk_hat[1] = sjk[1]/sjk_norm; sjk_hat[2] = sjk[2]/sjk_norm;
            
            n_hat_sign[0] = sjk_hat[1]*sji_hat[2] - sjk_hat[2]*sji_hat[1];
            n_hat_sign[1] = sjk_hat[2]*sji_hat[0] - sjk_hat[0]*sji_hat[2];
            n_hat_sign[2] = sjk_hat[0]*sji_hat[1] - sjk_hat[1]*sji_hat[0];
            
            // cos and sin of S
            
            cos_S = sji_hat[0]*sjk_hat[0] + sji_hat[1]*sjk_hat[1] + sji_hat[2]*sjk_hat[2];
            sin_S = ABS(rj_hat[0]*n_hat_sign[0] + rj_hat[1]*n_hat_sign[1] + rj_hat[2]*n_hat_sign[2]);
            
            // Matrix multiplication : cos_sin_wf = [-cos_S, sin_S; -sin_S, -cos_S]*cos_sin_wf;
            
            cos_sin_wf_temp[0] = - cos_S*cos_sin_wf[0] + sin_S*cos_sin_wf[1];
            cos_sin_wf_temp[1] = - sin_S*cos_sin_wf[0] - cos_S*cos_sin_wf[1];
            
            cos_sin_wf[0] = cos_sin_wf_temp[0];
            cos_sin_wf[1] = cos_sin_wf_temp[1];
            
        } // for gg
        
        sign_Omega = rj[0]*n_hat_sign[0] + rj[1]*n_hat_sign[1] + rj[2]*n_hat_sign[2];
        
        if (sign_Omega > 0.0) {
            
            sign_Omega = 1.0;
            
        } else if (sign_Omega < 0.0) {
            
            sign_Omega = -1.0;
            
        } else if (sign_Omega == 0.0) {
            
            sign_Omega = 0.0;
            
        } // for if
        
        Omega_temp = sign_Omega*2.0*atan2(1.0 - cos_sin_wf[0], cos_sin_wf[1]);
        
        Omega     += Omega_temp;
        
        // Edge 1 and 2 of the Surface
        
        Edge_1[0] = FacetVertex_2[ff][0] - FacetVertex_1[ff][0];  // [km]   Edge AB
        Edge_1[1] = FacetVertex_2[ff][1] - FacetVertex_1[ff][1];  // [km]   Edge AB
        Edge_1[2] = FacetVertex_2[ff][2] - FacetVertex_1[ff][2];  // [km]   Edge AB
        
        Edge_2[0] = FacetVertex_3[ff][0] - FacetVertex_1[ff][0];  // [km]   Edge AC
        Edge_2[1] = FacetVertex_3[ff][1] - FacetVertex_1[ff][1];  // [km]   Edge AC
        Edge_2[2] = FacetVertex_3[ff][2] - FacetVertex_1[ff][2];  // [km]   Edge AC
        
        // [n.d.] Surface normal
        
        cross_Edge_12[0] = Edge_1[1]*Edge_2[2] - Edge_1[2]*Edge_2[1];
        cross_Edge_12[1] = Edge_1[2]*Edge_2[0] - Edge_1[0]*Edge_2[2];
        cross_Edge_12[2] = Edge_1[0]*Edge_2[1] - Edge_1[1]*Edge_2[0];
        
        norm_cross_Edge_12 = sqrt(cross_Edge_12[0]*cross_Edge_12[0] + cross_Edge_12[1]*cross_Edge_12[1] + cross_Edge_12[2]*cross_Edge_12[2]);
        
        n_hat_f[0][ff] = cross_Edge_12[0]/norm_cross_Edge_12;
        n_hat_f[1][ff] = cross_Edge_12[1]/norm_cross_Edge_12;
        n_hat_f[2][ff] = cross_Edge_12[2]/norm_cross_Edge_12;
        
        // rf
        
        rf[0]  =  FacetVertex_1[ff][0] - SatPos[0];
        rf[1]  =  FacetVertex_1[ff][1] - SatPos[1];
        rf[2]  =  FacetVertex_1[ff][2] - SatPos[2];
        
        for (nn = 0; nn < 3; nn++) {
            
            F_temp_rf[nn] = 0.0;
            
            for (mm = 0; mm < 3; mm++) {
                
                F_temp[nn][mm]  = n_hat_f[nn][ff]*n_hat_f[mm][ff];
                F[nn][mm]      += F_temp[nn][mm];
                F_temp_rf[nn]  += F_temp[nn][mm]*rf[mm];
                
            } // for mm
            
        } // for nn
        
        rFrOmega += ( rf[0]*F_temp_rf[0] + rf[1]*F_temp_rf[1] + rf[2]*F_temp_rf[2] ) * Omega_temp;
        
        for (nn = 0; nn < 3; nn++) {
            
            FrOmega[nn] += F_temp_rf[nn]*Omega_temp;
            
            for (mm = 0; mm < 3; mm++) {
                
                FOmega[nn][mm] += F_temp[nn][mm]*Omega_temp;
                
            } // for mm
            
        } // for nn
        
    } // for ff
    
    //////////////////////////////////////////
    //// -- Compute the line integrals -- ////
    //////////////////////////////////////////
    
    for (ee = 0; ee < num_edges_total; ee++) {
        
        ///////////////////////////////////////////////////////////////
        //// -- Define the Edge/Vertices/Index for the Vertices -- ////
        ///////////////////////////////////////////////////////////////
        
        Vertex_1[0]   = Edges_1[ee][0]; Vertex_1[1] = Edges_1[ee][1]; Vertex_1[2] = Edges_1[ee][2]; // [km]   Position of the first vertex
        Vertex_2[0]   = Edges_2[ee][0]; Vertex_2[1] = Edges_2[ee][1]; Vertex_2[2] = Edges_2[ee][2]; // [km]   Position fo the second vertex
        IndexVertex_1 = (int) EdgesIndex[ee][0];     // [n.d.] Position of the first vertex
        IndexVertex_2 = (int) EdgesIndex[ee][1];     // [n.d.] Position of the second vertex
        
        //////////////////////////////
        //// -- Line Potential -- ////
        //////////////////////////////
        
        a_vec[0] = SatPos[0] - Vertex_1[0];   a_vec[1] = SatPos[1] - Vertex_1[1];   a_vec[2] = SatPos[2] - Vertex_1[2];         // [km]   Length between the vertex 1 and the satellite position
        b_vec[0] = SatPos[0] - Vertex_2[0];   b_vec[1] = SatPos[1] - Vertex_2[1];   b_vec[2] = SatPos[2] - Vertex_2[2];         // [km]   Length between the vertex 2 and the satellite position
        e_vec[0] = Vertex_2[0] - Vertex_1[0]; e_vec[1] = Vertex_2[1] - Vertex_1[1]; e_vec[2] = Vertex_2[2] - Vertex_1[2];   // [km]   Length of the edge
        
        a = sqrt(a_vec[0]*a_vec[0] + a_vec[1]*a_vec[1] + a_vec[2]*a_vec[2]);
        b = sqrt(b_vec[0]*b_vec[0] + b_vec[1]*b_vec[1] + b_vec[2]*b_vec[2]);
        e = sqrt(e_vec[0]*e_vec[0] + e_vec[1]*e_vec[1] + e_vec[2]*e_vec[2]);
        
        L_temp  = log((a + b + e)/(a + b - e));  // [n.d.] Line Potential
        L      += L_temp;
        
        /////////////////////////////////////////////////////////////////////////////////
        //// -- Define the face normal and the edge normal for the adjacent faces -- ////
        /////////////////////////////////////////////////////////////////////////////////
        
        // Facet and Edge direction
        
        FacetNumber_1     = (int) EdgeFacetNumberDirection_1[ee][0];
        VertexDirection_1 = (int) EdgeFacetNumberDirection_1[ee][1];
        FacetNumber_2     = (int) EdgeFacetNumberDirection_2[ee][0];
        VertexDirection_2 = (int) EdgeFacetNumberDirection_2[ee][1];
        
        // For the first Face
        
        n_hat_A[0] = n_hat_f[0][FacetNumber_1]; n_hat_A[1] = n_hat_f[1][FacetNumber_1]; n_hat_A[2] = n_hat_f[2][FacetNumber_1]; // [n.d.] Face normal for the first surface
                
        if (VertexDirection_1 == 1) {
            
            EdgeLine_A[0] = Vertex[IndexVertex_2][1] - Vertex[IndexVertex_1][1];
            EdgeLine_A[1] = Vertex[IndexVertex_2][2] - Vertex[IndexVertex_1][2];
            EdgeLine_A[2] = Vertex[IndexVertex_2][3] - Vertex[IndexVertex_1][3];
            
        } else if (VertexDirection_1 == -1 ) {
            
            EdgeLine_A[0] = Vertex[IndexVertex_1][1] - Vertex[IndexVertex_2][1];
            EdgeLine_A[1] = Vertex[IndexVertex_1][2] - Vertex[IndexVertex_2][2];
            EdgeLine_A[2] = Vertex[IndexVertex_1][3] - Vertex[IndexVertex_2][3];
            
        } // For if
        
        cross_hat_Edge_A[0] = n_hat_A[1]*EdgeLine_A[2] - n_hat_A[2]*EdgeLine_A[1];
        cross_hat_Edge_A[1] = n_hat_A[2]*EdgeLine_A[0] - n_hat_A[0]*EdgeLine_A[2];
        cross_hat_Edge_A[2] = n_hat_A[0]*EdgeLine_A[1] - n_hat_A[1]*EdgeLine_A[0];
        
        norm_cross_hat_Edge_A = sqrt(cross_hat_Edge_A[0]*cross_hat_Edge_A[0] + cross_hat_Edge_A[1]*cross_hat_Edge_A[1] + cross_hat_Edge_A[2]*cross_hat_Edge_A[2]);
        
        n_hat_A_12[0] = - cross_hat_Edge_A[0]/norm_cross_hat_Edge_A;
        n_hat_A_12[1] = - cross_hat_Edge_A[1]/norm_cross_hat_Edge_A;
        n_hat_A_12[2] = - cross_hat_Edge_A[2]/norm_cross_hat_Edge_A;
        
        // For the second Face
        
        n_hat_B[0] = n_hat_f[0][FacetNumber_2]; n_hat_B[1] = n_hat_f[1][FacetNumber_2]; n_hat_B[2] = n_hat_f[2][FacetNumber_2];
        
        
        
        if ( VertexDirection_2 == 1 ) {
            
            EdgeLine_B[0] = Vertex[IndexVertex_2][1] - Vertex[IndexVertex_1][1];
            EdgeLine_B[1] = Vertex[IndexVertex_2][2] - Vertex[IndexVertex_1][2];
            EdgeLine_B[2] = Vertex[IndexVertex_2][3] - Vertex[IndexVertex_1][3];
            
        } else if ( VertexDirection_2 == -1 ) {
            
            EdgeLine_B[0] = Vertex[IndexVertex_1][1] - Vertex[IndexVertex_2][1];
            EdgeLine_B[1] = Vertex[IndexVertex_1][2] - Vertex[IndexVertex_2][2];
            EdgeLine_B[2] = Vertex[IndexVertex_1][3] - Vertex[IndexVertex_2][3];
            
        } // For if
       
        cross_hat_Edge_B[0] = n_hat_B[1]*EdgeLine_B[2] - n_hat_B[2]*EdgeLine_B[1];
        cross_hat_Edge_B[1] = n_hat_B[2]*EdgeLine_B[0] - n_hat_B[0]*EdgeLine_B[2];
        cross_hat_Edge_B[2] = n_hat_B[0]*EdgeLine_B[1] - n_hat_B[1]*EdgeLine_B[0];
        
        norm_cross_hat_Edge_B = sqrt(cross_hat_Edge_B[0]*cross_hat_Edge_B[0] + cross_hat_Edge_B[1]*cross_hat_Edge_B[1] + cross_hat_Edge_B[2]*cross_hat_Edge_B[2]);
        
        n_hat_B_21[0] = - cross_hat_Edge_B[0]/norm_cross_hat_Edge_B;
        n_hat_B_21[1] = - cross_hat_Edge_B[1]/norm_cross_hat_Edge_B;
        n_hat_B_21[2] = - cross_hat_Edge_B[2]/norm_cross_hat_Edge_B;
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        //// -- Compute E and Other Terms to Compile -- ////
        ////////////////////////////////////////////////////////////////////////////////////////////////
        
        re[0] = Vertex_1[0] - SatPos[0]; re[1] = Vertex_1[1] - SatPos[1]; re[2] = Vertex_1[2] - SatPos[2]; // [n.d.] Vector from the field point to any point on edge e or its infinite extension
        
        for (nn = 0; nn < 3; nn++) {
            
            E_temp_re[nn] = 0.0;
            
            for (mm = 0; mm < 3; mm++) {
                
                E_temp[nn][mm]  = n_hat_A[nn]*n_hat_A_12[mm] + n_hat_B[nn]*n_hat_B_21[mm];
                E[nn][mm]      += E_temp[nn][mm];
                
                E_temp_re[nn]  += E_temp[nn][mm]*re[mm];
                
            } // for mm
            
        } // for nn
        
        rErL +=  ( re[0]*E_temp_re[0] + re[1]*E_temp_re[1] + re[2]*E_temp_re[2] ) * L_temp;
        
        for (nn = 0; nn < 3; nn++) {
            
            ErL[nn] += E_temp_re[nn]*L_temp;
            
            for (mm = 0; mm < 3; mm++) {
                
                EL[nn][mm] += E_temp[nn][mm]*L_temp;
                
            } // for mm
            
        } // for nn
        
    } // For each Line
    
    /////////////////////////////////////////////
    //// -- Compute the surface integrals -- ////
    /////////////////////////////////////////////
    
    U = 0.5*G*rho*rErL - 0.5*G*rho*rFrOmega; // [km^2/s^2]
    
    for (nn = 0; nn < 3; nn++) {
        
        del_U[nn]  += - G*rho*ErL[nn] + G*rho*FrOmega[nn]; // [km/s^2]
                
        for (mm = 0; mm < 3; mm++) {
            
            deldel_U[nn][mm] += G*rho*EL[nn][mm] - G*rho*FOmega[nn][mm]; // [1/s^2]
                        
        } // for mm
        
    } // for nn
    
    del2_U  = - G*rho*Omega; // [1/s^2]
        
} // For ComputePotentialPolygon

///////////////////////////////////////////////////////////////////////////