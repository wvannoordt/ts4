#ifndef RHS_HPP
#define RHS_HPP
#include "cmf.h"
#include "temp/MdArray.h"
#include "time_control.h"
#include "fluid_state.h"
using cmf::face_t;
using cmf::cell_t;

#define stencilIdx(v,j) ((v)+(5+3)*(j))
#define f_DivSplit(q,j,l,v1)         (0.500*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))]))
#define fg_QuadSplit(q,j,l,v1,v2)    (0.250*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))])*(q[stencilIdx((v2),(j))] + q[stencilIdx((v2),(j)+(l))]))
#define fg_CubeSplit(q,j,l,v1,v2,v3) (0.125*(q[stencilIdx((v1),(j))] + q[stencilIdx((v1),(j)+(l))])*(q[stencilIdx((v2),(j))] + q[stencilIdx((v2),(j)+(l))])*(q[stencilIdx((v3),(j))] + q[stencilIdx((v3),(j)+(l))]))
#define fg_DivSplit(q,j,l,v1,v2)     (0.500*((q[stencilIdx((v1),(j)+(l))]*q[stencilIdx((v2),(j))]) + (q[stencilIdx((v1),(j))]*q[stencilIdx((v2),(j)+(l))])))
void ComputeConv(cmf::CartesianMeshArray& prims, cmf::CartesianMeshArray& rhs, int centOrder, gas::perfect_gas_t<double>& gas, double alpha)
{
    double centerCoef[4] = {0.0};
    switch (centOrder)
    {
        case 2: {centerCoef[0] = 1.0/2.0; break;}
        case 4: {centerCoef[0] = 2.0/3.0; centerCoef[1] = -1.0/12.0; break;}
        case 6: {centerCoef[0] = 3.0/4.0; centerCoef[1] = -3.0/20.0; centerCoef[2] = 1.0/60.0; break;}
        case 8: {centerCoef[0] = 4.0/5.0; centerCoef[1] = -1.0/5.0 ; centerCoef[2] = 4.0/105 ; centerCoef[3] = -1.0/280.0; break;}
        default: {std::cout << "Bad central scheme order." << std::endl; abort();}
    }
    double Rgas = gas.R;
    for (auto lb: prims)
    {
        cmf::BlockArray<double, 1> primsLb = prims[lb];
        cmf::BlockArray<double, 1> rhsLb  = rhs[lb];
        cmf::BlockInfo info = rhs.Mesh()->GetBlockInfo(lb);
        int stencilWid = centOrder/2;
        for (int idir = 0; idir < 3; idir++)
        {
            int dijk[3] = {0};
            dijk[idir] = 1;
            for (cmf::cell_t k = primsLb.kmin; k < primsLb.kmax; k++)
            {
                for (cmf::cell_t j = primsLb.jmin; j < primsLb.jmax; j++)
                {
                    for (cmf::cell_t i = primsLb.imin; i < primsLb.imax; i++)
                    {
                        
                        double stencilData[9*(8)]; //ie,ke,T,P,rho,u,v,w
                        for (int n = 0; n < centOrder + 1; n++)
                        {
                            for (int v = 3; v < (5+3); v++)
                            {
                                int ii = i+dijk[0]*(n-stencilWid);
                                int jj = j+dijk[1]*(n-stencilWid);
                                int kk = k+dijk[2]*(n-stencilWid);
                                stencilData[stencilIdx(v,n)] = primsLb(v-3, ii, jj, kk);
                            }
                            // stencilData = ? ? ? P T u v w
                            // T
                            stencilData[stencilIdx(2,n)] = stencilData[stencilIdx(4,n)];
                            
                            // stencilData = ? ? T P T u v w
                            //rho
                            stencilData[stencilIdx(4,n)] = stencilData[stencilIdx(3,n)]/(Rgas*stencilData[stencilIdx(2,n)]);//p = rho r t -> rho = p/RT
                            
                            // IE = P/(rho*(gamma - 1))
                            stencilData[stencilIdx(0,n)] = stencilData[stencilIdx(3,n)]/(stencilData[stencilIdx(4,n)]*(gas.gamma - 1.0));
                            
                            // ke (don't care)
                            stencilData[stencilIdx(1,n)] = 0.0;

                            // Not needed per se starts
                            for (int vel_comp = 0; vel_comp < 3; vel_comp ++)
                            {
                                stencilData[stencilIdx(1,n)] += 0.5*stencilData[stencilIdx(5+vel_comp,n)]*stencilData[stencilIdx(5+vel_comp,n)];
                            }
                            // Not needed per se ends
                        }
                        
                        double C[2]     = {0.0};
                        double M[6]     = {0.0};
                        double PGRAD[2] = {0.0};
                        double KE[2]    = {0.0};
                        double IE[2]    = {0.0};
                        double PDIFF[2] = {0.0};
                        
                        for (int l = 1; l <= stencilWid; l++)
                        {
                            double al = centerCoef[l-1];
                            int jf = stencilWid;
                            for (int m = 0; m <= (l-1); m++)
                            {
                                C[1] += 2.0*centerCoef[l-1]*fg_QuadSplit(stencilData,jf-m, l,4,5+idir);
                                C[0] += 2.0*centerCoef[l-1]*fg_QuadSplit(stencilData,jf+m,-l,4,5+idir);
                                for (int idir_mom = 0; idir_mom < 3; idir_mom++)
                                {
                                    M[idir_mom    ] += 2.0*centerCoef[l-1]*fg_CubeSplit(stencilData,jf-m, l,4,5+idir,5+idir_mom);
                                    M[idir_mom + 3] += 2.0*centerCoef[l-1]*fg_CubeSplit(stencilData,jf+m,-l,4,5+idir,5+idir_mom);
                                }

                                PGRAD[1] += 2.0*centerCoef[l-1]*f_DivSplit(stencilData,jf-m, l,3);
                                PGRAD[0] += 2.0*centerCoef[l-1]*f_DivSplit(stencilData,jf+m,-l,3);

                                for (int vel_comp = 0;  vel_comp < 3; vel_comp ++)
                                {
                                    KE[1] += 2.0*centerCoef[l-1]*fg_QuadSplit(stencilData,jf-m, l,4,5+idir)*0.5*(stencilData[stencilIdx(5+vel_comp,jf-m)]*stencilData[stencilIdx(5+vel_comp,jf-m+l)]);
                                    KE[0] += 2.0*centerCoef[l-1]*fg_QuadSplit(stencilData,jf+m,-l,4,5+idir)*0.5*(stencilData[stencilIdx(5+vel_comp,jf+m)]*stencilData[stencilIdx(5+vel_comp,jf+m-l)]);
                                }

                                IE[1] += 2.0*centerCoef[l-1]*fg_CubeSplit(stencilData,jf-m, l,4,0,5+idir);
                                IE[0] += 2.0*centerCoef[l-1]*fg_CubeSplit(stencilData,jf+m,-l,4,0,5+idir);

                                PDIFF[1] += 2.0*centerCoef[l-1]*fg_DivSplit(stencilData,jf-m, l,5+idir,3);
                                PDIFF[0] += 2.0*centerCoef[l-1]*fg_DivSplit(stencilData,jf+m,-l,5+idir,3);
                            }
                        }
                        rhsLb(0, i, j, k)      -= (1.0-alpha)*info.dxInv[idir]*(C[1] - C[0]);
                        rhsLb(1, i, j, k)      -= (1.0-alpha)*info.dxInv[idir]*(IE[1] + KE[1] + PDIFF[1] - IE[0] - KE[0] - PDIFF[0]);
                        rhsLb(2, i, j, k)      -= (1.0-alpha)*info.dxInv[idir]*(M[0] - M[3]);
                        rhsLb(3, i, j, k)      -= (1.0-alpha)*info.dxInv[idir]*(M[1] - M[4]);
                        rhsLb(4, i, j, k)      -= (1.0-alpha)*info.dxInv[idir]*(M[2] - M[5]);
                        rhsLb(2+idir, i, j, k) -= (1.0-alpha)*info.dxInv[idir]*(PGRAD[1] - PGRAD[0]);
                    }
                }
            }
            dijk[idir] = 0;
        }
    }
}

void ComputeConvDiss(cmf::CartesianMeshArray& prims, cmf::CartesianMeshArray& rhs, gas::perfect_gas_t<double>& gas, double alpha)
{
    double Rgas = gas.R;
    for (auto lb: prims)
    {
        cmf::BlockArray<double, 1> primsLb = prims[lb];
        cmf::BlockArray<double, 1> rhsLb  = rhs[lb];
        cmf::BlockInfo info = rhs.Mesh()->GetBlockInfo(lb);
        DataView<double, Dims<5,5,2>> fluxStencil;
        double fluxR[5];
        double fluxL[5];
        for (int idir = 0; idir < 3; idir++)
        {
            int dijk[3] = {0};
            dijk[idir] = 1;
            for (cmf::cell_t k = primsLb.kmin; k < primsLb.kmax; k++)
            {
                for (cmf::cell_t j = primsLb.jmin; j < primsLb.jmax; j++)
                {
                    for (cmf::cell_t i = primsLb.imin; i < primsLb.imax; i++)
                    {
                        for (int n1 = 0; n1 < 5; n1++)
                        {
                            //P T U V W
                            fluid_state::prim_t<double> prim;
                            for (uint pp = 0; pp < 5; pp++) prim[pp] = primsLb(pp, i+(n1-2)*dijk[0], j+(n1-2)*dijk[1], k+(n1-2)*dijk[2]);
                            fluid_state::cons_t<double> cons;
                            
                            cons[0] = prim.p() / (gas.R*prim.T());
                            double rhoU2 = cons[0]*(prim.u()*prim.u()+prim.v()*prim.v()+prim.w()*prim.w());
                            cons[1] = 0.5*rhoU2 + (prim.p()/((gas.gamma - 1.0)));
                            cons[2] = cons[0]*prim.u();
                            cons[3] = cons[0]*prim.v();
                            cons[4] = cons[0]*prim.w();
                            
                            double rho   = cons.rho();
                            double rhoE  = cons.rho_H();
                            double rhoU  = cons.rho_u();
                            double rhoW  = cons.rho_v();
                            double rhoV  = cons.rho_w();
                            double Un    = cons[2+idir]/rho;
                            for (int pp = 0; pp < 2; pp++)
                            {
                                fluxStencil(n1, 0, pp) =  rho*Un*0.5;
                                fluxStencil(n1, 1, pp) = rhoE*Un*0.5;
                                fluxStencil(n1, 2, pp) = rhoU*Un*0.5;
                                fluxStencil(n1, 3, pp) = rhoV*Un*0.5;
                                fluxStencil(n1, 4, pp) = rhoW*Un*0.5;
                            }
                        }
                        
                        auto weno3 = [](double al, double bl, double cl, double dl, double ar, double br, double cr, double dr) -> double
                        {
                            //  a  b     c  d (left and right)
                            //  +  +  |  +  +
                            
                            //right
                            double rr1 = -0.5*ar+1.5*br;
                            double rr2 =  0.5*br+0.5*cr;
                            
                            double br1 = (ar-br)*(ar-br);
                            double br2 = (cr-br)*(cr-br);
                            double sumr = br1 + br2;
                            
                            br1/=sumr;
                            br2/=sumr;
                            
                            double wr1 = 0.3333333333/(1e-16+br1);
                            double wr2 = 0.6666666667/(1e-16+br2);
                            
                            //left
                            double rl1 = -0.5*dl+1.5*cl;
                            double rl2 =  0.5*bl+0.5*cl;
                            
                            double bl1 = (al-bl)*(al-bl);
                            double bl2 = (cl-bl)*(cl-bl);
                            double suml = bl1 + bl2;
                            
                            bl1/=suml;
                            bl2/=suml;
                            
                            double wl1 = 0.3333333333/(1e-16+bl1);
                            double wl2 = 0.6666666667/(1e-16+bl2);
                            
                            return rr1*wr1+rr2*wr2+rl1*wl1+rl2*wl2;
                        };
                        
                        for (int n1 = 0; n1 < 5; n1++)
                        {
                            // left
                            fluxR[n1] = weno3(fluxStencil(0, n1, 0), fluxStencil(1, n1, 0), fluxStencil(2, n1, 0), fluxStencil(3, n1, 0), fluxStencil(0, n1, 1), fluxStencil(1, n1, 1), fluxStencil(2, n1, 1), fluxStencil(3, n1, 1));
                            fluxL[n1] = weno3(fluxStencil(1, n1, 0), fluxStencil(2, n1, 0), fluxStencil(3, n1, 0), fluxStencil(4, n1, 0), fluxStencil(1, n1, 1), fluxStencil(2, n1, 1), fluxStencil(3, n1, 1), fluxStencil(4, n1, 1));
                            rhsLb(n1, i, j, k) -= alpha*(fluxR[n1]-fluxL[n1])*info.dxInv[idir];
                        }
                    }
                }
            }
            dijk[idir] = 0;
        }
    }
}

void ComputeVisc(cmf::CartesianMeshArray& prims, cmf::CartesianMeshArray& rhs, viscous_laws::constant_viscosity_t<double>& vlaw)
{
    double beta = 0.0;
    for (auto lb: prims)
    {
        cmf::BlockArray<double, 1> primsLb = prims[lb];
        cmf::BlockArray<double, 1> rhsLb  = rhs[lb];
        
        cmf::BlockInfo info = rhs.Mesh()->GetBlockInfo(lb);
        
        double fluxLRAr[10];
        cmf::MdArray<double, 1> fluxLeft (&fluxLRAr[0], 5);
        cmf::MdArray<double, 1> fluxRight(&fluxLRAr[5], 5);
        
        double velGradArr[18];
        cmf::MdArray<double, 2> velGradLeft (&velGradArr[0], 3, 3);
        cmf::MdArray<double, 2> velGradRight(&velGradArr[9], 3, 3);
        
        for (int idir = 0; idir < 3; idir++)
        {
            int dir1 = (idir+1) % 3;
            int dir2 = (idir+2) % 3;
            cell_t ijkCell[3] = {0};
            cell_t ijkCellR[3] = {0};
            cell_t ijkCellL[3] = {0};
            int dijk1[3] = {0};
            int dijk2[3] = {0};
            dijk1[dir1] = 1;
            dijk2[dir2] = 1;
            for (cell_t k = primsLb.kmin; k < primsLb.kmax; k++)
            {
                for (cell_t j = primsLb.jmin; j < primsLb.jmax; j++)
                {
                    for (cell_t i = primsLb.imin; i < primsLb.imax; i++)
                    {
                        ijkCell[0]  = i;
                        ijkCellR[0] = i;
                        ijkCellL[0] = i;
                        ijkCell[1]  = j;
                        ijkCellR[1] = j;
                        ijkCellL[1] = j;
                        ijkCell[2]  = k;
                        ijkCellR[2] = k;
                        ijkCellL[2] = k;
                        
                        ijkCellR[idir]++;
                        ijkCellL[idir]--;
                        
                        double u_UR, u_UL, u_LR, u_LL;
                        for (int uvw = 0; uvw < 3; uvw++)
                        {
                            velGradLeft (uvw, idir) = info.dxInv[idir]*(primsLb(2+uvw, ijkCell[0],  ijkCell[1],  ijkCell[2]) -primsLb(2+uvw, ijkCellL[0], ijkCellL[1], ijkCellL[2]));
                            velGradRight(uvw, idir) = info.dxInv[idir]*(primsLb(2+uvw, ijkCellR[0], ijkCellR[1], ijkCellR[2])-primsLb(2+uvw, ijkCell[0],  ijkCell[1],  ijkCell[2]));
                            
                            //left, dir1
                            u_UR = primsLb(2+uvw, ijkCell [0] + dijk1[0],  ijkCell [1] + dijk1[1],  ijkCell [2] + dijk1[2]);
                            u_UL = primsLb(2+uvw, ijkCellL[0] + dijk1[0],  ijkCellL[1] + dijk1[1],  ijkCellL[2] + dijk1[2]);
                            u_LR = primsLb(2+uvw, ijkCell [0] - dijk1[0],  ijkCell [1] - dijk1[1],  ijkCell [2] - dijk1[2]);
                            u_LL = primsLb(2+uvw, ijkCellL[0] - dijk1[0],  ijkCellL[1] - dijk1[1],  ijkCellL[2] - dijk1[2]);
                            velGradLeft (uvw, dir1) = 0.25*info.dxInv[dir1]*(u_UR + u_LR - u_UL - u_LL);
                            
                            //left, dir2
                            u_UR = primsLb(2+uvw, ijkCell [0] + dijk2[0],  ijkCell [1] + dijk2[1],  ijkCell [2] + dijk2[2]);
                            u_UL = primsLb(2+uvw, ijkCellL[0] + dijk2[0],  ijkCellL[1] + dijk2[1],  ijkCellL[2] + dijk2[2]);
                            u_LR = primsLb(2+uvw, ijkCell [0] - dijk2[0],  ijkCell [1] - dijk2[1],  ijkCell [2] - dijk2[2]);
                            u_LL = primsLb(2+uvw, ijkCellL[0] - dijk2[0],  ijkCellL[1] - dijk2[1],  ijkCellL[2] - dijk2[2]);
                            velGradLeft (uvw, dir2) = 0.25*info.dxInv[dir2]*(u_UR + u_LR - u_UL - u_LL);
                            
                            //right, dir1
                            u_UR = primsLb(2+uvw, ijkCellR[0] + dijk1[0],  ijkCellR[1] + dijk1[1],  ijkCellR[2] + dijk1[2]);
                            u_UL = primsLb(2+uvw, ijkCell [0] + dijk1[0],  ijkCell [1] + dijk1[1],  ijkCell [2] + dijk1[2]);
                            u_LR = primsLb(2+uvw, ijkCellR[0] - dijk1[0],  ijkCellR[1] - dijk1[1],  ijkCellR[2] - dijk1[2]);
                            u_LL = primsLb(2+uvw, ijkCell [0] - dijk1[0],  ijkCell [1] - dijk1[1],  ijkCell [2] - dijk1[2]);
                            velGradRight(uvw, dir1) = 0.25*info.dxInv[dir1]*(u_UR + u_LR - u_UL - u_LL);
                            
                            //right, dir2
                            u_UR = primsLb(2+uvw, ijkCellR[0] + dijk2[0],  ijkCellR[1] + dijk2[1],  ijkCellR[2] + dijk2[2]);
                            u_UL = primsLb(2+uvw, ijkCell [0] + dijk2[0],  ijkCell [1] + dijk2[1],  ijkCell [2] + dijk2[2]);
                            u_LR = primsLb(2+uvw, ijkCellR[0] - dijk2[0],  ijkCellR[1] - dijk2[1],  ijkCellR[2] - dijk2[2]);
                            u_LL = primsLb(2+uvw, ijkCell [0] - dijk2[0],  ijkCell [1] - dijk2[1],  ijkCell [2] - dijk2[2]);
                            velGradRight(uvw, dir2) = 0.25*info.dxInv[dir2]*(u_UR + u_LR - u_UL - u_LL);
                        }
                        
                        double tempGradLeft  = info.dxInv[idir]*(primsLb(1, ijkCell [0], ijkCell [1], ijkCell [2]) - primsLb(1, ijkCellL[0], ijkCellL[1], ijkCellL[2]));
                        double tempGradRight = info.dxInv[idir]*(primsLb(1, ijkCellR[0], ijkCellR[1], ijkCellR[2]) - primsLb(1, ijkCell [0], ijkCell [1], ijkCell[2]));
                        double uLeft[3] = {0};
                        double uRight[3] = {0};
                        double divLeft = 0;
                        double divRight = 0;
                        for (int d = 0; d < 3; d++)
                        {
                            uLeft[d]  = 0.5*(primsLb(2+d, ijkCell [0], ijkCell [1], ijkCell [2]) + primsLb(2+d, ijkCellL[0], ijkCellL[1], ijkCellL[2]));
                            uRight[d] = 0.5*(primsLb(2+d, ijkCellR[0], ijkCellR[1], ijkCellR[2]) + primsLb(2+d, ijkCell [0], ijkCell [1], ijkCell[2]));
                            divLeft  += velGradLeft (d, d);
                            divRight += velGradRight(d, d);
                        }
                        
                        fluxLeft (0) = 0.0;
                        fluxLeft (1) = (vlaw.visc / vlaw.prandtl) * tempGradLeft;
                        fluxLeft (2) = vlaw.visc*(velGradLeft(1, idir) + velGradLeft(idir, 1));
                        fluxLeft (3) = vlaw.visc*(velGradLeft(2, idir) + velGradLeft(idir, 2));
                        fluxLeft (4) = vlaw.visc*(velGradLeft(3, idir) + velGradLeft(idir, 3));
                        fluxLeft (2+idir) -= 0.666666666667*vlaw.visc*divLeft;
                        
                        
                        fluxRight(0) = 0.0;
                        fluxRight(1) = (vlaw.visc / vlaw.prandtl) * tempGradRight;
                        fluxRight(2) = vlaw.visc*(velGradRight(1, idir) + velGradRight(idir, 1));
                        fluxRight(3) = vlaw.visc*(velGradRight(2, idir) + velGradRight(idir, 2));
                        fluxRight(4) = vlaw.visc*(velGradRight(3, idir) + velGradRight(idir, 3));
                        fluxRight(2+idir) -= 0.666666666667*vlaw.visc*divRight;
                        
                        for (int d = 0; d < 3; d++)
                        {
                            fluxLeft (1) +=uLeft[d] *fluxLeft (2+d);
                            fluxRight(1) +=uRight[d]*fluxRight(2+d);
                        }
                        for (int v = 0; v < 5; v++)
                        {
                            rhsLb(v, i, j, k) += info.dxInv[idir]*(fluxRight(v) - fluxLeft(v));
                        }
                    }                    
                }
            }
        }
    }
}

void ComputeRhs_OLD(cmf::CartesianMeshArray& prims, cmf::CartesianMeshArray& rhs, gas::perfect_gas_t<double>& gas, viscous_laws::constant_viscosity_t<double>& vlaw)
{
    ComputeConv(prims, rhs, 2, gas, 0.5);
    ComputeConvDiss(prims, rhs, gas, 0.5);
    ComputeVisc(prims, rhs, vlaw);
}

void Advance(cmf::CartesianMeshArray& prims, cmf::CartesianMeshArray& rhs, double delta_t, gas::perfect_gas_t<double>& gas)
{
    for (auto lb: prims)
    {
        cmf::BlockArray<double, 1> primsLb= prims[lb];
        cmf::BlockArray<double, 1> rhsLb  = rhs[lb];
        for (cmf::cell_t k = primsLb.kmin; k < primsLb.kmax; k++)
        {
            for (cmf::cell_t j = primsLb.jmin; j < primsLb.jmax; j++)
            {
                for (cmf::cell_t i = primsLb.imin; i < primsLb.imax; i++)
                {
                    fluid_state::prim_t<double> prim;
                    for (uint pp = 0; pp < 5; pp++) prim[pp] = primsLb(pp, i, j, k);
                    fluid_state::cons_t<double> cons;
                    convert_state(prim, cons, gas);
                    for (uint pp = 0; pp < 5; pp++) cons[pp] += delta_t * rhsLb(pp, i, j, k);
                    convert_state(cons, prim, gas);
                    for (uint pp = 0; pp < 5; pp++) primsLb(pp, i, j, k) = prim[pp];
                }
            }
        }
    }
}

#endif