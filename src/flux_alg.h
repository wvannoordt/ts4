#pragma once

#include "cmf.h"
#include "ctrs.h"
#include "fluid_state.h"
namespace flux_alg
{
    
    template <typename dtype> struct weno3
    {
        constexpr static size_t  forward_extent() noexcept {return 2;}
        constexpr static size_t backward_extent() noexcept {return 2;}
        ctrs::array<dtype,5> compute_flux(ctrs::array<fluid_state::prim_t<dtype>*,4> states)
        {
            auto weno3lam = [](double al, double bl, double cl, double dl, double ar, double br, double cr, double dr) -> double
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
            ctrs::array<dtype,5> out(0.0);
            for (int n1 = 0; n1 < 5; n1++)
            {
                out(n1) = weno3lam((*(states[0]))[n1], (*(states[1]))[n1], (*(states[2]))[n1], (*(states[3]))[n1]);
            }
            return out;
        }
    };
    
    template <class state_t, const size_t extent> struct rotating_state_stencil
    {
        ctrs::array<state_t,extent> data;
        ctrs::array<state_t*,extent> permutation;
        size_t next = 0;
        rotating_state_stencil(void)
        {
            this->reset();
        }
        void reset(void)
        {
            for (size_t p = 0; p < extent; p++) permutation[p] = &data[p];
        }
    };
    
    // todo: write a concept for a flux calculation
    template <class line_flux_t> void line_flux_derivative(cmf::CartesianMeshArray& result, cmf::CartesianMeshArray& states, const line_flux_t& flux_calc)
    {
        for (auto lb: states)
        {
            cmf::BlockArray<double, 1> data = states[lb];
            cmf::BlockInfo info = states.GetBlockInfo(lb);
            cmf::cell_t i0 = info.imin;
            cmf::cell_t i1 = info.imax;
            cmf::cell_t j0 = info.jmin;
            cmf::cell_t j1 = info.jmax;
            cmf::cell_t k0 = info.kmin;
            cmf::cell_t k1 = info.kmax;
            constexpr size_t stencil_size = flux_calc.forward_extent()+flux_calc.backward_extent();
            rotating_state_stencil<fluid_state::prim_t<double>,stencil_size> stencil;
            ctrs::array<int,3> dijk(0);
            constexpr size_t numVars = 5;
            for (uint idir = 0; idir < cmf::Dim(); idir++)
            {
                dijk[idir] = 1;
                for (cmf::cell_t k = k0; k < k1; k++)
                {
                    for (cmf::cell_t j = j0; j < j1; j++)
                    {
                        for (cmf::cell_t i = i0; i < i1; i++)
                        {
                            for (int n = 0; n < stencil_size; n++)
                            {
                                for (int v = 0; v < numVars; v++)
                                {
                                    
                                }
                            }
                        }
                    }
                }
                dijk[idir] = 0;
            }
        }
    }
}