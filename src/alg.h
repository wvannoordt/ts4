#pragma once
#include <concepts>
#include <utility>
#include "cmf.h"
#include "guard_protocol.h"
#include "fluid_state.h"
using cmf::print;

template <class T> concept returns_array_from_coords = requires (T t)
{
    std::is_convertible<decltype(t(std::declval<cmf::Vec3<double>&>()).size()), size_t>::value;
    t(std::declval<cmf::Vec3<double>&>())[size_t()];
};

namespace alg
{
    //trying out some cool c++20 features
    template <
        returns_array_from_coords spatial_callable,
        const guard_protocol::guard_protocol protocol=guard_protocol::exclude_guards
        >
    void fill_array(cmf::CartesianMeshArray& array, const spatial_callable& func)
    {
        for (auto lb: array)
        {
            cmf::BlockArray<double, 1> data = array[lb];
            cmf::BlockInfo info = array.GetBlockInfo(lb);
            cmf::cell_t i0 = info.imin;
            cmf::cell_t i1 = info.imax;
            cmf::cell_t j0 = info.jmin;
            cmf::cell_t j1 = info.jmax;
            cmf::cell_t k0 = info.kmin;
            cmf::cell_t k1 = info.kmax;
            if constexpr (protocol==guard_protocol::include_guards)
            {
                i0 -= info.exchangeI;
                i1 += info.exchangeI;
                j0 -= info.exchangeJ;
                j1 += info.exchangeJ;
                k0 -= info.exchangeK;
                k1 += info.exchangeK;
            }
            cmf::Vec3<double> x(0,0,0);
            for (cmf::cell_t k = k0; k < k1; k++)
            {
                if constexpr(cmf::Dim()==3) x[2] = info.blockBounds[4]+((double)k+0.5)*info.dx[2];
                else x[2] = 0.0;
                for (cmf::cell_t j = j0; j < j1; j++)
                {
                    x[1] = info.blockBounds[2]+((double)j+0.5)*info.dx[1];
                    for (cmf::cell_t i = i0; i < i1; i++)
                    {
                        x[0] = info.blockBounds[0]+((double)i+0.5)*info.dx[0];
                        auto state = func(x);
                        for (uint n = 0; n < state.size(); n++)
                        {
                            data(n, i, j, k) = state[n];
                        }
                    }
                }
            }
        }
    }
}