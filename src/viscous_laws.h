#pragma once

#include "fluid_state.h"
#include <concepts>

//todo: figure out this concept
// template <class T> concept prim_state_type = std::is_floating_point_v<decltype(T::T())>;

namespace viscous_laws
{
    template <typename dtype> struct constant_viscosity
    {
        constant_viscosity(dtype visc_in) { this->visc = visc_in; }
        constant_viscosity(void) { this->visc = dtype(); }
        
        // template <prim_state_type state_type>dtype calc_visc(const state_type& q)
        template <class state_type>dtype calc_visc(const state_type& q)
        {
            return visc;
        }
        
        dtype visc;
        dtype prandtl;
    };
}