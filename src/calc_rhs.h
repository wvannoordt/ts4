#pragma once

#include "cmf.h"
#include "flux_alg.h"

template <class gas_type>
void calc_conv(cmf::CartesianMeshArray& rhs, cmf::CartesianMeshArray& prims, const gas_type& gas)
{
    //TODO: write this
    flux_alg::weno3<double> w3;
    flux_alg::line_flux_derivative(rhs, prims, w3);
}

template <class gas_type, class viscous_law>
void calc_visc(cmf::CartesianMeshArray& rhs, cmf::CartesianMeshArray& prims, const gas_type& gas, const viscous_law& visc_law)
{
    //TODO: write this
}

template <class gas_type, class viscous_law>
void calc_rhs(cmf::CartesianMeshArray& rhs, cmf::CartesianMeshArray& prims, const gas_type& gas, const viscous_law& visc_law)
{
    calc_conv(rhs, prims, gas);
    // calc_visc(rhs, prims, gas, visc_law);
}