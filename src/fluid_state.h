#pragma once

namespace fluid_state
{
    template <typename rtype> struct prim_t
    {
        const static size_t size = 5;
        rtype data[size];
        rtype& p() {return data[0];}
        rtype& T() {return data[1];}
        rtype& u() {return data[2];}
        rtype& v() {return data[3];}
        rtype& w() {return data[4];}
        rtype& operator [] (uint idx) {return data[idx];}
    };
    
    template <typename rtype> struct cons_t
    {
        const static size_t size = 5;
        rtype data[size];
        rtype& rho  () {return data[0];}
        rtype& rho_H() {return data[1];}
        rtype& rho_u() {return data[2];}
        rtype& rho_v() {return data[3];}
        rtype& rho_w() {return data[4];}
        rtype& operator [] (uint idx) {return data[idx];}
    };
}