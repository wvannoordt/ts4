#pragma once

namespace fluid_state
{
    template <typename rtype> struct prim_t
    {
        static constexpr size_t size(void) noexcept {return 5;}
        rtype data[size()];
        rtype& p() {return data[0];}
        rtype& T() {return data[1];}
        rtype& u() {return data[2];}
        rtype& v() {return data[3];}
        rtype& w() {return data[4];}
        rtype& operator [] (uint idx) {return data[idx];}
    };
    
    template <typename rtype> struct cons_t
    {
        static constexpr size_t size(void) noexcept {return 5;}
        rtype data[size()];
        rtype& rho  () {return data[0];}
        rtype& rho_H() {return data[1];}
        rtype& rho_u() {return data[2];}
        rtype& rho_v() {return data[3];}
        rtype& rho_w() {return data[4];}
        rtype& operator [] (uint idx) {return data[idx];}
    };
    
    template<typename ptype, typename ctype, class gastype> void convert_state(prim_t<ptype>& prim, cons_t<ctype>& cons, const gastype& gas)
    {
        double rho = prim.p() / (gas.R*prim.T());
        double rhoU2 = rho*(prim.u()*prim.u()+prim.v()*prim.v()+prim.w()*prim.w());
        double rhoE = 0.5*rhoU2 + (prim.p()/((gas.gamma - 1.0)));
        double rhoU = rho*prim.u();
        double rhoV = rho*prim.v();
        double rhoW = rho*prim.w();
        cons.rho()   = rho;
        cons.rho_H() = rhoE;
        cons.rho_u()  = rhoU;
        cons.rho_v()  = rhoV;
        cons.rho_w()  = rhoW;
    }
    
    template<typename ptype, typename ctype, class gastype> void convert_state(cons_t<ctype>& cons, prim_t<ptype>& prim, const gastype& gas)
    {
        double rho = cons.rho();
        double invrho = 1.0/rho;
        double u = invrho*cons.rho_u();
        double v = invrho*cons.rho_v();
        double w = invrho*cons.rho_w();
        double rhoU2 = rho*(u*u+v*v+w*w);
        double p = (gas.gamma - 1.0)*(cons.rho_H() - 0.5*rhoU2);
        double T = p/(gas.R*rho);
        prim.p() = p;
        prim.T() = T;
        prim.u() = u;
        prim.v() = v;
        prim.w() = w;
    }
}