#pragma once

#include "PTL.h"
#include "cmf.h"

#include "ctrs.h"
#include "gas.h"
#include "viscous_laws.h"
#include "alg.h"
#include "temp/rhs_temp.h"
namespace mms
{
    struct mms_settings_t
    {
        bool doMMS;
        mms_settings_t(PTL::PropertySection& section)
        {
            section["doMMS"].MapTo(&doMMS) = new PTL::PTLBoolean(false, "do method of manufactured solutions");
            section.StrictParse();
        }
    };
    template <typename dtype> struct cns_pergectgas_mms_t
    {
        gas::perfect_gas_t<dtype> gas;
        viscous_laws::constant_viscosity_t<dtype> vlaw;
        const double pi = 3.1415926535;
        const double two_pi = 2.0*pi;
        
        cns_pergectgas_mms_t(gas::perfect_gas_t<dtype>& gas_in, viscous_laws::constant_viscosity_t<dtype>& vlaw_in)
        {
            gas  = gas_in;
            vlaw = vlaw_in;
        }
        
        //P and derivatives
        dtype  P(const dtype& x, const dtype& y, const dtype& z)    const {return 5.0+2.0*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)+sin(4.0*two_pi*z);}
        dtype dP_dx(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*3.0*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y);}
        dtype dP_dy(const dtype& x, const dtype& y, const dtype& z) const {return 2.0*2.0*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y);}
        dtype dP_dz(const dtype& x, const dtype& y, const dtype& z) const {return 4.0*two_pi*cos(4.0*two_pi*z);}

        //T and derivatives
        dtype  T(const dtype& x, const dtype& y, const dtype& z)    const {return 10.0+2.0*cos(2.0*two_pi*x)*sin(3.0*two_pi*y)+sin(4.0*two_pi*z);}
        dtype dT_dx(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*2.0*two_pi*sin(2.0*two_pi*x)*sin(3.0*two_pi*y);}
        dtype dT_dy(const dtype& x, const dtype& y, const dtype& z) const {return  2.0*3.0*two_pi*cos(2.0*two_pi*x)*cos(3.0*two_pi*y);}
        dtype dT_dz(const dtype& x, const dtype& y, const dtype& z) const {return  4.0*two_pi*cos(4.0*two_pi*z);}
        dtype d2T_dx2(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*2.0*2.0*two_pi*two_pi*cos(2.0*two_pi*x)*sin(3.0*two_pi*y);}
        dtype d2T_dy2(const dtype& x, const dtype& y, const dtype& z) const {return  -2.0*3.0*3.0*two_pi*two_pi*cos(2.0*two_pi*x)*sin(3.0*two_pi*y);}
        dtype d2T_dz2(const dtype& x, const dtype& y, const dtype& z) const {return  -4.0*4.0*two_pi*two_pi*sin(4.0*two_pi*z);}
        
        //U and derivatives
        dtype  U(const dtype& x, const dtype& y, const dtype& z)    const {return sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype dU_dx(const dtype& x, const dtype& y, const dtype& z) const {return  3.0*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype dU_dy(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype dU_dz(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(2.0*two_pi*z);}
        
        dtype d2U_dxx(const dtype& x, const dtype& y, const dtype& z) const {return  -3.0*3.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype d2U_dxy(const dtype& x, const dtype& y, const dtype& z) const {return  -3.0*2.0*two_pi*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype d2U_dxz(const dtype& x, const dtype& y, const dtype& z) const {return  -3.0*2.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(2.0*two_pi*z);}
        
        dtype d2U_dyx(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype d2U_dyy(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        dtype d2U_dyz(const dtype& x, const dtype& y, const dtype& z) const {return  2.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(2.0*two_pi*z);}
        
        dtype d2U_dzx(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(2.0*two_pi*z);}
        dtype d2U_dzy(const dtype& x, const dtype& y, const dtype& z) const {return  2.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(2.0*two_pi*z);}
        dtype d2U_dzz(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(2.0*two_pi*z);}
        
        //V and derivatives
        dtype  V(const dtype& x, const dtype& y, const dtype& z)    const {return cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype dV_dx(const dtype& x, const dtype& y, const dtype& z) const {return -3.0*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype dV_dy(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype dV_dz(const dtype& x, const dtype& y, const dtype& z) const {return -3.0*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(3.0*two_pi*z);}
        
        dtype d2V_dxx(const dtype& x, const dtype& y, const dtype& z) const {return -3.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype d2V_dxy(const dtype& x, const dtype& y, const dtype& z) const {return  3.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype d2V_dxz(const dtype& x, const dtype& y, const dtype& z) const {return  3.0*3.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(3.0*two_pi*z);}
        
        dtype d2V_dyx(const dtype& x, const dtype& y, const dtype& z) const {return  2.0*3.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype d2V_dyy(const dtype& x, const dtype& y, const dtype& z) const {return -2.0*2.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        dtype d2V_dyz(const dtype& x, const dtype& y, const dtype& z) const {return  2.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(3.0*two_pi*z);}
        
        dtype d2V_dzx(const dtype& x, const dtype& y, const dtype& z) const {return  3.0*3.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(3.0*two_pi*z);}
        dtype d2V_dzy(const dtype& x, const dtype& y, const dtype& z) const {return  3.0*2.0*two_pi*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(3.0*two_pi*z);}
        dtype d2V_dzz(const dtype& x, const dtype& y, const dtype& z) const {return -3.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(3.0*two_pi*z);}
        
        //W and derivatives
        dtype  W(const dtype& x, const dtype& y, const dtype& z)    const {return sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype dW_dx(const dtype& x, const dtype& y, const dtype& z) const {return  3.0*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype dW_dy(const dtype& x, const dtype& y, const dtype& z) const {return  2.0*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype dW_dz(const dtype& x, const dtype& y, const dtype& z) const {return -4.0*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(4.0*two_pi*z);}
        
        dtype d2W_dxx(const dtype& x, const dtype& y, const dtype& z) const {return  -3.0*3.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype d2W_dxy(const dtype& x, const dtype& y, const dtype& z) const {return   3.0*2.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype d2W_dxz(const dtype& x, const dtype& y, const dtype& z) const {return  -3.0*4.0*two_pi*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(4.0*two_pi*z);}
        
        dtype d2W_dyx(const dtype& x, const dtype& y, const dtype& z) const {return   2.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*cos(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype d2W_dyy(const dtype& x, const dtype& y, const dtype& z) const {return  -2.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        dtype d2W_dyz(const dtype& x, const dtype& y, const dtype& z) const {return  -2.0*4.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(4.0*two_pi*z);}
        
        dtype d2W_dzx(const dtype& x, const dtype& y, const dtype& z) const {return -4.0*3.0*two_pi*two_pi*cos(3.0*two_pi*x)*sin(2.0*two_pi*y)*sin(4.0*two_pi*z);}
        dtype d2W_dzy(const dtype& x, const dtype& y, const dtype& z) const {return -4.0*2.0*two_pi*two_pi*sin(3.0*two_pi*x)*cos(2.0*two_pi*y)*sin(4.0*two_pi*z);}
        dtype d2W_dzz(const dtype& x, const dtype& y, const dtype& z) const {return -4.0*4.0*two_pi*two_pi*sin(3.0*two_pi*x)*sin(2.0*two_pi*y)*cos(4.0*two_pi*z);}
        
        dtype  Rho(const dtype& x, const dtype& y, const dtype& z) const {return P(x,y,z)/(gas.R*T(x,y,z));}
        dtype dRho_dx(const dtype& x, const dtype& y, const dtype& z) const
        {
            return (T(x,y,z)*dP_dx(x,y,z)-P(x,y,z)*dT_dx(x,y,z))/(T(x,y,z)*T(x,y,z)*gas.R);
        }
        
        dtype dRho_dy(const dtype& x, const dtype& y, const dtype& z) const
        {
            return (T(x,y,z)*dP_dy(x,y,z)-P(x,y,z)*dT_dy(x,y,z))/(T(x,y,z)*T(x,y,z)*gas.R);
        }
        
        dtype dRho_dz(const dtype& x, const dtype& y, const dtype& z) const
        {
            return (T(x,y,z)*dP_dz(x,y,z)-P(x,y,z)*dT_dz(x,y,z))/(T(x,y,z)*T(x,y,z)*gas.R);
        }
        
        dtype  H(const dtype& x, const dtype& y, const dtype& z) const
        {
            return 0.5*(U(x,y,z)*U(x,y,z) + V(x,y,z)*V(x,y,z) + W(x,y,z)*W(x,y,z)) + gas.gamma*P(x,y,z)/((gas.gamma-1.0)*Rho(x,y,z));
        }
        dtype dH_dx(const dtype& x, const dtype& y, const dtype& z) const
        {
            return U(x,y,z)*dU_dx(x,y,z)+V(x,y,z)*dV_dx(x,y,z)+W(x,y,z)*dW_dx(x,y,z)
                + gas.gamma*(Rho(x,y,z)*dP_dx(x,y,z)-dRho_dx(x,y,z)*P(x,y,z))/(Rho(x,y,z)*Rho(x,y,z)*(gas.gamma-1.0));
        }
        dtype dH_dy(const dtype& x, const dtype& y, const dtype& z) const
        {
            return U(x,y,z)*dU_dy(x,y,z)+V(x,y,z)*dV_dy(x,y,z)+W(x,y,z)*dW_dy(x,y,z)
                + gas.gamma*(Rho(x,y,z)*dP_dy(x,y,z)-dRho_dy(x,y,z)*P(x,y,z))/(Rho(x,y,z)*Rho(x,y,z)*(gas.gamma-1.0));
        }
        dtype dH_dz(const dtype& x, const dtype& y, const dtype& z) const
        {
            return U(x,y,z)*dU_dz(x,y,z)+V(x,y,z)*dV_dz(x,y,z)+W(x,y,z)*dW_dz(x,y,z)
                + gas.gamma*(Rho(x,y,z)*dP_dz(x,y,z)-dRho_dz(x,y,z)*P(x,y,z))/(Rho(x,y,z)*Rho(x,y,z)*(gas.gamma-1.0));
        }
        
        ctrs::array<dtype,5> conv_rhs(const dtype& x, const dtype& y, const dtype& z) const
        {
            ctrs::array<dtype,5> rhsAr;
            rhsAr[0] = 
                -U(x,y,z)*dRho_dx(x,y,z)-dU_dx(x,y,z)*Rho(x,y,z)
                -V(x,y,z)*dRho_dy(x,y,z)-dV_dy(x,y,z)*Rho(x,y,z)
                -W(x,y,z)*dRho_dz(x,y,z)-dW_dz(x,y,z)*Rho(x,y,z);
            rhsAr[1] = 
                -dRho_dx(x,y,z)*U    (x,y,z)*H    (x,y,z)
                -Rho    (x,y,z)*dU_dx(x,y,z)*H    (x,y,z)
                -Rho    (x,y,z)*U    (x,y,z)*dH_dx(x,y,z)
                -dRho_dy(x,y,z)*V    (x,y,z)*H    (x,y,z)
                -Rho    (x,y,z)*dV_dy(x,y,z)*H    (x,y,z)
                -Rho    (x,y,z)*V    (x,y,z)*dH_dy(x,y,z)
                -dRho_dz(x,y,z)*W    (x,y,z)*H    (x,y,z)
                -Rho    (x,y,z)*dW_dz(x,y,z)*H    (x,y,z)
                -Rho    (x,y,z)*W    (x,y,z)*dH_dz(x,y,z);
            rhsAr[2] = 
                -dRho_dx(x,y,z)*U    (x,y,z)*U    (x,y,z)
                -Rho    (x,y,z)*dU_dx(x,y,z)*U    (x,y,z)
                -Rho    (x,y,z)*U    (x,y,z)*dU_dx(x,y,z)
                -dRho_dy(x,y,z)*V    (x,y,z)*U    (x,y,z)
                -Rho    (x,y,z)*dV_dy(x,y,z)*U    (x,y,z)
                -Rho    (x,y,z)*V    (x,y,z)*dU_dy(x,y,z)
                -dRho_dz(x,y,z)*W    (x,y,z)*U    (x,y,z)
                -Rho    (x,y,z)*dW_dz(x,y,z)*U    (x,y,z)
                -Rho    (x,y,z)*W    (x,y,z)*dU_dz(x,y,z)
                -dP_dx  (x,y,z);
            rhsAr[3] = 
                -dRho_dx(x,y,z)*U    (x,y,z)*V    (x,y,z)
                -Rho    (x,y,z)*dU_dx(x,y,z)*V    (x,y,z)
                -Rho    (x,y,z)*U    (x,y,z)*dV_dx(x,y,z)
                -dRho_dy(x,y,z)*V    (x,y,z)*V    (x,y,z)
                -Rho    (x,y,z)*dV_dy(x,y,z)*V    (x,y,z)
                -Rho    (x,y,z)*V    (x,y,z)*dV_dy(x,y,z)
                -dRho_dz(x,y,z)*W    (x,y,z)*V    (x,y,z)
                -Rho    (x,y,z)*dW_dz(x,y,z)*V    (x,y,z)
                -Rho    (x,y,z)*W    (x,y,z)*dV_dz(x,y,z)
                -dP_dy  (x,y,z);
            rhsAr[4] = 
                -dRho_dx(x,y,z)*U    (x,y,z)*W    (x,y,z)
                -Rho    (x,y,z)*dU_dx(x,y,z)*W    (x,y,z)
                -Rho    (x,y,z)*U    (x,y,z)*dW_dx(x,y,z)
                -dRho_dy(x,y,z)*V    (x,y,z)*W    (x,y,z)
                -Rho    (x,y,z)*dV_dy(x,y,z)*W    (x,y,z)
                -Rho    (x,y,z)*V    (x,y,z)*dW_dy(x,y,z)
                -dRho_dz(x,y,z)*W    (x,y,z)*W    (x,y,z)
                -Rho    (x,y,z)*dW_dz(x,y,z)*W    (x,y,z)
                -Rho    (x,y,z)*W    (x,y,z)*dW_dz(x,y,z)
                -dP_dz  (x,y,z);
            return rhsAr;
        }
        
        dtype Tau_xx(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(dU_dx(x,y,z)+dU_dx(x,y,z))+vlaw.beta*(dU_dx(x,y,z)+dV_dy(x,y,z)+dW_dz(x,y,z)); }
        dtype Tau_xy(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(dU_dy(x,y,z)+dV_dx(x,y,z)); }
        dtype Tau_xz(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(dU_dz(x,y,z)+dW_dx(x,y,z)); }
        dtype Tau_yx(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(dV_dx(x,y,z)+dU_dy(x,y,z)); }
        dtype Tau_yy(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(dV_dy(x,y,z)+dV_dy(x,y,z))+vlaw.beta*(dU_dx(x,y,z)+dV_dy(x,y,z)+dW_dz(x,y,z)); }
        dtype Tau_yz(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(dV_dz(x,y,z)+dW_dy(x,y,z)); }
        dtype Tau_zx(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(dW_dx(x,y,z)+dU_dz(x,y,z)); }
        dtype Tau_zy(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(dW_dy(x,y,z)+dV_dz(x,y,z)); }
        dtype Tau_zz(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(dW_dz(x,y,z)+dW_dz(x,y,z))+vlaw.beta*(dU_dx(x,y,z)+dV_dy(x,y,z)+dW_dz(x,y,z)); }
        
        dtype dTau_xx_dx(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2U_dxx(x,y,z)+d2U_dxx(x,y,z))+vlaw.beta*(d2U_dxx(x,y,z)+d2V_dyx(x,y,z)+d2W_dzx(x,y,z)); }
        dtype dTau_xx_dy(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2U_dxy(x,y,z)+d2U_dxy(x,y,z))+vlaw.beta*(d2U_dxy(x,y,z)+d2V_dyy(x,y,z)+d2W_dzy(x,y,z)); }
        dtype dTau_xx_dz(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2U_dxz(x,y,z)+d2U_dxz(x,y,z))+vlaw.beta*(d2U_dxz(x,y,z)+d2V_dyz(x,y,z)+d2W_dzz(x,y,z)); }
        
        dtype dTau_xy_dx(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2U_dyx(x,y,z)+d2V_dxx(x,y,z)); }
        dtype dTau_xy_dy(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2U_dyy(x,y,z)+d2V_dxy(x,y,z)); }
        dtype dTau_xy_dz(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2U_dyz(x,y,z)+d2V_dxz(x,y,z)); }
        
        dtype dTau_xz_dx(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2U_dzx(x,y,z)+d2W_dxx(x,y,z)); }
        dtype dTau_xz_dy(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2U_dzy(x,y,z)+d2W_dxy(x,y,z)); }
        dtype dTau_xz_dz(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2U_dzz(x,y,z)+d2W_dxz(x,y,z)); }
        
        dtype dTau_yx_dx(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2V_dxx(x,y,z)+d2U_dyx(x,y,z)); }
        dtype dTau_yx_dy(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2V_dxy(x,y,z)+d2U_dyy(x,y,z)); }
        dtype dTau_yx_dz(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2V_dxz(x,y,z)+d2U_dyz(x,y,z)); }
        
        dtype dTau_yy_dx(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2V_dyx(x,y,z)+d2V_dyx(x,y,z))+vlaw.beta*(d2U_dxx(x,y,z)+d2V_dyx(x,y,z)+d2W_dzx(x,y,z)); }
        dtype dTau_yy_dy(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2V_dyy(x,y,z)+d2V_dyy(x,y,z))+vlaw.beta*(d2U_dxy(x,y,z)+d2V_dyy(x,y,z)+d2W_dzy(x,y,z)); }
        dtype dTau_yy_dz(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2V_dyz(x,y,z)+d2V_dyz(x,y,z))+vlaw.beta*(d2U_dxz(x,y,z)+d2V_dyz(x,y,z)+d2W_dzz(x,y,z)); }
        
        dtype dTau_yz_dx(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2V_dzx(x,y,z)+d2W_dyx(x,y,z)); }
        dtype dTau_yz_dy(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2V_dzy(x,y,z)+d2W_dyy(x,y,z)); }
        dtype dTau_yz_dz(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2V_dzz(x,y,z)+d2W_dyz(x,y,z)); }
        
        dtype dTau_zx_dx(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2W_dxx(x,y,z)+d2U_dzx(x,y,z)); }
        dtype dTau_zx_dy(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2W_dxy(x,y,z)+d2U_dzy(x,y,z)); }
        dtype dTau_zx_dz(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2W_dxz(x,y,z)+d2U_dzz(x,y,z)); }
        
        dtype dTau_zy_dx(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2W_dyx(x,y,z)+d2V_dzx(x,y,z)); }
        dtype dTau_zy_dy(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2W_dyy(x,y,z)+d2V_dzy(x,y,z)); }
        dtype dTau_zy_dz(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2W_dyz(x,y,z)+d2V_dzz(x,y,z)); }
        
        dtype dTau_zz_dx(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2W_dzx(x,y,z)+d2W_dzx(x,y,z))+vlaw.beta*(d2U_dxx(x,y,z)+d2V_dyx(x,y,z)+d2W_dzx(x,y,z)); }
        dtype dTau_zz_dy(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2W_dzy(x,y,z)+d2W_dzy(x,y,z))+vlaw.beta*(d2U_dxy(x,y,z)+d2V_dyy(x,y,z)+d2W_dzy(x,y,z)); }
        dtype dTau_zz_dz(const dtype& x, const dtype& y, const dtype& z) const { return vlaw.visc*(d2W_dzz(x,y,z)+d2W_dzz(x,y,z))+vlaw.beta*(d2U_dxz(x,y,z)+d2V_dyz(x,y,z)+d2W_dzz(x,y,z)); }
        
        dtype dQ_x_dx(const dtype& x, const dtype& y, const dtype& z) const { return (vlaw.visc/vlaw.prandtl)*d2T_dx2(x,y,z); }
        dtype dQ_y_dy(const dtype& x, const dtype& y, const dtype& z) const { return (vlaw.visc/vlaw.prandtl)*d2T_dy2(x,y,z); }
        dtype dQ_z_dz(const dtype& x, const dtype& y, const dtype& z) const { return (vlaw.visc/vlaw.prandtl)*d2T_dz2(x,y,z); }
        
        ctrs::array<dtype,5> visc_rhs(const dtype& x, const dtype& y, const dtype& z) const
        {
            ctrs::array<dtype,5> rhsAr;
            rhsAr[0] = 0.0;
            rhsAr[1] = 
                U(x,y,z)*dTau_xx_dx(x,y,z)+dU_dx(x,y,z)*Tau_xx(x,y,z)+
                V(x,y,z)*dTau_xy_dx(x,y,z)+dV_dx(x,y,z)*Tau_xy(x,y,z)+
                W(x,y,z)*dTau_xz_dx(x,y,z)+dW_dx(x,y,z)*Tau_xz(x,y,z)+
                U(x,y,z)*dTau_yx_dy(x,y,z)+dU_dy(x,y,z)*Tau_yx(x,y,z)+
                V(x,y,z)*dTau_yy_dy(x,y,z)+dV_dy(x,y,z)*Tau_yy(x,y,z)+
                W(x,y,z)*dTau_yz_dy(x,y,z)+dW_dy(x,y,z)*Tau_yz(x,y,z)+
                U(x,y,z)*dTau_zx_dz(x,y,z)+dU_dz(x,y,z)*Tau_zx(x,y,z)+
                V(x,y,z)*dTau_zy_dz(x,y,z)+dV_dz(x,y,z)*Tau_zy(x,y,z)+
                W(x,y,z)*dTau_zz_dz(x,y,z)+dW_dz(x,y,z)*Tau_zz(x,y,z)-
                dQ_x_dx(x,y,z)-dQ_y_dy(x,y,z)-dQ_z_dz(x,y,z);
            rhsAr[2] = dTau_xx_dx(x,y,z)+dTau_xy_dy(x,y,z)+dTau_xz_dz(x,y,z);
            rhsAr[3] = dTau_yx_dx(x,y,z)+dTau_yy_dy(x,y,z)+dTau_yz_dz(x,y,z);
            rhsAr[4] = dTau_zx_dx(x,y,z)+dTau_zy_dy(x,y,z)+dTau_zz_dz(x,y,z);
            return rhsAr;
        }
        
        fluid_state::prim_t<dtype> testFcn(const dtype& x, const dtype& y, const dtype& z) const
        {
            fluid_state::prim_t<dtype> primsAr;
            primsAr[0] = P(x, y, z);
            primsAr[1] = T(x, y, z);
            primsAr[2] = U(x, y, z);
            primsAr[3] = V(x, y, z);
            primsAr[4] = W(x, y, z);
            return primsAr;
        }
    };
    
    template <typename mms_type> void  run_mms_test(cmf::CartesianMeshArray& prims, cmf::CartesianMeshArray& rhs, mms_type& mms_obj)
    {
        using cmf::print;
        print("Running MMS test...");
        
        cmf::CartesianMesh& mesh = *(prims.Mesh());
        auto& rhs_mms = mesh.DefineVariable("rhs_mms", cmf::CmfArrayType::CmfDouble, {5});
        
        auto zero_all = [](cmf::Vec3<double>& x)  -> ctrs::array<double,5>       { return ctrs::array<double,5>(0.0);         };
        auto visc_rhs = [=](cmf::Vec3<double>& x) -> ctrs::array<double,5>       { return mms_obj.visc_rhs(x[0], x[1], x[2]); };
        auto conv_rhs = [=](cmf::Vec3<double>& x) -> ctrs::array<double,5>       { return mms_obj.conv_rhs(x[0], x[1], x[2]); };
        auto test_fcn = [=](cmf::Vec3<double>& x) -> fluid_state::prim_t<double> { return mms_obj.testFcn(x[0], x[1], x[2]);  };
        
        alg::fill_array(prims, test_fcn, guard_protocol::include_guards);
        
        print("Viscous...");
        alg::fill_array(rhs, zero_all, guard_protocol::include_guards);
        ComputeVisc(prims, rhs, mms_obj.vlaw);
        alg::fill_array(rhs_mms, visc_rhs, guard_protocol::include_guards);
        rhs.ExportFile("output", "rhs_visc_num");
        rhs_mms.ExportFile("output", "rhs_visc_exact");
        
        print("Convective (diss)...");
        alg::fill_array(rhs, zero_all, guard_protocol::include_guards);
        alg::fill_array(rhs_mms, conv_rhs, guard_protocol::include_guards);
        ComputeConvDiss(prims, rhs, mms_obj.gas, 1.0);
        rhs.ExportFile("output", "rhs_conv_diss_num");
        rhs_mms.ExportFile("output", "rhs_conv_diss_exact");
        
        print("Convective (cent)...");
        alg::fill_array(rhs, zero_all, guard_protocol::include_guards);
        ComputeConv(prims, rhs, 2, mms_obj.gas, 1.0);
        rhs.ExportFile("output", "rhs_conv_cent_num");
        rhs_mms.ExportFile("output", "rhs_conv_cent_exact");
        
        print("MMS done.");
    }
}