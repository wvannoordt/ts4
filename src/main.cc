#include <array>
#include "cmf.h"
#include "util.h"
#include "alg.h"
#include "gas.h"
#include "ctrs.h"
#include "fluid_state.h"
#include "timing.h"
#include "time_control.h"
#include "viscous_laws.h"
#include "temp/rhs_temp.h"

using cmf::print;
using cmf::strformat;
using cmf::strunformat;
using cmf::ZFill;


int main(int argc, char** argv)
{
	std::string inputFile = util::get_input_file(argc, argv);
	cmf::ReadInput(inputFile);
	cmf::globalSettings = cmf::GlobalSettings(cmf::mainInput["GlobalSettings"]);
	cmf::CreateParallelContext(&argc, &argv);
	cmf::CartesianMeshInputInfo inputInfo(cmf::mainInput["Domain"]);
	cmf::CartesianMesh domain(inputInfo);
	
	time_control control;
	control.Read(cmf::mainInput["Time"]);
	
	bool isRoot = cmf::globalGroup.IsRoot();
	
	auto& prims = domain.DefineVariable("prims", cmf::CmfArrayType::CmfDouble, {5});
	
	prims.ComponentName(0) = "P";
	prims.ComponentName(1) = "T";
	prims.ComponentName(2) = "U";
	prims.ComponentName(3) = "V";
	prims.ComponentName(4) = "W";
	
	auto& rhs = domain.DefineVariable("rhs", cmf::CmfArrayType::CmfDouble, {5});
	
	gas::perfect_gas_t<double> air{.R = 287.15, .gamma = 1.4};
	viscous_laws::constant_viscosity<double> vlaw{.visc = 1e-4};
	
	auto init_condition = [=](cmf::Vec3<double>& x) -> fluid_state::prim_t<double>
	{
		fluid_state::prim_t<double> out;
		out.p() = 1000;
		out.T() = 100;
		out.u() = cos(10*x[0]);
		out.v() = sin(15*x[1]);
		out.w() = 0.0;
		return out;
	};
	
	cmf::CreateDirectory("output");
	prims.ExportFile("output", "initialCondition");
	
	auto init_zero = [](cmf::Vec3<double>& x) -> ctrs::array<double,5> { return ctrs::array<double,5>(0.0); };
	
	alg::fill_array(prims, init_condition);
	
	prims.Exchange();
	
	for (; control.nt <= control.numSteps; control++)
	{
		// calc_rhs(rhs, prims, air, vlaw);
		print("===============================");
		{timing::scoped_tmr_t tmr("zero rhs"); alg::fill_array(rhs, init_zero);}
		{timing::scoped_tmr_t tmr("calc rhs"); ComputeRhs_OLD(prims, rhs, air, vlaw);} // check src/temp/rhs_temp.h for this function
		{timing::scoped_tmr_t tmr("advance");  Advance(prims, rhs, control.timestep, air);}
		{timing::scoped_tmr_t tmr("exchange"); prims.Exchange();}
	}
	
	
	
	return 0;
}
