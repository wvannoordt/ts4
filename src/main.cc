#include <array>
#include "cmf.h"
#include "util.h"
#include "alg.h"
#include "gas.h"
#include "ctrs.h"
#include "fluid_state.h"
#include "timing.h"
#include "time_control.h"

// #include "InputParams.h"
// #include "InitialCondition.hpp"
// #include "PrimsCons.hpp"
// #include "Rhs.hpp"
// #include "Integrate.hpp"
// #include "GetTimestep.hpp"

using cmf::print;
using cmf::strformat;
using cmf::strunformat;
using cmf::ZFill;

struct TimeControl
{
	double time;
	int nt;
	int minStep = 0;
	int maxStep = 0;
};

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
	
	gas::perfect_gas_t air{.R = 287.15, .gamma = 1.4};
	
	auto init_condition = [=](cmf::Vec3<double>& x) -> ctrs::array<double,5>
	{
		ctrs::array<double,5> out;
		out[0] = 1000;
		out[1] = 100;
		out[2] = cos(10*x[0]);
		out[3] = sin(15*x[1]);
		out[4] = 0.0;
		return out;
	};
	
	auto init_zero = [](cmf::Vec3<double>& x) -> ctrs::array<double,5>
	{
		return ctrs::array<double,5>(3.0);
	};
	
	alg::fill_array(prims, init_condition);
	alg::fill_array(rhs,   init_zero);
	
	prims.Exchange();
	
	cmf::CreateDirectory("output");
	rhs.ExportFile("output", "initialCondition");
	
	
	return 0;
}
