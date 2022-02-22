#pragma once

namespace util
{
	std::string get_input_file(int argc, char** argv)
	{
		if (argc<=1) return "input.ptl";
		std::string output(argv[1]);
		return output;
	}
}