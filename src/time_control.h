#pragma once

#include "PTL.h"
#include "cmf.h"
using cmf::strformat;
using cmf::strunformat;


struct time_control : public cmf::ICmfDataBaseReadWriteObject
{
    int numSteps;
    double timestep;
    double time = 0.0;
    int nt = 0;
    
    void Read(PTL::PropertySection& section)
    {
        section["numSteps"].MapTo(&numSteps) = new PTL::PTLInteger(1000, "Number of timesteps to take");
        section["timestep"].MapTo(&timestep) = new PTL::PTLDouble(1e-4, "Physical timestep");
        section.StrictParse();
    }
    
    inline time_control& operator ++ (int dummy)
    {
        nt++;
        time+=timestep;
        return *this;
    }
    
    virtual void ReadFromFile(cmf::ParallelFile& file) override
    {
        strunformat(file.Read(), "nt: {}", nt);
        strunformat(file.Read(), "time: {}", time);
    }
       
    virtual void WriteToFile(cmf::ParallelFile& file) override
    {
        file.Write(strformat("nt: {}", nt));
        file.Write(strformat("time: {}", time));
    }
};

static std::ostream & operator<<(std::ostream & os, const time_control & time)
{
    os << "Timestep: " << time.nt << ", time: " << time.time;
    return os;
}