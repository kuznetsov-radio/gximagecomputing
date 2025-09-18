#pragma once

#ifdef WINDOWS
extern "C" __declspec(dllexport) int ComputeEUV(int argc, void **argv);
#else
extern "C" int ComputeEUV(int argc, void **argv);
#endif