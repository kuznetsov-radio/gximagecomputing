#pragma once

#ifdef WINDOWS
extern "C" __declspec(dllexport) int ComputeMW(int argc, void **argv);
#else
extern "C" int ComputeMW(int argc, void **argv);
#endif