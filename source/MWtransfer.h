#pragma once

#ifdef WINDOWS
extern "C" __declspec(dllexport) int GET_MW(int argc, void **argv);
#else
extern "C" int GET_MW(int argc, void **argv);
#endif

#define InSize 17
#define OutSize 7
#define LpSize 5
#define RpSize 3

#define InSize_int 16