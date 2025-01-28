#pragma once

#define arrN 1000

#ifdef WINDOWS
extern "C" __declspec(dllexport) int RENDER(int argc, void **argv);
#else
extern "C" int RENDER(int argc, void **argv);
#endif