#pragma once

#define arrN 1000

#ifndef LINUX
extern "C" __declspec(dllexport) int RENDER(int argc, void **argv);
#else
extern "C" int RENDER(int argc, void **argv);
#endif
