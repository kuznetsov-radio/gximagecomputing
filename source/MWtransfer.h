#pragma once

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW(int argc, void **argv);
#else
extern "C" int GET_MW(int argc, void **argv);
#endif