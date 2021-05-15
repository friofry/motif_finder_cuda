#pragma once

#include "argo_cuda_params.h"


ArgoCudaParams read_ini_file(const char *inifile);

void print_argo_cuda_params(const ArgoCudaParams &params);