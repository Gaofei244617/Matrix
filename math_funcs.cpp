#include "math_funcs.h"

namespace mat
{
    // 符号函数
    double sgn(const double& num)
    {
        return num < 0 ? -1 : 1;
    }
}