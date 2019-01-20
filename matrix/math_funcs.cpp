#include "math_funcs.h"
#include <random>
#include <ctime>

namespace mat
{
    // 随机数引擎,用于生成一个随机的unsigned整数
    // 该引擎对象重载了()运算符
    static std::default_random_engine rand_eng((unsigned int)time(NULL));

    // 生成区间[0, 1)内的 double 类型的随机数
    double random_real()
    {
        return std::uniform_real_distribution<double>(0, 1.0)(rand_eng);
    }

    // 正态分布函数(均值为u,方差为t)
    double random_norm(const double& u, const double& t)
    {
        return std::normal_distribution<double>(u, t)(rand_eng);
    }

    // 符号函数
    double sgn(const double& num)
    {
        return num < 0 ? -1 : 1;
    }
}