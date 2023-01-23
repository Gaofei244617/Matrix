#ifndef _MATH_FUNCS_H_
#define _MATH_FUNCS_H_

namespace mat
{
	// 生成区间[0, 1)内的 double 类型的随机数
	double random_real();

	// 符号函数
	double sgn(const double& num);

	// 正态分布函数(均值为u,方差为t)
	double random_norm(const double& u, const double& t);
}

#endif
