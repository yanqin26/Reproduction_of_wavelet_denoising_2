function coefficent = improvedWaveletThreshold(x, miu, delta, lambda)
%论文中提出的参数可调阈值
%输入：1.miu用于控制新阈值输出系数在参数接近阈值时靠近硬阈值的程度（接近竖线），正比关系
%输入：2.delta用于控制新阈值输出系数在较大时时靠近硬阈值的程度（接近斜线，正比关系
%输出：经过阈值操作后的系数

%--------为了避免式子过长导致无法有效编程，将其中部分分式拆开单独编程-----------
%1.重复使用两次的指数项
e_H = exp(delta * (lambda - x)); %x > lambda
e_L = exp(delta * (lambda + x)); %x < lambda

%2.很长的一个分数式
f_long_H = (lambda ^ 2) / (sqrt(x^2 - (2 * x * (exp(miu)) * (exp(lambda-x) - 1)))); %x > lambda
f_long_L = (lambda ^ 2) / (sqrt(x^2 + (2 * x * (exp(miu)) * (exp(lambda+x) - 1)))); %x > lambda

%3.没那么长的一个分数式
f_short_H = (lambda ^ 2) / (x * (exp(delta * (x - lambda)))); %x > lambda
f_short_L = (lambda ^ 2) / ((-x) * (exp((-delta) * (x + lambda)))); %x > lambda

%-----------------------------阈值操作--------------------------------------

if x > lambda
    coefficent = x - (e_H * f_long_H) + ((1-e_H) * f_short_H);
elseif abs(x) <= lambda
    coefficent = 0;
else
    coefficent = x + (e_L * f_long_L) - ((1-e_L) * f_short_L);
end








