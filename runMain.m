%复现论文10.1186/s13634-023-01066-3 作者Chuanyun Ouyang1,2, Liming Cai2* , Bin Liu2 and Tianxiang Zhang1,2
%代码作者：陈羿乔 Yiqiao Chen
%日期：2024/4/27 - 2024/4/27

%论文中给的样本标准差计算公式好像不太适用，我自己调了个还可以的常数，后面会关注一下方差估计相关的文献

clc; clear; close all;
format compact; %命令行显示不换行

%==========================================================================
%-----------------------------常量设置与文件导入-----------------------------
%==========================================================================

%读取PPG文件
file1 = readmatrix("data.csv");
file1 = file1(1:length(file1),1);
rowPpg1 = file1; %补零延拓为2幂次的数据
clear file1;

%==========================================================================
%---------------------------------阈值作图对比------------------------------
%==========================================================================
%常量定义
miu_test = 3;
delta_test = 0.01;
lambda_test = 25;

%-----------调整delta与miu两个参数，绘制图像与硬阈值，软阈值对比---------------
x = -100 : 1 : 100;

%新阈值返回系数
coefficent1 = 1 : length(x);
for index = 1 : 1 : length(x)
    coefficent1(index) = improvedWaveletThreshold(x(index), miu_test, delta_test, lambda_test);
end

%-------------------------------硬阈值返回系数------------------------------
coefficent2 = 1 : length(x);
for index = 1 : 1 : length(x)
    if abs(x(index)) > lambda_test
        coefficent2(index) = x(index);
    else
        coefficent2(index) = 0;
    end
end

%-------------------------------软阈值返回系数------------------------------
coefficent3 = 1 : length(x);
for index = 1 : 1 : length(x)
    if x(index) > lambda_test
        coefficent3(index) = x(index) - lambda_test;
    elseif x(index) < -lambda_test
        coefficent3(index) = x(index) + lambda_test;
    else
        coefficent3(index) = 0;
    end
end

%---------------------------------作图-------------------------------------

figure(1);
grid on;
hold on;

plot(x, coefficent1, "r");
plot(x, coefficent2, "b");
plot(x, coefficent3, "k");

legend("new Thresholding", "hard Thresholding", "soft Thresholding", 'Location','northwest')
hold off;


%==========================================================================
%-----------------------------sym4小波分解----------------------------------
%==========================================================================
%论文中使用试错法选定为sym4小波分解
%论文中使用试错法选定分解层数，这里简化为分解为最低层数

%---------------------------------小波分解----------------------------------
[C, L] = wavedec(rowPpg1, 14, "sym4");

%-------------------------------第14层近似系数------------------------------
cA_14 = C(1 : L(1));

%---------------------------第1层到第14细节系数数组--------------------------
startIndex = L(1) + 1;
for level = 2 : 1 : 15
    endPoint = startIndex + (L(level) - 1); %向量拆分的终点：起始点加数据长度
    cD{16 - level, 1} = C(startIndex : endPoint); %拆到第16 - index层时细节系数
    startIndex = endPoint + 1; %起始点迭代为终点加1处
end


%==========================================================================
%-----------------------------新方法阈值去噪--------------------------------
%==========================================================================
%常量定义
miu = 6;
delta = 0.1;

%对每一层细节系数进行阈值去噪
for level = 1 : 1 : 14
    %每一层计算该层的新阈值(我感觉这个阈值计算方式有问题)
    %sigma = sqrt(median(cD{level, 1}) / 0.6745);
    sigma = 4;
    lambda_j = (sigma * sqrt(2 * log(length(rowPpg1)))) / log(level + 1);
    
    for index = 1 : 1 : length(cD{level, 1})
        %使用封装好的阈值操作去噪
         cD{level, 1}(index) = improvedWaveletThreshold(cD{level, 1}(index), miu, delta, lambda_j);
    end
end

%合成新系数向量new_C
new_C = [];
new_C = cat(1, new_C, cA_14);

for level = 14 : -1 : 1
    wj = cD{level, 1}; %分解到第j层时的细节系数向量
    new_C = cat(1, new_C, wj);
end

%==========================================================================
%------------------------------sym4小波重构---------------------------------
%==========================================================================

denoisiedPpg1 = waverec(new_C, L, "sym4");
figure(2);
plot(rowPpg1, "b");

grid on;
hold on;

plot(denoisiedPpg1, "r");

hold off;
legend("RowPPG", "FilteredPPG")


figure(3);
plot(denoisiedPpg1);

















