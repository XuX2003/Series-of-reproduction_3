E = [0 0 0 0 0 6.594;
    0.025 0 0 0 0 0;
    0.079 0.025 0.015 0 0 0.151;
    0.75 0 0.256 2.773 0 23.27;
    5.74 0 0 19.12 2.272 4.454;
    0 0 0 5.157 29.31 0];
% E 可以根据需求定义
% 节点ij表示j到i的转换
names = ['Sone  ';'Stwo  ';'Sthree';'Sfour ';'Sfive ';'Ssix  '];
L = loopset(E, names);
% Input
% E 生命周期图的邻接矩阵
% names 名称向量
% Output
% adj 生命周期图的邻接矩阵
% nms 输入名称向量
% fbl 一个矩阵，每行包含每个已识别循环的有序元素集（即节点）。矩阵中的循环按循环长度从短到长排序
% els 每个已识别循环的特征弹性系数
loop(L, index);
% Input
% L loopset函数的output
% index 欲查阅的loop循环索引，例：想查看Loop 4的特征弹性系数等，则输入loop(L, 4)
% Output
% 该loop的路径和特征弹性系数