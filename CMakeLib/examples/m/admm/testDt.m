clear;


% input1 = csvread('./iter_2nd/y3x_rho.xvx.csv');
% input2 = csvread('./iter_2nd/y3y_rho.xvy.csv');
input1 = csvread('y3xrhovx.csv');
input2 = csvread('y3yrhovy.csv');

[D, Dt]     = defDDt;

tmp = reshape([input1 input2], [size(input2), 2]);
result = Dt(tmp);
