function [randnum,muHat,sigmaHat,muCI,sigmaCI] = randnorm_wt(data, wt, dimsize)
% [randnum,muHat,sigmaHat,muCI,sigmaCI] = RANDGEN_NORMAL(data, wt, dimsize)
% This function is to generate a fitting to data which has weights
% not necessarily 1 with a normal (gaussian) distribution, then draw
% random numbers that comply with the estimated distribution 
%
% INPUT:
%   data:       raw observations 
%   wt:         weights of observations
%   dimsize:    dimension size of the output randoms
%   seed:          random seed state
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/09/30
% Last modified date:   2020/09/30


% returns estimates of normal distribution parameters, mu and sigma
% under 95% confidence interval
[muHat,sigmaHat,muCI,sigmaCI] = normfit(data, 0.05, [], wt);

% generate random numbers from fitted normal
% rng(seed);
% dimsize = [ceil(sum(wt)),1];  % size of desired randum numbers per sampling
randnum = normrnd(muHat, sigmaHat, dimsize);


