clc;clear;close all;

set(0,'DefaultFigureWindowStyle','normal')
format long;

addpath './frequency'
addpath './Colormaps'
addpath './distinguishable_colors/'
addpath './hex_and_rgb_v2'
addpath './fig'
addpath './legendflex'
addpath './setgetpos_V1.2'

cb = hex2rgb('#1f77b4');
cy = hex2rgb('#ff7f0e');
cg = hex2rgb('#2ca02c');
cr = hex2rgb('#d62728');
cp = hex2rgb('#9467bd');
ck = hex2rgb('#000000');
csil = hex2rgb('#7f7f7f');
colororder = [cb;cy;cg;cr;cp;ck;csil];
textwidth = 7;

lw = 0.5;
ms = 3;
fs = 8;