%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Demo of how to use mat2rgbCmpl.
% mat2rgbCmpl takes a complex matrix and converts it to a matrix of RGB
% values for plotting as an image.
% The intensity corresponds to the lightness in HSL
% The hue (color) correponds to the phase.
% Saturation is always set to 0.
% Any number with greater magnitude than the max magnitude chosen becomes
% black or white.
%
% Color representation on monitors is rather atrocious. On excellent
% monitors, there should be a smooth transition everywhere, but most
% monitors apply their own strange scalings making the plot appear
% different on every monitor.
%
% Mark Harfouche
% 2015 Sept 4th
%
%
% Copyright (c) 2015, Mark Harfouche
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright
%   notice, this list of conditions and the following disclaimer.
%2. Redistributions in binary form must reproduce the above copyright
%   notice, this list of conditions and the following disclaimer in the
%   documentation and/or other materials provided with the distribution.
%3. All advertising materials mentioning features or use of this software
%   must display the following acknowledgement:
%   This product includes software developed by the <organization>.
%4. Neither the name of the <organization> nor the
%   names of its contributors may be used to endorse or promote products
%   derived from this software without specific prior written permission.
%THIS SOFTWARE IS PROVIDED BY MARK HARFOUCHE ''AS IS'' AND ANY
%EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL MARK HARFOUCHE BE LIABLE FOR ANY
%DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
%close all

x = linspace(-2, 2, 1024);
y = linspace(-2, 2, 1024);

[x_mesh, y_mesh] = meshgrid(x, y);
z = x_mesh + 1i * y_mesh;

z_black_to_white = mat2rgbCplx(z, 1, 1);
z_white_to_black = mat2rgbCplx(z, 1, 0);

figure
imagesc(x, y, abs(z), 'CData', z_black_to_white)
title('Complex plane z = x + iy, black = 0')
xlabel('Real');
ylabel('Imaginary');
axis(gca, 'image')

figure
surf(x_mesh, y_mesh, abs(z), 'CData', z_white_to_black, 'EdgeColor', 'none')
title('Complex plane z = x + iy, white = 0')
xlabel('Real');
ylabel('Imaginary');
axis(gca, 'image')