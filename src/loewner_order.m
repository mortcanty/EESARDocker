function lo = loewner_order(X,Y)

% Loewner order of (1x1,) 2x2 or 3x3 Hermitian matrices X and Y: is X > Y?
%
% Loewner order: X > Y if (X-Y) is positive definite
%
% lo = loewner_order(X,Y)
%
% SAR input images are images of the (row-wise) upper-right elements of
% the covariance or coherency matrix, e.g. for full polarimetry:
% a vector with nine real (i.e., non-complex) elements, namely
% [ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*) ...
%  ShvShv* Re(ShvSvv*) Im(ShvSvv*) SvvSvv*]
%
% Input
%   X       - 3-D array of (1x1,) 2x2 or 3x3 complex covariance matrix
%             (polarimetric SAR image in covariance matrix form,
%             reading order is alphabetical (as read by freadenvisar):
%               hhhh (real), hhhv (complex), hhvv (complex),
%                            hvhv (real),    hvvv (complex),
%                                            vvvv (real)
%   Y       - 3-D array of (1x1,) 2x2 or 3x3 complex covariance matrix
%
% Output
%   lo - struct with fields
%           eigvals	- eigenvalues of X-Y
%           posdef	- X-Y is positive definite (logical)
%           negdef	- X-Y is negative definite (logical)
%
% Reference
%
% Allan A. Nielsen, Henning Skriver and Knut Conradsen (2019).
% "The Loewner Order Applied to Change Detection in Sentinel-1 and
% Radarsat-2 Data".
% Submitted

% (c) Copyright 2018-2019
% Allan Aasbjerg Nielsen, PhD, MSc
% alan@dtu.dk, people.compute.dtu.dk/alan
% 3-21 Sep 2018, 13 Apr 2019

if nargin<2
    help loewner_order
    error('loewner_order: wrong input');
end
if nargout>1
    help loewner_order
    error('loewner_order: wrong output');
end

if ischar(X), X = freadenvisar(X); end
if ischar(Y), Y = freadenvisar(Y); end

pX = size(X);
pY = size(Y);
if ~all(pX==pY)
    error('loewner_order: X and Y must have same dimensions');
end
nrows = pX(1);
ncols = pX(2);
if  ndims(X)==3
    nvars = pX(3);
elseif ismatrix(X)
    warning('scalar input, no reason to use loewner_order')
    nvars = 1;
else
    error('loewner_order: wrong sizes of input images');
end

if nvars==9 % full pol
%     fprintf('loewner_order: full pol\n')
    p = 3;
elseif nvars==4 % dual pol
%     fprintf('loewner_order: dual pol\n')
    p = 2;
elseif nvars==2 % dual pol, diagonal only
%     fprintf('loewner_order: dual pol, diagonal only\n')
    p = 2;
elseif nvars==1 % single channel
%     fprintf('loewner_order: single pol - makes no sense\n')
%     return
    p = 1;
elseif nvars==5 % full pol, azimuthal symmetry
%     fprintf('loewner_order: full pol, azimuthal symmetry\n')
    p = 3;
elseif nvars==3 % full pol, diagonal only
%     fprintf('loewner_order: full pol, diagonal only\n')    
    p = 3;
else
    error('loewner_order: wrong dimensionality of input images')
end

lo.eigvals = nan(nrows,ncols,p);

if nvars==2 || nvars==3 || nvars==1
    X = X-Y;
    clear Y
    lo.eigvals = sort(X,3,'descend');
elseif nvars==4 % dual pol
    X = X-Y;
    clear Y
    disc = sqrt((X(:,:,1)-X(:,:,4)).^2+4*(X(:,:,2).^2+X(:,:,3).^2));
    lo.eigvals(:,:,1) = 0.5*(X(:,:,1)+X(:,:,4) + disc);
    lo.eigvals(:,:,2) = 0.5*(X(:,:,1)+X(:,:,4) - disc);
    clear disc
elseif nvars==5 % azimuthal symmetry
    X = X-Y;
    clear Y
    disc = sqrt((X(:,:,1)-X(:,:,5)).^2+4*(X(:,:,2).^2+X(:,:,3).^2));
    lo.eigvals(:,:,1) = 0.5*(X(:,:,1)+X(:,:,5) + disc);
    lo.eigvals(:,:,2) = X(:,:,4);
    lo.eigvals(:,:,3) = 0.5*(X(:,:,1)+X(:,:,5) - disc);
    lo.eigvals = sort(lo.eigvals,3,'descend'); % to be on the safe side
    clear disc
else % full pol
    h = waitbar(0,'Implementation with for-loops over rows and cols...','name','Loewner order');
    for cc=1:ncols
        waitbar(cc/ncols,h);
        for rr=1:nrows
            x = squeeze(X(rr,cc,:));
            y = squeeze(Y(rr,cc,:));
% %             if sum(isnan(x))>0 || sum(isnan(y))>0 % are there missing values?
% % %             if ~(any(x) || any(y)) % are there missing values?
% %                 lo.eigvals(rr,cc,:) = nan(p,1); % zeros(p,1);
% %                 continue
% %             end
%             if nvars==9 % full pol
                xm = [x(1) x(2)+1i*x(3) x(4)+1i*x(5); x(2)-1i*x(3) x(6) x(7)+1i*x(8); x(4)-1i*x(5) x(7)-1i*x(8) x(9)];
                ym = [y(1) y(2)+1i*y(3) y(4)+1i*y(5); y(2)-1i*y(3) y(6) y(7)+1i*y(8); y(4)-1i*y(5) y(7)-1i*y(8) y(9)];
% %             elseif nvars==4 % dual pol
% %                 xm = [x(1) x(2)+1i*x(3); x(2)-1i*x(3) x(4)];
% %                 ym = [y(1) y(2)+1i*y(3); y(2)-1i*y(3) y(4)];
% %             elseif nvars==5 % full pol, azimuthal symmetry
% %                 xm = [x(1) 0 x(2)+1i*x(3); 0 x(4) 0; x(2)-1i*x(3) 0 x(5)];
% %                 ym = [y(1) 0 y(2)+1i*y(3); 0 y(4) 0; y(2)-1i*y(3) 0 y(5)];
% %             end
            lo.eigvals(rr,cc,:) = sort(eig(xm-ym),'descend'); % xm >= ym if (xm - ym) is positive semidefinite
%             lo.eigvals(rr,cc,:) = svd(xm-ym,'econ'); % xm >= ym if (xm - ym) is positive semidefinite
        end 
    end
    close(h)
end
clear X Y
lo.posdef = lo.eigvals(:,:,p) > 0;
lo.negdef = lo.eigvals(:,:,1) < 0;
% possible remaining obs are indefinite: ~(lo.posdef | lo.negdef)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while false
    
figure, imshow(255*cat(3,lo.posdef,lo.negdef,zeros(size(lo.posdef))))

figure, imshow(lo.posdef);

% To be combined with p-value from complex Wishart-based change detector, e.g.,
figure, imshow(and(wc(1).P>0.9999,lo.posdef));
figure, imshow(and(wc(1).P>0.9999,min(lo.eigvals,[],3)>0));

figure, imshow(255*(cat(3,wc(1).P>0.9999,lo.posdef,lo.posdef)))
figure, imshow(255*(cat(3,wc(1).P>0.9999,lo.posdef,zeros(size(lo.posdef)))))

% Fraport
cd fraportfloat
x1 = imread('sentinel1_VVVH_02_warp.tif');
x2 = imread('sentinel1_VVVH_03_warp.tif');
info = geotiffinfo('sentinel1_VVVH_02_warp.tif');
cd ..

wc = wcRjl(cat(4,x1,x2),4.4,'ddiag',false);
lo = loewner_order(x1,x2);
% lo.posdef = min(x1-x2,[],3) > 0;

pval = 0.99995;
figure
h1 = subplot_tight(2,2,1,[0.02,0.01]);
imshow(10*log10(x1(:,:,1)),[-24 6])
% title('VV t_1')
h2 = subplot_tight(2,2,2,[0.02,0.01]);
imshow(10*log10(x2(:,:,1)),[-24 6])
% title('VV t_2')
h3 = subplot_tight(2,2,3,[0.02,0.01]);
imshow(wc.P,[0 1])
% imshow(lo.posdef)
% title('t_2-t_1 pos def')
h4 = subplot_tight(2,2,4,[0.02,0.01]);
% imshow(cat(3,and(~lo.negdef,wc.P>pval),and(~lo.posdef,wc.P>pval),zeros(size(lo.posdef))))
imshow(cat(3,and(lo.posdef,wc.P>pval),and(lo.negdef,wc.P>pval),zeros(size(lo.posdef))))
% imshow(255*(catf change, G: t_2-t_1 pos def')
linkaxes([h1,h2,h3,h4],'xy')

% follow by "File | Open" in Google Earth Pro
outim = uint8(255*cat(3,and(lo.posdef,wc.P>0.9999),and(~lo.posdef,wc.P>0.9999),zeros(size(lo.posdef))));
geotiffwrite('wc_lo_change',outim,info.SpatialRef,...
    'CoordRefSysCode',info.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey)


% EMISAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd Foulum
L62 = flipud(freadenvisar('fl062_l/fl062_l','b'));
L63 = flipud(freadenvisar('fl063_l/fl063_l','b'));
L64 = flipud(freadenvisar('fl064_l/fl064_l','b'));
L65 = flipud(freadenvisar('fl065_l/fl065_l','b'));
L68 = flipud(freadenvisar('fl068_l/fl068_l','b'));
L74 = flipud(freadenvisar('fl074_l/fl074_l','b'));
cd ..

[nrows,ncols,~] = size(L65);

wc1 = wcRjl(cat(4,L62,L63),13,'full',false);
wc2 = wcRjl(cat(4,L63,L64),13,'full',false);
wc3 = wcRjl(cat(4,L64,L65),13,'full',false);
wc4 = wcRjl(cat(4,L65,L68),13,'full',false);
wc5 = wcRjl(cat(4,L68,L74),13,'full',false);

wc = wcRjl(cat(4,L62,L63,L64,L65,L68,L74),13,'full',false);

lo1 = loewner_order(L62,L63);
lo2 = loewner_order(L63,L64);
lo3 = loewner_order(L64,L65);
lo4 = loewner_order(L65,L68);
lo5 = loewner_order(L68,L74);

figure('Name','Loewner order, posdef (R), negdef (G), indef (K)')
h1 = subplot_tight(2,3,1,[0.02,0.01]);
imshow(255*cat(3,lo1.posdef,lo1.negdef,zeros(nrows,ncols)))
h2 = subplot_tight(2,3,2,[0.02,0.01]);
imshow(255*cat(3,lo2.posdef,lo2.negdef,zeros(nrows,ncols)))
h3 = subplot_tight(2,3,3,[0.02,0.01]);
imshow(255*cat(3,lo3.posdef,lo3.negdef,zeros(nrows,ncols)))
h4 = subplot_tight(2,3,4,[0.02,0.01]);
imshow(255*cat(3,lo4.posdef,lo4.negdef,zeros(nrows,ncols)))
h5 = subplot_tight(2,3,5,[0.02,0.01]);
imshow(255*cat(3,lo5.posdef,lo5.negdef,zeros(nrows,ncols)))
linkaxes([h1,h2,h3,h4,h5],'xy')

figure('Name','Loewner order, ~negdef (R), ~posdef (G), indef (Y)')
h1 = subplot_tight(2,3,1,[0.02,0.01]);
imshow(255*cat(3,~lo1.negdef,~lo1.posdef,zeros(nrows,ncols)))
h2 = subplot_tight(2,3,2,[0.02,0.01]);
imshow(255*cat(3,~lo2.negdef,~lo2.posdef,zeros(nrows,ncols)))
h3 = subplot_tight(2,3,3,[0.02,0.01]);
imshow(255*cat(3,~lo3.negdef,~lo3.posdef,zeros(nrows,ncols)))
h4 = subplot_tight(2,3,4,[0.02,0.01]);
imshow(255*cat(3,~lo4.negdef,~lo4.posdef,zeros(nrows,ncols)))
h5 = subplot_tight(2,3,5,[0.02,0.01]);
imshow(255*cat(3,~lo5.negdef,~lo5.posdef,zeros(nrows,ncols)))
linkaxes([h1,h2,h3,h4,h5],'xy')

% figure
% h1 = subplot_tight(2,3,1,[0.02,0.01]);
% imshow(255*cat(3,lo1.posdef,lo1.negdef,and(~lo1.posdef,~lo1.negdef)))
% h2 = subplot_tight(2,3,2,[0.02,0.01]);
% imshow(255*cat(3,lo2.posdef,lo2.negdef,and(~lo2.posdef,~lo2.negdef)))
% h3 = subplot_tight(2,3,3,[0.02,0.01]);
% imshow(255*cat(3,lo3.posdef,lo3.negdef,and(~lo3.posdef,~lo3.negdef)))
% h4 = subplot_tight(2,3,4,[0.02,0.01]);
% imshow(255*cat(3,lo4.posdef,lo4.negdef,and(~lo4.posdef,~lo4.negdef)))
% h5 = subplot_tight(2,3,5,[0.02,0.01]);
% imshow(255*cat(3,lo5.posdef,lo5.negdef,and(~lo5.posdef,~lo5.negdef)))
% linkaxes([h1,h2,h3,h4,h5],'xy')

pval = 0.9999; %pval = 0.99995
figure('Name','Wishart change signif, and posdef (R), and negdef (G)')
h1 = subplot_tight(2,3,1,[0.02,0.01]);
imshow(255*(cat(3,and(wc1.P>pval,lo1.posdef),and(wc1.P>pval,lo1.negdef),zeros(nrows,ncols))))
h2 = subplot_tight(2,3,2,[0.02,0.01]);
imshow(255*(cat(3,and(wc2.P>pval,lo2.posdef),and(wc2.P>pval,lo2.negdef),zeros(nrows,ncols))))
h3 = subplot_tight(2,3,3,[0.02,0.01]);
imshow(255*(cat(3,and(wc3.P>pval,lo3.posdef),and(wc3.P>pval,lo3.negdef),zeros(nrows,ncols))))
h4 = subplot_tight(2,3,4,[0.02,0.01]);
imshow(255*(cat(3,and(wc4.P>pval,lo4.posdef),and(wc4.P>pval,lo4.negdef),zeros(nrows,ncols))))
h5 = subplot_tight(2,3,5,[0.02,0.01]);
imshow(255*(cat(3,and(wc5.P>pval,lo5.posdef),and(wc5.P>pval,lo5.negdef),zeros(nrows,ncols))))
linkaxes([h1,h2,h3,h4,h5],'xy')

pval = 0.9999; %pval = 0.99995
figure('Name','Wishart change signif, and ~negdef (R), and ~posdef (G), and indef (Y)')
h1 = subplot_tight(2,3,1,[0.02,0.01]);
imshow(255*(cat(3,and(wc1.P>pval,~lo1.negdef),and(wc1.P>pval,~lo1.posdef),zeros(nrows,ncols))))
h2 = subplot_tight(2,3,2,[0.02,0.01]);
imshow(255*(cat(3,and(wc2.P>pval,~lo2.negdef),and(wc2.P>pval,~lo2.posdef),zeros(nrows,ncols))))
h3 = subplot_tight(2,3,3,[0.02,0.01]);
imshow(255*(cat(3,and(wc3.P>pval,~lo3.negdef),and(wc3.P>pval,~lo3.posdef),zeros(nrows,ncols))))
h4 = subplot_tight(2,3,4,[0.02,0.01]);
imshow(255*(cat(3,and(wc4.P>pval,~lo4.negdef),and(wc4.P>pval,~lo4.posdef),zeros(nrows,ncols))))
h5 = subplot_tight(2,3,5,[0.02,0.01]);
imshow(255*(cat(3,and(wc5.P>pval,~lo5.negdef),and(wc5.P>pval,~lo5.posdef),zeros(nrows,ncols))))
linkaxes([h1,h2,h3,h4,h5],'xy')

end % while false
