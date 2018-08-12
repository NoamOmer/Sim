% ------------------------------------------------------
% [nbe] implement an explicit, rathar than MATLAB's, FFT
% -- (capital) X is the input matrix
% -- 2D fft is done in two steps -- one for each dimension
% ------------------------------------------------------
% == 1st Dim ==
X1 = fftshift(X,1);
N  = size(X,1);
k  = (1:N)';
n  = (1:N);
E1 = exp(-1i*2*pi*(k-1)*(n-1)/N);
x1 = (E1 * X1) / sqrt(N);
x1 = fftshift(x1,1);

figure;
subplot(221); imagesc(abs(x));  subplot(223); imagesc(angle(x));
subplot(222); imagesc(abs(x1)); subplot(224); imagesc(angle(x1));


% == 2nd Dim ==
x2 = fftshift(x1,2);
N  = size(X,2);
k  = (1:N)';
n  = (1:N);
E2 = exp(-1i*2*pi*(k-1)*(n-1)/N);
% x2 = (transpose(E2) * x2) / sqrt(N);
x2 = transpose((E2 * transpose(x2)) / sqrt(N));
x2 = fftshift(x2,2);

figure;
subplot(221); imagesc(abs(x));  subplot(223); imagesc(angle(x));
subplot(222); imagesc(abs(x2)); subplot(224); imagesc(angle(x2));

