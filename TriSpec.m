% Example code which computes the spectrum of infinite periodic
% tri-diagonal matrices. Input is given in the form of three vectors
% representing one period for each diagonal.

clear;
M        = 100; 
theta    = linspace(0,2*pi,M);
L        = build_lattice(-2-3i, 4+3i, M);
npts     = length(L(:));
Spectrum = zeros(size(L));

% Lower, middle and upper diagonals:
Aminus = [2i; -3i; 2i; 0; 1i];
Adiag  = [1; 0; 1; 0; 2];
Aplus  = [-1; -2; 1; 3i; -5];

N      = length(Adiag); % Period

ctr = M/30;
disp('______________________________')
for k=1:length(theta)
    Aminus_theta = circshift(Aminus,1)*exp(-1i*theta(k)/N);
    Aplus_theta = circshift(Aplus,-1)*exp(1i*theta(k)/N);
    A = spdiags([Aminus_theta, Aplus_theta, Adiag, Aminus_theta, Aplus_theta], [-N+1, -1, 0, 1, N-1], N, N);
    S = zeros(size(L));
    p = charpoly(A).';
    parfor j=1:npts
        z = L(j);
        v = flip(z.^(0:N));
        S(j) = abs(v*p) < 10/sqrt(M);
    end
    Spectrum = Spectrum | S;
    if k >= ctr
        fprintf('=')
        ctr = ctr+M/30;
    end
end
disp(' ')

% Plot result:
figure
plot(L(Spectrum),'.','MarkerEdgeColor',[0.3 0.3 1], 'MarkerSize',5);
xlim([min(real(L(:))) max(real(L(:)))])
ylim([min(imag(L(:))) max(imag(L(:)))]);
title('Spectrum:');
