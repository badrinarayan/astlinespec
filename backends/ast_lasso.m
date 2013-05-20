function [x_debias,fs_lasso,c_debias,x,c,x_ls_debias,c_ls_debias]=ast_lasso(y, tau,nsmps)
n = length(y);

[c,c_debias] = SpaRSA(y,@(z) moment_pad(z,nsmps,n),tau/sqrt(nsmps),...
    'AT',@(z) fft(z,nsmps)/sqrt(nsmps),...
    'Psi',@(z,tau) max(abs(z)-tau,0).*exp(1i*angle(z)),... %shrinkage function
    'Phi',@(z) sum(abs(z)),... %penalty function (l1-norm)
    'MaxiterA',1000,...
    'ToleranceA',1e-3,...
    'Initialization',2,...
    'Verbose',0,...
    'Debias',1);

x = moment_pad(c,nsmps,n);
x_debias = moment_pad(c_debias,nsmps,n);

c_copy = c;
pole_regions = find(c_copy);
pole_est = [];

if isempty(pole_regions)
    c_ls_debias = zeros(size(c));
    x_ls_debias = zeros(size(x));
    return;
end
% find connected segments where the polynomial is near 1
break_points = [1;find(diff(pole_regions)>1)+1;length(pole_regions)];

% in each connected segment, find the point closest to 1 and call this a
% pole
for k=1:length(break_points)-1,
    lookup = pole_regions(break_points(k):(break_points(k+1)-1));
    %pole_est(k) = 2*pi*sum(lookup.*abs(polyval(lookup)))/grid_size/sum(abs(polyval(lookup)));

    if isempty(lookup)
        lookup = pole_regions(break_points(k));
    end

    [~,slot] = max(abs(c_copy(lookup)));
    pole_est(k) = lookup(slot)/nsmps;
end

fs_lasso = pole_est;

est_basis = exp( 1i*(0:(n-1))'*2*pi*fs_lasso(:)' )/sqrt(n);
c_ls_debias = est_basis\y;
x_ls_debias = est_basis*c_ls_debias;

% [U,S,V] = svd(est_basis);
% d = diag(S);
% I_th = find(diag(S)<1,1,'first');
% U = U(:,1:(I_th-1));
% c_ls_debias = U'*y;
% x_ls_debias = U*c_ls_debias;

end
