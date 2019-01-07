function [ w, flag ] = speciesAtElectrolyte ( params, ce, j_Li)

% N: number of transient sections
% h: spatial step ( in m )
% k: time step ( in s )
% ce is column vector
% j_Li is column vector

%% equation
% eps * d(ce)/dt = d(D_eff*d(ce)/dx)/dx + a * (1-t_+) * j
% D_eff*d(ce)/dx = 0 @ x = +0
% D_eff*d(ce)/dx = 0 @ x = -0

%% list of parameters
    
% F = params.F;  % in C/mol

tran_num_plus = params.tran_num_plus;

% De = params.De;

as_a = params.as_a; % in m^2/m^3 
as_c = params.as_c; % in m^2/m^3 

%%
h = params.h; % spatial step 
k = params.k; % time step

m_c = params.m_c; m_s = params.m_s; m_a = params.m_a;
m = params.m; % number of "interior" spatial points 
 
eps = params.eps;
% brug = params.brug;

% De_eff = De * eps .^ brug;

%%
ib1 = params.ib1;
ib2 = params.ib2;
a = params.speciesAtElectrolyte_a;
b = params.speciesAtElectrolyte_b;
S = params.speciesAtElectrolyte_S;

%%

% % NOT considering convection!!
% % ve=zeros(m+2,1);

% a: size (m+2)*(m+2)

% % sharp jump??
% % % a = diag([-3; 1/k*eps(2:(end-1)) + 2/h^2*De_eff(2:(end-1)); -3])...
% % %     +diag([4; -1/h^2*De_eff(3:end)],1)...
% % %     +diag([-1/h^2*De_eff(1:(end-2)); 4],-1);
% a = diag([-3; 1/k*eps(2:(end-1)) + 2/h^2*De_eff(2:(end-1)); -3])...
%     +diag([4; -1/h^2*De_eff(2:(end-1))],1)...
%     +diag([-1/h^2*De_eff(2:(end-1)); 4],-1);
% 
% a(1,3) = -1; 
% a(m+2,m) = -1;
% %
% ib1 = m_c+2;
% a(ib1,ib1-2) = -1*De_eff(ib1-1); a(ib1,ib1-1) =  4*De_eff(ib1-1); 
% a(ib1,ib1) = -3*De_eff(ib1-1) -3*De_eff(ib1+1);
% a(ib1,ib1+1) =  4*De_eff(ib1+1); a(ib1,ib1+2) = -1*De_eff(ib1+1); 
% %
% ib2 = m_c+m_s+3;
% a(ib2,ib2-2) = -1*De_eff(ib2-1); a(ib2,ib2-1) =  4*De_eff(ib2-1); 
% a(ib2,ib2) = -3*De_eff(ib2-1) -3*De_eff(ib2+1);
% a(ib2,ib2+1) =  4*De_eff(ib2+1); a(ib2,ib2+2) = -1*De_eff(ib2+1); 

%%
% b: size (m+2)*1
% b = zeros(m+2,1);
b(2:(end-1)) = 1/k * eps(2:(end-1)) .* ce(2:(end-1));
%
b(ib1) = 0; 
b(ib2) = 0;

% S: size (m+2)*1
% S = zeros(m+2,1);
S(2:(m_c+2)) = as_c * (1-tran_num_plus) * j_Li(2:(m_c+2));
S((m_c+m_s+3):(m+1)) = as_a * (1-tran_num_plus) * j_Li((m_c+m_s+3):(m+1));
%
S(ib1) = 0; 
S(ib2) = 0;

w=a \ (b+S); % do not miss the parenthesis

%%
flag = 0;
% -1
check = find(w >= 0);
if numel(check) ~= numel(w)
    flag = -1;
end


%% 8.1.3 Backward Difference Method, Numerical Analysis, Tim Sauer
% % function w = heatbd( xl,xr,yb,yt,M,N )
% % 
% % f=@(x) sin(2*pi*x).^2;
% % l=@(t) 0*t;
% % r=@(t) 0*t;
% % 
% % D=1;
% % 
% % h=(xr-xl)/M;k=(yt-yb)/N;
% % m=M-1;n=N;
% % 
% % sigma=D*k/(h*h); % must be less than 1/2
% % 
% % a=diag(1+2*sigma*ones(m,1))+diag(-sigma*ones(m-1,1),1)+diag(-sigma*ones(m-1,1),-1);
% % lside=l(yb+(0:n)*k);rside=r(yb+(0:n)*k);
% % w(:,1)=f(xl+(1:m)*h)';
% % 
% % for j=1:n
% %     w(:,j+1)=a\(w(:,j)+sigma*[lside(j);zeros(m-2,1);rside(j)]);
% % end
% % 
% % w=[lside;w;rside];
% % x=(0:m+1)*h;t=(0:n)*k;
% % mesh(x,t,w');
% % view(60,30);axis ([xl xr yb yt -1 1])
% % 
% % 
% % 
% % end

