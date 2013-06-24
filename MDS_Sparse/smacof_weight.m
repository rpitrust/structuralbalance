function [X,hist] = smacof_weight(Adj,X0,iter,verbose,xhistory,rtol,atol, tol, wp, wn, wu, dp, dn, du, wwp, wwn, dwp, dwn)

% check input correctness
if nargin < 13,
    error('Incorrect number of arguments, exiting.')     
end


if size(Adj,1) ~= size(Adj,2),
    error('Matrix D must be square, exiting.')     
end

if size(Adj,1) ~= size(X0,1),
    error('X0 and D dimensions mismatch, exiting. X0 must be a size(D,1)*dim matrix.')     
end

if rtol < 0
    error('rtol must be non-negative, exiting.')     
end


% check input and output flags
if strcmp(lower(verbose),'iter'),
    VERBOSE = 1;
else
    VERBOSE = 0;
end

if nargout == 2,
    HISTORY = 1;
else
    HISTORY = 0;
end

if strcmp(lower(xhistory),'on'),
    XHISTORY = 1;
else
    XHISTORY = 0;
end



% initialize
iii = 1;

X = X0;
[n,d]=size(X);

Lw=calc_Lw (Adj, wp, wn, wu, wwp, wwn);

%dd=diag(Lw);
%find(dd(:)==0)


[i, W, D]=calc_i(Adj, dp,dn, du, wp, wn, wu, wwp, wwn, dwp, dwn);
ij=calc_ij(i,n);

[norm_X, norms]=calc_normX(ij, X);
fprintf('Fatorization begins!\n');

%full(Adj)
%full(Lw)
%make Lw strict diagonally dominant matrix
delta=tol*speye(size(Lw,1));
Lw=Lw+delta;

%Lw=Lw(2:n,2:n);
%full(Lw)
%[R,p] = chol(Lw);

F=factorize(Lw,'chol',1) ;
%S=inverse(Lw); 
%[U,S,V] = svds(Lw,kk);
%Lwp=V*pinv(S)*U';

fprintf('initialized!\n');

% initialize history
if HISTORY || VERBOSE,
    %hist.time = zeros(iter);
    %hist.s = zeros(iter);
    hist.s(1) = calc_stress (W,D,norms);
end

if XHISTORY
   hist.X{1} = X0; 
end

if HISTORY && VERBOSE,
    fprintf(1,'iter         stress   time (sec)\n') 
    fprintf(1,'INIT   %12.3g   ----------\n', hist.s(1)); 
end


while (iii <= iter),
    t = cputime;       
    
    B=calc_B (norm_X, W, D, ij, n, d);
    %B=B(2:n,:);
    %fprintf('B solved!\n', cputime);
    
    %rX = R\(R'\B);
    X=F\B;

    %rX=F\B;
    %X(2:n,:)=rX;
    %X=F\B;
    
    %fprintf('X solved!\n', cputime);

    [norm_X, norms]=calc_normX(ij, X);
    S = calc_stress (W,D,norms);
    

    % add history
    if HISTORY || VERBOSE,    
        hist.time(iii) = cputime-t;
        hist.s(iii) = S;
    end
    
    if XHISTORY
        hist.X{iii} = X; 
    end

    if HISTORY && VERBOSE,
        fprintf(1,'%4d   %12.3g   %10.3g\n', iii,hist.s(iii),hist.time(iii)) 
    end
    
    
    % check stopping conditions
    if S < atol,
        fprintf(1,'atol=%g reached, exiting\n',atol)
        return 
    end

    if (iii > 1) && (HISTORY || VERBOSE),
        if (hist.s(iii-1)/hist.s(iii)-1) < rtol,
            fprintf(1,'rtol=%g reached, exiting\n',rtol)
            return 
        end
    end

    iii = iii+1;
        
end    

% SERVICE FUNCTIONS
function [Lw] = calc_Lw (Adj, wp, wn, wu, wwp ,wwn)

Lw = sparse(Adj);


i = find(Adj(:) == 5);
Lw(i) = - wp;
i = find(Adj(:) == 4);
Lw(i) = - wwp;
i = find(Adj(:) == 3);
Lw(i) = - wu;
i = find(Adj(:) == 2);
Lw(i) = - wwn;
i = find(Adj(:) == 1);
Lw(i) = - wn;

d = sum(Lw);
Lw = Lw - diag(d);
return

%compute non-zero single-digit index from sparse Adj
function [i, W, D]=calc_i(Adj, dp,dn, du, wp, wn, wu, wwp, wwn, dwp, dwn)
i=find(Adj(:)~=0);
W=zeros(size(i,1),1);
D=zeros(size(i,1),1);
for k=1: size(i,1),
    if (Adj(i(k))==5),
	W(k)=wp;
	D(k)=dp;
    elseif (Adj(i(k))==4),
	W(k)=wwp;
	D(k)=dwp;
    elseif (Adj(i(k))==3),
	W(k)=wu;
	D(k)=du;
    elseif (Adj(i(k))==2), 	
	W(k)=wwn;
	D(k)=dwn;
    elseif (Adj(i(k))==1), 	
	W(k)=wn;
	D(k)=dn;
    else,
	fprintf('type error!\n');
    end
end	
return

%compute non-zero i,j index from single-digit index
function [ij]=calc_ij(i, n)
ij=zeros(size(i,1),2);
ij(:,1)=floor((i-0.01)/n)+1;
ij(:,2)=i-(ij(:,1)-1)*n;
return

%compute norms and norm_X list based on X and ij
function [norm_X, norms]=calc_normX(ij, X)
Xi=X(ij(:,1),:);
Xj=X(ij(:,2),:);
deltaX=Xi-Xj;
norms= sqrt(sum((deltaX).^2,2));
if (size(X,2)==2),
    norm_X=deltaX./[norms,norms];
elseif (size(X,2)==4),
    norm_X=deltaX./[norms,norms,norms,norms];
elseif (size(X,2)==3),
    norm_X=deltaX./[norms,norms,norms];
else,
    fprintf('dim=?\n');
end

% compute the stress by sparse Adjj
function [S] = calc_stress (W,D, norms)
S=sum(W.*((norms-D).^2));

return

%compute LxX by sparse index
function [B] = calc_B (norm_X, W, D, ij, n, d)
B = zeros(n,d);
wd=W.*D;
if (d==2),
    WD=[wd,wd];
elseif (d==4),
    WD=[wd,wd,wd,wd];
elseif (d==3),
    nWD=[wd,wd,wd];
else,
    fprintf('dim=?\n');
end
WD=WD.*norm_X;
for k=1:size(ij,1),
    B(ij(k,1),:)=B(ij(k,1),:)+WD(k,:);
end
   
return

