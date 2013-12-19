function [X,hist,total_time] = smacof_sparse(Adj,X0,iter,verbose,xhistory,rtol,atol, tol, wp, wn, wu, dp, dn, du)

% check input correctness
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

if nargout == 3,
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
iii = 1;  %iteration counter
X = X0;   %initial guess of X
[n,d]=size(X);

Lw=calc_Lw (Adj, wp, wn, wu);  %calculate L^w before the recursive step

[it, jt, W, D]=calc_i(Adj, dp,dn,du, wp, wn, wu); %find out the non-zeros terms in the adjacency matrix

[norm_X, norms]=calc_normX(it,jt, X);  %compute Xi-Xj/|Xi-Xj| with non-zero weights

fprintf('Fatorization begins!\n');

%Approximate L^w by L^w+ delta*I, where delta is small. See
%Leevnberg-Marquardt Method for reference.
delta=tol*speye(size(Lw,1));
Lw=Lw+delta;

%1. Use Tim Davis's object oriented factorize package to factorize sparse, symmetric and positive definite Lw
%F=factorize(Lw,'chol',1) ;

%2. construct function handler for matrix Lw.
Ldiag=diag(Lw);
afun=@(x)lwx(it, jt, W, Ldiag, x);

fprintf('initialized!\n');

% initialize history
if HISTORY || VERBOSE,
    hist.s(1) = calc_stress (W,D,norms);
end

if XHISTORY
   hist.X{1} = X0; 
end

if HISTORY && VERBOSE,
    fprintf(1,'iter         stress   time (sec)\n') 
    fprintf(1,'INIT   %12.3g   ----------\n', hist.s(1)); 
end

clock=cputime;
while (iii <= iter),   
    t = cputime;       
    %compute Lx^X based on updated X
    B=calc_B (norm_X, W, D, it, n, d);
    
    %1. solve X by cholesky factorization
    %X=F\B;
    
    %2. solve X by Preconditioned Conjugate Gradients Method with function
    %handler input @afun
    
    for j=1:d
            X(:,j)=pcg(afun,B(:,j),0.1,100);
    end
    
    
    %compute Xi-Xj/|Xi-Xj| with non-zero weights based on updated X
    [norm_X, norms]=calc_normX(it,jt, X);
    
    %update current stress function
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
        fprintf(1,'atol=%g reached, exiting\n',atol) %absolute stopping condition
        return 
    end

    if (iii > 1) && (HISTORY || VERBOSE),
        if (hist.s(iii-1)/hist.s(iii)-1) < rtol, %relative stopping condition
            total_time=cputime-clock;
            fprintf(1,'rtol=%g reached, exiting\n',rtol)
            fprintf(1,'total running time(sec): \n')
            fprintf(1,' %10.4g', total_time) 
            return 
        end
    end

    iii = iii+1;
        
end    




%SERVICE FUNCTIONS
%------------------------------------------------------------------------------
%compute L^w matrix
function [Lw] = calc_Lw (Adj, wp, wn, wu)
Lw = sparse(Adj);

i = find(Adj(:) == 3);
Lw(i) = - wp;
i = find(Adj(:) == 2);
Lw(i) = - wu;
i = find(Adj(:) == 1);
Lw(i) = - wn;

d = sum(Lw);
Lw = Lw - diag(d);
return

%compute the non-zero terms from sparse adjacency matrix Adj
%W,D are the weight/distance values of these non-zero terms
function [it, jt, W, D]=calc_i(Adj, dp,dn,du, wp, wn, wu)
[it,jt,s]=find(Adj);
W=zeros(size(s));
D=zeros(size(s));
for k=1: size(s,1),
    if (s(k)==3),
	W(k)=wp;
	D(k)=dp;
    elseif (s(k)==2),
	W(k)=wu;
	D(k)=du;
    elseif (s(k)==1), 	
	W(k)=wn;
	D(k)=dn;
    else
	fprintf('type error!\n');
    end
end	

return

%compute the norm Xi-Xj/|Xi-Xj| for the non-zero elements
function [norm_X, norms]=calc_normX(it,jt, X)
Xi=X(it,:);
Xj=X(jt,:);
deltaX=Xi-Xj;
norms= sqrt(sum((deltaX).^2,2));
if (size(X,2)==2),
    norm_X=deltaX./[norms,norms];
elseif (size(X,2)==4),
    norm_X=deltaX./[norms,norms,norms,norms];
elseif (size(X,2)==3),
    norm_X=deltaX./[norms,norms,norms];
else
    fprintf('dim=?\n');
end

% compute the current stress
function [S] = calc_stress (W,D, norms)
S=sum(W.*((norms-D).^2));

return

%compute b=LxX by the summation of all non-zero wijdij Xi-Xj/|Xi-Xj|
function [B] = calc_B (norm_X, W, D, it, n, d)
B = zeros(d,n);
wd=W.*D;
if (d==2),
    WD=[wd,wd];
elseif (d==4),
    WD=[wd,wd,wd,wd];
elseif (d==3),
    WD=[wd,wd,wd];
else
    fprintf('dim=?\n');
end

WD=(WD.*norm_X)';
for k=1:size(it,1),
    B(:,it(k))=B(:,it(k))+WD(:,k);
end

B=B';
return

%define the function handler of L^wX
%In particular, we do LwX=L'X+ diag(Lw)X, where diag(Lw) is
%the diagonal vector of Lw.
function y=lwx(it, jt, W, Ldiag, X)

y=zeros(size(X));
for k=1:size(W,1),
    y(it(k))=y(it(k))-W(k)*X(jt(k));
end
y=y+Ldiag.*X;
    


