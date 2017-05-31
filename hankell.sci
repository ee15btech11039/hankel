function [nk,W]=hankell(g,tol)
    
    [lhs,rhs]=argn(0),
    if rhs <1 then
        error(msprintf(gettext("%s: Wrong number of input arguments: At least %d expected.\n"),"hankel",1))
    end,
    if (typeof(g)=='rational') then
        if(degree(g.num)>degree(g.den))
            error(msprintf(gettext("The %s command cannot be used for models with more zeroes than poles","hankel")))
        end;
        
    end,

    if and(typeof(g)<>['rational','state-space']) then
        error(msprintf(gettext("%s: Wrong type for input argument #%d: Linear state space or a transfer function expected.\n"),"hankel",1))
    end
    if rhs==1 then 
    tol=100*%eps,
    else
        if type(tol)<>1|size(tol,'*')<>1|~isreal(tol)|tol<=0 then
      error(msprintf(gettext("%s: Wrong value for input argument #%d: Must be a positive scalar.\n"),"hankel",2))
        end
    end,
    
    if typeof(g)=='rational' then
        g=tf2ss(g); 
    end;
    lf=spec(g(2)),
    if min(abs(lf))<=tol then
        error(msprintf(gettext("%s: Wrong value for input argument #%d: Pure imaginary poles unexpected.\n"),"hankel",1)),
    end
    
    if g.dt==[] then
        warning(msprintf(gettext("%s: Input argument %d is assumed continuous time.\n"),"hankel",1));
        g.dt='c'
    end;

    if(g.dt=='c') then
        [sla,sls,d]=dtsi(g);
        lc=ctr_gram(sls),lo=obs_gram(sls),W=lc*lo;
        nk=gsort(real(spec(W)));
        return;
    end
    
    [a,b,c,d]=abcd(g),gi=d,
    [n1,n2,t]=size(g),
    [a,u0]=balanc(a);b=u0\b;c=c*u0;
    [u,n]=schur(a,'d'),
    a=u'*a*u,
    flag=0;
    if n==t then
        ga=0,
        flag='stable';
        gs=g,
    end,
    if n==0 then
        gs=0,
        ga=g,
    end,
    //if k==0 then 
      //  flag='stable';
    //else
      //  flag='unstable';
    //end,
        
    if(g.dt~='c') then
        if flag=='stable' then
            lc=ctr_gram(g),lo=obs_gram(g),W=lc*lo;
            nk=gsort(real(spec(W)));
            return;
        
        else
    //      [ab,w,bs]=bdiag(a);
            a1=a(1:n,1:n),a4=a(n+1:t,n+1:t),x=a(1:n,n+1:t),
            z=sylv(a1,-a4,-x,'c'),
            w=[eye(n,n),z;0*ones(t-n,n),eye(t-n,t-n)],
            wi=[eye(n,n),-z;0*ones(t-n,n),eye(t-n,t-n)],
            tr=u*w,tri=wi*u';
            bb=tri*b,b1=bb(1:n,:),b2=bb(n+1:t,:),
            cc=c*tr,c1=cc(:,1:n),c2=cc(:,n+1:t),
            ga=syslin('d',a4,b2,c2),
            gs=syslin('d',a1,b1,c1);
            lc=ctr_gram(gs),lo=obs_gram(gs),W=lc*lo;
            nk=gsort(real(spec(W)));
            return;
        end
        
    end;
endfunction

