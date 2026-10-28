#!/usr/bin/env python3
"""Independent PGLib sparse-moment/SOS LOWER-BOUND checker (G0 repaired).
Input = a MATPOWER .m decimal-data file plus lbsparse_*.json certificate.
Uses ONLY standard Python, numpy for shift proposals, mpmath.iv for all proof decisions.
No imports from research-hub. Does not silently accept invalid certificates.

This version accepts any case/JSON with the schema documented by lb_sparse.save_cert.
All certificate contributions loc, thermal, eqs use their POSITIVE hex scale s.
It has NOT been replayed on original private repository cases in this environment.
"""
import hashlib, json, math, re, sys
from decimal import Decimal
import numpy as np
from mpmath import iv
iv.dps=70

def I(v):return iv.mpf(str(v))
def H(s):
    """Enclose the exact binary64 value, not its shortest decimal repr."""
    x=float.fromhex(s)
    if not math.isfinite(x):raise ValueError('nonfinite hex float')
    n,d=x.as_integer_ratio()
    return iv.mpf(n)/iv.mpf(d)

def F64(x):
    """Exact dyadic embedding of a parsed JSON binary64 into interval numbers."""
    v=float(x)
    if not math.isfinite(v):raise ValueError('NONFINITE BINARY64')
    n,d=v.as_integer_ratio()
    return iv.mpf(n)/iv.mpf(d)
def pos_scale(entry):
    """CERTIFICATE SCALE: loc, thermal and eqs terms are multiplied by s."""
    q=H(entry.get("s","0x1.0p+0"))
    if not q.a>0: raise ValueError("nonpositive or ambiguous certificate scale")
    return q
def thin_upper(x):return x.b
def Babs(x):return abs(x).b

def parse(b):
    s=b.decode('utf-8'); m=re.search(r'mpc\.baseMVA\s*=\s*([^;]+);',s)
    if m is None:raise ValueError('No baseMVA')
    def tab(k):
        m=re.search(r'mpc\.'+k+r'\s*=\s*\[([\s\S]*?)\]\s*;',s)
        if m is None:raise ValueError('Missing '+k)
        return [x.split('%')[0].replace(';','').split() for x in m.group(1).splitlines() if x.split('%')[0].replace(';','').strip()]
    return I(m.group(1).strip()),tab('bus'),tab('gen'),tab('branch'),tab('gencost')

Z=()
def mulmon(a,b):
    d=dict(a)
    for k,p in b:d[k]=d.get(k,0)+p
    return tuple(sorted((k,p) for k,p in d.items() if p))

def P(c=0):return {Z:c if hasattr(c,'ctx') else I(c)}
def add(x,y,sign=1):
    z=dict(x)
    for k,v in y.items():z[k]=z.get(k,I(0))+sign*v
    return z

def scale(x,c):return {k:v*c for k,v in x.items()}
def mul(x,y):
    out={}
    for a,p in x.items():
        for b,q in y.items():
            key=mulmon(a,b);out[key]=out.get(key,I(0))+p*q
    return out

def square_poly(p):
    # Algebraic polynomial square, unlike naive interval product x*x of an interval value.
    return mul(p,p)

def var(k):return {((k,1),):I(1)}

def mincost(a2,b,c,lo,hi):
    f=lambda x:(a2*x+b)*x+c
    x,y=f(lo),f(hi)
    if a2==I(0):return x.a if x.a<y.a else y.a
    if not a2.a>0:raise ValueError('nonconvex generator cost')
    # Unconstrained minimum always gives a valid lower bound, even if outside box.
    # Taking max(unconstrained, safe endpoint) would be invalid without proof.
    v=-b/(2*a2)
    if v.b<lo.a:return x.a
    if v.a>hi.b:return y.a
    return (c-b*b/(4*a2)).a


def strict_ldl(X, shift):
    """Sufficient exact interval LDL pivot positivity test."""
    n=len(X); L=[[I(0) for _ in range(n)] for _ in range(n)];D=[I(0)]*n
    for j in range(n):
        v=X[j][j]+shift
        for k in range(j):v-=L[j][k]**2*D[k]
        if not (v.a>0):return False
        D[j]=v
        for i in range(j+1,n):
            z=X[i][j]
            for k in range(j):z-=L[i][k]*L[j][k]*D[k]
            L[i][j]=z/D[j]
    return True


def psd_shift_and_bound(entries,basis,vertices):
    N=len(basis);X=[[I(0) for _ in range(N)] for _ in range(N)]
    for i,j,h in entries:
        assert 0<=i<=j<N
        X[i][j]=X[j][i]=H(h)
    mid=np.array([[float(x.mid) for x in row] for row in X],dtype=float)
    ev=float(np.linalg.eigvalsh((mid+mid.T)/2)[0]); sc=max(1.,float(np.max(np.abs(mid))))
    lo=0.; hi=max(sc*1e-16, -ev + sc*1e-14)
    if strict_ldl(X,I(0)):hi=0.
    else:
        for _ in range(64):
            if strict_ldl(X,F64(float(hi))):break
            hi*=2
        else:raise ArithmeticError('interval PSD shift cannot be bracketed')
        # Minimize positive shift while retaining strictly positive LDL pivots.
        for _ in range(26):
            m=lo+(hi-lo)/2
            if strict_ldl(X,F64(float(m))):hi=m
            else:lo=m
    delta=F64(float(hi))
    assert strict_ldl(X,delta)
    pol={}
    for i in range(N):
        for j in range(N):
            key=mulmon(basis[i],basis[j]);pol[key]=pol.get(key,I(0))+X[i][j]
    norm=sum((bmon(m,vertices)**2 for m in basis),I(0))
    return pol,delta,norm

def bmon(m,verts):
    result=I(1)
    for k,p in m:result*=verts[k]**p
    return result.b



def reject_nonfinite_certificate_values(obj, path="$"):
    """Fail closed on NaN/Inf even in currently unused certificate fields."""
    import math, re
    if isinstance(obj, dict):
        for k,v in obj.items():reject_nonfinite_certificate_values(v, path+"."+str(k))
    elif isinstance(obj, list):
        for i,v in enumerate(obj):reject_nonfinite_certificate_values(v, path+"["+str(i)+"]")
    elif isinstance(obj, float):
        if not math.isfinite(obj):raise ValueError("NONFINITE CERTIFICATE FIELD "+path)
    elif isinstance(obj, str):
        t=obj.strip().lower()
        if t in ("nan","+nan","-nan","inf","+inf","-inf","infinity","+infinity","-infinity"):
            raise ValueError("NONFINITE CERTIFICATE HEX/STRING "+path)
        if re.match(r"^[+-]?0x",t) and not math.isfinite(float.fromhex(t)):
            raise ValueError("NONFINITE HEX FLOAT "+path)



def main(mpath,jpath):
    raw=open(mpath,'rb').read();cert=json.load(open(jpath));reject_nonfinite_certificate_values(cert);assert hashlib.sha256(raw).hexdigest()==cert['sha256'],'SHA mismatch'
    B,bus,gen,br,gc=parse(raw);n=len(bus);ref=next(i for i,r in enumerate(bus) if int(r[1])==3)
    ids={int(r[0]):i for i,r in enumerate(bus)}
    gb=[ids[int(r[0])] for r in gen if Decimal(r[7])>0]
    gen2=[(r,c) for r,c in zip(gen,gc) if Decimal(r[7])>0]
    ng=len(gen2)
    variant=cert.get('variant','pglib')
    assert variant in ('pglib','current')
    assert len(set(ids))==n and len([r for r in bus if int(r[1])==3])==1
    assert len(gen2)==len(gb)
    active=[row for row in br if Decimal(row[10])>0]
    assert all(int(row[0]) in ids and int(row[1]) in ids for row in active)
    assert all(int(r[0]) in ids for r,_ in gen2)
    idx={('e',i):i for i in range(n)};idx.update({('f',i):n+i-(i>ref) for i in range(n) if i!=ref})
    e=[var(idx[('e',i)]) for i in range(n)];f=[var(idx[('f',i)]) if i!=ref else P(0) for i in range(n)]
    def w(i,j):return add(mul(e[i],e[j]),mul(f[i],f[j])),add(mul(f[i],e[j]),mul(e[i],f[j]),-1)
    diag=[w(i,i)[0] for i in range(n)]
    Pr=[add(P(I(r[2])/B),scale(diag[i],I(r[4])/B)) for i,r in enumerate(bus)]
    Qr=[add(P(I(r[3])/B),scale(diag[i],-I(r[5])/B)) for i,r in enumerate(bus)]
    desc={};put=lambda k,v:desc.__setitem__(tuple(k),v)
    for k,z in enumerate(active):
        i,j=ids[int(z[0])],ids[int(z[1])];r,x,bc=I(z[2]),I(z[3]),I(z[4]);G=r/(r*r+x*x);b=-x/(r*r+x*x)
        tau=I(z[8]) if Decimal(z[8]) else I(1);phi=I(z[9])*iv.pi/180
        cs,sn=iv.cos(phi),iv.sin(phi);R,Im=w(i,j)
        ar=(G*cs-b*sn)/tau;ai=-(G*sn+b*cs)/tau
        brre=(G*cs+b*sn)/tau;brim=(G*sn-b*cs)/tau
        Pf=add(scale(diag[i],G/(tau*tau)),add(scale(R,ar),scale(Im,ai),-1),-1)
        Qf=add(scale(diag[i],-(b+bc/2)/(tau*tau)),add(scale(Im,ar),scale(R,ai)),-1)
        Pt=add(scale(diag[j],G),add(scale(R,brre),scale(Im,brim)),-1)
        Qt=add(scale(diag[j],-(b+bc/2)),add(scale(R,brim),scale(Im,brre),-1),-1)
        Pr[i]=add(Pr[i],Pf);Pr[j]=add(Pr[j],Pt);Qr[i]=add(Qr[i],Qf);Qr[j]=add(Qr[j],Qt)
        rr=I(z[5])/B
        if variant=='current':
            if Decimal(z[5])>0:
                # I_f = (Y+jbc/2)/tau² * V_i - Y/conj(T) * V_j
                # I_t = -Y/T * V_i + (Y+jbc/2) * V_j
                def cmul(a,b,k):
                    ar,ai=a;R_=add(scale(e[k],ar),scale(f[k],ai),-1)
                    I_=add(scale(f[k],ar),scale(e[k],ai))
                    return R_,I_
                def cadd(x,y):return add(x[0],y[0]),add(x[1],y[1])
                aa=(G/(tau*tau),(b+bc/2)/(tau*tau))
                beta=(-(G*cs-b*sn)/tau,-(G*sn+b*cs)/tau)
                gamma=(-(G*cs+b*sn)/tau,-(b*cs-G*sn)/tau)
                delta=(G,b+bc/2)
                for side,(A_,B_) in [('from',(cadd(cmul(aa,None,i),cmul(beta,None,j)),None)),
                                      ('to',(cadd(cmul(gamma,None,i),cmul(delta,None,j)),None))]:
                    ir,ii=A_;put(('I',k,side),('g',add(P(rr*rr),add(mul(ir,ir),mul(ii,ii)),-1)))
        else:
            if Decimal(z[5])>0:
                put(('S',k,'from'),('S',(rr,Pf,Qf)));put(('S',k,'to'),('S',(rr,Pt,Qt)))
            assert abs(Decimal(z[11]))<90 and abs(Decimal(z[12]))<90
            put(('angmax',k),('g',add(scale(R,iv.tan(I(z[12])*iv.pi/180)),Im,-1)))
            put(('angmin',k),('g',add(Im,scale(R,iv.tan(I(z[11])*iv.pi/180)),-1)))
    put(('eref',),('g',add(e[ref],P(-I(bus[ref][12])))))
    for i,r in enumerate(bus):
        put(('vmax',i),('g',add(P(I(r[11])**2),diag[i],-1)))
        put(('vmin',i),('g',add(diag[i],P(-I(r[12])**2))))
        ks=[k for k,(g,c) in enumerate(gen2) if ids[int(g[0])]==i]
        if ks:
            summer=lambda col:sum((I(gen2[k][0][col])/B for k in ks),I(0))
            put(('pmax',i),('g',add(P(summer(8)),Pr[i],-1)))
            put(('pmin',i),('g',add(Pr[i],P(-summer(9)))))
            put(('qmax',i),('g',add(P(summer(3)),Qr[i],-1)))
            put(('qmin',i),('g',add(Qr[i],P(-summer(4)))))
        else:
            put(('hP',i),('h',Pr[i]));put(('hQ',i),('h',Qr[i]))
    varsbound={v:I(bus[i][11]) for (kind,i),v in idx.items()}
    def bs(m):return bmon(m,varsbound)
    def basis_entry(e):return tuple(sorted((idx[(k,int(i))],int(d)) for k,i,d in e))
    lhs=P(H(cert['t']));pen=I(0)
    proofcounts={'gram':0,'loc':0,'thermal':0,'eqs':0,'pieces':0}
    for gr in cert['gram']:
        basis=[basis_entry(v) for v in gr['basis']]
        pol,d,norm=psd_shift_and_bound(gr['X'],basis,varsbound);pen+=d*norm;lhs=add(lhs,pol);proofcounts['gram']+=1
    for loc in cert['loc']:
        kind,cons=desc[tuple(cert['desc'][loc['j']])];assert kind=='g'
        cons=scale(cons,pos_scale(loc))
        if 'scalar' in loc:
            c=H(loc['scalar']);lhs=add(lhs,scale(cons,c if c.a>=0 else I(0)))
        else:
            basis=[basis_entry(v) for v in loc['basis']]
            pol,d,norm=psd_shift_and_bound(loc['X'],basis,varsbound)
            lhs=add(lhs,mul(pol,cons));gmax=sum((Babs(v)*bs(k) for k,v in cons.items()),I(0));pen+=d*norm*gmax
        proofcounts['loc']+=1
    for th in cert['thermal']:
        kind,(r,Pf,Qf)=desc[tuple(cert['desc'][th['j']])];assert kind=='S'
        sf=pos_scale(th)
        if th['kind']=='quartic':
            sc=H(th['v'][0]);sc=(sc if sc.a>=0 else I(0));sc*=sf;lhs=add(lhs,scale(add(P(r*r),add(square_poly(Pf),square_poly(Qf)),-1),sc))
        else:
            ss,aa,bb=[H(v) for v in th['v']];sig=iv.sqrt(aa**2+bb**2);sigma=ss if ss.a>=sig.b else sig.b
            lhs=add(lhs,scale(add(P(sigma*r),add(scale(Pf,aa),scale(Qf,bb))),sf))
        proofcounts['thermal']+=1
    for eq in cert['eqs']:
        kind,h=desc[tuple(cert['desc'][eq['j']])];assert kind=='h'
        poly={basis_entry(b):H(t) for b,t in zip(eq['basis'],eq['p'])}
        lhs=add(lhs,mul(poly,scale(h,pos_scale(eq))));proofcounts['eqs']+=1
    objective=P(0);covered=[]
    for i in sorted(set(gb)):
        ks=[k for k,b in enumerate(gb) if b==i]
        if len(ks)==1:
            k=ks[0];row=gen2[k][1];mw=scale(Pr[i],B)
            objective=add(objective,add(scale(mul(mw,mw),I(row[4])),add(scale(mw,I(row[5])),P(I(row[6])))))
            covered.append(i)
    for piece in cert['pieces']:
        i=piece['bus'];ks=piece['gens'];assert len(ks)>1 and i not in covered
        rho=[(x if x.a>=0 else I(0)) for x in (H(r) for r in piece['rho'])]
        total=sum(rho,I(0));assert total.a>0
        for lam,rr in zip(piece['lam'],rho):
            la=H(lam);psi=I(0)
            for k in ks:
                g,c=gen2[k];a2=I(c[4])*B*B; bb=(I(c[5])-la)*B
                psi+=mincost(a2,bb,I(c[6]),I(g[9])/B,I(g[8])/B)
            objective=add(objective,scale(add(scale(Pr[i],B*la),P(psi)),rr/total))
        covered.append(i);proofcounts['pieces']+=1
    assert sorted(covered)==sorted(set(gb)), 'not all generator buses covered'
    residual=add(objective,lhs,-1)
    rbound=sum((Babs(v)*bs(m) for m,v in residual.items()),I(0))
    L=H(cert['t'])-pen-rbound
    Lout=math.nextafter(float(L.a),-math.inf)
    print('source_sha256',hashlib.sha256(raw).hexdigest())
    print('case',cert['case'],'variant',variant,'n',n,'active_branches',len(active))
    print('certificate_scales', {k:[str(pos_scale(t)) for t in cert[k][:3]] for k in ('loc','thermal','eqs')})
    print('cert_parts',proofcounts)
    print('PSD_SHIFT_AND_RESIDUAL_PROVED',True)
    print('PSD_penalty_upper',str(pen.b),'residual_upper',str(rbound.b))
    print('L_strict_unrounded',str(L),'L_outward',repr(Lout),'cert_L',repr(cert['L']))
    if not F64(cert['L']).b<=L.a:raise AssertionError('declared L exceeds strictly proved bound')
    print('ACCEPT')

if __name__=='__main__':
    if len(sys.argv)!=3:raise SystemExit('usage: independent_sparse.py case.m lbsparse.json')
    main(*sys.argv[1:])
