#!/usr/bin/env python3
"""Standalone generic PGLib upper-bound Krawczyk checker; no repo imports.

Usage: python independent_ub_krawczyk.py cases/pglib_opf_case118_ieee.m ubcert_case118_ieee.json
Requires numpy, mpmath. Checks SHA of complete raw .m input, computes its own
interval-Jacobian Krawczyk enclosure, checks all original inequalities, and
reports an outward-rounded binary64 cost U if all checks pass. Fail closed.
Usage: python independent_ub_krawczyk_r12.py case.m ubcert.json
Exit 0 = ACCEPT, 2 = failed Krawczyk, 3 = operational limits,
1 = other REJECT/ERROR. Automatically uses scipy.sparse inverse rows for N>100.
N=4000 3GB/2h performance target is UNVERIFIED: test on Claude host before relying on it.
No generator or project implementation is imported.
Exit policy: 0 ACCEPT; 2 REJECT_MATH; 4 UNSUPPORTED_FORMAT; 5 INPUT_MISMATCH; 6 CHECKER_ERROR.
Decimal-exact fixed bounds are treated as identities, not ambiguous interval inequalities.
"""
import hashlib, json, math, re, sys
from decimal import Decimal
import numpy as np
from mpmath import iv
iv.dps = 75

def I(s): return iv.mpf(str(s))
def H(s):
    """Exact binary64 represented as dyadic rational, then bounded interval."""
    x = float.fromhex(s)
    if not math.isfinite(x): raise ValueError('nonfinite hex float')
    n, d = x.as_integer_ratio()
    return iv.mpf(n) / iv.mpf(d)

def F64(x):
    """Exact dyadic embedding of a parsed JSON binary64 into interval numbers."""
    v=float(x)
    if not math.isfinite(v):raise ValueError('NONFINITE BINARY64')
    n,d=v.as_integer_ratio()
    return iv.mpf(n)/iv.mpf(d)
def lower(x): return float(x.a)
def upper(x): return float(x.b)
def lo(x, z): return x.a > I(z).b
def hi(x, z): return x.b < I(z).a

def parse(raw):
    s = raw.decode('utf-8'); base = re.search(r'mpc\.baseMVA\s*=\s*([^;]+);',s).group(1).strip()
    def tab(k):
        m=re.search(r'mpc\.'+k+r'\s*=\s*\[([\s\S]*?)\]\s*;',s)
        if m is None: raise ValueError('missing '+k)
        return [row.split('%')[0].replace(';','').split() for row in m.group(1).splitlines()
                if row.split('%')[0].replace(';','').strip()]
    return I(base),tab('bus'),tab('gen'),tab('branch'),tab('gencost')

class Dual:
    def __init__(self, v, g=None): self.v=v; self.g={} if g is None else g
    @staticmethod
    def cv(o): return o if isinstance(o,Dual) else Dual(I(o))
    def __add__(a,b):
        b=Dual.cv(b); g=a.g.copy()
        for k,v in b.g.items(): g[k]=g.get(k,I(0))+v
        return Dual(a.v+b.v,g)
    __radd__=__add__
    def __neg__(a):return Dual(-a.v,{k:-v for k,v in a.g.items()})
    def __sub__(a,b):return a+-Dual.cv(b)
    def __rsub__(a,b):return Dual.cv(b)+-a
    def __mul__(a,b):
        b=Dual.cv(b); g={k:v*b.v for k,v in a.g.items()}
        for k,v in b.g.items(): g[k]=g.get(k,I(0))+v*a.v
        return Dual(a.v*b.v,g)
    __rmul__=__mul__
    def __truediv__(a,b):
        b=Dual.cv(b); return a*Dual(1/b.v,{k:-v/(b.v*b.v) for k,v in b.g.items()})
    def __rtruediv__(a,b):return Dual.cv(b)/a
    def __pow__(a,n):
        if n != 2:raise NotImplementedError('only interval square supported')
        return Dual(a.v**2, {k:2*a.v*der for k,der in a.g.items()})
    def sin(a):return Dual(iv.sin(a.v),{k:v*iv.cos(a.v) for k,v in a.g.items()})
    def cos(a):return Dual(iv.cos(a.v),{k:-v*iv.sin(a.v) for k,v in a.g.items()})


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


def main(casefile,certfile):
    raw=open(casefile,'rb').read(); z=json.load(open(certfile))
    reject_nonfinite_certificate_values(z)
    assert hashlib.sha256(raw).hexdigest()==z['sha256'], 'SHA mismatch'
    B,bus,generator,branch,cost=parse(raw)
    bus=[row for row in bus if Decimal(row[1])!=4]
    n=len(bus);refs=[i for i,row in enumerate(bus) if Decimal(row[1])==3];assert len(refs)==1
    ref=refs[0]; ids={int(row[0]):i for i,row in enumerate(bus)}
    branch=[row for row in branch if Decimal(row[10])>0]
    variant=z.get('premises',{}).get('model_variant',z.get('variant','pglib'))
    assert variant in ('pglib','current')
    gens=[(k,row) for k,row in enumerate(generator) if Decimal(row[7])>0]; ng=len(gens)
    gb=[ids[int(row[0])] for _,row in gens]
    qfixed=[Decimal(row[3])==Decimal(row[4]) for _,row in gens]
    vbuses=sorted({gb[k] for k in range(ng) if not qfixed[k]}); slack=int(z['slack_gen'])
    assert 0<=slack<ng and Decimal(gens[slack][1][8]) != Decimal(gens[slack][1][9])
    coords=[H(s) for s in z['y_tilde']]; pgfix=[H(s) for s in z['pg_fix']]
    vmfix=[H(s) for s in z['vm_fix']]; qbase=[H(s) for s in z['qg_cand']]
    layout=[tuple(u) for u in z['layout']]; N=len(layout); assert N==2*n
    assert len(set(layout))==N and len(coords)==N
    assert set(x for x in layout if x[0]=='va')=={('va',i) for i in range(n) if i!=ref}
    assert set(x for x in layout if x[0]=='qb')=={('qb',i) for i in vbuses}
    assert set(x for x in layout if x[0]=='vm')=={('vm',i) for i in range(n) if i not in vbuses}
    assert [x for x in layout if x[0]=='ps']==[('ps',slack)]
    rad=F64(z['radius']); assert rad.a>0
    def fields(yy,grad=False):
        vv=[Dual(vmfix[i]) for i in range(n)]; aa=[Dual(I(0)) for _ in range(n)]
        pp=[Dual(pgfix[k]) for k in range(ng)]; qq=[Dual(qbase[k]) for k in range(ng)]
        qb={}; pgsl=None
        for k,(kind,j) in enumerate(layout):
            a=Dual(yy[k],{k:I(1)} if grad else {})
            if kind=='vm':vv[j]=a
            elif kind=='va':aa[j]=a
            elif kind=='qb':qb[j]=a
            elif kind=='ps':pgsl=a
        pp[slack]=pgsl
        for b in vbuses:
            free=[k for k in range(ng) if gb[k]==b and not qfixed[k]]
            fix=[k for k in range(ng) if gb[k]==b and qfixed[k]]
            fixed_total=sum((I(gens[k][1][3])/B for k in fix),I(0))
            rest=qb[b]-fixed_total
            if len(free)==1:qq[free[0]]=rest
            else:
                mean=sum((qbase[k] for k in free),I(0))/len(free)
                for k in free:qq[k]=rest/len(free)+(qbase[k]-mean)
        for i,row in enumerate(bus):
            if Decimal(row[11])==Decimal(row[12]):vv[i]=Dual(I(row[11]))
        for k,(_,row) in enumerate(gens):
            if Decimal(row[8])==Decimal(row[9]):pp[k]=Dual(I(row[8])/B)
            if Decimal(row[3])==Decimal(row[4]):qq[k]=Dual(I(row[3])/B)
        p=[Dual(-I(row[2])/B)-vv[i]**2*(I(row[4])/B) for i,row in enumerate(bus)]
        q=[Dual(-I(row[3])/B)+vv[i]**2*(I(row[5])/B) for i,row in enumerate(bus)]
        for k,j in enumerate(gb):p[j]=p[j]+pp[k];q[j]=q[j]+qq[k]
        flows=[]
        for edge in branch:
            i,j=ids[int(edge[0])],ids[int(edge[1])];r,x,bc=I(edge[2]),I(edge[3]),I(edge[4]); tau=I(edge[8]) if Decimal(edge[8]) else I(1)
            g=r/(r**2+x**2);b=-x/(r**2+x**2);phase=I(edge[9])*iv.pi/180
            delta=aa[i]-aa[j]-phase; cs=delta.cos(); sn=delta.sin(); vprod=vv[i]*vv[j]/tau
            pf=vv[i]**2*(g/(tau**2))-vprod*(cs*g+sn*b)
            qf=vv[i]**2*(-(b+bc/2)/(tau**2))-vprod*(sn*g-cs*b)
            pt=vv[j]**2*g-vprod*(cs*g-sn*b)
            qt=vv[j]**2*(-(b+bc/2))+vprod*(sn*g+cs*b)
            p[i]=p[i]-pf; q[i]=q[i]-qf;p[j]=p[j]-pt;q[j]=q[j]-qt
            flows.append((edge,i,j,pf,qf,pt,qt))
        return p+q,(vv,aa,pp,qq,flows)
    center,_=fields(coords,False); Fx=[u.v for u in center]
    Y=[iv.mpf([(u-rad).a,(u+rad).b]) for u in coords]
    Fbox,_=fields(Y,True)
    # Never rely on floating BLAS for an acceptance decision. Float LU constructs
    # only a candidate preconditioner, immediately embedded as exact dyadics.
    # The interval summations below certify the full Krawczyk operator.
    if N > 100:
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla
        from collections import defaultdict
        entries=[(i,j,x) for i,row in enumerate(Fbox) for j,x in row.g.items()]
        midsp=sp.csc_matrix(([float(x.mid) for i,j,x in entries],
                             ([i for i,j,x in entries],[j for i,j,x in entries])),shape=(N,N))
        try: lu=spla.splu(midsp.T.tocsc())
        except Exception as exc: raise UnsupportedFormat('sparse LU candidate unavailable') from exc
        K=[]
        block_size=2 # Bound inverse-row residence; no dense N x N inverse stored.
        for start in range(0,N,block_size):
            loc=list(range(start,min(start+block_size,N)))
            E=np.zeros((N,len(loc)), dtype=np.float64); E[loc,np.arange(len(loc))]=1.0
            cb=lu.solve(E).T   # each row is C_i of the candidate inverse
            for rnum,ii in enumerate(loc):
                c=[F64(float(v)) for v in cb[rnum]]
                cif=sum((c[k]*Fx[k] for k in range(N)),I(0))
                accum={}     # rigorous sparse column sums C_i * J(Y)
                for k in range(N):
                    if cb[rnum,k]==0.0: continue
                    for j, d in Fbox[k].g.items():
                        accum[j]=accum.get(j,I(0))+c[k]*d
                # Include EVERY omitted column through its exact Kronecker delta.
                # Unrepresented J entries are exactly zero by the AD derivative.
                radbound=I(0)
                for j in range(N):
                    rj=(I(1) if ii==j else I(0))-accum.get(j,I(0))
                    radbound+=abs(rj).b
                radbound*=rad.b
                cent=coords[ii]-cif
                K.append(iv.mpf([(cent-radbound).a,(cent+radbound).b]))
        print('KRAWCZYK_BACKEND','sparse-interval-rows','dim',N,'jac_nnz',len(entries),flush=True)
    else:
        J=[[v.g.get(j,I(0)) for j in range(N)] for v in Fbox]
        mid=np.array([[float(x.mid) for x in row] for row in J],float)
        C=np.linalg.inv(mid)  # candidate inverse ONLY; proof uses intervals of C
        Cb=[[F64(float(v)) for v in row] for row in C]
        K=[]
        for i in range(N):
            zval=coords[i]-sum((Cb[i][k]*Fx[k] for k in range(N)),I(0))
            for j in range(N):
                bij=(I(1) if i==j else I(0))-sum((Cb[i][k]*J[k][j] for k in range(N)),I(0))
                zval+=bij*iv.mpf([-rad.b,rad.b])
            K.append(zval)
    fail=[i for i in range(N) if not (K[i].a>Y[i].a and K[i].b<Y[i].b)]
    print('SHA256',hashlib.sha256(raw).hexdigest())
    print('case',z['case'],'variant',variant,'buses',n,'active_lines',len(branch),'free_q_buses',len(vbuses),'slack_gen',slack)
    print('Krawczyk', 'PASS' if not fail else 'FAIL','dimensions',N,'failures',fail[:8])
    if fail: raise SystemExit(2)
    _,(vv,aa,pp,qq,flows)=fields(K,False)
    violations=[];passed=0
    def checked(name,t):
        nonlocal passed
        if t is True:passed+=1
        else:violations.append(name)
    for i,row in enumerate(bus):
        if Decimal(row[11])==Decimal(row[12]):
            # Fixed binary-decimal input represents ONE exact model parameter.
            checked('vm'+str(i),True)
        else: checked('vm'+str(i),vv[i].v>=I(row[12]) and vv[i].v<=I(row[11]))
    for k,(_,row) in enumerate(gens):
        if Decimal(row[8])==Decimal(row[9]):
            checked('pg'+str(k),True)  # constructed directly from exact fixed P
        else: checked('pg'+str(k),pp[k].v>=I(row[9])/B and pp[k].v<=I(row[8])/B)
        if Decimal(row[3])==Decimal(row[4]):
            checked('qg'+str(k),True)  # constructed directly from exact fixed Q
        else: checked('qg'+str(k),qq[k].v>=I(row[4])/B and qq[k].v<=I(row[3])/B)
    for ed,i,j,pf,qf,pt,qt in flows:
        if Decimal(ed[5])>0:
            R=(I(ed[5])/B)**2
            if variant=='pglib':
                checked('thermal_f_'+ed[0]+'_'+ed[1],pf.v**2+qf.v**2<=R)
                checked('thermal_t_'+ed[0]+'_'+ed[1],pt.v**2+qt.v**2<=R)
            else:
                r,x,bc=I(ed[2]),I(ed[3]),I(ed[4]);YYr=r/(r**2+x**2);YYi=-x/(r**2+x**2)
                tau=I(ed[8]) if Decimal(ed[8]) else I(1);phi=I(ed[9])*iv.pi/180
                c,s=iv.cos(phi),iv.sin(phi);a,b=YYr,YYi+bc/2
                def cm(ra,ia,Vr,Vi):return ra*Vr-ia*Vi,ra*Vi+ia*Vr
                vi=(vv[i].v*iv.cos(aa[i].v),vv[i].v*iv.sin(aa[i].v))
                vj=(vv[j].v*iv.cos(aa[j].v),vv[j].v*iv.sin(aa[j].v))
                af=cm(a/(tau**2),b/(tau**2),*vi)
                bf=cm(-(YYr*c-YYi*s)/tau,-(YYr*s+YYi*c)/tau,*vj)
                at=cm(-(YYr*c+YYi*s)/tau,-(YYi*c-YYr*s)/tau,*vi)
                bt=cm(a,b,*vj)
                for side, u,v in [('f',af,bf),('t',at,bt)]:
                    rr_=u[0]+v[0];ii_=u[1]+v[1]
                    checked('current_'+side+'_'+ed[0]+'_'+ed[1],rr_**2+ii_**2<=R)
        if variant=='pglib':
            checked('angle_'+ed[0]+'_'+ed[1],aa[i].v-aa[j].v>=I(ed[11])*iv.pi/180 and aa[i].v-aa[j].v<=I(ed[12])*iv.pi/180)
    print('limits',passed,'/',passed+len(violations),'failures',violations[:10])
    if violations:raise SystemExit(3)
    val=I(0)
    for k,(gidx,row) in enumerate(gens):
        crow=cost[gidx]; m=pp[k].v*B
        val+=I(crow[4])*m**2+I(crow[5])*m+I(crow[6])
    # Convert output to an outward binary64 number. For a robust implementation
    # compare I(U) >= interval upper endpoint after nextafter and fail otherwise.
    U=math.nextafter(float(val.b),math.inf)
    if not (F64(U).a>=val.b):raise ArithmeticError('upper rounding failure')
    if not (F64(z['U']).a>=val.b):raise AssertionError('declared U not safe')
    print('U',repr(U),'declared',repr(z['U']),'ball',str(val))

class UnsupportedFormat(Exception): pass
class InputMismatch(Exception): pass

def classify(exc):
    if isinstance(exc, InputMismatch): return 5, 'INPUT_MISMATCH'
    if isinstance(exc, UnsupportedFormat): return 4, 'UNSUPPORTED_FORMAT'
    if isinstance(exc, (ValueError, KeyError, IndexError, TypeError)):
        if 'nonfinite' in str(exc).lower() or 'missing' in str(exc).lower():return 4,'UNSUPPORTED_FORMAT'
        return 6,'CHECKER_ERROR'
    if isinstance(exc, AssertionError):
        msg=str(exc)
        if 'SHA mismatch' in msg:return 5,'INPUT_MISMATCH'
        if 'declared U not safe' in msg:return 2,'REJECT_MATH'
        return 4,'UNSUPPORTED_FORMAT'
    if isinstance(exc, ArithmeticError):
        if 'rounding failure' in str(exc):return 6,'CHECKER_ERROR'
        return 2,'REJECT_MATH'
    return 6,'CHECKER_ERROR'

if __name__=='__main__':
    if len(sys.argv)!=3: raise SystemExit('usage: independent_ub_krawczyk_r12.py case.m ubcert.json')
    try:main(*sys.argv[1:])
    except SystemExit as e:
        if e.code==2:print('REJECT_MATH: Krawczyk containment did not hold',flush=True)
        elif e.code==3:print('REJECT_MATH: physical limits not verified',flush=True);raise SystemExit(2)
        raise
    except Exception as exc:
        code,category=classify(exc)
        print(category,type(exc).__name__,str(exc),flush=True)
        raise SystemExit(code)
