#!/usr/bin/env python3
"""Independent chordal ACOPF lower-bound certificate checker.

Usage: python3 independent_chordal_r13.py cases/pglib_opf_case14_ieee.m lbchordal_case14_ieee.json
Exit 0 = ACCEPT; 1 = REJECT/ERROR. Uses no repository imports.
Inputs ONLY from command-line files. No imports from research-hub.
Interval checks use mpmath.iv; double eigenvalues only propose PSD shifts.
For cases above 200 buses use branch-local sparse polynomial reconstruction and sparse coefficient accumulation.
For small-case fast path comparison set TOP035_FORCE_SPARSE_CHORDAL=1.
No generator or source-verifier imports.
NEVER accept a certificate unless all interval PSD and residual tests complete.
"""
import sys, re, json, hashlib, math
from decimal import Decimal
import numpy as np
from mpmath import iv
iv.dps=80
I=lambda x: x if hasattr(x, "ctx") else iv.mpf(str(x))
def H(s):
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
J=lambda x=1: (I(0),I(x))

def C(a=0,b=0):return (I(a),I(b))
def ca(x,y):return (x[0]+y[0],x[1]+y[1])
def cn(x):return (-x[0],-x[1])
def cs(x,r):return (x[0]*r,x[1]*r)
def cm(x,y):return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def cc(x):return (x[0],-x[1])
def cd(x,y):return cm(x,cs(cc(y),I(1)/(y[0]*y[0]+y[1]*y[1])))
def abup(v):return abs(v).b
def cop(x):return (x[0],x[1])

def parse(raw):
    st=raw.decode('utf-8');bm=re.search(r'mpc\.baseMVA\s*=\s*([^;]+);',st)
    assert bm,'baseMVA missing'
    def tab(k):
        m=re.search(r'mpc\.'+k+r'\s*=\s*\[([\s\S]*?)\]\s*;',st)
        assert m,k
        return [line.split('%')[0].replace(';',' ').split() for line in m.group(1).splitlines()
                if line.split('%')[0].replace(';',' ').strip()]
    return I(bm.group(1).strip()),tab('bus'),tab('gen'),tab('branch'),tab('gencost')

def strict_psd(A):
    """Interval LDL with positive interval pivots: sufficient proof A >= 0."""
    n=len(A);D=[I(0)]*n;L=[[I(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        p=A[i][i]
        for k in range(i):p-=L[i][k]*L[i][k]*D[k]
        if not p.a>0:return False,i
        D[i]=p
        for j in range(i+1,n):
            y=A[j][i]
            for k in range(i):y-=L[j][k]*L[i][k]*D[k]
            L[j][i]=y/D[i]
    return True,n

def shift_psd(A):
    """Dyadic-ish bisection of positive PSD correction; every decision via interval LDL."""
    n=len(A)
    mid=np.array([[float(x.mid) for x in row] for row in A]); mid=(mid+mid.T)/2
    sc=max(1.,float(np.max(np.abs(mid)))); eig=float(np.linalg.eigvalsh(mid)[0])
    low=0.; high=max(sc*1e-16, -eig+sc*1e-14)
    if strict_psd(A)[0]:high=0.
    else:
        for _ in range(64):
            shifted=[[A[i][j]+(F64(float(high)) if i==j else I(0)) for j in range(n)] for i in range(n)]
            if strict_psd(shifted)[0]:break
            high*=2
        else:raise ArithmeticError('PSD bracketing failed')
        for _ in range(26):
            med=(low+high)/2
            dd=F64(float(med))
            shifted=[[A[i][j]+(dd if i==j else I(0)) for j in range(n)] for i in range(n)]
            if strict_psd(shifted)[0]:high=med
            else:low=med
    delta=F64(float(high))
    if not strict_psd([[A[i][j]+(delta if i==j else I(0)) for j in range(n)] for i in range(n)])[0]:
        raise ArithmeticError('PSD proof failed')
    return delta

def min_quad(a2,b,c,lo,hi):
    f=lambda x:(a2*x+b)*x+c
    if a2==I(0):return min(f(lo).a,f(hi).a)
    if not a2.a>0:raise ValueError('nonconvex/ambiguous quadratic cost')
    xs=-b/(2*a2)
    if xs.b<lo.a:return f(lo).a
    if xs.a>hi.b:return f(hi).a
    return (c-b*b/(4*a2)).a


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


def main(mpath,cpath):
    raw=open(mpath,'rb').read();cert=json.load(open(cpath));sha=hashlib.sha256(raw).hexdigest()
    reject_nonfinite_certificate_values(cert)
    assert sha==cert['sha256'],'sha256 mismatch'
    base,bus0,gen0,branch0,gc0=parse(raw)
    bus=[row for row in bus0 if Decimal(row[1])!=4]
    pos={int(row[0]):i for i,row in enumerate(bus)};n=len(bus)
    gen=[(r,c) for r,c in zip(gen0,gc0) if Decimal(r[7])>0]
    ng=len(gen);gb=[pos[int(r[0])] for r,_ in gen]
    lines=[r for r in branch0 if Decimal(r[10])>0]
    variant=cert.get('variant','pglib');assert variant in ('pglib','current')
    offs=cert['offs'];yy=[H(x) for x in cert['y']]
    def piece(name,nonneg=False):
        a,b=offs[name];result=yy[a:b]
        return [(x if x.a>=0 else I(0)) if nonneg else x for x in result]
    lam=piece('lam');mu=piece('mu');nup=piece('nup',True);nlo=piece('nlo',True)
    aaa=piece('a');bbb=piece('b');ep=piece('etp',True);em=piece('etm',True)
    cur=piece('cur',True) if 'cur' in offs else []
    assert len(lam)==len(mu)==len(nup)==len(nlo)==n
    active=[k for k,r in enumerate(lines) if Decimal(r[5])>0]
    if variant=='pglib':assert len(aaa)==len(bbb)==2*len(active) and len(ep)==len(em)==len(lines) and not cur
    else:assert len(cur)==2*len(active) and not aaa and not ep
    vmin=[I(r[12]) for r in bus];vmax=[I(r[11]) for r in bus]
    # Generator-only constant terms after minimising P_g,Q_g.
    Lconst=I(0)
    for k,(r,c) in enumerate(gen):
        i=gb[k]
        a2=I(c[4])*base*base;b1=I(c[5])*base-lam[i];c0=I(c[6])
        Lconst+=min_quad(a2,b1,c0,I(r[9])/base,I(r[8])/base)
        Lconst+=min_quad(I(0),-mu[i],I(0),I(r[4])/base,I(r[3])/base)
    for i,r in enumerate(bus):
        Lconst+=lam[i]*I(r[2])/base+mu[i]*I(r[3])/base-nup[i]*vmax[i]**2+nlo[i]*vmin[i]**2
    for k,e in enumerate(active):
        if variant=='pglib':
            for t in (2*k,2*k+1):
                Lconst-=iv.sqrt(aaa[t]**2+bbb[t]**2)*(I(lines[e][5])/base)
        else:
            for t in (2*k,2*k+1):Lconst-=cur[t]*(I(lines[e][5])/base)**2
    # Define interval quadratic form q(V) from MODEL.tex directly.
    if n>200 or __import__('os').environ.get('TOP035_FORCE_SPARSE_CHORDAL')=='1':
        # Fast path for medium/large cases.  The physical OPF Lagrangian contains
        # bus-diagonal terms and branch-local couplings only.  Reconstruct those
        # interval coefficients *from the original model*, not from certificate Z.
        # This eliminates O((n+|E|)*|E|) repeated qform calculations and all n*n
        # interval storage.  FLOAT is used ONLY inside shift_psd as a candidate.
        active_index={edge:t for t,edge in enumerate(active)}
        diag=[(lam[i]*I(row[4])/base-mu[i]*I(row[5])/base+nup[i]-nlo[i])
              for i,row in enumerate(bus)]
        Q_sparse={}
        def addq(i,j,v):
            old=Q_sparse.get((i,j))
            Q_sparse[i,j]=ca(old,v) if old is not None else v
        def branch_qform(k,r,i,j,Vi,Vj):
            # Identical physical terms to the R12 qform, but just on one branch.
            r0,x0=I(r[2]),I(r[3]);Y=cd(C(1),C(r0,x0));yc=cc(Y)
            bc=I(r[4]);tau=I(r[8]) if Decimal(r[8]) else I(1)
            phi=I(r[9])*iv.pi/180;T=C(tau*iv.cos(phi),tau*iv.sin(phi));Tc=cc(T)
            A=ca(yc,cn(J(bc/2)))
            Wi=cm(Vi,cc(Vi));Wj=cm(Vj,cc(Vj))
            Wij=cm(Vi,cc(Vj));Wji=cc(Wij)
            sf=ca(cd(cm(A,Wi),cm(T,Tc)),cn(cd(cm(yc,Wij),T)))
            st=ca(cm(A,Wj),cn(cd(cm(yc,Wji),Tc)))
            out=lam[i]*sf[0]+mu[i]*sf[1]+lam[j]*st[0]+mu[j]*st[1]
            if variant=='pglib':
                mx=I(r[12])*iv.pi/180;mn=I(r[11])*iv.pi/180
                assert mx.b<iv.pi.a/2 and mn.a>-(iv.pi.b/2)
                out+=ep[k]*(Wij[1]-iv.tan(mx)*Wij[0])+em[k]*(iv.tan(mn)*Wij[0]-Wij[1])
                if Decimal(r[5])>0:
                    t=active_index[k]
                    out-=aaa[2*t]*sf[0]+bbb[2*t]*sf[1]+aaa[2*t+1]*st[0]+bbb[2*t+1]*st[1]
            elif Decimal(r[5])>0:
                t=active_index[k]
                a=ca(Y,J(bc/2))
                If=ca(cd(cm(a,Vi),cm(T,Tc)),cn(cd(cm(Y,Vj),Tc)))
                It=ca(cn(cd(cm(Y,Vi),T)),cm(a,Vj))
                out+=cur[2*t]*cm(If,cc(If))[0]+cur[2*t+1]*cm(It,cc(It))[0]
            return out
        one=C(1);zero=C();imag=J(1)
        # Preserve hermiticity at the coefficient level.  Each edge contribution
        # is reconstructed by exact interval polarization of its two-bus form.
        for k,r in enumerate(lines):
            i=pos[int(r[0])];j=pos[int(r[1])]
            di=branch_qform(k,r,i,j,one,zero)
            dj=branch_qform(k,r,i,j,zero,one)
            # For the off-diagonal evaluate both polarization points.
            re=(branch_qform(k,r,i,j,one,one)-di-dj)/2
            im=(di+dj-branch_qform(k,r,i,j,one,imag))/2
            diag[i]+=di;diag[j]+=dj
            addq(i,j,C(re,im));addq(j,i,C(re,-im))
        for i,v in enumerate(diag):addq(i,i,C(v))
        # PSD check from R12 is unchanged, including interval compressed Z,
        # shift, and the strictly outward penalty accounting.
        sumZ={};penalty=I(0)
        assert len(cert['X'])==len(cert['cliques'])
        for clique,rawX in zip(cert['cliques'],cert['X']):
            s=len(clique);assert len(rawX)==4*s*s
            x=[H(t) for t in rawX]
            A=[[(x[i*2*s+j]+x[j*2*s+i])/2 for j in range(2*s)] for i in range(2*s)]
            Rm=[[I(0) for _ in range(2*s)] for _ in range(2*s)]
            reZ=[[I(0) for _ in range(s)] for _ in range(s)]
            imZ=[[I(0) for _ in range(s)] for _ in range(s)]
            for p in range(s):
                for q in range(s):
                    re=(A[p][q]+A[s+p][s+q])/2
                    im=(A[s+p][q]-A[p][s+q])/2
                    reZ[p][q]=re;imZ[p][q]=im
                    Rm[p][q]=re;Rm[s+p][s+q]=re
                    Rm[p][s+q]=-im;Rm[s+p][q]=im
            delta=shift_psd(Rm)
            clique_penalty=delta*sum((vmax[i]**2 for i in clique),I(0))
            penalty+=clique_penalty
            if float(clique_penalty.b)>0:
                print('PSD_BLOCK',len(clique),'COMPRESSED_DELTA',str(delta.b),
                      'PENALTY_UPPER',str(clique_penalty.b),flush=True)
            for p,i in enumerate(clique):
                for q,j in enumerate(clique):
                    val=C(reZ[p][q],imZ[p][q]);old=sumZ.get((i,j))
                    sumZ[i,j]=ca(old,val) if old is not None else val
        residual=I(0);maxv=I(0)
        for i,j in (Q_sparse.keys() | sumZ.keys()):
            r=ca(Q_sparse.get((i,j),C()),cn(sumZ.get((i,j),C())))
            m=iv.sqrt(r[0]**2+r[1]**2).b
            residual+=m*vmax[i]*vmax[j]
            if m.b>maxv.b:maxv=m
        print('H1_SPARSE_PROOF', 'model_branch_edges',len(lines),
              'nonzero_coeff_slots',len(Q_sparse),'nonzero_Z_slots',len(sumZ),
              'residual_slots',len(Q_sparse.keys() | sumZ.keys()),flush=True)
    else:
        bybus=[[] for _ in range(n)]
        for k,r in enumerate(lines):
            i=pos[int(r[0])];j=pos[int(r[1])]
            r0,x0=I(r[2]),I(r[3]);Y=cd(C(1),C(r0,x0));yc=cc(Y)
            bc=I(r[4]);tau=I(r[8]) if Decimal(r[8]) else I(1)
            phi=I(r[9])*iv.pi/180;T=C(tau*iv.cos(phi),tau*iv.sin(phi))
            A=ca(yc,cn(J(bc/2)))
            bybus[i].append((k,0,i,j,Y,yc,A,T,bc));bybus[j].append((k,1,i,j,Y,yc,A,T,bc))
        def qform(V):
            out=I(0)
            for i,row in enumerate(bus):
                wi=cm(V[i],cc(V[i]))[0]
                out+=(lam[i]*I(row[4])/base-mu[i]*I(row[5])/base+nup[i]-nlo[i])*wi
            for k,r in enumerate(lines):
                i=pos[int(r[0])];j=pos[int(r[1])]
                r0,x0=I(r[2]),I(r[3]);Y=cd(C(1),C(r0,x0));yc=cc(Y)
                bc=I(r[4]);tau=I(r[8]) if Decimal(r[8]) else I(1);phi=I(r[9])*iv.pi/180
                T=C(tau*iv.cos(phi),tau*iv.sin(phi));Tc=cc(T)
                A=ca(yc,cn(J(bc/2)))
                Wi=cm(V[i],cc(V[i]));Wj=cm(V[j],cc(V[j]));Wij=cm(V[i],cc(V[j]));Wji=cc(Wij)
                sf=ca(cd(cm(A,Wi),cm(T,Tc)),cn(cd(cm(yc,Wij),T)))
                st=ca(cm(A,Wj),cn(cd(cm(yc,Wji),Tc)))
                out+=lam[i]*sf[0]+mu[i]*sf[1]+lam[j]*st[0]+mu[j]*st[1]
                if variant=='pglib':
                    mx=I(r[12])*iv.pi/180;mn=I(r[11])*iv.pi/180
                    assert mx.b<iv.pi.a/2 and mn.a>-(iv.pi.b/2)
                    out+=ep[k]*(Wij[1]-iv.tan(mx)*Wij[0])+em[k]*(iv.tan(mn)*Wij[0]-Wij[1])
                    if Decimal(r[5])>0:
                        t=active.index(k)
                        out-=aaa[2*t]*sf[0]+bbb[2*t]*sf[1]+aaa[2*t+1]*st[0]+bbb[2*t+1]*st[1]
                elif Decimal(r[5])>0:
                    t=active.index(k)
                    a=ca(Y,J(bc/2))
                    If=ca(cd(cm(a,V[i]),cm(T,Tc)),cn(cd(cm(Y,V[j]),Tc)))
                    It=ca(cn(cd(cm(Y,V[i]),T)),cm(a,V[j]))
                    out+=cur[2*t]*cm(If,cc(If))[0]+cur[2*t+1]*cm(It,cc(It))[0]
            return out
        zero=[C() for _ in range(n)];diag=[]
        for i in range(n):
            V=zero.copy();V[i]=C(1);diag.append(qform(V))
        Q=[[C() for _ in range(n)] for _ in range(n)]
        for i in range(n):Q[i][i]=C(diag[i])
        pair=set()
        for line in lines:
            i=pos[int(line[0])];j=pos[int(line[1])];pair.add((min(i,j),max(i,j)))
        for clique in cert['cliques']:
            for i in clique:
                for j in clique:
                    if i<j:pair.add((i,j))
        for i,j in pair:
            V=zero.copy();V[i]=C(1);V[j]=C(1)
            re=(qform(V)-diag[i]-diag[j])/2
            V[j]=J(1)
            im=(diag[i]+diag[j]-qform(V))/2
            Q[i][j]=C(re,im);Q[j][i]=C(re,-im)
        sumZ=[[C() for _ in range(n)] for _ in range(n)]
        penalty=I(0)
        assert len(cert['X'])==len(cert['cliques'])
        for clique,rawX in zip(cert['cliques'],cert['X']):
            s=len(clique);assert len(rawX)==4*s*s
            x=[H(t) for t in rawX]
            A=[[ (x[i*2*s+j]+x[j*2*s+i])/2 for j in range(2*s)] for i in range(2*s)]
            # IMPORTANT: the dual certificate contains a REAL lifted X, but only its
            # COMPLEX COMPRESSION Z contributes to M. X can have negative directions
            # in the kernel of the compression without invalidating the bound.
            # Requiring X+delta*I PSD is mathematically sufficient but can waste
            # significant delta (case30_as/case30_ieee/case89). Certify Z via the
            # real 2s-by-2s representation [Re Z, -Im Z; Im Z, Re Z].
            Rm=[[I(0) for _ in range(2*s)] for _ in range(2*s)]
            reZ=[[I(0) for _ in range(s)] for _ in range(s)]
            imZ=[[I(0) for _ in range(s)] for _ in range(s)]
            for p in range(s):
                for q in range(s):
                    re=(A[p][q]+A[s+p][s+q])/2
                    im=(A[s+p][q]-A[p][s+q])/2
                    reZ[p][q]=re; imZ[p][q]=im
                    Rm[p][q]=re;Rm[s+p][s+q]=re
                    Rm[p][s+q]=-im;Rm[s+p][q]=im
            delta=shift_psd(Rm)
            clique_penalty=delta*sum((vmax[i]**2 for i in clique),I(0))
            penalty+=clique_penalty
            if float(clique_penalty.b)>0:
                print('PSD_BLOCK',len(clique),'COMPRESSED_DELTA',str(delta.b),
                      'PENALTY_UPPER',str(clique_penalty.b),flush=True)
            for p,i in enumerate(clique):
                for q,j in enumerate(clique):
                    sumZ[i][j]=ca(sumZ[i][j],C(reZ[p][q],imZ[p][q]))
        residual=I(0);maxv=I(0)
        for i in range(n):
            for j in range(n):
                r=ca(Q[i][j],cn(sumZ[i][j]))
                m=iv.sqrt(r[0]**2+r[1]**2).b
                residual+=m*vmax[i]*vmax[j]
                if m.b>maxv.b:maxv=m
    lb=Lconst-penalty-residual
    L=float(lb.a);L=math.nextafter(L,-math.inf)
    print('L_DIAGNOSTIC_CONSTANT',str(Lconst),flush=True)
    print('L_DIAGNOSTIC_PSD_PENALTY',str(penalty.b),flush=True)
    print('L_DIAGNOSTIC_CHORDAL_RESIDUAL',str(residual.b),flush=True)
    print('L_DIAGNOSTIC_VERIFIED_UNROUNDED',str(lb),flush=True)
    print('L_DIAGNOSTIC_DECLARED',repr(cert['L']),flush=True)
    if not F64(cert['L']).b<=lb.a:raise ValueError('declared L exceeds verified lower bound')
    print('INPUT_SHA256',sha)
    print('CASE',cert['case'],'variant',variant,'n',n,'cliques',len(cert['cliques']))
    print('CONSTANT',Lconst,'PSD_PENALTY_UPPER',penalty.b,'CHORDAL_RESIDUAL_UPPER',residual.b)
    print('RIGOROUS_L',repr(L),'DECLARED_L',repr(cert['L']))
    print('ACCEPT exit=0')

if __name__=='__main__':
    if len(sys.argv)!=3:raise SystemExit('usage: independent_chordal_r13.py case.m lbchordal.json')
    try: main(*sys.argv[1:])
    except Exception as exc:
        msg=str(exc)
        if 'sha256 mismatch' in msg or 'SHA mismatch' in msg:code=5;kind='INPUT_MISMATCH'
        elif 'declared L exceeds verified lower bound' in msg:code=2;kind='REJECT_MATH'
        elif isinstance(exc,(KeyError,ValueError,AssertionError,IndexError)):code=4;kind='UNSUPPORTED_FORMAT'
        elif isinstance(exc,ArithmeticError):code=2;kind='REJECT_MATH'
        else:code=6;kind='CHECKER_ERROR'
        print(kind,type(exc).__name__,msg,flush=True)
        raise SystemExit(code)
