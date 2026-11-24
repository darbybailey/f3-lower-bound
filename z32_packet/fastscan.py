from fractions import Fraction
import time,sys

def weights(p):
 h=p//2
 return [((k-h)*(k-h+1))//2 for k in range(p)]

def tables(p,q):
 w=weights(p); wmax=max(w); R=q*wmax
 cnt=[0]*(R+1); mn=[None]*(R+1); mx=[None]*(R+1);cnt[0]=1;mn[0]=mx[0]=0
 b=2*p-1
 for i in range(q):
  pw=b**i; nr=(i+1)*wmax
  nc=[0]*(R+1); nmn=[None]*(R+1); nmx=[None]*(R+1)
  oldmax=i*wmax
  for s in range(oldmax+1):
   c=cnt[s]
   if not c: continue
   a=mn[s];x=mx[s]
   for d,wd in enumerate(w):
    ns=s+wd; val=d*pw
    nc[ns]+=c
    va=a+val; vx=x+val
    if nmn[ns] is None or va<nmn[ns]: nmn[ns]=va
    if nmx[ns] is None or vx>nmx[ns]: nmx[ns]=vx
  cnt,mn,mx=nc,nmn,nmx
 return cnt,mn,mx

def best_for(z,p,q,window=None):
 cnt,mn,mx=tables(p,q); b=2*p-1; mean=q*sum(weights(p))/p
 lo=max(0,int(mean-(window or 150))); hi=min(len(cnt)-1,int(mean+(window or 50)))
 best=None
 for r in range(lo,hi+1):
  n=cnt[r]
  if not n: continue
  ell=mn[r]; t=mx[r]-ell+1; m=max(z,t)
  s=(m+z+1,3*m+2*t+z+3,3*m+4*t+z+4)
  mu=Fraction(b**q-1,4)-(ell-1)
  g=n*sum(Fraction(1,sj+mu) for sj in s)
  if best is None or g>best[0]: best=(g,r,n,t,ell,3*m+5*t+z+4)
 return best

if __name__=='__main__':
 z=int(sys.argv[1]); q=int(sys.argv[2]); ps=list(map(int,sys.argv[3:]))
 for p in ps:
  t=time.time(); best=best_for(z,p,q)
  print(p,'r',best[1],'gain',float(best[0]),'t',best[3],'znew',best[5],'sec',round(time.time()-t,2),flush=True)
