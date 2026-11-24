import sys,json
from fractions import Fraction
sys.path.insert(0,'/mnt/data/f3_lower_frontier')
from fastscan import tables,weights
z=int(sys.argv[1]);p=int(sys.argv[2]);q=int(sys.argv[3]);r=int(sys.argv[4])
cnt,mn,mx=tables(p,q);n=cnt[r];ell=mn[r];rawmax=mx[r];off=ell-1;t=rawmax-off;m=max(z,t);b=2*p-1
s=(m+z+1,3*m+2*t+z+3,3*m+4*t+z+4);mu=Fraction(b**q-1,4)-off
g=n*sum(Fraction(1,sj+mu) for sj in s);zn=3*m+5*t+z+4
print(json.dumps({'p':p,'q':q,'r':r,'count':str(n),'block_min':str(ell),'block_max_raw':str(rawmax),'offset':str(off),'block_max_norm':str(t),'z_prev':str(z),'m':str(m),'s1':str(s[0]),'s2':str(s[1]),'s3':str(s[2]),'mu_num':str(mu.numerator),'mu_den':str(mu.denominator),'gain_num':str(g.numerator),'gain_den':str(g.denominator),'gain_decimal':format(float(g),'.17g'),'z_new':str(zn)},indent=2))
