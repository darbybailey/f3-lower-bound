#!/usr/bin/env python3
from fractions import Fraction
import json
from fastscan import tables
SC=10**40
# Exact Z24 lower endpoint and terminal state from the outward-rounded moment certificate.
lo24=30085007914684292896021976970576819728782
z=151034783145739756303329159115181
H=Fraction(lo24,SC)
steps=[(25,80,15,3910),(26,98,15,5878),(27,122,15,9151),(28,108,16,7633),
       (29,132,16,11464),(30,164,16,17779),(31,144,17,14536),(32,172,17,20803)]
for idx,p,q,r in steps:
    assert p%2==0
    cnt,mn,mx=tables(p,q)
    n=cnt[r]; ell=mn[r]; rawmax=mx[r]
    off=ell-1; t=rawmax-off; m=max(z,t); b=2*p-1
    s=(m+z+1,3*m+2*t+z+3,3*m+4*t+z+4)
    mu=Fraction(b**q-1,4)-off
    # Jensen: sum_{v in normalized B} 1/(S+v) >= |B|/(S+E[v]);
    # even-p complement symmetry gives E[raw v]=(b^q-1)/4 exactly.
    gain=n*sum(Fraction(1,sj+mu) for sj in s)
    stored=json.load(open(f'Z{idx}_jensen.json'))
    assert str(n)==stored['count']
    assert str(ell)==stored['block_min']
    assert str(rawmax)==stored['block_max_raw']
    assert str(gain.numerator)==stored['gain_num'] and str(gain.denominator)==stored['gain_den']
    z=3*m+5*t+z+4
    assert str(z)==stored['z_new']
    H += gain

def trunc(fr,places):
    a=fr.numerator*10**places//fr.denominator
    return f'{a//10**places}.{a%10**places:0{places}d}'
print('PASS')
print('f(3) >=',trunc(H,18))
assert trunc(H,18)=='3.008500858304852623'
