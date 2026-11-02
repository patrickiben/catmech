import csv,sys,statistics as st,math
rows=list(csv.DictReader(open(sys.argv[1])))
for r in rows:
    r['ok']=r['ok']=='TRUE'
    for k in 'fraction r estimate lower upper df b target icc'.split():
        r[k]=float(r[k]) if r[k] not in ('NA','') else None
key=lambda r:(r['dgp'],r['fraction'],r['r'])
orc={key(r):r['estimate'] for r in rows if r['variant']=='oracle'}
order=[]
for r in rows:
    if r['variant'] not in order: order.append(r['variant'])
print('dgp      frac variant        n   cover  medwidth wide rmse    recov  ratio  meanB   B/recov^2')
for d in sorted({r['dgp'] for r in rows}):
  for f in sorted({r['fraction'] for r in rows}):
    base=None
    for v in order:
      g=[r for r in rows if r['dgp']==d and r['fraction']==f and r['variant']==v and r['ok']]
      n=len(g); cov=100*sum(r['lower']<=r['target']<=r['upper'] for r in g)/n
      w=[r['upper']-r['lower'] for r in g]; mw=st.median(w); wide=sum(x>5*mw for x in w)
      rm=math.sqrt(sum((r['estimate']-r['target'])**2 for r in g)/n)
      rec=math.sqrt(sum((r['estimate']-orc[key(r)])**2 for r in g)/n)
      if v=='existing': base=rec
      ratio=(rec/base) if base and v!='oracle' else float('nan')
      mb=st.mean(r['b'] for r in g); mid=st.median(r['b'] for r in g)
      print(f"{d:8} {f:4} {v:13} {n:4} {cov:5.1f} {mw:8.4g} {wide:3d} {rm:7.4g} {rec:7.4g} {ratio:6.3f} {mb:9.3g} {(mid/rec**2 if rec else float('nan')):6.2f}")
    print('   median residual ICC', round(st.median(r['icc'] for r in rows if r['dgp']==d and r['fraction']==f),2))
