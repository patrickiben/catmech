#!/usr/bin/env python3
"""Compare dependent numerical replay to preserved v0.10.0 evidence."""
import csv,json,math,sys
from pathlib import Path
old,new=map(Path,sys.argv[1:3])
def records(p):
 r=list(csv.DictReader((p/"records.csv").open()))
 d={(z["cell"],z["replicate"],z["method"]):z for z in r}
 assert len(d)==len(r)==27600
 return d
a,b=records(old),records(new);assert a.keys()==b.keys()
assert (new/"COMPLETE.json").exists()
comparisons=0; recovered=[]
for k,x in a.items():
 y=b[k];assert y["status"]=="ok",k
 if x["status"]!="ok":recovered.append(k);continue
 for field,v in x.items():
  if field in ("source_fingerprint","analysis_seconds"):continue
  w=y[field]
  if v in ("NA","NaN",""):
   assert v==w,(k,field);continue
  try: u,z=float(v),float(w)
  except ValueError:assert v==w,(k,field)
  else:assert math.isclose(u,z,rel_tol=1e-10,abs_tol=1e-11),(k,field,u,z)
  comparisons+=1
assert len(recovered)==4,recovered
r=dict(status="PASS",scope="dependent numerical replay; not independent confirmation",unchanged_successful_analyses=27596,field_comparisons=comparisons,recovered_analyses=recovered)
(new/"comparison-with-v0.10.0.json").write_text(json.dumps(r,indent=2)+"\n")
print(json.dumps(r,indent=2))
