from pathlib import Path
import csv, json, hashlib, sys, collections, math
import openpyxl
import pandas as pd

archive=Path(sys.argv[1]); repo=Path(sys.argv[2]); out=Path(sys.argv[3]); out.mkdir(parents=True,exist_ok=True)
roots=['Quercus','Nahuel','Salix','Gottin','Garraf-Montseny-Olot']
files=[]
for name in roots:
 p=archive/name
 files+=list((p/'raw-data').rglob('*'))
 files+=list(p.glob('*metadata*'))
 files+=list(p.glob('*processing*'))
files+=list((archive/'Salix/lib').glob('*'))
files+=[repo/'scripts/prepare_zero_presence.py',repo/'data/README.md']
files=sorted(set(p for p in files if p.is_file() and not p.name.startswith('~$')))
manifest=[];books=[];tables=[]
for p in files:
 manifest.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 if p.suffix=='.xlsx':
  wb=openpyxl.load_workbook(p,read_only=False,data_only=True)
  sheets=[]
  for ws in wb:
   rows=list(ws.values); nonempty=[(i+1,list(r)) for i,r in enumerate(rows) if any(x is not None for x in r)]
   nums=[float(x) for _,r in nonempty for x in r[1:] if isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)]
   notes=[dict(cell=c.coordinate,text=c.comment.text) for row in ws for c in row if c.comment]
   sheets.append(dict(sheet=ws.title,dimensions=ws.calculate_dimension(),nonempty_rows=len(nonempty),sample_rows=nonempty[:8],last_rows=nonempty[-2:],numeric_min=min(nums) if nums else None,numeric_max=max(nums) if nums else None,numeric_gt1=sum(x>1 for x in nums),fractional_numeric=sum(x!=int(x) for x in nums),comments=notes))
  books.append(dict(path=str(p),sheets=sheets));wb.close()
 elif p.suffix=='.csv':
  if p.name=='Salix_webs.csv':
   raw=list(csv.reader(p.open(encoding='utf-8-sig')));head=raw[10];rows=raw[11:]
   df=pd.DataFrame(rows,columns=[f'{i+1}:{v}' for i,v in enumerate(head)])
  else:
   df=pd.read_csv(p,sep=';' if 'all_interactions' in p.name else ',',dtype=str,keep_default_na=False).dropna(how='all')
  cols=[]
  for c in df:
   v=df[c];non=v[~v.isin(['','NA','NaN','nan'])];num=pd.to_numeric(non,errors='coerce')
   cols.append(dict(column=c,nonmissing=len(non),unique=non.nunique(),examples=list(non.unique()[:8]),numeric_min=float(num.min()) if num.notna().any() else None,numeric_max=float(num.max()) if num.notna().any() else None))
  tables.append(dict(path=str(p),rows=len(df),columns=cols))
 elif p.suffix=='.txt':
  df=pd.read_csv(p,sep='\t',header=None)
  values=df.apply(pd.to_numeric,errors='coerce').to_numpy().flatten();values=[float(x) for x in values if math.isfinite(x)]
  tables.append(dict(path=str(p),rows=len(df),cols=len(df.columns),numeric_max=max(values),positive_cells=sum(x>0 for x in values),values_gt1=sum(x>1 for x in values)))
pd.DataFrame(manifest).to_csv(out/'file_inventory.csv',index=False)
(out/'workbook_inspection.json').write_text(json.dumps(books,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
(out/'tabular_inspection.json').write_text(json.dumps(tables,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
print('Inventoried',len(files),'files;',len(books),'workbooks;',len(tables),'text tables')
for b in books:
 if 'Quercus' not in b['path'] or 'Web1 ' in b['path']:
  print(Path(b['path']).name,[(s['sheet'],s['nonempty_rows'],s['numeric_max']) for s in b['sheets']])
for t in tables:
 print(Path(t['path']).name,t['rows'],[(c['column'],c['nonmissing'],c['unique'],c['examples'][:3]) for c in t.get('columns',[])[:26] if c['nonmissing']])
