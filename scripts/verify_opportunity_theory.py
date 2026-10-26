from pathlib import Path
import csv, math, hashlib, json, sys
import numpy as np
import pandas as pd

root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
out=Path(sys.argv[2]) if len(sys.argv)>2 else root/'outputs/statistical description of local support/opportunity_theory'
env=pd.read_csv(out/'envelope_curves.csv')
vectors=pd.read_csv(out/'extremal_support_vectors.csv')
qc=pd.read_csv(out/'quality_control.csv')
worst=0.0
for (ds,e), g in env.groupby(['dataset','envelope']):
    q=qc[qc.dataset==ds]; n=int(q.N.iloc[0]); v=vectors[(vectors.dataset==ds)&(vectors.envelope==e)]
    for col,ks in [('empirical',q.K),('minimum',v.min_gap_support),('maximum',v.max_gap_support)]:
        for row in g.itertuples():
            exact=row.m/n-np.mean([math.comb(row.m,int(k))/math.comb(n,int(k)) if row.m>=k else 0 for k in ks])
            worst=max(worst,abs(exact-getattr(row,col)))
    assert np.all(g.minimum-1e-10<=g.empirical) and np.all(g.empirical<=g.maximum+1e-10)
assert worst<1e-12
pred=pd.read_csv(out/'prediction_curves.csv');met=pd.read_csv(out/'prediction_metrics.csv')
for row in met.itertuples():
    v=pred[(pred.dataset==row.dataset)&(pred.scope==row.scope)&(pred.model==row.model)&(pred.removal_fraction<=.8)]
    err=100*(v.predicted-v.empirical)
    assert np.isclose(abs(err).mean(),row.mae_pp)
    assert np.isclose(err.mean(),row.signed_error_pp)
    assert np.isclose(100*abs(err).mean()/(100*v.empirical.mean()),row.relative_error_percent)
fits=pd.read_csv(out/'fitted_parameters.csv')
assert all(z.dataset not in z.training_datasets.split(';') for z in fits[fits.scope=='LODO'].itertuples())
assert all(len(z.training_datasets.split(';'))==9 for z in fits[fits.scope=='LODO'].itertuples())
summary=pd.read_csv(out/'envelope_summary.csv').pivot(index='dataset',columns='envelope')
summary['narrowing_percent']=100*(1-summary['mean_width_pp']['B']/summary['mean_width_pp']['A'])
summary['narrowing_percent'].to_csv(out/'opportunity_narrowing.csv')
sc=pd.read_csv(out/'controlled_scenarios.csv')
sc['activity_error_percent']=100*(sc.realised_M/sc.target_M-1)
ag=sc.groupby(['dataset','model','h']).agg(mean_M=('realised_M','mean'),sd_M=('realised_M','std'),target_M=('target_M','first'),mean_activity_error_percent=('activity_error_percent','mean'),sd_activity_error_percent=('activity_error_percent','std'),mean_singleton=('singleton','mean'),mean_gap_pp=('integrated_gap_pp','mean')).reset_index()
ag.to_csv(out/'controlled_scenario_summary.csv',index=False)
print('Envelope error',worst)
print(summary['narrowing_percent'].to_string())
print('Max mean activity deviation %',ag.mean_activity_error_percent.abs().max())
print('Pair within improvements',sum(met[(met.scope=='Within')&(met.model=='Pair')].mae_pp.to_numpy()<met[(met.scope=='Within')&(met.model=='Common')].mae_pp.to_numpy()))
print('DIAGNOSTIC COVERAGE')
dg=pd.read_csv(out/'prediction_diagnostics.csv')
print(dg[dg.scope=='LODO'].groupby(['model','statistic']).covered.sum().to_string())
print('WITHIN ERRORS')
print(met[met.scope=='Within'][['dataset','model','mae_pp','brier']].to_string(index=False))
print('DIAGNOSTIC FAILURES')
print(dg[(dg.scope=='LODO')&(dg.statistic.isin(['M','singleton','tail_ge5']))][['dataset','model','statistic','observed','predicted','covered']].to_string(index=False))
(out/'verification.json').write_text(json.dumps({'exact_combination_max_difference':worst,'curve_metrics_verified':len(met),'fold_exclusion_verified':True},indent=2))
