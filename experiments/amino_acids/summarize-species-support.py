"""Distinguish new species assignments from crossing a descriptive support threshold."""
import collections,csv,json,sys
from pathlib import Path
d=Path(sys.argv[1]);rows=list(csv.DictReader(open(d/'species-changes.tsv'),delimiter='\t'));counts=collections.defaultdict(collections.Counter)
for x in rows:
 if x['minimum_exemplars']!='3' or x['change']!='gained':continue
 before=int(x['baseline_exemplars']);assert before<3<=int(x['arm_exemplars'])
 for group in ['all_species']+(['human_annotated'] if x['human_annotation']=='1' else []):
  counts[x['arm']+'/'+group].update({'gained_at_three':1,'baseline_support_'+str(before):1})
summary={'minimum_exemplars':3,'counts':{k:dict(v) for k,v in counts.items()},'interpretation':'Descriptive breakdown of gained sample-species calls. Prior support of one or two is a threshold crossing, not an entirely new species assignment. Prior support of zero is an additional assignment, not known-truth biological detection.'}
(d/'support-transitions.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
