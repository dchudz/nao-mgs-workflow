"""Experimental DIAMOND-to-production-hit adapter; no fabricated SAM or DNA scores.

Input hits: qseqid sseqid bitscore evalue length pident qlen qstart qend
sstart send qframe slen qseq_gapped sseq_gapped (DIAMOND outfmt 6).
Taxonomy uses the unchanged production LCA implementation. Near-best hits are
selected per mate, then both mates contribute. Primary geometry is selected
independently per mate, deterministically. Original FASTQ records supply qualities.
"""
import argparse
import csv
import gzip
import json
import math
import sys
from itertools import groupby
from dataclasses import dataclass
from pathlib import Path

from diamond_geometry import Hsp, geometry
from union_pairs import pairs

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'modules/local/lcaTsv/resources/usr/bin'))
import lca_tsv as lca

HIT_FIELDS = 'qseqid sseqid bitscore evalue length pident qlen qstart qend sstart send qframe slen qseq_gapped sseq_gapped'.split()
LEGACY = 'aligner_length_normalized_score_mean aligner_length_normalized_score_mean_combined aligner_length_normalized_score_mean_artificial prim_align_fragment_length prim_align_best_alignment_score prim_align_best_alignment_score_rev prim_align_edit_distance prim_align_edit_distance_rev'.split()
MATE_FIELDS = 'query_len query_seq query_qual prim_align_ref_start prim_align_ref_start_unclipped prim_align_ref_end_unclipped prim_align_query_rc protein_bitscore protein_subject protein_query_frame'.split()
FIELDS = (['seq_id', 'sample'] + lca.get_output_header('aligner','seq_id','taxid','protein_bitscore')[1:]
          + LEGACY + ['prim_align_genome_id_all','prim_align_taxid_all','prim_align_pair_status','protein_aligned_mates']
          + MATE_FIELDS + [s+'_rev' for s in MATE_FIELDS])

def opener(path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)

@dataclass(frozen=True)
class Hit:
    query: str
    taxid: int
    score: float
    evalue: float
    length: int
    hsp: Hsp

    def order(self):
        return (-self.score, self.evalue, -self.length, self.hsp.subject,
                self.hsp.subject_start, self.hsp.subject_end, self.hsp.query_start, self.hsp.query_end)

def read_hits(path, genomes):
    with opener(path) as f:
        for values in csv.reader(f, delimiter='\t'):
            if len(values) != len(HIT_FIELDS):
                raise ValueError('Unexpected DIAMOND output fields')
            row = dict(zip(HIT_FIELDS, values))
            gid = row['sseqid'].rsplit('~', 2)[0]
            taxid, genome_length = genomes[gid]
            if int(row['length']) < 25:
                continue
            score, evalue = float(row['bitscore']), float(row['evalue'])
            if not math.isfinite(score) or score <= 0 or not math.isfinite(evalue) or evalue < 0:
                raise ValueError('Invalid DIAMOND score')
            h = Hsp(row['sseqid'], genome_length, int(row['slen']), int(row['qlen']),
                    int(row['qframe']), int(row['qstart']), int(row['qend']),
                    int(row['sstart']), int(row['send']), row['qseq_gapped'], row['sseq_gapped'])
            geometry(h)  # Validate every eligible HSP, including filtered-out queries.
            yield Hit(row['qseqid'], taxid, score, evalue, int(row['length']), h)


def query_pair(hit):
    pair, mate = hit.query.rsplit('/', 1)
    if mate not in ('1', '2'):
        raise ValueError('Protein query must identify mate 1 or 2')
    return pair


def near_best(hits):
    threshold = .9 * max((h.score for h in hits), default=0)
    return sorted((h for h in hits if h.score >= threshold), key=Hit.order)


def grouped_hits(path, genomes, wanted):
    previous = None
    for pair, records in groupby(read_hits(path, genomes), key=query_pair):
        if previous is not None and pair <= previous:
            raise ValueError('DIAMOND evidence must be sorted by query ID')
        previous = pair
        if pair not in wanted:
            for _ in records:
                pass
            continue
        mates = [[], []]
        for hit in records:
            mates[int(hit.query[-1])-1].append(hit)
        yield pair, [near_best(hs) for hs in mates]


def group_from_hits(seq_id, retained, artificial, unclassified):
    # Populate the exact production Subgroup representation in one pass; its
    # existing summarizer still decides classification, ties, scores and LCA.
    def subgroup(hits):
        return lca.Subgroup([h.taxid for h in hits],
                            [h.taxid not in unclassified for h in hits],
                            [h.score for h in hits], None)
    natural = [h for h in retained if h.taxid not in artificial]
    constructs = [h for h in retained if h.taxid in artificial]
    return lca.Group(seq_id, subgroup(retained), subgroup(natural), subgroup(constructs))


def paired_row(seq_id, sample, sequences, qualities, mate_hits, parents, artificial, unclassified):
    retained = [h for hs in mate_hits for h in sorted(hs,key=Hit.order)]
    if not retained:
        raise ValueError('Cannot create a hit row without an eligible alignment')
    group = group_from_hits(seq_id, retained, artificial, unclassified)
    group, _ = lca.summarize_group(group, parents, {})
    row = dict.fromkeys(FIELDS, 'NA')
    row.update(zip(lca.get_output_header('aligner','seq_id','taxid','protein_bitscore'), lca.output_group(group)))
    row['sample'] = sample
    primary = [min(hs, key=Hit.order) if hs else None for hs in mate_hits]
    gids = [h.hsp.subject.rsplit('~',2)[0] for h in primary if h]
    taxa = [str(h.taxid) for h in primary if h]
    row['prim_align_genome_id_all'] = '/'.join(dict.fromkeys(gids))
    row['prim_align_taxid_all'] = '/'.join(dict.fromkeys(taxa))
    # Protein HSPs do not establish Bowtie2 pair concordance or template length.
    row['protein_aligned_mates'] = len(gids)
    for i, suffix in enumerate(('', '_rev')):
        seq, qual, h = sequences[i], qualities[i], primary[i]
        if len(seq) != len(qual):
            raise ValueError('Sequence/quality length mismatch')
        row['query_seq'+suffix], row['query_qual'+suffix] = seq, qual
        row['query_len'+suffix] = len(seq)
        if h:
            if h.hsp.query_length != len(seq):
                raise ValueError('HSP query length does not match original read')
            g = geometry(h.hsp)
            for field in ('ref_start','ref_start_unclipped','ref_end_unclipped','query_rc'):
                row['prim_align_'+field+suffix] = getattr(g,field)
            row['protein_bitscore'+suffix] = h.score
            row['protein_subject'+suffix] = h.hsp.subject
            row['protein_query_frame'+suffix] = h.hsp.query_frame
    return row

def members(path):
    mapping = {}
    with opener(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['seq_id'] in mapping:
                raise ValueError('Repeated original read ID cannot safely enter production downstream')
            mapping[row['seq_id']] = row['query_pair']
    return mapping

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--hits', required=True); p.add_argument('--members', required=True)
    p.add_argument('--query-metadata', required=True)
    p.add_argument('--reads', required=True); p.add_argument('--sample', required=True)
    p.add_argument('--reference', required=True); p.add_argument('--genome-lengths', required=True)
    p.add_argument('--output', required=True); p.add_argument('--extract-only', action='store_true')
    a = p.parse_args(); ref = Path(a.reference)/'results'
    lengths = json.loads(Path(a.genome_lengths).read_text())
    with opener(ref/'total-virus-db-annotated.tsv.gz') as f:
        viral = {row['taxid'] for row in csv.DictReader(f,delimiter='\t')}
    with opener(ref/'virus-genome-metadata-gid.tsv.gz') as f:
        genomes = {r['genome_id']: (int(r['taxid'] if r['taxid'] in viral or r['species_taxid'] not in viral
                                         else r['species_taxid']), lengths[r['genome_id']])
                   for r in csv.DictReader(f,delimiter='\t')}
    membership = members(a.members)
    with opener(a.query_metadata) as f:
        query_sequences = {r['query_pair']:(r['query_seq'],r['query_seq_rev'])
                           for r in csv.DictReader(f,delimiter='\t')}
    if a.extract_only:
        # Eligibility needs existence, not all HSP objects or score sorting.
        eligible = {query_pair(h) for h in read_hits(a.hits, genomes)}
    else:
        parents, children = lca.parse_nodes_db(str(ref/'taxonomy-nodes.dmp'),81077)
        artificial = lca.get_descendants({81077},children)
        unclassified = lca.get_unclassified_taxids(lca.parse_names_db(str(ref/'taxonomy-names.dmp')))
        representatives = {}
        for key, raw in pairs(a.reads):
            sid,s1,q1,s2,q2 = [x.decode() for x in key]; sid=sid.removeprefix('@')
            query = membership[sid]
            if (s1,s2) != query_sequences[query]:
                raise ValueError('Original read sequence differs from the searched query pair')
            representatives.setdefault(query, ((s1,s2),(q1,q2)))
        templates = {}
        for query, hs in grouped_hits(a.hits, genomes, representatives):
            sequences, qualities = representatives[query]
            templates[query] = paired_row(query,a.sample,sequences,qualities,hs,parents,artificial,unclassified)
        eligible = templates.keys()
    count = 0; seen = set()
    with opener(a.output, 'wb' if a.extract_only else 'wt') as out:
        if not a.extract_only:
            writer = csv.DictWriter(out, fieldnames=FIELDS, delimiter='\t'); writer.writeheader()
        for key, raw in pairs(a.reads):
            sid,s1,q1,s2,q2 = [x.decode() for x in key]; sid=sid.removeprefix('@')
            if sid in seen:
                raise ValueError('Repeated original read ID')
            seen.add(sid)
            query=membership[sid]
            if (s1,s2) != query_sequences[query]:
                raise ValueError('Original read sequence differs from the searched query pair')
            if query not in eligible:continue
            if a.extract_only:out.write(raw)
            else:
                row = templates[query].copy()
                row.update(seq_id=sid,query_qual=q1,query_qual_rev=q2)
                if len(s1)!=len(q1) or len(s2)!=len(q2):
                    raise ValueError('Sequence/quality length mismatch')
                writer.writerow(row)
            count+=1
    Path(a.output+'.json').write_text(json.dumps({'pairs':count,'extract_only':a.extract_only,
        'score_semantics':'Protein bit scores are labeled separately. Legacy nucleotide scores and Bowtie2 pair status are NA.',
        'taxonomy':'Unchanged production natural/artificial/unclassified LCA rules over per-mate near-best eligible protein hits.',
        'primary':'Independent best HSP per mate; score, evalue, length, subject, then coordinate tie breaking.',
        'near_best_fraction':.9,'minimum_alignment_aa':25,
        'memory_bound':'At most one unique query pair of HSPs; compact per-query output templates are reused for exact sequence duplicates.'},indent=2)+'\n')

if __name__ == '__main__':main()
