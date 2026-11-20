from pathlib import Path
import re
import sys

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
PATH = ROOT / 'manuscript/PEDS_Research_Article_English_draft_v2_2026-11-17.md'
text = PATH.read_text(encoding='utf-8')
WORD = re.compile(r"\b[A-Za-z0-9]+(?:[’'\-–][A-Za-z0-9]+)*\b")
checks = []

def check(name, condition, detail=''):
    checks.append({'name': name, 'passed': bool(condition), 'detail': str(detail)})


def wc(value):
    value = re.sub(r'<!--.*?-->', '', value, flags=re.S)
    value = re.sub(r'`[^`]*`', '', value)
    return len(WORD.findall(value))


def section(name, level=2):
    match = re.search(rf'^{"#" * level} {re.escape(name)}\s*$', text, re.M)
    if not match:
        raise ValueError(name)
    following = re.search(rf'^#{{1,{level}}} ', text[match.end():], re.M)
    end = match.end() + following.start() if following else len(text)
    return text[match.end():end]

abstract = section('Abstract')
intro = section('Introduction')
results = section('Results')
discussion = section('Discussion')
intro_paragraphs = [p for p in re.split(r'\n\s*\n', intro.strip()) if p.strip()]
counts = {name: wc(section(name)) for name in ['Abstract', 'Introduction', 'Methods', 'Results', 'Discussion']}
counts['Introduction + Results + Discussion'] = sum(counts[name] for name in ['Introduction', 'Results', 'Discussion'])
for name, value in counts.items():
    print(f'{name}: {value}')
check('abstract_max_250_words', counts['Abstract'] <= 250, counts['Abstract'])
check('introduction_four_paragraphs', len(intro_paragraphs) == 4, len(intro_paragraphs))
check('main_text_3500_to_5000_words', 3500 <= counts['Introduction + Results + Discussion'] <= 5000,
      counts['Introduction + Results + Discussion'])
check('title_focused_on_probability_across_targets',
      text.startswith('# Evaluating Boltz-2 ipTM as an experimental-success probability across binder-design targets'))

anchors = [
    '2,084', '748 primary and 678 strict', '0.675 (0.642–0.705)', '0.629 (0.593–0.667)',
    '0.620 (0.582–0.655)', '0.416', '−1.417', '0.3086', '0.8643', '−0.0809', '−1.1205',
    '−0.1393', '0.0153', 'four of six targets improving and two worsening',
    'three targets improved and three worsened', '+0.0014', '0.2590', '0.7120', '−0.0687',
    '−0.3211', '−0.2184', '+0.0809', '0.1918', '0.1861', 'seven of eight targets',
    'only in PD-L1 and TREM2', '+18.08', '−3.39', '−0.81', '−38.64', '2,074', '0.6366',
    '0.6770', '0.2949', '0.3139', '0.4561 to 0.6146', '0.6429–0.7182', '0.6341–0.7239',
    'Mythos Preview (480', 'Opus 4.8 (240', 'PXDesign (178)', 'RFdiffusion3 (138)',
    '0.1693', '45.1%', '0.3140', '0.1793', '42.9%', '6.5% to 71.7%',
    '225.08', '207.00', '+19.87', '+22.66'
]
missing = [anchor for anchor in anchors if anchor not in text]
check('all_required_numerical_anchors_present', not missing, missing)

check('results_contains_no_prior_work_citations',
      not re.search(r'Overath|Badonyi|Sagadinov|Smorodina|Garcia et al\.|Passaro et al\.|Lee et al\.|Van Calster|Steyerberg|Watson et al\.|Pacesa et al\.', results))
check('no_internal_rq_labels_in_results_or_discussion',
      not re.search(r'\bRQ[1-6][A-Z]?\b', results + '\n' + discussion))
check('plain_language_definitions_in_introduction',
      'We use “application to a new target” to mean' in intro and '“Local adaptation” then means' in intro)
check('rq4e_causal_boundary_present',
      'they do not estimate the causal contribution of baseline prevalence' in results
      and 'does not establish baseline prevalence as the largest causal mechanism' in discussion)
check('rq4c_cost_boundary_present',
      'after charging pilot measurements to the same budget' in results
      and 'not evidence of prospective monetary savings' in discussion)

required_authors = ['Abramson', 'Cao', 'Dauparas', 'Garcia', 'Jumper', 'Pacesa', 'Watson', 'Yang',
                    'Overath', 'Badonyi', 'Passaro', 'Sagadinov', 'Lee', 'Smorodina', 'Van Calster',
                    'Steyerberg', 'Roberts', 'Varma', 'Adaptyv Bio']
missing_authors = [author for author in required_authors if author not in text]
check('all_expected_reference_authors_present', not missing_authors, missing_authors)
for author in ['Abramson', 'Cao', 'Dauparas', 'Garcia', 'Jumper', 'Pacesa', 'Watson', 'Yang']:
    check(f'new_reference_cited_outside_bibliography_{author}', text.count(author) >= 2, text.count(author))

for phrase in ['single design-generation process', 'one design generator', 'six independent cohorts',
               'Pinc is universally inferior', 'cross-predictor transport was validated']:
    check(f'forbidden_claim_absent_{phrase}', phrase.lower() not in text.lower())

paths = sorted(set(re.findall(r'`((?:analysis|data)/[^`]+)`', text)))
missing_paths = [rel for rel in paths if not (ROOT / rel).exists()]
check('all_local_source_paths_exist', not missing_paths, missing_paths)
figure_links = re.findall(r'!\[[^\]]*\]\((figures/[^)]+\.png)\)', text)
missing_figures = [rel for rel in figure_links if not (ROOT / 'manuscript' / rel).is_file()]
check('seven_unique_figures_exist', len(figure_links) == 7 and len(set(figure_links)) == 7 and not missing_figures,
      {'count': len(figure_links), 'missing': missing_figures})

check('main_cohort_counts', (sum([133,268,1028,277,268,110]), sum([19,158,101,198,226,46]),
                             sum([19,155,100,159,204,41])) == (2084,748,678))
check('external_cohort_counts', (8*90, sum([20,47,23,19,19,39,65,55]),
                                 sum([18,36,20,19,17,35,62,47])) == (720,287,254))
check('original_v1_preserved', (ROOT / 'manuscript/PEDS_Research_Article_English_draft_v1_2026-11-15.md').is_file())

n_passed = sum(item['passed'] for item in checks)
print(f'\nQC: {n_passed}/{len(checks)} passed')
for item in checks:
    if not item['passed']:
        print(f"FAIL: {item['name']}: {item['detail']}")
if n_passed != len(checks):
    sys.exit(1)
print('All manuscript v2 checks passed.')
