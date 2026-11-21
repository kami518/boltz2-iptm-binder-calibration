from pathlib import Path
import json
import re
from collections import Counter
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[2]
PDF = ROOT / 'publication/manuscript/PEDS_Main_and_Supplement_Review_Draft_v1.1.0.pdf'
BUILD = ROOT / 'publication/manuscript/review_pdf_build.json'
RESULT = ROOT / 'publication/manuscript/review_pdf_qc_results.json'
MAIN = ROOT / 'publication/manuscript/main_article.md'
SUPP = ROOT / 'publication/manuscript/supplementary_materials.md'

reader = PdfReader(str(PDF))
page_texts = [page.extract_text() or '' for page in reader.pages]
normalized_pages = [re.sub(r'\s+', ' ', text).strip() for text in page_texts]
full_text = ' '.join(normalized_pages)
full_text = full_text.replace('\u200b', '').replace('\x00', '')
build = json.loads(BUILD.read_text(encoding='utf-8'))
main_source = MAIN.read_text(encoding='utf-8')
supp_source = SUPP.read_text(encoding='utf-8')

main_tables = re.findall(r'^### (Table [1-4]\. .+)$', main_source, re.M)
main_figures = re.findall(r'^### (Figure [1-7]\. .+)$', main_source, re.M)
supp_tables = re.findall(r'^### (Table S(?:[1-9]|1[0-7])\. .+)$', supp_source, re.M)
supp_figures = re.findall(r'^### (Supplementary Figure S(?:[1-9]|1[01])\. .+)$', supp_source, re.M)
required_titles = main_tables + main_figures + supp_tables + supp_figures

checks = {}
def check(name, condition, detail):
    checks[name] = {'passed': bool(condition), 'detail': detail}

check('page_count_matches_build', len(reader.pages) == build['pages'],
      {'pdf_pages': len(reader.pages), 'build_pages': build['pages']})
check('file_size_matches_build', PDF.stat().st_size == build['bytes'],
      {'pdf_bytes': PDF.stat().st_size, 'build_bytes': build['bytes']})
wrong_footers = [n for n, text in enumerate(page_texts[1:], start=2)
                 if len(re.findall(rf'^Page {n}$', text, re.M)) != 1]
check('continuous_page_numbers', not wrong_footers,
      {'expected': 'Page 2 through Page ' + str(len(reader.pages)), 'wrong_pages': wrong_footers,
       'cover_unnumbered': not re.search(r'^Page 1$', page_texts[0], re.M)})
check('a4_pages', all(abs(float(p.mediabox.width) - 595.2756) < 1 and
                      abs(float(p.mediabox.height) - 841.8898) < 1 for p in reader.pages),
      'All pages are A4-sized.')

check('source_heading_counts', [len(main_tables), len(main_figures), len(supp_tables), len(supp_figures)] == [4, 7, 17, 11],
      {'main_tables': len(main_tables), 'main_figures': len(main_figures),
       'supplementary_tables': len(supp_tables), 'supplementary_figures': len(supp_figures)})
missing_titles = [title for title in required_titles if re.sub(r'\s+', ' ', title) not in full_text]
check('all_39_figure_and_table_titles_present', not missing_titles,
      {'expected': len(required_titles), 'missing': missing_titles})

image_counts = {n: len(page.images) for n, page in enumerate(reader.pages, start=1) if len(page.images)}
figure_title_pages = {}
for title in main_figures + supp_figures:
    # The contents pages repeat every title; resolve the unique body occurrence.
    matches = [n for n, text in enumerate(normalized_pages, start=1) if n > 4 and title in text]
    if len(matches) == 1:
        figure_title_pages[title] = matches[0]
missing_figure_images = [title for title, page in figure_title_pages.items() if image_counts.get(page, 0) != 1]
check('all_18_figures_embedded_on_title_pages',
      sum(image_counts.values()) == 18 and len(figure_title_pages) == 18 and not missing_figure_images,
      {'embedded_images': sum(image_counts.values()), 'image_pages': image_counts,
       'resolved_figure_title_pages': len(figure_title_pages),
       'titles_without_exactly_one_image': missing_figure_images})


def flatten_outline(items):
    for item in items:
        if isinstance(item, list):
            yield from flatten_outline(item)
        else:
            yield item

outline_items = list(flatten_outline(reader.outline))
outline_destinations = []
invalid_outline_pages = []
for item in outline_items:
    try:
        page_index = reader.get_destination_page_number(item)
        title = getattr(item, 'title', str(item))
        outline_destinations.append({'title': title, 'page': page_index + 1})
        if not 0 <= page_index < len(reader.pages):
            invalid_outline_pages.append(title)
    except Exception:
        invalid_outline_pages.append(getattr(item, 'title', str(item)))
outline_titles = {entry['title'] for entry in outline_destinations}
missing_outline_titles = [title for title in required_titles if title not in outline_titles]
check('bookmarks_and_destinations',
      len(outline_destinations) >= len(required_titles) and not invalid_outline_pages and not missing_outline_titles,
      {'destination_count': len(outline_destinations), 'invalid_destinations': invalid_outline_pages,
       'required_titles_missing_from_bookmarks': missing_outline_titles})

link_counts = Counter()
link_pages = {}
invalid_internal_links = []
invalid_external_uris = []
page_reference_ids = {page.indirect_reference.idnum: n for n, page in enumerate(reader.pages, start=1)}
for page_number, page in enumerate(reader.pages, start=1):
    for ref in page.get('/Annots', []):
        annotation = ref.get_object()
        if annotation.get('/Subtype') != '/Link':
            continue
        link_counts['total'] += 1
        link_pages[page_number] = link_pages.get(page_number, 0) + 1
        action = annotation.get('/A')
        action = action.get_object() if action else None
        uri = action.get('/URI') if action and action.get('/S') == '/URI' else None
        if uri:
            uri = str(uri)
            if uri.startswith(('https://', 'http://')):
                link_counts['web'] += 1
            else:
                link_counts['other_uri'] += 1
                invalid_external_uris.append({'page': page_number, 'uri': uri[:300]})
            continue
        destination = annotation.get('/Dest')
        if destination is None and action and action.get('/S') == '/GoTo':
            destination = action.get('/D')
        if destination is not None:
            link_counts['internal'] += 1
            if isinstance(destination, list) and destination:
                target = destination[0]
                target_id = getattr(target, 'idnum', None)
                if target_id not in page_reference_ids:
                    invalid_internal_links.append({'page': page_number, 'destination': str(destination)[:160]})
            elif str(destination) not in reader.named_destinations:
                invalid_internal_links.append({'page': page_number, 'destination': str(destination)[:160]})
        else:
            link_counts['unclassified'] += 1
check('working_link_annotations',
      link_counts['internal'] > 0 and link_counts['web'] > 0 and link_counts['local_file'] == 0
      and link_counts['other_uri'] == 0 and link_counts['unclassified'] == 0
      and not invalid_internal_links and not invalid_external_uris,
      {'counts': dict(link_counts), 'invalid_internal_links': invalid_internal_links,
       'invalid_external_uris': invalid_external_uris,
       'contents_pages_link_counts': {n: link_pages.get(n, 0) for n in (2, 3, 4)}})

anchors = [
    '2,084', '748', '678', '0.3086', '−0.0809', '45.1%', '42.9%',
    '25.9%', '39.5% to 49.8%', '38.4% to 46.8%', '19.3% to 31.9%',
    '0.1693', '0.1793', '+18.08', '−38.64', '2,074', '3,669', '3,766',
    '0.2624', '0.7254', '−0.7152', '+0.0466',
    '1,048', '821', '0.0011', '0.0211'
]
missing_anchors = [anchor for anchor in anchors if anchor not in full_text]
check('key_numerical_anchors_present', not missing_anchors,
      {'expected': anchors, 'missing': missing_anchors})
confirmed_author_metadata = ['KangMing Chen', '0009-0006-9634-5397', 'LeoChen070809@outlook.com',
                             'College of Biological Science and Medical Engineering']
missing_author_metadata = [marker for marker in confirmed_author_metadata if marker not in full_text]
check('confirmed_author_metadata_present', not missing_author_metadata,
      {'expected': confirmed_author_metadata, 'missing': missing_author_metadata})
declarations = [
    'Funding: This research received no specific grant from any funding agency.',
    'Conflicts of interest: The author declares no conflicts of interest.',
    'Acknowledgements: The author thanks Kai Li for supervision, guidance, and support during this work.',
    'Author contributions: KangMing Chen: Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Visualization, Writing - original draft, Writing - review & editing.'
]
missing_declarations = [marker for marker in declarations if marker not in full_text]
check('author_declarations_complete', not missing_declarations,
      {'expected': declarations, 'missing': missing_declarations})
unexpected_release_markers = ['REPOSITORY AND DOI TO BE ADDED', 'Not for submission', 'SCIENTIFIC REVIEW DRAFT']
check('publication_repository_and_reserved_doi_present',
      'github.com/kami518/boltz2-iptm-binder-calibration' in full_text
      and '10.5281/zenodo.24400262' in full_text
      and not any(marker in full_text for marker in unexpected_release_markers),
      {'repository_present': 'github.com/kami518/boltz2-iptm-binder-calibration' in full_text,
       'reserved_doi_present': '10.5281/zenodo.24400262' in full_text,
       'unexpected': [marker for marker in unexpected_release_markers if marker in full_text]})
resolved_source_markers = ['VERIFIED ANTHROPIC DATASET CITATION REQUIRED',
                           'EXTERNAL DATASET CITATION REQUIRED',
                           'VERIFIED ANTHROPIC EXTERNAL-DATASET REFERENCE TO BE INSERTED']
unexpected_source_markers = [marker for marker in resolved_source_markers if marker in full_text]
check('verified_source_placeholders_removed', not unexpected_source_markers,
      {'resolved_markers_absent': resolved_source_markers, 'unexpected': unexpected_source_markers})
raw_markers = {marker: full_text.count(marker) for marker in ['![', '**', '##', '�']}
check('no_raw_markdown_or_replacement_glyphs', not any(raw_markers.values()), raw_markers)
old_evidence_phrases = [
    'pilot value',
    'Pilot measurements can therefore improve allocation',
    'conditional decision boundary',
]
present_old_evidence_phrases = [phrase for phrase in old_evidence_phrases if phrase.casefold() in full_text.casefold()]
check('retired_pilot_wording_absent', not present_old_evidence_phrases,
      {'checked': old_evidence_phrases, 'present': present_old_evidence_phrases})
local_path_pattern = re.compile(r'(?i)(?:\bfile:|(?<![A-Za-z0-9])[A-Z]:[\\/]|<USERS>/|\\Users\\|<HOME>/)')
local_text_findings = sorted(set(local_path_pattern.findall(full_text)))
local_uri_findings = [item for item in invalid_external_uris if str(item['uri']).startswith('file:')]
check('no_local_paths_or_file_uris', not local_text_findings and not local_uri_findings,
      {'text_findings': local_text_findings, 'link_findings': local_uri_findings})
cohort_flow_anchors = ['2,185', '101', '2,758', '809', '1,949', '79', '60 Latent GDF-8']
missing_flow_anchors = [anchor for anchor in cohort_flow_anchors if anchor not in full_text]
check('cohort_flow_anchors_present', not missing_flow_anchors,
      {'expected': cohort_flow_anchors, 'missing': missing_flow_anchors})
nonempty_body_pages = [n for n, text in enumerate(normalized_pages, start=1)
                       if len(re.sub(r'SCIENTIFIC REVIEW DRAFT.*?Page \d+', '', text)) < 20]
check('no_text_empty_pages', not nonempty_body_pages,
      {'pages_with_less_than_20_extracted_characters_after_header_footer_removal': nonempty_body_pages})

result = {
    'status': 'passed' if all(entry['passed'] for entry in checks.values()) else 'failed',
    'pdf': PDF.relative_to(ROOT).as_posix(),
    'pages': len(reader.pages),
    'bytes': PDF.stat().st_size,
    'main_tables': len(main_tables),
    'main_figures': len(main_figures),
    'supplementary_table_guides': len(supp_tables),
    'supplementary_figures': len(supp_figures),
    'detailed_supplementary_panels': 'The original 75 detailed panels remain in the companion Excel workbook; post-review Table S17 is documented in the guide and linked to frozen CSV outputs.',
    'checks': checks,
}
RESULT.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(result, ensure_ascii=False, indent=2))
if result['status'] != 'passed':
    raise SystemExit(1)
