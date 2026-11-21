from pathlib import Path
import html
import json
import re

import matplotlib
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate, Frame, HRFlowable, Image, PageBreak, PageTemplate,
    Paragraph, Spacer, Table, TableStyle
)
from reportlab.platypus.tableofcontents import TableOfContents

ROOT = Path(__file__).resolve().parents[2]
MAIN = ROOT / 'publication/manuscript/main_article.md'
SUPP = ROOT / 'publication/manuscript/supplementary_materials.md'
OUTDIR = ROOT / 'publication/manuscript'
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTPUT = OUTDIR / 'PEDS_Main_and_Supplement_Review_Draft_v1.1.0.pdf'
BUILD_INFO = OUTDIR / 'review_pdf_build.json'
REPOSITORY_URL = 'https://github.com/kami518/boltz2-iptm-binder-calibration'
RELEASE_BASE_URL = f'{REPOSITORY_URL}/releases/download/v1.1.0'
BLOB_BASE_URL = f'{REPOSITORY_URL}/blob/v1.1.0'
PAGE_W, PAGE_H = A4
LEFT = RIGHT = 1.75 * cm
TOP = 2.15 * cm
BOTTOM = 1.75 * cm
FRAME_W = PAGE_W - LEFT - RIGHT
NAVY = HexColor('#183B56')
TEXT = HexColor('#28343D')
GREY = HexColor('#5C6B76')
LIGHT = HexColor('#DCE5EB')
PENDING = HexColor('#8A4B08')
PENDING_BG = HexColor('#FFF6E5')

FONT_DIR = Path(matplotlib.get_data_path()) / 'fonts/ttf'
font_files = {
    'ReviewSerif': 'DejaVuSerif.ttf',
    'ReviewSerif-Bold': 'DejaVuSerif-Bold.ttf',
    'ReviewSerif-Italic': 'DejaVuSerif-Italic.ttf',
    'ReviewSerif-BoldItalic': 'DejaVuSerif-BoldItalic.ttf',
    'ReviewSans': 'DejaVuSans.ttf',
    'ReviewSans-Bold': 'DejaVuSans-Bold.ttf',
    'ReviewSans-Italic': 'DejaVuSans-Oblique.ttf',
    'ReviewMono': 'DejaVuSansMono.ttf',
}
for name, filename in font_files.items():
    path = FONT_DIR / filename
    if not path.is_file():
        raise FileNotFoundError(path)
    pdfmetrics.registerFont(TTFont(name, str(path)))
pdfmetrics.registerFontFamily('ReviewSerif', normal='ReviewSerif', bold='ReviewSerif-Bold', italic='ReviewSerif-Italic', boldItalic='ReviewSerif-BoldItalic')
pdfmetrics.registerFontFamily('ReviewSans', normal='ReviewSans', bold='ReviewSans-Bold', italic='ReviewSans-Italic', boldItalic='ReviewSans-Bold')

styles = getSampleStyleSheet()
BODY = ParagraphStyle('ReviewBody', fontName='ReviewSerif', fontSize=9.25, leading=12.9, textColor=TEXT, spaceAfter=6.5, splitLongWords=True, allowWidows=0, allowOrphans=0)
META = ParagraphStyle('ReviewMeta', parent=BODY, fontName='ReviewSans', fontSize=8.3, leading=11.5, textColor=GREY, spaceAfter=4)
PENDING_STYLE = ParagraphStyle('ReviewPending', parent=BODY, fontName='ReviewSans', fontSize=8.5, leading=12, textColor=PENDING, backColor=PENDING_BG, borderPadding=5, spaceBefore=3, spaceAfter=7)
H1 = ParagraphStyle('ReviewH1', fontName='ReviewSans-Bold', fontSize=17, leading=21, textColor=NAVY, spaceBefore=4, spaceAfter=12, keepWithNext=True)
H2 = ParagraphStyle('ReviewH2', fontName='ReviewSans-Bold', fontSize=12.5, leading=16, textColor=NAVY, spaceBefore=11, spaceAfter=7, keepWithNext=True)
H3 = ParagraphStyle('ReviewH3', fontName='ReviewSans-Bold', fontSize=10.4, leading=13.5, textColor=NAVY, spaceBefore=8, spaceAfter=5, keepWithNext=True)
CAPTION = ParagraphStyle('ReviewCaption', parent=BODY, fontName='ReviewSerif', fontSize=8.2, leading=11.5, textColor=TEXT, spaceBefore=5, spaceAfter=6)
SOURCE = ParagraphStyle('ReviewSource', parent=BODY, fontName='ReviewSans', fontSize=7.1, leading=10.1, textColor=GREY, spaceAfter=5)
QUOTE = ParagraphStyle('ReviewQuote', parent=BODY, fontName='ReviewSerif-Italic', fontSize=8.4, leading=12, leftIndent=13, rightIndent=8, textColor=GREY, borderColor=LIGHT, borderWidth=0, borderPadding=5)
BULLET = ParagraphStyle('ReviewBullet', parent=BODY, leftIndent=16, firstLineIndent=0, bulletIndent=3, spaceAfter=4)
REFERENCE = ParagraphStyle('ReviewReference', parent=BODY, fontSize=8.4, leading=12, leftIndent=20, firstLineIndent=-20, spaceAfter=7)
COVER_KICKER = ParagraphStyle('CoverKicker', fontName='ReviewSans-Bold', fontSize=10, leading=14, alignment=TA_CENTER, textColor=NAVY, spaceAfter=18)
COVER_TITLE = ParagraphStyle('CoverTitle', fontName='ReviewSans-Bold', fontSize=22, leading=29, alignment=TA_CENTER, textColor=NAVY, spaceAfter=18)
COVER_SUB = ParagraphStyle('CoverSub', fontName='ReviewSerif', fontSize=11, leading=16, alignment=TA_CENTER, textColor=TEXT, spaceAfter=8)
COVER_NOTE = ParagraphStyle('CoverNote', fontName='ReviewSans', fontSize=8.2, leading=12, alignment=TA_LEFT, textColor=GREY, spaceAfter=6)
TABLE_BODY = ParagraphStyle('ReviewTableBody', fontName='ReviewSans', fontSize=7.1, leading=9.1, textColor=TEXT, splitLongWords=True)
TABLE_HEAD = ParagraphStyle('ReviewTableHead', fontName='ReviewSans-Bold', fontSize=7.1, leading=9.1, textColor=NAVY, splitLongWords=True)


def plain_math(text):
    replacements = {
        r'\operatorname{logit}': 'logit', r'\operatorname{clip}': 'clip',
        r'\alpha': 'α', r'\beta': 'β', r'\sum': 'Σ',
        r'\hat p_i': 'p̂ᵢ', r'^{-1}': '⁻¹', '_i': 'ᵢ',
        r'\{': '{', r'\}': '}', r'\(': '', r'\)': '',
        r'\,': ' ', r'\;': ' ', '$': '',
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def code_text(value):
    escaped = html.escape(value, quote=False)
    return f'<font name="ReviewMono" size="7.1" color="#425563">{escaped}</font>'


CITATIONS = {
    'Abramson et al., 2024': 'ref-abramson-2024', 'Abramson et al. (2024)': 'ref-abramson-2024',
    'Adaptyv Bio, 2026': 'ref-adaptyv-2026', 'Adaptyv Bio (2026)': 'ref-adaptyv-2026',
    'Baek et al., 2021': 'ref-baek-2021', 'Baek et al. (2021)': 'ref-baek-2021',
    'Bennett et al., 2023': 'ref-bennett-2023', 'Bennett et al. (2023)': 'ref-bennett-2023',
    'Breiman, 2001': 'ref-breiman-2001', 'Breiman (2001)': 'ref-breiman-2001',
    'Brier, 1950': 'ref-brier-1950', 'Brier (1950)': 'ref-brier-1950',
    'Cao et al., 2020': 'ref-cao-2020', 'Cao et al. (2020)': 'ref-cao-2020',
    'Cao et al., 2022': 'ref-cao-2022', 'Cao et al. (2022)': 'ref-cao-2022',
    'Claude Science and Shanehsazzadeh, 2026': 'ref-claude-science-2026',
    'Claude Science and Shanehsazzadeh (2026)': 'ref-claude-science-2026',
    'Daily, 2016': 'ref-daily-2016', 'Daily (2016)': 'ref-daily-2016',
    'Dauparas et al., 2022': 'ref-dauparas-2022', 'Dauparas et al. (2022)': 'ref-dauparas-2022',
    'Dunbrack, 2025': 'ref-dunbrack-2025', 'Dunbrack (2025)': 'ref-dunbrack-2025',
    'Field and Welsh, 2007': 'ref-field-welsh-2007', 'Field and Welsh (2007)': 'ref-field-welsh-2007',
    'Garcia et al., 2026': 'ref-garcia-2026', 'Garcia et al. (2026)': 'ref-garcia-2026',
    'Huang et al., 2016': 'ref-huang-2016', 'Huang et al. (2016)': 'ref-huang-2016',
    'Jumper et al., 2021': 'ref-jumper-2021', 'Jumper et al. (2021)': 'ref-jumper-2021',
    'Kaufman et al., 2012': 'ref-kaufman-2012', 'Kaufman et al. (2012)': 'ref-kaufman-2012',
    'Lee et al., 2026': 'ref-lee-2026', 'Lee et al. (2026)': 'ref-lee-2026',
    'Li et al., 2026': 'ref-li-2026', 'Li et al. (2026)': 'ref-li-2026',
    'Lin et al., 2023': 'ref-lin-2023', 'Lin et al. (2023)': 'ref-lin-2023',
    'Pacesa et al., 2025': 'ref-pacesa-2025', 'Pacesa et al. (2025)': 'ref-pacesa-2025',
    'Sagadinov et al., 2026': 'ref-sagadinov-2026', 'Sagadinov et al. (2026)': 'ref-sagadinov-2026',
    'Watson et al., 2023': 'ref-watson-2023', 'Watson et al. (2023)': 'ref-watson-2023',
    'Yang et al., 2026': 'ref-yang-2026', 'Yang et al. (2026)': 'ref-yang-2026',
    'Badonyi and Toth-Petroczy, 2026': 'ref-badonyi-2026', 'Badonyi and Toth-Petroczy (2026)': 'ref-badonyi-2026',
    'Overath et al., 2025': 'ref-overath-2025', 'Overath et al. (2025)': 'ref-overath-2025',
    'Passaro et al., 2025': 'ref-passaro-2025', 'Passaro et al. (2025)': 'ref-passaro-2025',
    'Pedregosa et al., 2011': 'ref-pedregosa-2011', 'Pedregosa et al. (2011)': 'ref-pedregosa-2011',
    'Qian et al., 2026': 'ref-qian-2026', 'Qian et al. (2026)': 'ref-qian-2026',
    'Roberts et al., 2017': 'ref-roberts-2017', 'Roberts et al. (2017)': 'ref-roberts-2017',
    'Saito and Rehmsmeier, 2015': 'ref-saito-rehmsmeier-2015', 'Saito and Rehmsmeier (2015)': 'ref-saito-rehmsmeier-2015',
    'Smorodina et al., 2026': 'ref-smorodina-2026', 'Smorodina et al. (2026)': 'ref-smorodina-2026',
    'Steinegger and Söding, 2017': 'ref-steinegger-soding-2017', 'Steinegger and Söding (2017)': 'ref-steinegger-soding-2017',
    'Steyerberg and Harrell, 2016': 'ref-steyerberg-harrell-2016', 'Steyerberg and Harrell (2016)': 'ref-steyerberg-harrell-2016',
    'Steyerberg et al., 2010': 'ref-steyerberg-2010', 'Steyerberg et al. (2010)': 'ref-steyerberg-2010',
    'Van Calster et al., 2019': 'ref-vancalster-2019', 'Van Calster et al. (2019)': 'ref-vancalster-2019',
    'Van Calster et al., 2016': 'ref-vancalster-2016', 'Van Calster et al. (2016)': 'ref-vancalster-2016',
    'Varma and Simon, 2006': 'ref-varma-2006', 'Varma and Simon (2006)': 'ref-varma-2006',
    'Vickers and Elkin, 2006': 'ref-vickers-2006', 'Vickers and Elkin (2006)': 'ref-vickers-2006',
    'Zambaldi et al., 2024': 'ref-zambaldi-2024', 'Zambaldi et al. (2024)': 'ref-zambaldi-2024',
    'The protein binder design data release authors, 2026': 'ref-protein-binder-data-2026',
    'The protein binder design data release authors (2026)': 'ref-protein-binder-data-2026',
}
REF_STARTS = {
    'Abramson, J.': 'ref-abramson-2024', 'Adaptyv Bio (2026)': 'ref-adaptyv-2026',
    'Badonyi, M. and Toth-Petroczy, A. (2026)': 'ref-badonyi-2026',
    'Baek, M.': 'ref-baek-2021', 'Bennett, N.R.': 'ref-bennett-2023',
    'Breiman, L.': 'ref-breiman-2001', 'Brier, G.W.': 'ref-brier-1950',
    'Cao, L. et al. (2020)': 'ref-cao-2020', 'Cao, L. et al. (2022)': 'ref-cao-2022',
    'Claude Science and Shanehsazzadeh': 'ref-claude-science-2026',
    'Daily, J.': 'ref-daily-2016', 'Dauparas, J.': 'ref-dauparas-2022',
    'Dunbrack, R.L.': 'ref-dunbrack-2025', 'Field, C.A. and Welsh': 'ref-field-welsh-2007',
    'Garcia, M.': 'ref-garcia-2026', 'Huang, P.-S.': 'ref-huang-2016',
    'Jumper, J.': 'ref-jumper-2021', 'Kaufman, S.': 'ref-kaufman-2012',
    'Lee, G.': 'ref-lee-2026',
    'Li, W.': 'ref-li-2026', 'Lin, Z.': 'ref-lin-2023',
    'Overath, M.D.': 'ref-overath-2025', 'Pacesa, M.': 'ref-pacesa-2025',
    'Passaro, S.': 'ref-passaro-2025', 'Pedregosa, F.': 'ref-pedregosa-2011',
    'Qian, J.': 'ref-qian-2026',
    'Roberts, D.R.': 'ref-roberts-2017', 'Sagadinov, S.': 'ref-sagadinov-2026',
    'Saito, T. and Rehmsmeier': 'ref-saito-rehmsmeier-2015',
    'Smorodina, E.': 'ref-smorodina-2026', 'Steinegger, M. and Söding': 'ref-steinegger-soding-2017',
    'Steyerberg, E.W. and Harrell': 'ref-steyerberg-harrell-2016',
    'Steyerberg, E.W. et al. (2010)': 'ref-steyerberg-2010',
    'The protein binder design data release authors': 'ref-protein-binder-data-2026',
    'Van Calster, B. et al. (2016)': 'ref-vancalster-2016',
    'Van Calster, B. et al. (2019)': 'ref-vancalster-2019',
    'Varma, S. and Simon': 'ref-varma-2006',
    'Vickers, A.J. and Elkin': 'ref-vickers-2006',
    'Watson, J.L.': 'ref-watson-2023', 'Yang, W.': 'ref-yang-2026',
    'Zambaldi, V.': 'ref-zambaldi-2024',
}


def add_crossrefs(text):
    text = text.replace('Supplementary Figures 3–4', 'Supplementary Figures <a href="#figure-s3" color="#183B56">3</a>–<a href="#figure-s4" color="#183B56">4</a>')
    text = re.sub(r'Supplementary Tables S1–S16', 'Supplementary Tables <a href="#table-s1" color="#183B56">S1</a>–<a href="#table-s16" color="#183B56">S16</a>', text)
    text = re.sub(r'Supplementary Table S(\d+)', r'<a href="#table-s\1" color="#183B56">Supplementary Table S\1</a>', text)
    text = re.sub(r'Supplementary Figure S(\d+)', r'<a href="#figure-s\1" color="#183B56">Supplementary Figure S\1</a>', text)
    text = re.sub(r'(?<!Supplementary )Table ([1-4])\b', r'<a href="#table-\1" color="#183B56">Table \1</a>', text)
    text = re.sub(r'(?<!Supplementary )Fig\. ([1-7])\b', r'<a href="#figure-\1" color="#183B56">Fig. \1</a>', text)
    text = re.sub(r'(?<!Supplementary )Figure ([1-7])\b', r'<a href="#figure-\1" color="#183B56">Figure \1</a>', text)
    for phrase, key in sorted(CITATIONS.items(), key=lambda item: -len(item[0])):
        text = text.replace(phrase, f'<a href="#{key}" color="#28343D">{phrase}</a>')
    return text


def inline(text, source_path, crossrefs=True):
    text = plain_math(text.strip())
    # Preserve Markdown links as placeholders before escaping and inline styling.
    links = []
    def stash_link(match):
        label, target = match.group(1), html.unescape(match.group(2))
        token = f'@@LINK{len(links)}@@'
        links.append((token, label, target))
        return token
    text = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', stash_link, text)
    text = html.escape(text, quote=False)
    text = re.sub(r'`([^`]+)`', lambda m: code_text(html.unescape(m.group(1))), text)
    text = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', text)
    text = re.sub(r'(?<!\*)\*([^*\n]+)\*(?!\*)', r'<i>\1</i>', text)
    if crossrefs:
        text = add_crossrefs(text)
    # Plain web addresses in the bibliography remain clickable.
    def url_link(match):
        url = match.group(0)
        trailing = ''
        while url and url[-1] in '.,;':
            trailing = url[-1] + trailing
            url = url[:-1]
        return f'<a href="{url}" color="#183B56">{url}</a>{trailing}'
    text = re.sub(r'https?://[^\s<]+', url_link, text)
    for token, label, target in links:
        label_html = html.escape(label, quote=False)
        label_html = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', label_html)
        if target.startswith(('http://', 'https://')):
            href = target
        elif target.startswith('#'):
            href = target
        else:
            resolved = (source_path.parent / target).resolve()
            try:
                published = resolved.relative_to(ROOT).as_posix()
            except ValueError:
                href = ''
            else:
                base = RELEASE_BASE_URL if resolved.suffix.lower() in {'.xlsx', '.csv', '.pdf'} else BLOB_BASE_URL
                href = f'{base}/{published}'
        if not href:
            text = text.replace(token, label_html)
            continue
        text = text.replace(token, f'<a href="{html.escape(href, quote=True)}" color="#183B56">{label_html}</a>')
    return text


class ReviewDocTemplate(BaseDocTemplate):
    def afterFlowable(self, flowable):
        key = getattr(flowable, '_bookmark', None)
        if key:
            self.canv.bookmarkPage(key)
        if isinstance(flowable, Paragraph):
            level = getattr(flowable, '_toc_level', None)
            title = getattr(flowable, '_plain_text', '')
            if level is not None and key:
                self.canv.addOutlineEntry(title, key, level=level, closed=(level == 0))
                self.notify('TOCEntry', (level, title, self.page, key))


def header_footer(canvas, doc):
    canvas.saveState()
    canvas.setTitle('Context-dependent calibration and experimental decision value of Boltz-2 ipTM in binder design')
    canvas.setAuthor('KangMing Chen')
    canvas.setSubject('Main article and supplementary materials')
    if doc.page > 1:
        canvas.setStrokeColor(LIGHT)
        canvas.setLineWidth(0.6)
        canvas.line(LEFT, PAGE_H - 1.35 * cm, PAGE_W - RIGHT, PAGE_H - 1.35 * cm)
        canvas.setFont('ReviewSans', 7.3)
        canvas.setFillColor(GREY)
        canvas.drawString(LEFT, PAGE_H - 1.12 * cm, 'MAIN ARTICLE + SUPPLEMENTARY MATERIALS')
        canvas.drawRightString(PAGE_W - RIGHT, PAGE_H - 1.12 * cm, '21 November 2026')
        canvas.setStrokeColor(LIGHT)
        canvas.line(LEFT, 1.17 * cm, PAGE_W - RIGHT, 1.17 * cm)
        canvas.setFont('ReviewSans', 7.2)
        canvas.drawString(LEFT, 0.78 * cm, 'v1.1.0 submission candidate')
        canvas.drawRightString(PAGE_W - RIGHT, 0.78 * cm, f'Page {doc.page}')
    canvas.restoreState()


def heading(text, level, key, toc_level):
    style = {1: H1, 2: H2, 3: H3}[level]
    p = Paragraph(inline(text, MAIN, crossrefs=False), style)
    p._bookmark = key
    p._toc_level = toc_level
    p._plain_text = re.sub(r'[*`]', '', text)
    return p


def add_page_break(story):
    if not story or not isinstance(story[-1], PageBreak):
        story.append(PageBreak())


def image_flowable(path):
    with PILImage.open(path) as im:
        px_w, px_h = im.size
    max_w = FRAME_W
    max_h = 17.8 * cm
    scale = min(max_w / px_w, max_h / px_h)
    img = Image(str(path), width=px_w * scale, height=px_h * scale)
    img.hAlign = 'CENTER'
    return img


def table_widths(rows, ncols):
    if ncols == 2:
        return [FRAME_W * 0.24, FRAME_W * 0.76]
    if ncols == 3:
        return [FRAME_W * 0.19, FRAME_W * 0.405, FRAME_W * 0.405]
    if ncols == 8:
        fractions = [0.075, 0.07, 0.115, 0.115, 0.115, 0.105, 0.105, 0.30]
        return [FRAME_W * f for f in fractions]
    lengths = []
    for col in range(ncols):
        values = [len(re.sub(r'[*`]', '', row[col])) if col < len(row) else 0 for row in rows]
        lengths.append(max(7, min(32, max(values))))
    # Numeric columns stay compact while descriptive columns receive room.
    weights = [max(2.8, length ** 0.55) for length in lengths]
    total = sum(weights)
    return [FRAME_W * w / total for w in weights]


def make_table(rows, source_path):
    ncols = len(rows[0])
    assert all(len(row) == ncols for row in rows), rows
    font_size = 7.1 if ncols <= 5 else 6.2
    body_style = ParagraphStyle('TableBodyDynamic', parent=TABLE_BODY, fontSize=font_size, leading=font_size * 1.30)
    head_style = ParagraphStyle('TableHeadDynamic', parent=TABLE_HEAD, fontSize=font_size, leading=font_size * 1.30)
    data = []
    for ridx, row in enumerate(rows):
        style = head_style if ridx == 0 else body_style
        data.append([Paragraph(inline(cell, source_path), style) for cell in row])
    table = Table(data, colWidths=table_widths(rows, ncols), repeatRows=1, splitByRow=1, hAlign='LEFT')
    table.setStyle(TableStyle([
        ('LINEABOVE', (0, 0), (-1, 0), 1.1, NAVY),
        ('LINEBELOW', (0, 0), (-1, 0), 0.65, NAVY),
        ('LINEBELOW', (0, -1), (-1, -1), 1.1, NAVY),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING', (0, 0), (-1, -1), 3),
        ('RIGHTPADDING', (0, 0), (-1, -1), 3),
        ('TOPPADDING', (0, 0), (-1, -1), 4),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, HexColor('#F7F9FA')]),
    ]))
    return table


def is_pending(text):
    markers = ['To be supplied', 'TO BE ADDED', 'CITATION REQUIRED', 'to be verified', 'to be confirmed', 'VERIFIED ANTHROPIC']
    return any(marker in text for marker in markers)


def paragraph_for(text, source_path, current_major, after_image=False):
    if after_image:
        style = CAPTION
    elif text.startswith(('*Frozen panel sources:*', '*Numerical source:*', '*Accepted panel sources:*', '*Evidence boundary:*', '*Numerical or accepted panel source:*')):
        style = SOURCE
    elif is_pending(text):
        style = PENDING_STYLE
    elif current_major == 'References':
        style = REFERENCE
    elif text.startswith('**') and (':**' in text or text.startswith('**Panel')):
        style = META
    else:
        style = BODY
    p = Paragraph(inline(text, source_path), style)
    if current_major == 'References':
        clean = re.sub(r'[*`]', '', text)
        for prefix, key in REF_STARTS.items():
            if clean.startswith(prefix):
                p._bookmark = key
                break
    return p


def parse_markdown(story, path, source_kind, start_line=0):
    lines = path.read_text(encoding='utf-8').splitlines()
    i = start_line
    current_major = ''
    after_image = False
    first_major = True
    counters = {'heading': 0}
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if not stripped:
            i += 1
            continue
        heading_match = re.match(r'^(#{1,4})\s+(.+)$', stripped)
        if heading_match:
            marks, text = heading_match.groups()
            level = len(marks)
            if level == 1:
                i += 1
                continue
            if level == 2:
                current_major = re.sub(r'[*`]', '', text)
                if not first_major:
                    add_page_break(story)
                first_major = False
            if level == 3 and re.match(r'(?:Supplementary )?Figure (?:S)?\d+\.', text):
                add_page_break(story)
            key = None
            toc_level = None
            table_match = re.match(r'Table (S?\d+)\.', text)
            figure_match = re.match(r'(?:Supplementary )?Figure (S?\d+)\.', text)
            if table_match:
                key = 'table-' + table_match.group(1).lower()
            elif figure_match:
                key = 'figure-' + figure_match.group(1).lower()
            else:
                counters['heading'] += 1
                key = f'{source_kind}-heading-{counters["heading"]}'
            if level == 2:
                toc_level = 0
            elif level == 3:
                toc_level = 1
            story.append(heading(text, min(level - 1, 3), key, toc_level))
            after_image = False
            i += 1
            continue
        image_match = re.match(r'^!\[([^\]]*)\]\(([^)]+)\)$', stripped)
        if image_match:
            target = image_match.group(2)
            image_path = (path.parent / target).resolve()
            if not image_path.is_file():
                raise FileNotFoundError(image_path)
            story.append(image_flowable(image_path))
            after_image = True
            i += 1
            continue
        if stripped.startswith('|'):
            table_lines = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                table_lines.append(lines[i].strip())
                i += 1
            rows = []
            for idx, table_line in enumerate(table_lines):
                cells = [cell.strip() for cell in table_line.strip('|').split('|')]
                if idx == 1 and all(re.fullmatch(r':?-{3,}:?', cell.replace(' ', '')) for cell in cells):
                    continue
                rows.append(cells)
            story.append(make_table(rows, path))
            story.append(Spacer(1, 7))
            after_image = False
            continue
        bullet_match = re.match(r'^[-*]\s+(.+)$', stripped)
        numbered_match = re.match(r'^(\d+)\.\s+(.+)$', stripped)
        if bullet_match or numbered_match:
            if bullet_match:
                content = bullet_match.group(1)
                bullet = '•'
            else:
                content = numbered_match.group(2)
                bullet = numbered_match.group(1) + '.'
            story.append(Paragraph(inline(content, path), BULLET, bulletText=bullet))
            after_image = False
            i += 1
            continue
        if stripped.startswith('>'):
            content = stripped.lstrip('>').strip()
            story.append(Paragraph(inline(content, path), QUOTE))
            after_image = False
            i += 1
            continue
        # Gather a normal paragraph until a structural Markdown marker.
        paragraph_lines = [stripped]
        i += 1
        while i < len(lines) and lines[i].strip():
            nxt = lines[i].strip()
            if re.match(r'^(#{1,4})\s+', nxt) or nxt.startswith('|') or nxt.startswith('![') or re.match(r'^[-*]\s+', nxt) or re.match(r'^\d+\.\s+', nxt) or nxt.startswith('>'):
                break
            paragraph_lines.append(nxt)
            i += 1
        text = ' '.join(part[:-2] if part.endswith('  ') else part for part in paragraph_lines)
        story.append(paragraph_for(text, path, current_major, after_image=after_image))
        after_image = False
    return story


story = []
# Cover.
story += [
    Spacer(1, 2.4 * cm),
    Paragraph('ORIGINAL ARTICLE  |  VERSION 1.1.0', COVER_KICKER),
    Paragraph('Context-dependent calibration and experimental decision value of Boltz-2 ipTM in binder design', COVER_TITLE),
    HRFlowable(width='70%', thickness=1.2, color=NAVY, spaceBefore=4, spaceAfter=20, hAlign='CENTER'),
    Paragraph('Main article and supplementary materials', COVER_SUB),
    Paragraph('Complete submission candidate with continuous pagination, bookmarks, and linked contents', COVER_SUB),
    Spacer(1, 1.2 * cm),
]
cover_data = [
    ['Article type', 'Original Article'],
    ['Authors', 'KangMing Chen (ORCID 0009-0006-9634-5397)'],
    ['Prepared', '21 November 2026'],
    ['Main article', 'Abstract, main text, four main tables, seven main figures, references, and declarations'],
    ['Supplement in this PDF', 'Supplementary guide, Tables S1–S17 documentation, and Supplementary Figures S1–S11'],
    ['Required companion', 'Supplementary_Tables_S1_S16_v1.1.0.xlsx contains all 75 detailed Excel tables at source precision'],
]
cover_table = Table([[Paragraph(f'<b>{html.escape(a)}</b>', COVER_NOTE), Paragraph(html.escape(b), COVER_NOTE)] for a, b in cover_data], colWidths=[3.4 * cm, 11.8 * cm], hAlign='CENTER')
cover_table.setStyle(TableStyle([
    ('LINEABOVE', (0, 0), (-1, 0), 0.8, NAVY), ('LINEBELOW', (0, -1), (-1, -1), 0.8, NAVY),
    ('VALIGN', (0, 0), (-1, -1), 'TOP'), ('TOPPADDING', (0, 0), (-1, -1), 6), ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
]))
story += [
    cover_table,
    Spacer(1, 0.9 * cm),
    Paragraph('<b>Release status.</b> Source versions, dataset citations, post-review sensitivity analyses, vector figure sources, redistribution terms, author metadata, and declarations have been checked. The public tag, final archival record, DOI activation, and journal-specific file conversion remain subject to mainline release approval.', PENDING_STYLE),
    Spacer(1, 0.4 * cm),
    Paragraph('Source documents: publication/manuscript/main_article.md and publication/manuscript/supplementary_materials.md', COVER_NOTE),
    PageBreak(),
]
# TOC.
story.append(heading('Contents', 1, 'contents', None))
toc = TableOfContents()
toc.levelStyles = [
    ParagraphStyle('TOCLevel0', fontName='ReviewSans-Bold', fontSize=9, leading=13, leftIndent=0, firstLineIndent=0, textColor=NAVY, spaceBefore=5),
    ParagraphStyle('TOCLevel1', fontName='ReviewSerif', fontSize=8.2, leading=11.5, leftIndent=15, firstLineIndent=0, textColor=TEXT, spaceBefore=1),
]
toc.dotsMinLevel = 0
story += [toc, PageBreak()]
# Main article starts with its metadata and abstract.
story.append(heading('Main Article', 1, 'main-article', 0))
parse_markdown(story, MAIN, 'main')
# Supplement guide starts on a fresh page and retains its complete content.
add_page_break(story)
story.append(heading('Supplementary Materials: Guide and Figures', 1, 'supplement-guide', 0))
parse_markdown(story, SUPP, 'supplement', start_line=4)

frame = Frame(LEFT, BOTTOM, FRAME_W, PAGE_H - TOP - BOTTOM, id='review-frame', leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
doc = ReviewDocTemplate(str(OUTPUT), pagesize=A4, leftMargin=LEFT, rightMargin=RIGHT, topMargin=TOP, bottomMargin=BOTTOM, title='Context-dependent calibration and experimental decision value of Boltz-2 ipTM in binder design', author='KangMing Chen', subject='Post-review main article and supplementary materials')
doc.addPageTemplates([PageTemplate(id='review', frames=[frame], onPage=header_footer)])
doc.multiBuild(story)

from pypdf import PdfReader
reader = PdfReader(str(OUTPUT))
info = {
    'output': OUTPUT.relative_to(ROOT).as_posix(),
    'pages': len(reader.pages),
    'bytes': OUTPUT.stat().st_size,
    'main_source': MAIN.relative_to(ROOT).as_posix(),
    'supplement_source': SUPP.relative_to(ROOT).as_posix(),
    'scope': 'Complete post-review main article plus supplementary guide and all 18 figures; the original 75 detailed Excel tables in S1-S16 remain in the companion workbook, and S17 is documented in the guide.',
    'unresolved_placeholders_retained': False,
}
BUILD_INFO.write_text(json.dumps(info, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(info, ensure_ascii=False, indent=2))
