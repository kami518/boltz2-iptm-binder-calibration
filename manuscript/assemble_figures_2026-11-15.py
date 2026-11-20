from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
import csv
import hashlib
import textwrap

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
OUT = ROOT / 'manuscript' / 'figures'
OUT.mkdir(parents=True, exist_ok=True)
WIDTH = 3600
MARGIN = 60
CONTENT_WIDTH = WIDTH - 2 * MARGIN
PANEL_HEADER = 100
GAP = 45
DPI = (300, 300)
FONT_REGULAR = r'C:/Windows/Fonts/arial.ttf'
FONT_BOLD = r'C:/Windows/Fonts/arialbd.ttf'


def font(size, bold=False):
    return ImageFont.truetype(FONT_BOLD if bold else FONT_REGULAR, size)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fitted_text(draw, text, max_width, start_size=48, bold=False):
    for size in range(start_size, 25, -2):
        candidate = font(size, bold)
        if draw.textbbox((0, 0), text, font=candidate)[2] <= max_width:
            return candidate
    return font(26, bold)


def draw_centered(draw, box, text, size=44, bold=False, fill='#162b3a', spacing=12):
    x0, y0, x1, y1 = box
    max_width = x1 - x0 - 60
    lines = []
    for paragraph in text.split('\n'):
        if not paragraph:
            lines.append('')
            continue
        words = paragraph.split()
        line = ''
        for word in words:
            trial = word if not line else line + ' ' + word
            f = fitted_text(draw, trial, max_width, size, bold)
            if draw.textbbox((0, 0), trial, font=f)[2] <= max_width:
                line = trial
            else:
                lines.append(line)
                line = word
        if line:
            lines.append(line)
    f = font(size, bold)
    while True:
        bbox = draw.multiline_textbbox((0, 0), '\n'.join(lines), font=f, spacing=spacing, align='center')
        if bbox[2] - bbox[0] <= max_width and bbox[3] - bbox[1] <= y1 - y0 - 30:
            break
        size -= 2
        f = font(size, bold)
    bbox = draw.multiline_textbbox((0, 0), '\n'.join(lines), font=f, spacing=spacing, align='center')
    x = (x0 + x1 - (bbox[2] - bbox[0])) / 2
    y = (y0 + y1 - (bbox[3] - bbox[1])) / 2 - bbox[1]
    draw.multiline_text((x, y), '\n'.join(lines), font=f, fill=fill, spacing=spacing, align='center')


def box(draw, coords, text, fill, outline, size=44, bold=False):
    draw.rounded_rectangle(coords, radius=30, fill=fill, outline=outline, width=4)
    draw_centered(draw, coords, text, size=size, bold=bold)


def arrow(draw, start, end, color='#527284', width=10):
    draw.line([start, end], fill=color, width=width)
    x, y = end
    if abs(end[1] - start[1]) >= abs(end[0] - start[0]):
        sign = 1 if end[1] > start[1] else -1
        draw.polygon([(x, y), (x - 20, y - sign * 38), (x + 20, y - sign * 38)], fill=color)
    else:
        sign = 1 if end[0] > start[0] else -1
        draw.polygon([(x, y), (x - sign * 38, y - 20), (x - sign * 38, y + 20)], fill=color)


def make_figure1():
    image = Image.new('RGB', (WIDTH, 2400), 'white')
    draw = ImageDraw.Draw(image)
    draw.text((MARGIN, 45), 'Study design and evaluation layers', font=font(70, True), fill='#102f40')
    draw.text((MARGIN, 140), 'Frozen retrospective analyses; evidence layers are not interchangeable', font=font(38), fill='#527284')

    left = (100, 250, 1710, 2020)
    right = (1890, 250, 3500, 2020)
    draw.rounded_rectangle(left, radius=35, fill='#f5f9fc', outline='#a8c5d5', width=4)
    draw.rounded_rectangle(right, radius=35, fill='#f7fbf8', outline='#a8cbb5', width=4)
    draw.text((160, 290), 'MAIN COHORT', font=font(48, True), fill='#21617b')
    draw.text((1950, 290), 'EXTERNAL EVALUATION', font=font(48, True), fill='#2d6d50')

    box(draw, (180, 390, 1630, 700),
        'Proteinbase frozen snapshot\n4,280 designs | 4,943 design–target pairs\n2,893 pairs with released Boltz-2 ipTM',
        '#ffffff', '#74a8bf', size=43)
    arrow(draw, (905, 700), (905, 790))
    box(draw, (180, 790, 1630, 1080),
        'Primary analysis cohort\n2,084 pairs | 6 targets\n748 Weak-inclusive positives | 678 Medium/Strong positives',
        '#e2f0f6', '#4f8da7', size=43, bold=True)
    arrow(draw, (905, 1080), (905, 1170))
    box(draw, (180, 1170, 1630, 1930),
        'Evaluation layers\n• Within-target discrimination and raw probability error\n• Target-blocked zero-label U\n• Out-of-fold local intercept/slope adaptation (I/F)\n• Fixed-budget and shared-budget decision analyses\n• Prespecified exploratory Pinc increment',
        '#ffffff', '#74a8bf', size=42)

    box(draw, (1970, 390, 3420, 700),
        'External release\n1,440 designs | 16 targets\n1,320 with interpretable integrated outcomes',
        '#ffffff', '#76a98b', size=43)
    arrow(draw, (2695, 700), (2695, 790), color='#527c64')
    box(draw, (1970, 790, 3420, 1130),
        'Confirmatory zero-label U\n720 pairs | 8 targets\n4 new targets | 4 same-target new campaigns\nDifferent five-seed score-generation protocol',
        '#e3f1e8', '#5d9373', size=41, bold=True)
    arrow(draw, (2695, 1130), (2695, 1210), color='#527c64')
    box(draw, (1970, 1210, 3420, 1490),
        'Descriptive low-event stress test\n521 pairs | 6 targets\nNo refitting or mapping updates',
        '#ffffff', '#76a98b', size=41)
    box(draw, (1970, 1550, 3420, 1930),
        'Post hoc / exploratory extensions\n• Out-of-fold external I/F (label-dependent)\n• Paired description of 10 predictor pipelines\n• No cross-predictor U/I/F transport analysis',
        '#ffffff', '#76a98b', size=40)

    draw.rounded_rectangle((100, 2110, 3500, 2320), radius=30, fill='#fff4e5', outline='#d9a85c', width=4)
    draw_centered(draw, (100, 2110, 3500, 2320),
                  'Evidence kept separate: ranking  |  probability calibration  |  zero-label transport  |  experimental-selection value',
                  size=43, bold=True, fill='#6c471b')
    path = OUT / 'Figure_1_study_design.png'
    image.save(path, dpi=DPI, optimize=True)
    return path


PANELS = {
    2: [
        ('A', 'ipTM distributions by experimental outcome and target', 'analysis/rq12_outputs/RQ1_ipTM_distributions.png'),
        ('B', 'Within-target discrimination: AUC (left) and AP (right)', 'analysis/rq12_outputs/RQ1_target_discrimination.png'),
    ],
    3: [
        ('A', 'Raw-ipTM reliability curves', 'analysis/rq12_outputs/RQ2_raw_reliability.png'),
        ('B', 'Target-specific intercepts and slopes before and after partial pooling', 'analysis/rq12_outputs/RQ2_partial_pooling.png'),
    ],
    4: [
        ('A', 'Main-cohort leave-one-target-out zero-label transport', 'analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/RQ3_U_LOTO.png'),
        ('B', 'Main-cohort predictive out-of-fold local adaptation', 'analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/RQ3_adaptation.png'),
        ('C', 'Frozen external zero-label performance by target and transport stratum', 'analysis/rq3e_outputs/run_20261114_164950/RQ3E_U_external_target_metrics.png'),
        ('D', 'Frozen external zero-label reliability curves', 'analysis/rq3e_outputs/run_20261114_164950/RQ3E_U_external_reliability.png'),
    ],
    5: [
        ('A', 'Fixed-proportion within-target selection', 'analysis/rq4_outputs/run_20261114_231925/RQ4A_fixed_budget.png'),
        ('B', 'Raw-ipTM operating points', 'analysis/rq4_outputs/run_20261114_231925/RQ4B_raw_thresholds.png'),
        ('C', 'Post hoc exploratory retrospective shared-budget simulation', 'analysis/rq4_outputs/run_20261114_231925/RQ4C_shared_budget.png'),
    ],
    6: [
        ('A', 'Same-row raw-score discrimination for ipTM and Pinc', 'analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720/RQ5_score_discrimination.png'),
        ('B', 'Target-blocked zero-label transport on common complete cases', 'analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720/RQ5_U_transport.png'),
        ('C', 'Local-adaptation increments on common complete cases', 'analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720/RQ5_adaptation_increment.png'),
    ],
    7: [
        ('A', 'Paired score shifts and within-target Spearman agreement', 'analysis/audit_outputs/RQ6_acceptance_run_20261115_151147/RQ6_figure17_score_shift.png'),
        ('B', 'Target-specific AUC and raw Brier under Weak+ and Medium+ outcomes', 'analysis/audit_outputs/RQ6_acceptance_run_20261115_151147/RQ6_figure18_discrimination_error.png'),
    ],
}

FIGURE_NAMES = {
    2: 'Figure_2_ipTM_discrimination.png',
    3: 'Figure_3_raw_calibration.png',
    4: 'Figure_4_transport_adaptation.png',
    5: 'Figure_5_experimental_selection.png',
    6: 'Figure_6_Pinc_increment.png',
    7: 'Figure_7_predictor_pipelines.png',
}


def compose(number, panels):
    loaded = []
    for label, title, rel in panels:
        source = Image.open(ROOT / rel).convert('RGB')
        height = round(source.height * CONTENT_WIDTH / source.width)
        source = source.resize((CONTENT_WIDTH, height), Image.Resampling.LANCZOS)
        loaded.append((label, title, rel, source))
    total_height = MARGIN + sum(PANEL_HEADER + im.height for _, _, _, im in loaded) + GAP * (len(loaded) - 1) + MARGIN
    canvas = Image.new('RGB', (WIDTH, total_height), 'white')
    draw = ImageDraw.Draw(canvas)
    y = MARGIN
    for label, title, rel, source in loaded:
        draw.rounded_rectangle((MARGIN, y, WIDTH - MARGIN, y + PANEL_HEADER - 12), radius=18, fill='#edf3f6')
        draw.text((MARGIN + 24, y + 13), label, font=font(58, True), fill='#174f68')
        title_font = fitted_text(draw, title, CONTENT_WIDTH - 150, start_size=43, bold=True)
        draw.text((MARGIN + 115, y + 24), title, font=title_font, fill='#183b4c')
        y += PANEL_HEADER
        canvas.paste(source, (MARGIN, y))
        y += source.height + GAP
    path = OUT / FIGURE_NAMES[number]
    canvas.save(path, dpi=DPI, optimize=True)
    return path


outputs = {1: make_figure1()}
for number, panels in PANELS.items():
    outputs[number] = compose(number, panels)

manifest_rows = []
for number, path in outputs.items():
    with Image.open(path) as im:
        dimensions = f'{im.width}×{im.height}'
        dpi = im.info.get('dpi', DPI)
    manifest_rows.append({
        'figure': f'Figure {number}',
        'file': str(path.relative_to(ROOT)).replace('\\', '/'),
        'dimensions_px': dimensions,
        'dpi': f'{dpi[0]:.0f}',
        'bytes': path.stat().st_size,
        'sha256': sha256(path),
    })
with (OUT / 'figure_file_manifest.csv').open('w', encoding='utf-8-sig', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]))
    writer.writeheader()
    writer.writerows(manifest_rows)

lines = [
    '# Figure 1–7 来源与拼版清单',
    '',
    '生成日期：2026-11-15。拼版仅包括统一分面标题、等比例缩放与上下排列；没有修改统计数值、裁剪数据区域或重新拟合模型。原始图内标题和内部标记保留，以便追溯冻结输出。',
    '',
    '## 成品文件',
    '',
    '| 图 | 文件 | 像素 | DPI | SHA-256 |',
    '|---|---|---:|---:|---|',
]
for row in manifest_rows:
    lines.append(f"| {row['figure']} | `{row['file']}` | {row['dimensions_px']} | {row['dpi']} | `{row['sha256']}` |")
lines.extend(['', '## 分面来源', ''])
lines.append('### Figure 1. Study design and cohort flow')
lines.append('')
lines.append('- 新绘流程示意图；数字与队列定义来自 `manuscript/PEDS_Research_Article_English_draft_v1_2026-11-15.md` 的 Table 1 和 `analysis/paper_results_RQ1_RQ6_draft_2026-11-15.md`。')
lines.append('- 示意图不改变外部队列、低事件压力测试或 post hoc I/F 的证据等级。')
for number, panels in PANELS.items():
    lines.extend(['', f'### Figure {number}', ''])
    for label, title, rel in panels:
        source_path = ROOT / rel
        lines.append(f'- **{label} — {title}**')
        lines.append(f'  - 来源：`{rel}`')
        lines.append(f'  - SHA-256：`{sha256(source_path)}`')
        lines.append('  - 操作：完整保留源图，等比例缩放，无数据裁剪。')
lines.extend([
    '',
    '## 使用边界',
    '',
    '- Figure 4 的外部 I/F 不进入主图；主图仅包含确认性的 frozen external U，post hoc external I/F 仍留给补充材料。',
    '- Figure 5C 是 post hoc exploratory retrospective shared-budget simulation，不能解释为前瞻实验节省。',
    '- Figure 6 是相同 complete cases 上的探索性 Pinc 增量分析，不证明 Pinc 普遍优于或劣于 ipTM。',
    '- Figure 7 是十个完整预测器管线的配对描述，不是主队列跨模型验证，也不是跨预测器 U/I/F 迁移分析。',
    '- 原始低事件旧图未用于 Figure 1–7。',
])
(OUT / 'figure_source_manifest_2026-11-15.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')

for number, path in outputs.items():
    with Image.open(path) as im:
        print(f'Figure {number}: {path} | {im.width}x{im.height} | {path.stat().st_size} bytes')
print(f'Manifest: {OUT / "figure_source_manifest_2026-11-15.md"}')
