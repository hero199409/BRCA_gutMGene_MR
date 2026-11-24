"""从本地分章、图和 CSV 生成 Frontiers 投稿包。用 --help 查看参数。"""

import argparse
import csv
from decimal import Decimal
from html.parser import HTMLParser
from io import BytesIO
import hashlib
import os
import json
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import tempfile
from xml.etree import ElementTree as ET
import zipfile

from PIL import Image
import pymupdf
import xlsxwriter

ROOT = Path(os.environ.get("BRCA_PROJECT_ROOT", Path(__file__).resolve().parent.parent))
FORMAT = Path(__file__).resolve().parent / 'frontiers_format'
CHAPTERS = (
    '01_Title_Abstract_Intro_Methods.md',
    '02_Results.md',
    '03_Discussion_Declarations_Refs.md',
    '04_Figure_Legends_and_Supp_Inventory.md',
)
MARKDOWN = 'markdown+superscript+subscript+pipe_tables+grid_tables+autolink_bare_uris-raw_html-raw_tex-raw_attribute-yaml_metadata_block-implicit_figures'
INVENTORY = re.compile(r'^# 8\. Supplementary Tables Inventory[^\n]*$', re.M)
FIGURE = re.compile(r'^\*\*Figure (S?[1-9][0-9]*)\.', re.M)
FIGURE_FILE = re.compile(r'^Fig(S?[1-9][0-9]*)_.+\.png$')
TABLE = re.compile(r'^\|\s*\*\*(Table S[1-9][0-9]*[a-z]?)\*\*\s*\|\s*`([^`/\\|]+\.csv)`\s*\|\s*([^|]+?)\s*\|\s*(.*?)\s*\|\s*$', re.M)
NUMBER = re.compile(r'[+-]?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?\Z|[+-]?\.[0-9]+(?:[eE][+-]?[0-9]+)?\Z')
W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
R = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
REL = 'http://schemas.openxmlformats.org/package/2006/relationships'
CT = 'http://schemas.openxmlformats.org/package/2006/content-types'
ET.register_namespace('w', W)
ET.register_namespace('r', R)


def read(path):
    return path.read_text(encoding='utf-8-sig')


def sha256(path):
    with path.open('rb') as source:
        digest = hashlib.sha256()
        for chunk in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def pandoc(*args, text=None, cwd=None):
    result = subprocess.run(['pandoc', *map(str, args)], input=text, cwd=cwd,
                            encoding='utf-8', capture_output=True)
    if result.returncode or result.stderr.strip():
        raise ValueError(f'Pandoc {' '.join(map(str, args))}: {result.stderr.strip()}')
    return result.stdout


def polish(text):
    # 只在输出副本中修正早期草稿的上标与比较符号，保留代码跨度。
    text = re.sub(r'\$\^\{([^}]+)\}\$',
                  lambda match: '^' + match[1].replace(r'\dagger', '†').replace('*', r'\*') + '^', text)
    spans = re.split(r'(`[^`\n]*`)', text)
    for index in range(0, len(spans), 2):
        span = spans[index]
        span = re.sub(r'10\^(-\d+)(?![\d^])', r'10^\1^', span)
        span = re.sub(r'(\*[Ir]\*|\))\^2(?![\d^])', r'\1^2^', span)
        span = re.sub(r'2\^log2\(RSEM\+1\)(?!\^)', r'2^log2(RSEM+1)^', span)
        span = span.replace('<=', '≤').replace('>=', '≥')
        spans[index] = span.replace(' -> ', ' → ')
    return ''.join(spans)


def walk(node):
    if isinstance(node, dict):
        if 't' in node:
            yield node
        for child in node.values():
            yield from walk(child)
    elif isinstance(node, list):
        for child in node:
            yield from walk(child)


def document(markdown, images=()):
    parsed = json.loads(pandoc('--from', MARKDOWN, '--to', 'json', text=polish(markdown)))
    actual_images = [node['c'][2][0] for node in walk(parsed) if node['t'] == 'Image']
    if actual_images != list(images):
        raise ValueError('Unexpected image reference in Markdown. Only explicitly staged local images are allowed.')
    if any(node['t'] in {'RawBlock', 'RawInline', 'Math'} for node in walk(parsed)):
        raise ValueError('Unrendered raw/TeX content in Markdown. Convert math to supported Unicode or Pandoc inlines.')
    return json.dumps(parsed, ensure_ascii=False)


def w(name):
    return f'{{{W}}}{name}'


def xml_bytes(element):
    return ET.tostring(element, encoding='utf-8', xml_declaration=True)


def configure_docx(path, numbered):
    with zipfile.ZipFile(path) as source:
        parts = {name: source.read(name) for name in source.namelist()}
    doc = ET.fromstring(parts['word/document.xml'])
    sections = list(doc.iter(w('sectPr')))
    if len(sections) != 1:
        raise ValueError(f'Expected a single Word section in {path.name}, found {len(sections)}')
    section = sections[0]
    for child in list(section):
        if child.tag in {w('pgSz'), w('pgMar')}:
            section.remove(child)
    geometry_position = next((i for i, child in enumerate(section)
                              if child.tag not in {w('headerReference'), w('footerReference'),
                                                   w('footnotePr'), w('endnotePr'), w('type')}), len(section))
    section.insert(geometry_position, ET.Element(w('pgSz'), {w('w'): '11906', w('h'): '16838'}))
    section.insert(geometry_position + 1, ET.Element(w('pgMar'), {
        w('left'): '1224', w('right'): '1224', w('top'): '1260', w('bottom'): '1260',
        w('header'): '450', w('footer'): '450', w('gutter'): '0',
    }))
    if numbered:
        footer_path = 'word/frontiers-footer.xml'
        rel_id = 'rIdFrontiersPageFooter'
        relationships = ET.fromstring(parts['word/_rels/document.xml.rels'])
        if any(rel.get('Id') == rel_id for rel in relationships) or footer_path in parts:
            raise ValueError('Unexpected collision with the page-number footer.')
        ET.SubElement(relationships, f'{{{REL}}}Relationship', {
            'Id': rel_id, 'Type': f'{R}/footer', 'Target': Path(footer_path).name,
        })
        types = ET.fromstring(parts['[Content_Types].xml'])
        ET.SubElement(types, f'{{{CT}}}Override', {
            'PartName': '/' + footer_path,
            'ContentType': 'application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml',
        })
        section.insert(0, ET.Element(w('footerReference'), {w('type'): 'default', f'{{{R}}}id': rel_id}))
        # lnNumType 在页码/分栏设置之前，遵循 OOXML CT_SectPr 顺序。
        after_numbering = {'pgNumType', 'cols', 'formProt', 'vAlign', 'noEndnote', 'titlePg',
                           'textDirection', 'bidi', 'rtlGutter', 'docGrid', 'printerSettings', 'sectPrChange'}
        position = next((i for i, child in enumerate(section)
                         if child.tag.split('}')[-1] in after_numbering), len(section))
        section.insert(position, ET.Element(w('lnNumType'), {
            w('countBy'): '1', w('start'): '1', w('distance'): '180', w('restart'): 'continuous',
        }))
        footer = ET.Element(w('ftr'))
        paragraph = ET.SubElement(footer, w('p'))
        ET.SubElement(ET.SubElement(paragraph, w('pPr')), w('jc'), {w('val'): 'center'})
        field = ET.SubElement(paragraph, w('fldSimple'), {w('instr'): ' PAGE '})
        ET.SubElement(ET.SubElement(field, w('r')), w('t')).text = '1'
        parts[footer_path] = xml_bytes(footer)
        parts['word/_rels/document.xml.rels'] = xml_bytes(relationships)
        parts['[Content_Types].xml'] = xml_bytes(types)
    parts['word/document.xml'] = xml_bytes(doc)
    with zipfile.ZipFile(path, mode='w', compression=zipfile.ZIP_DEFLATED) as dest:
        for name, data in parts.items():
            dest.writestr(name, data)


def make_docx(ast, output, work, numbered=False):
    options = ['--from=json', '--to=docx', '--output', output]
    pandoc(*options, text=ast, cwd=work)
    configure_docx(output, numbered)


class MainHTML(HTMLParser):
    """收集页面主体的可见文本，用于检查 PDF 排版覆盖情况。"""
    def __init__(self):
        super().__init__()
        self.text = []

    def handle_data(self, data):
        self.text.append(data)


def pdf_tokens(text):
    # 排版可以在换行处拆分连字符和下标，逐页核验另在输出后进行。
    return re.findall(r'\w+', text.casefold())


def make_pdf(ast, output):
    html = pandoc('--from=json', '--to=html5', text=ast)
    story = pymupdf.Story(html=html, user_css=read(FORMAT / 'manuscript.css'))
    a4 = pymupdf.paper_rect('a4')
    body = a4 + (68, 54, -54, -54)
    raw = BytesIO()
    writer = pymupdf.DocumentWriter(raw)
    try:
        for _ in range(500):
            more, filled_bounds = story.place(body)
            filled = pymupdf.Rect(filled_bounds)
            if filled.is_empty:
                raise ValueError(f'PDF layout stuck or overflowed: {filled}')
            device = writer.begin_page(a4)
            story.draw(device)
            writer.end_page()
            if not more:
                break
        else:
            raise ValueError('PDF layout did not complete within 500 pages.')
    finally:
        writer.close()

    line_number = 0
    pages = []
    source = MainHTML()
    source.feed(html)
    expected_tokens = set(pdf_tokens(' '.join(source.text)))
    actual_tokens = set()
    with pymupdf.open(stream=raw.getvalue(), filetype='pdf') as pdf:
        pdf.set_metadata({'title': 'Frontiers in Immunology manuscript'})
        for page in pdf:
            data = page.get_text('dict')
            actual_tokens.update(pdf_tokens(page.get_text()))
            baselines = []
            for block in data['blocks']:
                for line in block.get('lines', []):
                    spans = [span for span in line['spans'] if span['text'].strip()]
                    if not spans:
                        continue
                    x0, y0, x1, y1 = line['bbox']
                    if x0 < body.x0 - 1 or x1 > body.x1 + 2 or y1 > a4.height - 41 or y0 < body.y0 - 6:
                        raise ValueError(f'Page {page.number + 1} text outside PDF body: {line["bbox"]}')
                    nominal = max(span['size'] for span in spans)
                    baselines.append(statistics.median(span['origin'][1] for span in spans
                                                       if span['size'] >= nominal * 0.9))
            merged = []
            for baseline in sorted(baselines):
                if not merged or baseline - merged[-1] > 2.8:
                    merged.append(baseline)
            if not merged:
                raise ValueError(f'Page {page.number + 1} has no lines.')
            previous = line_number
            for baseline in merged:
                line_number += 1
                label = str(line_number)
                width = pymupdf.get_text_length(label, fontname='helv', fontsize=6)
                page.insert_text((body.x0 - 11 - width, baseline), label, fontsize=6, color=(0.25,) * 3)
            label = str(page.number + 1)
            width = pymupdf.get_text_length(label, fontname='helv', fontsize=9)
            page.insert_text(((a4.width - width) / 2, a4.height - 25), label, fontsize=9)
            pages.append({'page': page.number + 1, 'first_line': previous + 1, 'last_line': line_number})
        missing = expected_tokens - actual_tokens
        if len(missing) / max(len(expected_tokens), 1) > 0.04:
            raise ValueError(f'PDF text coverage low ({len(missing)} unmatched tokens): {sorted(missing)[:36]}')
        pdf.save(output, garbage=4, deflate=True)

    with pymupdf.open(output) as pdf:
        for page, record in zip(pdf, pages):
            margin = pymupdf.Rect(0, 0, body.x0 - 4, a4.height)
            footer = pymupdf.Rect(0, a4.height - 40, a4.width, a4.height)
            ink = page.get_text(clip=margin).split()
            expected = list(map(str, range(record['first_line'], record['last_line'] + 1)))
            if ink != expected or page.get_text(clip=footer).strip() != str(record['page']):
                raise ValueError(f'PDF line/page label verification failed on page {record["page"]}.')
    return {'pages': len(pages), 'lines': line_number, 'page_line_ranges': pages}


def source_images(directory, legends):
    by_id = {}
    for path in directory.glob('*.png'):
        match = FIGURE_FILE.fullmatch(path.name)
        if match:
            if match[1] in by_id:
                raise ValueError(f'Multiple input PNGs match Figure {match[1]} in {directory}')
            by_id[match[1]] = path
    figure_ids = FIGURE.findall(legends)
    if not figure_ids or len(set(figure_ids)) != len(figure_ids) or set(figure_ids) != set(by_id):
        raise ValueError(f'Legend/PNG ID mismatch: legends={figure_ids}, images={sorted(by_id)}')
    return {key: by_id[key] for key in figure_ids}


def make_tiff(source, target):
    with Image.open(source) as image:
        rgba = image.convert('RGBA')
        rgb = Image.new('RGB', image.size, 'white')
        rgb.paste(rgba, mask=rgba.getchannel('A'))
        width = round(180 * 350 / 25.4)
        height = round(rgb.height * width / rgb.width)
        if rgb.size != (width, height):
            rgb = rgb.resize((width, height), Image.Resampling.LANCZOS)
        rgb.save(target, format='TIFF', dpi=(350, 350), compression='tiff_lzw')
    with Image.open(target) as image:
        image.load()
        if (image.format, image.mode, image.info.get('dpi'), image.tag_v2[259], image.tag_v2[296]) != (
            'TIFF', 'RGB', (350.0, 350.0), 5, 2
        ) or image.size != rgb.size:
            raise ValueError(f'TIFF format verification failed for {target.name}')
        return {'source': source.name, 'source_sha256': sha256(source), 'file': target.name,
                'width_px': image.width, 'height_px': image.height, 'dpi': 350,
                'mode': image.mode, 'compression': 'LZW'}


def supplement_markdown(legends, inventory, images, work):
    start = re.search(r'^## Supplementary Figures[ \t]*$', legends, re.M)
    if not start:
        raise ValueError('Cannot find the Supplementary Figures heading.')
    block = legends[start.end():]
    positions = list(FIGURE.finditer(block))
    references = []
    parts = ['# Supplementary Material', '## Supplementary Figures']
    for index, match in enumerate(positions):
        key = match[1]
        if not key.startswith('S'):
            raise ValueError(f'Unexpected main figure in the supplementary legends: {key}')
        end = positions[index + 1].start() if index + 1 < len(positions) else len(block)
        legend = block[match.start():end].strip().removesuffix('---').rstrip()
        name = f'{key}.png'
        shutil.copy2(images[key], work / name)
        references.append(name)
        # 图片不另生成 Pandoc 图注，采用分章中完整的原图注。
        parts += [f'![Supplementary Figure {key}]({name}){{width=6.4in}}', legend]
    return '\n\n'.join(parts + [inventory]), references


def excel_value(value):
    if not NUMBER.fullmatch(value):
        return value
    decimal = Decimal(value)
    if (decimal and not -307 <= decimal.adjusted() <= 307) or len(decimal.normalize().as_tuple().digits) > 15:
        return value
    return float(decimal)


def read_tables(inventory, directory):
    indexed = TABLE.findall(inventory)
    disk = {path.name for path in directory.iterdir() if path.is_file() and path.suffix.lower() == '.csv'}
    names = [name for _, name, _, _ in indexed]
    ids = [identifier.casefold() for identifier, *_ in indexed]
    if not indexed or set(names) != disk or len(names) != len(set(names)) or len(ids) != len(set(ids)):
        raise ValueError(f'Inventory/CSV mapping mismatch. Unlisted={sorted(disk - set(names))}, missing={sorted(set(names) - disk)}.')
    tables = []
    for identifier, filename, dimensions, description in indexed:
        path = directory / filename
        digest = sha256(path)
        with path.open(encoding='utf-8-sig', newline='') as stream:
            rows = list(csv.reader(stream, strict=True))
        if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
            raise ValueError(f'Empty or non-rectangular CSV: {filename}')
        actual = f'{len(rows) - 1} × {len(rows[0])}'
        tables.append((identifier, filename, dimensions, actual, description, rows, digest))
    return tables


def write_workbook(tables, output):
    records = []
    with xlsxwriter.Workbook(str(output), {'strings_to_formulas': False, 'strings_to_urls': False}) as workbook:
        heading = workbook.add_format({'bold': True, 'bg_color': '#1f4e78', 'font_color': '#ffffff'})
        wrapped = workbook.add_format({'text_wrap': True, 'valign': 'top'})
        index = workbook.add_worksheet('00_Table_Index')
        index.write_row(0, 0, ['Table ID', 'Worksheet', 'CSV file', 'Data rows', 'Columns', 'Content and scope'], heading)
        index.freeze_panes(1, 0)
        index.set_column(0, 1, 15)
        index.set_column(2, 2, 61)
        index.set_column(3, 4, 11)
        index.set_column(5, 5, 94)
        for index_row, (identifier, filename, stated, actual, description, rows, digest) in enumerate(tables, start=1):
            tab = identifier.removeprefix('Table ')
            sheet = workbook.add_worksheet(tab)
            if len(rows) > 1048576 or len(rows[0]) > 16384:
                raise ValueError(f'CSV exceeds Excel dimensions: {filename}')
            sheet.write_row(0, 0, rows[0], heading)
            sheet.freeze_panes(1, 0)
            sheet.autofilter(0, 0, len(rows) - 1, len(rows[0]) - 1)
            for col, header in enumerate(rows[0]):
                sheet.set_column(col, col, min(42, max(11, len(header) + 2)))
            for row_no, row in enumerate(rows[1:], start=1):
                for col, value in enumerate(row):
                    if value and sheet.write(row_no, col, excel_value(value)) != 0:
                        raise ValueError(f'Excel cannot store {filename} row {row_no + 1}, col {col + 1}.')
            index.write(index_row, 0, identifier)
            index.write_url(index_row, 1, f"internal:'{tab}'!A1", string=tab)
            index.write_row(index_row, 2, [filename, len(rows) - 1, len(rows[0])])
            index.write(index_row, 5, description, wrapped)
            records.append({'table': identifier, 'sheet': tab, 'csv': filename, 'sha256': digest,
                            'rows': len(rows) - 1, 'columns': len(rows[0]),
                            'source_inventory_dimensions': stated, 'actual_dimensions': actual})
    return records


def build(args, destination, work):
    ms = args.root / '04_manuscript'
    sources = [ms / name for name in CHAPTERS]
    input_hashes = {path: sha256(path) for path in (*sources, args.scope, FORMAT / 'manuscript.css')}
    first, results, discussion, figures = (read(path) for path in sources)
    inventory_match = INVENTORY.search(figures)
    if not inventory_match:
        raise ValueError('Cannot find the supplementary table inventory heading in chapter 04.')
    legends = figures[:inventory_match.start()].strip().removesuffix('---').rstrip()
    inventory = '##' + figures[inventory_match.start() + 1:]
    image_dir = args.figures or args.root / '03_figures' / 'Publication_PlanA_CNS'
    original_images = source_images(image_dir, legends)
    snapshot = work / 'figure-sources'
    snapshot.mkdir()
    images = {key: Path(shutil.copyfile(path, snapshot / path.name)) for key, path in original_images.items()}
    input_hashes.update({path: sha256(images[key]) for key, path in original_images.items()})
    csv_dir = args.csv_dir or args.root / '02_results' / 'Supplementary_Tables_PlanA'
    tables = read_tables(inventory, csv_dir)
    input_hashes.update({csv_dir / row[1]: row[-1] for row in tables})
    mismatches = []
    for identifier, filename, stated, actual, *_ in tables:
        normalized = re.sub(r'\s+', '', stated.replace(',', '')).replace('x', '×')
        if normalized != actual.replace(' ', ''):
            mismatches.append({'table': identifier, 'csv': filename, 'source': stated, 'actual': actual})
            pattern = rf'(^\|\s*\*\*{re.escape(identifier)}\*\*\s*\|[^|]+\|\s*){re.escape(stated)}(?=\s*\|)'
            inventory, count = re.subn(pattern, lambda match: match[1] + actual, inventory, flags=re.M)
            if count != 1:
                raise ValueError(f'Could not refresh the dimensions in the generated inventory for {identifier}.')

    scope = read(args.scope).strip()
    scope_words = max(len(scope.split()), len(re.findall(r"\b\w+(?:['’]\w+)*\b", scope)))
    if not 0 < scope_words <= 200:
        raise ValueError(f'Scope must contain 1 to 200 words, found {scope_words}.')
    (destination / '02_Scope_Statement.txt').write_text(scope + '\n', encoding='utf-8')
    main = document('\n\n'.join(chapter.strip() for chapter in (first, results, discussion, legends)))
    make_docx(main, destination / '01_Manuscript_Main_Text_Frontiers.docx', work, numbered=True)
    pdf = make_pdf(main, destination / '01_Manuscript_Main_Text_Frontiers.pdf')

    tiff_dir = destination / '03_Main_Figures_TIFF'
    tiff_dir.mkdir()
    tiffs = [make_tiff(path, tiff_dir / f'Figure_{key}.tif') for key, path in images.items() if not key.startswith('S')]
    supplement, references = supplement_markdown(legends, inventory, images, work)
    make_docx(document(supplement, references), destination / '04_Supplementary_Material.docx', work)
    sheets = write_workbook(tables, destination / '05_Supplementary_Tables_Combined_PlanA.xlsx')
    manifest = {'scope_words': scope_words, 'pdf': pdf, 'main_tiffs': tiffs,
                'supplementary_figures': [{'source': p.name, 'sha256': sha256(p)} for key, p in images.items()
                                          if key.startswith('S')],
                'xlsx': {'worksheets': 1 + len(sheets), 'data_rows': sum(sheet['rows'] for sheet in sheets),
                         'sheets': sheets, 'inventory_dimension_mismatches': mismatches},
                'markdown_sources': {path.name: input_hashes[path] for path in sources}}
    changed = [str(path) for path, digest in input_hashes.items() if not path.is_file() or sha256(path) != digest]
    if changed:
        raise ValueError(f'Source files changed during the build. Retry after inputs are stable: {changed}')
    if source_images(image_dir, legends) != original_images:
        raise ValueError('Figure inputs changed during the build. Retry after inputs are stable.')
    (destination / 'packaging_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n',
                                                         encoding='utf-8')
    print(f'PDF: {pdf["pages"]} pages, {pdf["lines"]} consecutive lines. Scope: {scope_words} words.')
    print(f'Figures: {len(tiffs)} TIFFs. Supplementary DOCX: {len(references)} embedded figures.')
    print(f'Excel: {1 + len(sheets)} worksheets, {manifest["xlsx"]["data_rows"]} data rows, '
          f'{len(mismatches)} dimension discrepancies in the source inventory.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT, help='Project root (default: this script’s parent project).')
    parser.add_argument('--output', type=Path, required=True, help='A new output directory. Existing paths are refused.')
    parser.add_argument('--figures', type=Path, help='Directory with the local Figure PNGs.')
    parser.add_argument('--csv-dir', type=Path, help='Directory with the local supplementary CSVs.')
    parser.add_argument('--scope', type=Path, default=FORMAT / 'scope.txt', help='Local UTF-8 scope statement.')
    args = parser.parse_args()
    output = args.output.absolute()
    if output.exists() or output.is_symlink():
        parser.error(f'Output already exists: {output}')
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.frontiers-', dir=output.parent) as temp:
            work = Path(temp)
            stage = work / 'package'
            stage.mkdir()
            build(args, stage, work)
            if output.exists() or output.is_symlink():
                raise ValueError(f'Output appeared during the build: {output}')
            stage.rename(output)
    except (OSError, ValueError, RuntimeError, csv.Error) as error:
        parser.exit(1, f'Packaging failed: {error}\n')
    print(f'Package written to {output}')


if __name__ == '__main__':
    main()
