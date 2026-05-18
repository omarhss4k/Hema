#!/usr/bin/env python3
"""
Replace "Simulation de population virtuelle" section in corrige.docx
with extended §2.4 content as specified.
"""

from docx import Document
from docx.shared import Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
import os

COLOR_TITLE = RGBColor(44, 62, 80)
COLOR_EQUATION = RGBColor(60, 60, 120)
COLOR_HEADER_BG = "2C3E50"
COLOR_ROW_ALT = "F2F2F2"
COLOR_ROW_WHITE = "FFFFFF"


def set_cell_background(cell, hex_color):
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    shd = OxmlElement('w:shd')
    shd.set(qn('w:val'), 'clear')
    shd.set(qn('w:color'), 'auto')
    shd.set(qn('w:fill'), hex_color)
    tcPr.append(shd)


def set_cell_borders(cell):
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    tcBorders = OxmlElement('w:tcBorders')
    for border_name in ['top', 'left', 'bottom', 'right']:
        border = OxmlElement(f'w:{border_name}')
        border.set(qn('w:val'), 'single')
        border.set(qn('w:sz'), '4')
        border.set(qn('w:space'), '0')
        border.set(qn('w:color'), '000000')
        tcBorders.append(border)
    tcPr.append(tcBorders)


def set_run_font(run, font_name='Times New Roman'):
    rPr = run._element.get_or_add_rPr()
    rFonts = rPr.find(qn('w:rFonts'))
    if rFonts is None:
        rFonts = OxmlElement('w:rFonts')
        rPr.insert(0, rFonts)
    rFonts.set(qn('w:ascii'), font_name)
    rFonts.set(qn('w:hAnsi'), font_name)


def add_spacing(p_element, before_twips=None, after_twips=None):
    pPr = p_element.get_or_add_pPr()
    spacing = OxmlElement('w:spacing')
    if before_twips is not None:
        spacing.set(qn('w:before'), str(before_twips))
    if after_twips is not None:
        spacing.set(qn('w:after'), str(after_twips))
    pPr.append(spacing)


def make_heading2(doc, text):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    add_spacing(p._element, before_twips=240, after_twips=120)
    run = p.add_run(text)
    run.bold = True
    run.font.name = 'Times New Roman'
    run.font.size = Pt(14)
    run.font.color.rgb = COLOR_TITLE
    set_run_font(run, 'Times New Roman')
    return p


def make_heading3(doc, text):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    add_spacing(p._element, before_twips=180, after_twips=90)
    run = p.add_run(text)
    run.bold = True
    run.font.name = 'Times New Roman'
    run.font.size = Pt(12)
    run.font.color.rgb = COLOR_TITLE
    set_run_font(run, 'Times New Roman')
    return p


def make_body_para(doc, text, first_line_indent=True):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    pPr = p._element.get_or_add_pPr()
    if first_line_indent:
        ind = OxmlElement('w:ind')
        ind.set(qn('w:firstLine'), '567')
        pPr.append(ind)
    add_spacing(p._element, after_twips=120)
    run = p.add_run(text)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(12)
    set_run_font(run, 'Times New Roman')
    return p


def make_equation_para(doc, text):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    add_spacing(p._element, before_twips=60, after_twips=60)
    run = p.add_run(text)
    run.italic = True
    run.font.name = 'Courier New'
    run.font.size = Pt(11)
    run.font.color.rgb = COLOR_EQUATION
    set_run_font(run, 'Courier New')
    return p


def make_bullet_para(doc, text):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    pPr = p._element.get_or_add_pPr()
    ind = OxmlElement('w:ind')
    ind.set(qn('w:left'), '284')
    ind.set(qn('w:hanging'), '284')
    pPr.append(ind)
    add_spacing(p._element, after_twips=60)
    run = p.add_run(text)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(12)
    set_run_font(run, 'Times New Roman')
    return p


def make_table(doc):
    headers = ['Paramètre', 'Valeur typique (θ_pop)', 'CV inter-individuel (ω)', 'Source']
    rows_data = [
        ['CL (L/h)', '0,50', '30%', 'FDA BLA 761139'],
        ['V1 (L)', '3,1', '25%', 'FDA BLA 761139'],
        ['Q (L/h)', '0,80', '30%', 'FDA BLA 761139'],
        ['V2 (L)', '2,5', '25%', 'FDA BLA 761139'],
        ['Slope_MPP (µM⁻¹)', '—', '20%', 'Calibration rat'],
        ['Slope_CMP (µM⁻¹)', '—', '20%', 'Calibration rat'],
        ['Slope_MEP (µM⁻¹)', '—', '20%', 'Calibration rat'],
        ['Poids corporel (kg)', '70,0', '15%', 'Littérature clinique'],
    ]

    table = doc.add_table(rows=1 + len(rows_data), cols=4)
    table.style = 'Table Grid'

    hrow = table.rows[0]
    for j, hdr in enumerate(headers):
        cell = hrow.cells[j]
        set_cell_background(cell, COLOR_HEADER_BG)
        set_cell_borders(cell)
        p = cell.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run(hdr)
        run.bold = True
        run.font.name = 'Times New Roman'
        run.font.size = Pt(11)
        run.font.color.rgb = RGBColor(255, 255, 255)
        set_run_font(run, 'Times New Roman')

    for i, row_data in enumerate(rows_data):
        row = table.rows[i + 1]
        bg = COLOR_ROW_ALT if i % 2 == 0 else COLOR_ROW_WHITE
        for j, cell_text in enumerate(row_data):
            cell = row.cells[j]
            set_cell_background(cell, bg)
            set_cell_borders(cell)
            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT if j == 0 else WD_ALIGN_PARAGRAPH.CENTER
            run = p.add_run(cell_text)
            run.font.name = 'Times New Roman'
            run.font.size = Pt(11)
            run.font.color.rgb = RGBColor(0, 0, 0)
            set_run_font(run, 'Times New Roman')

    return table


def main():
    path = '/home/user/Hema/memoire_stage_M2_PKPD_hematotoxicite_corrige.docx'
    doc = Document(path)
    body = doc.element.body

    # 1. Find start and end paragraphs
    start_para = None
    end_para = None
    for p in doc.paragraphs:
        text = p.text.strip()
        if 'Simulation de population virtuelle' in text and start_para is None:
            start_para = p
            print(f"Start: {repr(text[:80])}")
        elif start_para is not None and 'Grading CTCAE' in text and end_para is None:
            end_para = p
            print(f"End:   {repr(text[:80])}")
            break

    if start_para is None:
        raise ValueError("Start not found")
    if end_para is None:
        raise ValueError("End not found")

    # 2. Determine body indices and remove old elements
    body_children = list(body)
    start_idx = next(i for i, c in enumerate(body_children) if c is start_para._element)
    end_idx = next(i for i, c in enumerate(body_children) if c is end_para._element)

    print(f"Body indices: start={start_idx}, end={end_idx}, removing {end_idx - start_idx} elements")
    for el in body_children[start_idx:end_idx]:
        el.getparent().remove(el)

    # 3. Build new content by appending to doc (at end of body temporarily)
    #    We'll collect the elements and then move them before end_el.

    # Marker: record current last child of body before we start adding
    # We'll add to doc normally, then move elements into position.

    staging_elements = []

    def capture(x):
        """After adding element to doc, capture its element reference."""
        staging_elements.append(x._element)
        return x

    # Build content
    capture(make_heading2(doc, "2.4 Simulation de population virtuelle"))
    capture(make_body_para(doc,
        "Afin d'évaluer la distribution des grades de toxicité hématologique attendus en population, "
        "une approche de simulation Monte-Carlo a été mise en œuvre. Une cohorte virtuelle de N=300 "
        "patients a été générée pour le T-DXd (5,4 mg/kg Q3W × 6 cycles), en intégrant la variabilité "
        "inter-individuelle sur les paramètres PK et PD."))
    capture(make_heading3(doc, "Modélisation de la variabilité inter-individuelle"))
    capture(make_body_para(doc,
        "Les paramètres individuels sont supposés distribués selon une loi log-normale, ce qui garantit "
        "leur positivité et est cohérent avec la distribution observée des paramètres PK en population clinique :",
        first_line_indent=False))
    capture(make_equation_para(doc, "θᵢ = θ_pop × exp(ηᵢ)    avec ηᵢ ~ N(0, ω²)"))
    capture(make_body_para(doc,
        "où θ_pop est la valeur typique de population, ηᵢ l'effet aléatoire individuel et ω² la variance "
        "inter-individuelle. Le coefficient de variation (CV%) associé est approximé par CV% ≈ ω × 100 pour "
        "des valeurs de ω < 0,5. Les valeurs retenues, issues de l'analyse de population FDA (BLA 761139) "
        "et de la littérature, sont :",
        first_line_indent=False))
    tbl = make_table(doc)
    staging_elements.append(tbl._element)
    capture(make_body_para(doc, "", first_line_indent=False))  # spacer
    capture(make_heading3(doc, "Procédure de simulation Monte-Carlo"))
    capture(make_body_para(doc,
        "Pour chaque patient simulé i (i = 1, …, 300), la procédure suit les étapes suivantes :",
        first_line_indent=False))
    for b in [
        "• Tirage aléatoire des paramètres individuels : θᵢ = θ_pop × exp(ηᵢ), avec ηᵢ ~ N(0, ω²) pour chaque paramètre PK et PD",
        "• Calcul de la dose individuelle : Dose_mg = 5,4 mg/kg × BWᵢ, arrondie à la dizaine de mg (pratique clinique standard)",
        "• Résolution numérique du système ODE complet (PK + Damage + PD) via rxode2 (solveur LSODA, pas adaptatif) sur 126 jours, avec administration IV aux jours 1, 22, 43, 64, 85, 106",
        "• Extraction du nadir pour chaque lignée : min(Neut(t)), min(Plt(t)), min(RBC(t)) sur t ∈ [0, 126 jours]",
        "• Attribution du grade CTCAE v5 par comparaison du nadir aux seuils (cf. §2.5)",
    ]:
        capture(make_bullet_para(doc, b))
    capture(make_body_para(doc,
        "Le générateur pseudo-aléatoire est initialisé avec une graine fixe (set.seed(42)) garantissant "
        "la reproductibilité exacte des résultats. La taille de N=300 patients a été choisie pour assurer "
        "une estimation stable des proportions de grades rares (G4 < 5%) avec une erreur standard inférieure "
        "à 1,5 point de pourcentage (intervalle de confiance à 95% : ±1,5%).",
        first_line_indent=False))
    capture(make_heading3(doc, "Résumé statistique des grades simulés"))
    capture(make_body_para(doc,
        "Pour chaque toxicité, les résultats sont résumés par la proportion de patients atteignant chaque "
        "grade (G0 à G4), le taux de tout grade (G≥1 = 100% − %G0) et le taux de grade sévère "
        "(G3-4 = %G3 + %G4). Ces métriques sont directement comparables aux données de fréquence rapportées "
        "dans les notices médicamenteuses et les publications d'essais cliniques.",
        first_line_indent=False))

    # 4. Move all staged elements from their current (end-of-body) position
    #    to just before end_el
    end_el = end_para._element
    for el in staging_elements:
        # Remove from current position (end of body)
        parent = el.getparent()
        if parent is not None:
            parent.remove(el)
        # Insert before end_el
        end_el.addprevious(el)

    doc.save(path)
    print(f"\nSaved: {path}")
    size = os.path.getsize(path)
    print(f"File size: {size:,} bytes ({size/1024:.1f} KB)")
    if size < 500 * 1024:
        print("WARNING: File < 500 KB!")
    else:
        print("OK: File > 500 KB.")

    # Verify
    doc2 = Document(path)
    found_24 = any('2.4 Simulation de population virtuelle' in p.text for p in doc2.paragraphs)
    found_grading = any('Grading CTCAE' in p.text for p in doc2.paragraphs)
    found_h3 = any('Modélisation de la variabilité' in p.text for p in doc2.paragraphs)
    found_montecarlo = any('Procédure de simulation Monte-Carlo' in p.text for p in doc2.paragraphs)
    found_resume = any('Résumé statistique' in p.text for p in doc2.paragraphs)
    print(f"\nVerification:")
    print(f"  §2.4 title: {found_24}")
    print(f"  Grading CTCAE: {found_grading}")
    print(f"  H3 variabilité: {found_h3}")
    print(f"  H3 Monte-Carlo: {found_montecarlo}")
    print(f"  H3 Résumé: {found_resume}")
    # Check table
    tbl_count = len(doc2.tables)
    print(f"  Total tables: {tbl_count}")


if __name__ == '__main__':
    main()
