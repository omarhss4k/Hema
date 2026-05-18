#!/usr/bin/env python3
"""
Update section §2.4 in memoire_stage_M2_PKPD_hematotoxicite_corrige.docx
Replace the short §2.5 "Simulation de population virtuelle" with extended §2.4 content.
"""

import copy
from docx import Document
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.shared import Pt, RGBColor, Cm, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from lxml import etree

INPUT_PATH = "/home/user/Hema/memoire_stage_M2_PKPD_hematotoxicite_corrige.docx"
OUTPUT_PATH = "/home/user/Hema/memoire_stage_M2_PKPD_hematotoxicite_corrige.docx"

# Colors
COLOR_TITLE = RGBColor(44, 62, 80)       # #2C3E50
COLOR_HEADER_BG = "2C3E50"               # hex for XML
COLOR_EQUATION = RGBColor(60, 60, 120)   # #3C3C78
COLOR_ROW_ALT = "F2F2F2"
COLOR_WHITE = "FFFFFF"

def set_cell_background(cell, color_hex):
    """Set cell background color."""
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    shd = OxmlElement('w:shd')
    shd.set(qn('w:val'), 'clear')
    shd.set(qn('w:color'), 'auto')
    shd.set(qn('w:fill'), color_hex)
    # Remove existing shd if any
    existing = tcPr.find(qn('w:shd'))
    if existing is not None:
        tcPr.remove(existing)
    tcPr.append(shd)

def make_paragraph_element(text, font_name="Times New Roman", font_size_pt=12,
                            bold=False, italic=False, color=None,
                            alignment=WD_ALIGN_PARAGRAPH.JUSTIFY,
                            indent_cm=None, space_after_pt=6,
                            bullet=False, bullet_indent_cm=0.5):
    """Create a paragraph XML element with specified formatting."""
    p = OxmlElement('w:p')
    pPr = OxmlElement('w:pPr')
    p.append(pPr)

    # Alignment
    jc = OxmlElement('w:jc')
    align_map = {
        WD_ALIGN_PARAGRAPH.JUSTIFY: 'both',
        WD_ALIGN_PARAGRAPH.CENTER: 'center',
        WD_ALIGN_PARAGRAPH.LEFT: 'left',
    }
    jc.set(qn('w:val'), align_map.get(alignment, 'both'))
    pPr.append(jc)

    # Spacing
    spacing = OxmlElement('w:spacing')
    spacing.set(qn('w:after'), str(int(space_after_pt * 20)))  # twips
    pPr.append(spacing)

    # Indent
    if indent_cm is not None or bullet:
        ind = OxmlElement('w:ind')
        if bullet:
            ind.set(qn('w:left'), str(int(bullet_indent_cm * 360)))  # 360 twips per cm
            ind.set(qn('w:hanging'), '360')
        elif indent_cm is not None:
            ind.set(qn('w:firstLine'), str(int(indent_cm * 360)))
        pPr.append(ind)

    # Run
    r = OxmlElement('w:r')
    rPr = OxmlElement('w:rPr')
    r.append(rPr)

    # Font
    rFonts = OxmlElement('w:rFonts')
    rFonts.set(qn('w:ascii'), font_name)
    rFonts.set(qn('w:hAnsi'), font_name)
    rPr.append(rFonts)

    # Size
    sz = OxmlElement('w:sz')
    sz.set(qn('w:val'), str(int(font_size_pt * 2)))
    rPr.append(sz)
    szCs = OxmlElement('w:szCs')
    szCs.set(qn('w:val'), str(int(font_size_pt * 2)))
    rPr.append(szCs)

    # Bold
    if bold:
        b = OxmlElement('w:b')
        rPr.append(b)
        bCs = OxmlElement('w:bCs')
        rPr.append(bCs)

    # Italic
    if italic:
        i_elem = OxmlElement('w:i')
        rPr.append(i_elem)
        iCs = OxmlElement('w:iCs')
        rPr.append(iCs)

    # Color
    if color is not None:
        color_elem = OxmlElement('w:color')
        if isinstance(color, RGBColor):
            hex_color = '{:02X}{:02X}{:02X}'.format(color[0], color[1], color[2])
        else:
            hex_color = color
        color_elem.set(qn('w:val'), hex_color)
        rPr.append(color_elem)

    # Text
    t = OxmlElement('w:t')
    t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    t.text = text
    r.append(t)
    p.append(r)

    return p


def make_table_element(doc, headers, rows):
    """Create a formatted table XML element."""
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    table = doc.add_table(rows=1 + len(rows), cols=len(headers))
    table.style = 'Table Grid'

    # Header row
    hdr_row = table.rows[0]
    for i, (cell, hdr_text) in enumerate(zip(hdr_row.cells, headers)):
        set_cell_background(cell, COLOR_HEADER_BG)
        cell.text = ''
        p = cell.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run(hdr_text)
        run.bold = True
        run.font.name = 'Times New Roman'
        run.font.size = Pt(11)
        run.font.color.rgb = RGBColor(255, 255, 255)

    # Data rows
    for r_idx, row_data in enumerate(rows):
        row = table.rows[r_idx + 1]
        bg = COLOR_ROW_ALT if r_idx % 2 == 0 else COLOR_WHITE
        for c_idx, (cell, cell_text) in enumerate(zip(row.cells, row_data)):
            set_cell_background(cell, bg)
            cell.text = ''
            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER if c_idx > 0 else WD_ALIGN_PARAGRAPH.LEFT
            run = p.add_run(cell_text)
            run.font.name = 'Times New Roman'
            run.font.size = Pt(11)

    return table._tbl


def main():
    doc = Document(INPUT_PATH)
    body = doc.element.body

    # ---- Step 1: Find the range to replace ----
    # Find paragraph with "2.5 Simulation de population virtuelle" -> start
    # Find paragraph with "2.6 Grading CTCAE" -> end (keep this one)

    start_idx = None
    end_idx = None

    body_children = list(body)
    for idx, child in enumerate(body_children):
        if child.tag.endswith('}p'):
            texts = [r.text for r in child.iter() if r.tag.endswith('}t') and r.text]
            full_text = ''.join(texts).strip()
            if start_idx is None and "Simulation de population virtuelle" in full_text:
                start_idx = idx
                print(f"Found start at body[{idx}]: {full_text[:80]}")
            elif start_idx is not None and end_idx is None and "Grading CTCAE" in full_text:
                end_idx = idx
                print(f"Found end at body[{idx}]: {full_text[:80]}")
                break

    if start_idx is None or end_idx is None:
        raise ValueError(f"Could not find section boundaries. start_idx={start_idx}, end_idx={end_idx}")

    print(f"Will remove body[{start_idx}] to body[{end_idx-1}] (inclusive)")

    # ---- Step 2: Remove elements from start_idx to end_idx-1 ----
    # Collect elements to remove (re-index after each removal)
    elements_to_remove = body_children[start_idx:end_idx]
    for elem in elements_to_remove:
        body.remove(elem)

    # The anchor is now the element that was at end_idx (now at start_idx position)
    # We need to insert new content before this element
    anchor_element = list(body)[start_idx]
    print(f"Anchor element (§2.6 Grading CTCAE): {anchor_element.tag}")

    # ---- Step 3: Build new content ----
    new_elements = []

    # 1. Title H2: "2.4 Simulation de population virtuelle"
    p1 = make_paragraph_element(
        "2.4 Simulation de population virtuelle",
        font_name="Times New Roman", font_size_pt=14,
        bold=True, color=COLOR_TITLE,
        alignment=WD_ALIGN_PARAGRAPH.LEFT,
        space_after_pt=6
    )
    new_elements.append(p1)

    # 2. Body paragraph
    p2 = make_paragraph_element(
        "Afin d'évaluer la distribution des grades de toxicité hématologique attendus en population, "
        "une approche de simulation Monte-Carlo a été mise en œuvre. Une cohorte virtuelle de N=300 patients "
        "a été générée pour le T-DXd (5,4 mg/kg Q3W × 6 cycles), en intégrant la variabilité inter-individuelle "
        "sur les paramètres PK et PD.",
        font_name="Times New Roman", font_size_pt=12,
        alignment=WD_ALIGN_PARAGRAPH.JUSTIFY,
        indent_cm=1.0, space_after_pt=6
    )
    new_elements.append(p2)

    # 3. Title H3: "Modélisation de la variabilité inter-individuelle"
    p3 = make_paragraph_element(
        "Modélisation de la variabilité inter-individuelle",
        font_name="Times New Roman", font_size_pt=12,
        bold=True, color=COLOR_TITLE,
        alignment=WD_ALIGN_PARAGRAPH.LEFT,
        space_after_pt=6
    )
    new_elements.append(p3)

    # 4. Body paragraph
    p4 = make_paragraph_element(
        "Les paramètres individuels sont supposés distribués selon une loi log-normale, ce qui garantit leur "
        "positivité et est cohérent avec la distribution observée des paramètres PK en population clinique :",
        font_name="Times New Roman", font_size_pt=12,
        alignment=WD_ALIGN_PARAGRAPH.JUSTIFY,
        space_after_pt=6
    )
    new_elements.append(p4)

    # 5. Equation (centered, Courier New 11pt italic)
    p5 = make_paragraph_element(
        "θᵢ = θ_pop × exp(ηᵢ)    avec ηᵢ ~ N(0, ω²)",
        font_name="Courier New", font_size_pt=11,
        italic=True, color=COLOR_EQUATION,
        alignment=WD_ALIGN_PARAGRAPH.CENTER,
        space_after_pt=6
    )
    new_elements.append(p5)

    # 6. Body paragraph
    p6 = make_paragraph_element(
        "où θ_pop est la valeur typique de population, ηᵢ l'effet aléatoire individuel et ω² la variance "
        "inter-individuelle. Le coefficient de variation (CV%) associé est approximé par CV% ≈ ω × 100 pour "
        "des valeurs de ω < 0,5. Les valeurs retenues, issues de l'analyse de population FDA (BLA 761139) "
        "et de la littérature, sont :",
        font_name="Times New Roman", font_size_pt=12,
        alignment=WD_ALIGN_PARAGRAPH.JUSTIFY,
        space_after_pt=6
    )
    new_elements.append(p6)

    # 7. Table: PK parameters
    table_headers = ["Paramètre", "Valeur typique (θ_pop)", "CV inter-individuel (ω)", "Source"]
    table_rows = [
        ["CL (L/h)", "0,50", "30%", "FDA BLA 761139"],
        ["V1 (L)", "3,1", "25%", "FDA BLA 761139"],
        ["Q (L/h)", "0,80", "30%", "FDA BLA 761139"],
        ["V2 (L)", "2,5", "25%", "FDA BLA 761139"],
        ["Slope_MPP (µM⁻¹)", "—", "20%", "Calibration rat"],
        ["Slope_CMP (µM⁻¹)", "—", "20%", "Calibration rat"],
        ["Slope_MEP (µM⁻¹)", "—", "20%", "Calibration rat"],
        ["Poids corporel (kg)", "70,0", "15%", "Littérature clinique"],
    ]
    tbl_elem = make_table_element(doc, table_headers, table_rows)
    new_elements.append(tbl_elem)

    # Add a spacing paragraph after table
    p_space = make_paragraph_element("", font_size_pt=6, space_after_pt=0)
    new_elements.append(p_space)

    # 8. Title H3: "Procédure de simulation Monte-Carlo"
    p8 = make_paragraph_element(
        "Procédure de simulation Monte-Carlo",
        font_name="Times New Roman", font_size_pt=12,
        bold=True, color=COLOR_TITLE,
        alignment=WD_ALIGN_PARAGRAPH.LEFT,
        space_after_pt=6
    )
    new_elements.append(p8)

    # 9. Body paragraph
    p9 = make_paragraph_element(
        "Pour chaque patient simulé i (i = 1, …, 300), la procédure suit les étapes suivantes :",
        font_name="Times New Roman", font_size_pt=12,
        alignment=WD_ALIGN_PARAGRAPH.JUSTIFY,
        space_after_pt=6
    )
    new_elements.append(p9)

    # 10. Bullet points
    bullets = [
        "Tirage aléatoire des paramètres individuels : θᵢ = θ_pop × exp(ηᵢ), avec ηᵢ ~ N(0, ω²) pour chaque paramètre PK et PD",
        "Calcul de la dose individuelle : Dose_mg = 5,4 mg/kg × BWᵢ, arrondie à la dizaine de mg (pratique clinique standard)",
        "Résolution numérique du système ODE complet (PK + Damage + PD) via rxode2 (solveur LSODA, pas adaptatif) sur 126 jours, avec administration IV aux jours 1, 22, 43, 64, 85, 106",
        "Extraction du nadir pour chaque lignée : min(Neut(t)), min(Plt(t)), min(RBC(t)) sur t ∈ [0, 126 jours]",
        "Attribution du grade CTCAE v5 par comparaison du nadir aux seuils (cf. §2.5)",
    ]
    for bullet_text in bullets:
        pb = make_paragraph_element(
            "• " + bullet_text,
            font_name="Times New Roman", font_size_pt=12,
            alignment=WD_ALIGN_PARAGRAPH.JUSTIFY,
            bullet=True, bullet_indent_cm=0.5,
            space_after_pt=3
        )
        new_elements.append(pb)

    # 11. Body paragraph (seed/reproducibility)
    p11 = make_paragraph_element(
        "Le générateur pseudo-aléatoire est initialisé avec une graine fixe (set.seed(42)) garantissant la "
        "reproductibilité exacte des résultats. La taille de N=300 patients a été choisie pour assurer une "
        "estimation stable des proportions de grades rares (G4 < 5%) avec une erreur standard inférieure à "
        "1,5 point de pourcentage (intervalle de confiance à 95% : ±1,5%).",
        font_name="Times New Roman", font_size_pt=12,
        alignment=WD_ALIGN_PARAGRAPH.JUSTIFY,
        space_after_pt=6
    )
    new_elements.append(p11)

    # 12. Title H3: "Résumé statistique des grades simulés"
    p12 = make_paragraph_element(
        "Résumé statistique des grades simulés",
        font_name="Times New Roman", font_size_pt=12,
        bold=True, color=COLOR_TITLE,
        alignment=WD_ALIGN_PARAGRAPH.LEFT,
        space_after_pt=6
    )
    new_elements.append(p12)

    # 13. Body paragraph
    p13 = make_paragraph_element(
        "Pour chaque toxicité, les résultats sont résumés par la proportion de patients atteignant chaque "
        "grade (G0 à G4), le taux de tout grade (G≥1 = 100% − %G0) et le taux de grade sévère (G3-4 = %G3 + %G4). "
        "Ces métriques sont directement comparables aux données de fréquence rapportées dans les notices "
        "médicamenteuses et les publications d'essais cliniques.",
        font_name="Times New Roman", font_size_pt=12,
        alignment=WD_ALIGN_PARAGRAPH.JUSTIFY,
        space_after_pt=6
    )
    new_elements.append(p13)

    # ---- Step 4: Insert all new elements before the anchor ----
    for elem in new_elements:
        anchor_element.addprevious(elem)

    print(f"Inserted {len(new_elements)} new elements before anchor")

    # ---- Step 5: Save ----
    doc.save(OUTPUT_PATH)
    print(f"Saved to {OUTPUT_PATH}")

    # ---- Step 6: Verify file size ----
    import os
    size = os.path.getsize(OUTPUT_PATH)
    print(f"File size: {size:,} bytes ({size/1024:.1f} KB)")
    if size < 500 * 1024:
        print("WARNING: File is smaller than 500 KB!")
    else:
        print("OK: File size is above 500 KB")


if __name__ == "__main__":
    main()
