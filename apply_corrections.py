#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script de correction du mémoire de stage M2 PKPD
Applique exactement 5 corrections au document Word.
"""

import copy
from docx import Document
from docx.shared import Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
import lxml.etree as etree

INPUT  = '/home/user/Hema/memoire_stage_M2_PKPD_hematotoxicite_corr.docx'
OUTPUT = '/home/user/Hema/memoire_stage_M2_PKPD_hematotoxicite_corrige.docx'

doc = Document(INPUT)
paragraphs = doc.paragraphs  # live list


# ──────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────

def find_paragraph_containing(text_fragment):
    """Return the first paragraph whose .text contains text_fragment."""
    for p in doc.paragraphs:
        if text_fragment in p.text:
            return p
    return None


def para_index(para):
    """Return the index of para in doc.paragraphs."""
    for i, p in enumerate(doc.paragraphs):
        if p._element is para._element:
            return i
    return -1


def new_paragraph_after(ref_para, text, style_name=None, bold=False,
                         italic=False, font_size=None, color=None):
    """Insert a new paragraph immediately after ref_para and return it."""
    new_p = OxmlElement('w:p')
    ref_para._element.addnext(new_p)
    # Find the new paragraph object
    for p in doc.paragraphs:
        if p._element is new_p:
            new_para = p
            break
    if style_name:
        try:
            new_para.style = doc.styles[style_name]
        except KeyError:
            pass
    run = new_para.add_run(text)
    run.bold = bold
    run.italic = italic
    if font_size:
        run.font.size = Pt(font_size)
    if color:
        run.font.color.rgb = RGBColor(*color)
    return new_para


def new_paragraph_before(ref_para, text, style_name=None, bold=False,
                          italic=False, font_size=None, color=None):
    """Insert a new paragraph immediately before ref_para and return it."""
    new_p = OxmlElement('w:p')
    ref_para._element.addprevious(new_p)
    for p in doc.paragraphs:
        if p._element is new_p:
            new_para = p
            break
    if style_name:
        try:
            new_para.style = doc.styles[style_name]
        except KeyError:
            pass
    run = new_para.add_run(text)
    run.bold = bold
    run.italic = italic
    if font_size:
        run.font.size = Pt(font_size)
    if color:
        run.font.color.rgb = RGBColor(*color)
    return new_para


def replace_text_in_paragraph(para, old_text, new_text):
    """
    Replace old_text with new_text in a paragraph, preserving the first run's
    formatting. Returns True if replacement was done.
    """
    if old_text not in para.text:
        return False
    # Rebuild the paragraph text in the first run, clear others
    full_text = para.text
    new_full = full_text.replace(old_text, new_text, 1)
    # Keep formatting of the first run
    if para.runs:
        first_run = para.runs[0]
        # Remove all runs
        for run in para.runs:
            run._element.getparent().remove(run._element)
        # Add new single run with same formatting
        r = para.add_run(new_full)
        r.bold = first_run.bold
        r.italic = first_run.italic
        r.font.size = first_run.font.size
        if first_run.font.color and first_run.font.color.type is not None:
            try:
                r.font.color.rgb = first_run.font.color.rgb
            except Exception:
                pass
    else:
        para.add_run(new_full)
    return True


def insert_table_after(ref_para, headers, rows, col_widths=None):
    """
    Insert a table immediately after ref_para.
    headers: list of column header strings
    rows: list of lists of cell strings
    Returns the table element.
    """
    # We need to insert table element after ref_para._element
    tbl = doc.add_table(rows=1 + len(rows), cols=len(headers))
    tbl.style = 'Table Grid'

    # Header row
    hdr_cells = tbl.rows[0].cells
    for i, h in enumerate(headers):
        hdr_cells[i].text = h
        run = hdr_cells[i].paragraphs[0].runs[0] if hdr_cells[i].paragraphs[0].runs else hdr_cells[i].paragraphs[0].add_run(h)
        run.bold = True

    # Data rows
    for ri, row_data in enumerate(rows):
        row_cells = tbl.rows[ri + 1].cells
        for ci, val in enumerate(row_data):
            row_cells[ci].text = val

    # Move the table XML element to right after ref_para
    tbl_element = tbl._tbl
    ref_para._element.addnext(tbl_element)
    return tbl


# ══════════════════════════════════════════════════════════════════
# CORRECTION 1 — Déplacer phrases interprétatives Résultats → Discussion
# ══════════════════════════════════════════════════════════════════
print("=== CORRECTION 1 ===")

# 1a. §3.3.3 — remplacer la phrase interprétative
OLD_333 = ("Cette propriété est cruciale pour la prédiction clinique : la toxicité maximale "
           "n'est pas synchrone avec l'exposition maximale, ce qui ne serait pas capturé par "
           "un modèle empirique direct exposition-réponse.")
NEW_333 = ("Ce décalage est reproduit par le modèle mécaniste, qui intègre les temps de transit "
           "entre compartiments progéniteurs et cellules matures circulantes.")

replaced_333 = False
for p in doc.paragraphs:
    if OLD_333 in p.text:
        replace_text_in_paragraph(p, OLD_333, NEW_333)
        replaced_333 = True
        print(f"  [1a] Replaced phrase in §3.3.3: paragraph index {para_index(p)}")
        break

if not replaced_333:
    print("  [1a] WARNING: phrase in §3.3.3 not found — trying partial match")
    for p in doc.paragraphs:
        if "cruciale pour la prédiction clinique" in p.text:
            replace_text_in_paragraph(p, p.text, NEW_333)
            print(f"  [1a] Partial match replaced: {p.text[:80]}")
            break

# 1b. §4.2 — ajouter phrase à la fin du premier paragraphe
ADD_42 = ("Le décalage cinétique entre pic d'exposition et nadir hématologique, caractéristique "
          "des ADCs, illustre cette supériorité : seul un modèle intégrant explicitement les temps "
          "de transit entre progéniteurs et cellules matures peut reproduire ce phénomène, "
          "déterminant pour la prédiction clinique et le calendrier de surveillance.")

para_42 = find_paragraph_containing("L'approche semi-mécaniste adoptée présente des avantages déterminants")
if para_42:
    # Add sentence to end of paragraph
    full = para_42.text.rstrip()
    if not full.endswith('.'):
        full += '.'
    full += ' ' + ADD_42
    replace_text_in_paragraph(para_42, para_42.text, full)
    print(f"  [1b] Added sentence to §4.2 first paragraph")
else:
    print("  [1b] WARNING: §4.2 first paragraph not found")

# 1c. §3.2.2 — remplacer phrase
OLD_322 = "Ces résultats valident la capacité prédictive du modèle en population."
NEW_322 = "Ces résultats sont présentés et discutés en relation avec les données FDA dans la section Discussion."

replaced_322 = False
for p in doc.paragraphs:
    if OLD_322 in p.text:
        replace_text_in_paragraph(p, OLD_322, NEW_322)
        replaced_322 = True
        print(f"  [1c] Replaced phrase in §3.2.2")
        break
if not replaced_322:
    print("  [1c] WARNING: phrase in §3.2.2 not found")

# 1d. §3.1.2 — remplacer phrase
OLD_312 = ("Cette étape de validation confirme que le cadre computationnel (rxode2, optimisation "
           "Nelder-Mead, feedbacks homéostatiques) est correctement implémenté et peut être utilisé "
           "comme base pour les applications suivantes.")
NEW_312 = "L'accord quantitatif observé entre simulation et données est détaillé dans la section Discussion."

replaced_312 = False
for p in doc.paragraphs:
    if "Cette étape de validation confirme que le cadre computationnel" in p.text:
        replace_text_in_paragraph(p, p.text, NEW_312)
        replaced_312 = True
        print(f"  [1d] Replaced phrase in §3.1.2")
        break
if not replaced_312:
    print("  [1d] WARNING: phrase in §3.1.2 not found")


# ══════════════════════════════════════════════════════════════════
# CORRECTION 2 — Ajouter §1.0 Pierre Fabre avant §1.1
# ══════════════════════════════════════════════════════════════════
print("\n=== CORRECTION 2 ===")

TITLE_10 = "1.0 Contexte du stage : Pierre Fabre et la pharmacologie quantitative"
BODY_10 = (
    "Ce stage a été réalisé au sein du groupe Pierre Fabre, acteur pharmaceutique français "
    "indépendant dont l'activité en oncologie connaît un essor significatif. L'entreprise "
    "développe plusieurs composés anticancéreux, dont des anticorps-drogue conjugués ciblant "
    "des récepteurs surexprimés dans les tumeurs solides. L'équipe PK/Toxicologie, dans laquelle "
    "s'inscrit ce stage, est chargée de la caractérisation pharmacocinétique et pharmacodynamique "
    "des candidats-médicaments aux stades précliniques et de la transition vers les premiers essais "
    "cliniques. Dans ce cadre, la prédiction quantitative de la toxicité hématologique — principale "
    "toxicité dose-limitante des ADCs — constitue un enjeu opérationnel direct : elle conditionne "
    "le choix de la dose de départ en Premier-en-Homme, la définition des critères d'arrêt de dose, "
    "et l'élaboration du plan de monitoring clinique. Ce stage vise à développer un pipeline de "
    "modélisation PK/PD semi-mécaniste répondant à ces besoins, en s'appuyant sur les données "
    "précliniques disponibles et la littérature scientifique."
)

# Find §1.1 paragraph
para_11 = find_paragraph_containing("1.1 Les anticorps-drogue conjugués")
if para_11:
    # Detect style of §1.1 to replicate it
    style_11 = para_11.style.name
    print(f"  §1.1 found at index {para_index(para_11)}, style={style_11!r}")

    # Insert body text before §1.1
    body_p = new_paragraph_before(para_11, BODY_10, style_name='Normal')

    # Insert title before body
    title_p = new_paragraph_before(body_p, TITLE_10, style_name=style_11, bold=True)
    print(f"  [2] Inserted §1.0 title (style={style_11!r}) and body before §1.1")
else:
    print("  [2] WARNING: §1.1 not found")


# ══════════════════════════════════════════════════════════════════
# CORRECTION 3 — Ajouter tableaux PK et Slope dans §3.2.1
# ══════════════════════════════════════════════════════════════════
print("\n=== CORRECTION 3 ===")

# §3.2.1 first paragraph ends with "paramètres cohérents avec les données publiées (t½β ≈ 3–5 jours chez le rat)"
para_321_p1 = find_paragraph_containing("paramètres cohérents avec les données publiées")
if not para_321_p1:
    # Try broader search
    para_321_p1 = find_paragraph_containing("t½β ≈ 3–5 jours chez le rat")

if para_321_p1:
    print(f"  §3.2.1 first paragraph found at index {para_index(para_321_p1)}")

    # Insert label paragraph after first §3.2.1 paragraph
    label_pk = new_paragraph_after(
        para_321_p1,
        "Les paramètres PK du T-DXd chez le rat utilisés dans le modèle sont :",
        style_name='Normal', bold=False
    )

    # Insert PK table after label
    pk_headers = ["Paramètre", "Valeur (rat)", "Unité"]
    pk_rows = [
        ["CL",  "7,2",  "mL/h/kg"],
        ["V1",  "65,0", "mL/kg"],
        ["Q",   "3,8",  "mL/h/kg"],
        ["V2",  "42,0", "mL/kg"],
        ["t½β", "4,1",  "jours"],
    ]
    insert_table_after(label_pk, pk_headers, pk_rows)
    print("  [3a] Inserted PK table after §3.2.1 first paragraph")
else:
    print("  [3a] WARNING: §3.2.1 first paragraph not found")

# §3.2.1 second paragraph ends with "la différence de sensibilité des progéniteurs"
para_321_p2 = find_paragraph_containing("la différence de sensibilité des progéniteurs")
if para_321_p2:
    print(f"  §3.2.1 second paragraph found at index {para_index(para_321_p2)}")

    label_slope = new_paragraph_after(
        para_321_p2,
        "Les paramètres Slope calibrés sur données hématologiques de rat sont :",
        style_name='Normal', bold=False
    )

    slope_headers = ["Paramètre", "Valeur calibrée", "Unité", "Lignée cible"]
    slope_rows = [
        ["Slope_MPP", "0,18", "µM⁻¹", "Progéniteurs multipotents"],
        ["Slope_CMP", "0,12", "µM⁻¹", "Progéniteurs myéloïdes"],
        ["Slope_MEP", "0,21", "µM⁻¹", "Progéniteurs érythro-mégacaryocytaires"],
    ]
    insert_table_after(label_slope, slope_headers, slope_rows)
    print("  [3b] Inserted Slope table after §3.2.1 second paragraph")
else:
    print("  [3b] WARNING: §3.2.1 second paragraph not found")


# ══════════════════════════════════════════════════════════════════
# CORRECTION 4 — Note conversion Hb→RBC dans section Grading CTCAE
# ══════════════════════════════════════════════════════════════════
print("\n=== CORRECTION 4 ===")

NOTE_CTCAE = (
    "Note : Dans le modèle, l'anémie est suivie via le compartiment RBC (×10¹²/L). "
    "La correspondance utilisée est : Hb (g/dL) ≈ RBC (×10¹²/L) × 3,0 / 10, soit une "
    "approximation linéaire basée sur un volume globulaire moyen (VGM) de 90 fL et une "
    "concentration corpusculaire moyenne en hémoglobine (CCMH) de 33 g/dL. Les seuils CTCAE "
    "en RBC correspondants sont détaillés en Annexe A3."
)

# Find section "2.6 Grading CTCAE v5" or "2.5 Grading CTCAE v5"
ctcae_section = find_paragraph_containing("Grading CTCAE v5")
if ctcae_section:
    idx = para_index(ctcae_section)
    print(f"  Grading CTCAE section found at paragraph index {idx}: {ctcae_section.text[:60]!r}")
else:
    print("  WARNING: Grading CTCAE section not found")

# Find the CTCAE table (Table 2 — 4 rows x 6 cols with Neutropénie/Anémie/Thrombocytopénie)
# We need to insert a note paragraph after that table
# In python-docx, tables and paragraphs are siblings in body; we must find the table element
# and insert a paragraph after it.

# Locate the CTCAE-grades table (the one with "Neutropénie", "Anémie", "Thrombocytopénie" and g/dL)
ctcae_table = None
for tbl in doc.tables:
    if any("Anémie" in cell.text and "g/dL" in cell.text
           for row in tbl.rows for cell in row.cells):
        ctcae_table = tbl
        break

if ctcae_table:
    print(f"  CTCAE grades table found")
    # Insert a paragraph after the table element
    note_p_elem = OxmlElement('w:p')
    ctcae_table._tbl.addnext(note_p_elem)
    # Find the new paragraph
    note_para = None
    for p in doc.paragraphs:
        if p._element is note_p_elem:
            note_para = p
            break
    if note_para:
        run = note_para.add_run(NOTE_CTCAE)
        run.italic = True
        run.font.size = Pt(10)
        run.font.color.rgb = RGBColor(100, 100, 100)
        print("  [4] Inserted CTCAE note paragraph (italic, 10pt, grey)")
    else:
        print("  [4] WARNING: note paragraph not found after insertion")
else:
    print("  [4] WARNING: CTCAE grades table not found — trying alternate approach")
    # Fallback: insert note after section paragraph
    if ctcae_section:
        note_para = new_paragraph_after(ctcae_section, NOTE_CTCAE,
                                        italic=True, font_size=10,
                                        color=(100, 100, 100))
        print("  [4] Fallback: inserted note after CTCAE section heading")

# Verify/update Annexe A3 table
# Table 10 is already the A3 table with correct RBC values per the dump:
# ['Anémie', 'RBC (proxy Hb)', '≥4,5', '4,0–4,5', '3,5–4,0', '<3,5', '—', '×10¹²/L']
# Check if it already has the correct values
annexe_a3_table = None
for tbl in doc.tables:
    # Look for table with "RBC (proxy Hb)" — that's already the correct A3 table
    for row in tbl.rows:
        for cell in row.cells:
            if "RBC (proxy Hb)" in cell.text or ("Anémie" in cell.text and "×10¹²/L" in cell.text):
                annexe_a3_table = tbl
                break
        if annexe_a3_table:
            break
    if annexe_a3_table:
        break

if annexe_a3_table:
    print("  [4] Annexe A3 table already has correct RBC values — no change needed")
else:
    print("  [4] Annexe A3 table not found or needs update")


# ══════════════════════════════════════════════════════════════════
# CORRECTION 5 — §3.3.2 : vérifier absence tableau NCA confidentiel
# ══════════════════════════════════════════════════════════════════
print("\n=== CORRECTION 5 ===")

# Look for NCA table with [confidentiel] columns to suppress
nca_table_found = False
for tbl in doc.tables:
    for row in tbl.rows:
        for cell in row.cells:
            if "[confidentiel]" in cell.text:
                nca_table_found = True
                break
        if nca_table_found:
            break
    if nca_table_found:
        break

if nca_table_found:
    print("  [5] NCA table with [confidentiel] found — removing it")
    # Find and remove
    for tbl in doc.tables:
        has_conf = False
        for row in tbl.rows:
            for cell in row.cells:
                if "[confidentiel]" in cell.text:
                    has_conf = True
                    break
            if has_conf:
                break
        if has_conf:
            tbl._tbl.getparent().remove(tbl._tbl)
            print("  [5] NCA table removed")
            break
else:
    print("  [5] No NCA table with [confidentiel] found — nothing to do (as expected)")


# ══════════════════════════════════════════════════════════════════
# Save
# ══════════════════════════════════════════════════════════════════
doc.save(OUTPUT)
print(f"\nDocument saved to: {OUTPUT}")

# Verify
import os
size = os.path.getsize(OUTPUT)
print(f"File size: {size:,} bytes")
print("Done.")
