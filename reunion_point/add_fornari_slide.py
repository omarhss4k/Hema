#!/usr/bin/env python3
"""
Add a slide explaining what to keep and what to replace from the Fornari
carboplatin PK/PD model when adapting it for an ADC.
"""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Cm

PPTX_PATH = "/home/user/Hema/reunion_point/reunion_prochaine.pptx"

# Colors
GREEN_HEADER = RGBColor(0x00, 0xB0, 0x50)   # #00B050
RED_HEADER   = RGBColor(0xC0, 0x00, 0x00)   # #C00000 dark red/orange
WHITE        = RGBColor(0xFF, 0xFF, 0xFF)
BLACK        = RGBColor(0x00, 0x00, 0x00)
LIGHT_GREEN_BG = RGBColor(0xE2, 0xEF, 0xDA)  # soft green tint
LIGHT_RED_BG   = RGBColor(0xFF, 0xE0, 0xE0)  # soft red tint
NOTE_BG      = RGBColor(0xFF, 0xFF, 0xCC)    # light yellow for note box


def set_cell_background(cell, color: RGBColor):
    """Fill a table cell background with a solid color."""
    from pptx.oxml.ns import qn
    from lxml import etree
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    # Remove existing solidFill if present
    for existing in tcPr.findall(qn('a:solidFill')):
        tcPr.remove(existing)
    solidFill = etree.SubElement(tcPr, qn('a:solidFill'))
    srgbClr = etree.SubElement(solidFill, qn('a:srgbClr'))
    srgbClr.set('val', f'{color.rgb:06X}')


def add_text_box(slide, text, left, top, width, height,
                 font_size=14, bold=False, color=BLACK,
                 bg_color=None, align=PP_ALIGN.LEFT,
                 italic=False):
    """Add a styled text box to a slide."""
    from pptx.oxml.ns import qn
    from lxml import etree

    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True

    # Background fill for the text box shape
    if bg_color is not None:
        fill = txBox.fill
        fill.solid()
        fill.fore_color.rgb = bg_color

    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(font_size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color

    return txBox


def add_bullet_box(slide, lines, left, top, width, height,
                   font_size=13, bg_color=None, text_color=BLACK):
    """Add a text box with multiple bullet lines."""
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True

    if bg_color is not None:
        txBox.fill.solid()
        txBox.fill.fore_color.rgb = bg_color

    for i, line in enumerate(lines):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.alignment = PP_ALIGN.LEFT
        run = p.add_run()
        run.text = line
        run.font.size = Pt(font_size)
        run.font.color.rgb = text_color
        # small spacing between bullets
        p.space_before = Pt(2)
        p.space_after = Pt(2)

    return txBox


def add_rounded_rect_with_text(slide, text, left, top, width, height,
                                bg_color, text_color=WHITE,
                                font_size=15, bold=True):
    """Add a rectangle shape with centered text (used for section headers)."""
    from pptx.util import Pt
    shape = slide.shapes.add_shape(
        1,  # MSO_SHAPE_TYPE.RECTANGLE = 1
        left, top, width, height
    )
    shape.fill.solid()
    shape.fill.fore_color.rgb = bg_color
    shape.line.color.rgb = bg_color  # no visible border

    tf = shape.text_frame
    tf.word_wrap = False
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    run = p.add_run()
    run.text = text
    run.font.size = Pt(font_size)
    run.font.bold = bold
    run.font.color.rgb = text_color

    return shape


prs = Presentation(PPTX_PATH)
slide_layout = prs.slide_layouts[6]  # blank layout
slide = prs.slides.add_slide(slide_layout)

W = prs.slide_width    # ~12,192,000 EMU  (33.9 cm)
H = prs.slide_height   # ~6,858,000 EMU  (19.1 cm)

margin = Cm(0.5)
col_gap = Cm(0.4)

# ------------------------------------------------------------------
# TITLE BAR
# ------------------------------------------------------------------
title_height = Cm(1.5)
title_box = slide.shapes.add_shape(1, 0, 0, W, title_height)
title_box.fill.solid()
title_box.fill.fore_color.rgb = RGBColor(0x1F, 0x49, 0x7D)  # dark blue
title_box.line.color.rgb = RGBColor(0x1F, 0x49, 0x7D)

tf = title_box.text_frame
tf.word_wrap = False
p = tf.paragraphs[0]
p.alignment = PP_ALIGN.CENTER
run = p.add_run()
run.text = "Adaptation du modèle Fornari : ce qu'on garde / ce qu'on remplace"
run.font.size = Pt(20)
run.font.bold = True
run.font.color.rgb = WHITE

# ------------------------------------------------------------------
# LAYOUT: two columns
# ------------------------------------------------------------------
body_top    = title_height + Cm(0.3)
note_height = Cm(1.6)
note_top    = H - margin - note_height
body_height = note_top - body_top - Cm(0.3)

col_width = (W - 2 * margin - col_gap) / 2
left_col_x  = margin
right_col_x = margin + col_width + col_gap

header_height = Cm(1.1)
content_top   = body_top + header_height + Cm(0.15)
content_height = body_height - header_height - Cm(0.15)

# ------------------------------------------------------------------
# LEFT COLUMN HEADER  – green
# ------------------------------------------------------------------
add_rounded_rect_with_text(
    slide,
    text="✅  Ce qu'on garde (Fornari)",
    left=left_col_x, top=body_top,
    width=col_width, height=header_height,
    bg_color=GREEN_HEADER,
    text_color=WHITE,
    font_size=15, bold=True
)

# LEFT COLUMN CONTENT
left_bullets = [
    "• ODEs de l'hématopoïèse (compartiments 1–9)",
    "",
    "• Baselines :  Plt₀,  Neut₀,  Ret₀,  RBC₀",
    "",
    "• MTT_Plt = 168 h,  MTT_Neut = 210 h",
    "",
    "• Feedbacks γ  (mécanisme de régulation)",
    "",
    "• Grades NCI-CTCAE (toxicité hématologique)",
]

left_content = slide.shapes.add_textbox(
    left_col_x + Cm(0.3), content_top,
    col_width - Cm(0.6), content_height
)
left_content.fill.solid()
left_content.fill.fore_color.rgb = LIGHT_GREEN_BG
tf_l = left_content.text_frame
tf_l.word_wrap = True

for i, line in enumerate(left_bullets):
    p = tf_l.paragraphs[0] if i == 0 else tf_l.add_paragraph()
    p.alignment = PP_ALIGN.LEFT
    run = p.add_run()
    run.text = line
    run.font.size = Pt(13)
    run.font.color.rgb = RGBColor(0x1A, 0x5C, 0x1A) if line.startswith("•") else BLACK
    run.font.bold = line.startswith("•")
    p.space_before = Pt(1)

# add a subtle border / background rectangle behind left content
left_bg = slide.shapes.add_shape(
    1, left_col_x, content_top, col_width, content_height
)
left_bg.fill.solid()
left_bg.fill.fore_color.rgb = LIGHT_GREEN_BG
left_bg.line.color.rgb = GREEN_HEADER
left_bg.line.width = Pt(1)
# move behind the text box
left_bg.shape_id  # just to reference; z-order handled below

# We'll re-add text directly on this shape
tf_bg = left_bg.text_frame
tf_bg.word_wrap = True

# clear and re-populate
from pptx.util import Pt as PtU
for i, line in enumerate(left_bullets):
    p = tf_bg.paragraphs[0] if i == 0 else tf_bg.add_paragraph()
    p.alignment = PP_ALIGN.LEFT
    run = p.add_run()
    run.text = line
    run.font.size = PtU(13.5)
    run.font.color.rgb = RGBColor(0x1A, 0x5C, 0x1A) if line.startswith("•") else BLACK
    run.font.bold = line.startswith("•")
    p.space_before = PtU(2)
    p.space_after  = PtU(2)

# Remove the plain textbox (we'll use the shape instead)
sp = left_content._element
sp.getparent().remove(sp)

# ------------------------------------------------------------------
# RIGHT COLUMN HEADER  – dark red
# ------------------------------------------------------------------
add_rounded_rect_with_text(
    slide,
    text="❌  Ce qu'on remplace (ADC)",
    left=right_col_x, top=body_top,
    width=col_width, height=header_height,
    bg_color=RED_HEADER,
    text_color=WHITE,
    font_size=15, bold=True
)

# RIGHT COLUMN CONTENT
right_bullets = [
    "• PK carboplatine  →  Cavg fixe  (pas de modèle PK)",
    "",
    "• Damage (adducts ADN)  →  Emax :",
    "   Effect = Cavg / (IC50 + Cavg)",
    "",
    "• Slope / δ_Plt  →  à recalibrer avec données ADC",
]

right_bg = slide.shapes.add_shape(
    1, right_col_x, content_top, col_width, content_height
)
right_bg.fill.solid()
right_bg.fill.fore_color.rgb = LIGHT_RED_BG
right_bg.line.color.rgb = RED_HEADER
right_bg.line.width = Pt(1)

tf_rbg = right_bg.text_frame
tf_rbg.word_wrap = True

for i, line in enumerate(right_bullets):
    p = tf_rbg.paragraphs[0] if i == 0 else tf_rbg.add_paragraph()
    p.alignment = PP_ALIGN.LEFT
    run = p.add_run()
    run.text = line
    is_bullet = line.startswith("•")
    is_formula = line.strip().startswith("Effect")
    run.font.size = PtU(13.5)
    run.font.color.rgb = (
        RGBColor(0x8B, 0x00, 0x00) if is_bullet else
        RGBColor(0x44, 0x44, 0x44)
    )
    run.font.bold = is_bullet
    run.font.italic = is_formula
    p.space_before = PtU(2)
    p.space_after  = PtU(2)

# ------------------------------------------------------------------
# NOTE BOX at the bottom
# ------------------------------------------------------------------
note_shape = slide.shapes.add_shape(
    1, margin, note_top, W - 2 * margin, note_height
)
note_shape.fill.solid()
note_shape.fill.fore_color.rgb = NOTE_BG
note_shape.line.color.rgb = RGBColor(0xCC, 0xAA, 0x00)
note_shape.line.width = Pt(1.5)

tf_note = note_shape.text_frame
tf_note.word_wrap = True
p_note = tf_note.paragraphs[0]
p_note.alignment = PP_ALIGN.CENTER
run_note = p_note.add_run()
run_note.text = (
    "Seul inconnu : Slope / Emax  →  estimation par analyse de sensibilité"
)
run_note.font.size = PtU(14)
run_note.font.bold = True
run_note.font.italic = True
run_note.font.color.rgb = RGBColor(0x7F, 0x60, 0x00)

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------
prs.save(PPTX_PATH)
print(f"Done. Total slides now: {len(prs.slides)}")
print(f"New slide index: {len(prs.slides) - 1} (0-based), slide #{len(prs.slides)} (1-based)")
