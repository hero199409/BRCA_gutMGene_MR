import os
"""Build same-width Old-vs-New comparison previews for Fig 6 and Fig 8, plus Fig 1 draft preview."""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

PROJ_ROOT = Path(os.environ.get("BRCA_PROJECT_ROOT", Path(__file__).resolve().parents[2]))
old_dir = PROJ_ROOT / "03_figures" / "Publication_PlanA_CNS"
new_dir = PROJ_ROOT / "05_revision_v2" / "03_figures"
prv_dir = new_dir / "previews"
prv_dir.mkdir(parents=True, exist_ok=True)


def make_comparison(old_name, new_name, out_name, title_old, title_new, target_col_w=1750):
    with Image.open(old_dir / old_name) as im_old_raw, Image.open(new_dir / new_name) as im_new_raw:
        im_old = im_old_raw.convert("RGB")
        im_new = im_new_raw.convert("RGB")

    # Normalize both to identical physical display width (both were authored for 180 mm)
    old_h = round(im_old.height * target_col_w / im_old.width)
    new_h = round(im_new.height * target_col_w / im_new.width)
    im_old = im_old.resize((target_col_w, old_h), Image.Resampling.LANCZOS)
    im_new = im_new.resize((target_col_w, new_h), Image.Resampling.LANCZOS)

    pad = 40
    header_h = 72
    canvas_w = target_col_w * 2 + pad * 3
    canvas_h = max(old_h, new_h) + header_h + pad * 2
    canvas = Image.new("RGB", (canvas_w, canvas_h), "white")
    draw = ImageDraw.Draw(canvas)

    try:
        font_hdr = ImageFont.truetype("arialbd.ttf", 34)
    except OSError:
        font_hdr = ImageFont.load_default()

    # Headers
    draw.rectangle([pad, 16, pad + target_col_w, header_h], fill=(241, 245, 249), outline=(203, 213, 225), width=2)
    draw.text((pad + 18, 28), title_old, fill=(51, 65, 85), font=font_hdr)

    x_new = pad * 2 + target_col_w
    draw.rectangle([x_new, 16, x_new + target_col_w, header_h], fill=(254, 242, 242), outline=(200, 62, 43), width=2)
    draw.text((x_new + 18, 28), title_new, fill=(142, 31, 19), font=font_hdr)

    # Paste both figures aligned at top
    y_img = header_h + 16
    canvas.paste(im_old, (pad, y_img))
    canvas.paste(im_new, (x_new, y_img))
    draw.rectangle([pad - 1, y_img - 1, pad + target_col_w, y_img + old_h], outline=(203, 213, 225), width=2)
    draw.rectangle([x_new - 1, y_img - 1, x_new + target_col_w, y_img + new_h], outline=(200, 62, 43), width=2)

    out_path = prv_dir / out_name
    canvas.save(out_path, dpi=(300, 300))
    print(f"Saved {out_path} ({canvas_w}x{canvas_h} px)")


if __name__ == "__main__":
    make_comparison(
        "Fig6_Tumor_scRNA_and_HPA_Pathology.png",
        "Fig6_Tumor_scRNA_and_HPA_Pathology.png",
        "Fig6_Old_vs_New_Comparison.png",
        "OLD Fig 6 (180 x 232 mm | Publication_PlanA_CNS)",
        "NEW Fig 6 (180 x 224 mm | 05_revision_v2 Clean)",
    )
    make_comparison(
        "Fig8_Therapy_Response_ISPY2_and_GDSC2.png",
        "Fig8_Therapy_Response_ISPY2_and_GDSC2.png",
        "Fig8_Old_vs_New_Comparison.png",
        "OLD Fig 8 (180 x 253 mm | Publication_PlanA_CNS)",
        "NEW Fig 8 (180 x 216 mm | 05_revision_v2 Clean)",
    )
