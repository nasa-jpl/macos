-- figure_captions.lua -- the manual's figures the 3.2 way (Dave 2026-10-19): the image centered,
-- the "**FIGURE n** caption" paragraph centered beneath it.  Panel labels are drawn INTO the
-- engine-made images by tools/regen_engine_figs.py (under their panels), so the text carries none.
--   a paragraph that is one image          -> Div custom-style "Figure"        (centered)
--   a paragraph starting with "**FIGURE n**" -> Div custom-style "FigureCaption" (centered)
-- The two styles are in reference.docx (tools/make_reference_docx.py); the HTML gets a class and
-- tools/manual.css centers it.
local function wrap(block, style)
  return pandoc.Div({block}, pandoc.Attr("", {style}, {["custom-style"] = style}))
end

function Para(el)
  local c = el.content
  if #c == 1 and c[1].t == "Image" then
    return wrap(el, "Figure")
  end
  if #c >= 1 and c[1].t == "Strong" then
    local s = pandoc.utils.stringify(c[1])
    if s:match("^FIGURE %d+$") then
      return wrap(el, "FigureCaption")
    end
  end
  return nil
end
