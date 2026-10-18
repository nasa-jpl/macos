-- number_sections.lua -- the manual's heading numbers in the 3.2 pattern, and a page break
-- before every Section / Appendix (2026-10-18, Dave's end-product edits 1 and 2).
--
--   "## SECTION 4: Title"   -> unchanged text, sets the chapter number 4, page break before
--   "## APPENDIX A: Title"  -> unchanged, chapter "A", page break before
--   "## REFERENCES" etc.    -> unchanged, no chapter (nothing below it is numbered), page break
--   "### Title"             -> "4.1 Title", "4.2 Title", ...
--   "#### Title"            -> "4.1.1 Title", ...  (a #### with no ### above it in the chapter
--                              is numbered "4.1", "4.2", ... as if it were the second level)
--   level-1 headings (the front matter) and anything outside a chapter are left alone.
-- The numbers are put into the heading text, so the table of contents carries them too.

local chap, s2, s3 = nil, 0, 0

local function pagebreak()
  if FORMAT:match("docx") then
    return pandoc.RawBlock("openxml", '<w:p><w:r><w:br w:type="page"/></w:r></w:p>')
  elseif FORMAT:match("html") then
    return pandoc.RawBlock("html", '<div style="page-break-before:always"></div>')
  end
  return nil
end

function Header(el)
  if el.level == 2 then
    local txt = pandoc.utils.stringify(el.content)
    local n = txt:match("^SECTION%s+(%d+)")
    local a = txt:match("^APPENDIX%s+(%a)")
    chap = n or a   -- nil for REFERENCES and the like
    s2, s3 = 0, 0
    local pb = pagebreak()
    if pb then return { pb, el } end
    return el
  elseif el.level == 3 and chap then
    s2 = s2 + 1; s3 = 0
    el.content:insert(1, pandoc.Str(string.format("%s.%d", chap, s2)))
    el.content:insert(2, pandoc.Space())
    return el
  elseif el.level == 4 and chap then
    if s2 == 0 then
      s3 = s3 + 1
      el.content:insert(1, pandoc.Str(string.format("%s.%d", chap, s3)))
    else
      s3 = s3 + 1
      el.content:insert(1, pandoc.Str(string.format("%s.%d.%d", chap, s2, s3)))
    end
    el.content:insert(2, pandoc.Space())
    return el
  end
  return el
end
