-- number_sections.lua -- the manual's heading numbers in the 3.2 pattern, a page break before every
-- Section / Appendix, and CROSS-REFERENCES resolved at build time (2026-10-18, Dave's edits 1, 2 and
-- the numbering question).
--
--   "## SECTION 4: Title"   -> unchanged text, sets the chapter number 4, page break before
--   "## APPENDIX A: Title"  -> unchanged, chapter "A", page break before
--   "## REFERENCES" etc.    -> unchanged, no chapter (nothing below it is numbered), page break
--   "### Title"             -> "4.1 Title", "4.2 Title", ...
--   "#### Title"            -> "4.1.1 Title", ...  (a #### with no ### above it in the chapter
--                              is numbered "4.1", "4.2", ... as if it were the second level)
--   level-1 headings (the front matter) and anything outside a chapter are left alone.
--
-- A cross-reference is an ordinary internal link whose text holds a "?":
--   "see [Section ?](#conic-surfaces)"  ->  "see Section 4.5.2", still a link to the heading
--   "[?](#conic-surfaces)"              ->  "4.5.2"
-- The target is the heading's pandoc id (its text lower-cased, spaces to hyphens, punctuation
-- dropped -- `pandoc -t native` shows them).  A "?" whose target is not a numbered heading is left
-- as "?" and reported on stderr, so a stale reference is visible in every build.
-- The numbers are generated at build time, so the sources carry none and nothing goes stale.

local numbers = {}   -- heading id -> "4.5.2"

local function pagebreak()
  if FORMAT:match("docx") then
    return pandoc.RawBlock("openxml", '<w:p><w:r><w:br w:type="page"/></w:r></w:p>')
  elseif FORMAT:match("html") then
    return pandoc.RawBlock("html", '<div style="page-break-before:always"></div>')
  end
  return nil
end

-- pass 1: assign the numbers (in document order) and put them into the headings
local function number_headers(doc)
  local chap, s2, s3 = nil, 0, 0
  local out = {}
  for _, el in ipairs(doc.blocks) do
    if el.t == "Header" then
      if el.level == 2 then
        local txt = pandoc.utils.stringify(el.content)
        chap = txt:match("^SECTION%s+(%d+)") or txt:match("^APPENDIX%s+(%a)")
        s2, s3 = 0, 0
        local pb = pagebreak()
        if pb then table.insert(out, pb) end
      elseif el.level == 3 and chap then
        s2 = s2 + 1; s3 = 0
        local n = string.format("%s.%d", chap, s2)
        numbers[el.identifier] = n
        el.content:insert(1, pandoc.Str(n)); el.content:insert(2, pandoc.Space())
      elseif el.level == 4 and chap then
        s3 = s3 + 1
        local n = (s2 == 0) and string.format("%s.%d", chap, s3) or string.format("%s.%d.%d", chap, s2, s3)
        numbers[el.identifier] = n
        el.content:insert(1, pandoc.Str(n)); el.content:insert(2, pandoc.Space())
      end
    end
    table.insert(out, el)
  end
  doc.blocks = out
  return doc
end

-- pass 2: resolve "?" in internal links
local function resolve(link)
  local target = link.target
  if not target:match("^#") then return nil end
  local id = target:sub(2)
  local n = numbers[id]
  local out = pandoc.List()
  for _, il in ipairs(link.content) do
    if il.t == "Str" and il.text:find("%?") then
      if n then
        out:insert(pandoc.Str((il.text:gsub("%?", n))))
      else
        io.stderr:write("number_sections.lua: no numbered heading for #" .. id .. "\n")
        out:insert(il)
      end
    else
      out:insert(il)
    end
  end
  link.content = out
  return link
end

function Pandoc(doc)
  doc = number_headers(doc)
  return doc:walk({ Link = resolve })
end
