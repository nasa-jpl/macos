-- style_map.lua — used with `pandoc -f docx+styles` when converting the
-- styled manual .docx to Markdown (one-time migration; kept for provenance).
--
-- The styled docx tags every paragraph with a semantic style
-- (Body Text, Code Block, List Paragraph, Figure Caption, ...).  The
-- docx reader wraps each such paragraph in a Div with
-- custom-style="<name>", and indentation makes most of them arrive as
-- BlockQuotes.  This filter maps them to native Markdown constructs:
--
--   Code Block      -> fenced CodeBlock (adjacent blocks merged)
--   everything else -> unwrapped plain paragraphs (BlockQuote removed)
--   empty Headers   -> dropped
--   Header attrs    -> stripped (pandoc regenerates identical auto-ids)

local function inlines_to_text(inls)
  local parts = {}
  for _, il in ipairs(inls) do
    if il.t == "Str" then
      table.insert(parts, il.text)
    elseif il.t == "Space" then
      table.insert(parts, " ")
    elseif il.t == "LineBreak" or il.t == "SoftBreak" then
      table.insert(parts, "\n")
    elseif il.content ~= nil then
      -- Strong/Emph/Span/etc: keep only the text
      table.insert(parts, inlines_to_text(il.content))
    end
  end
  return table.concat(parts)
end

local function gather_code_lines(blocks, out)
  for _, b in ipairs(blocks) do
    if b.t == "Para" or b.t == "Plain" then
      table.insert(out, inlines_to_text(b.content))
    elseif b.t == "BlockQuote" or b.t == "Div" then
      gather_code_lines(b.content, out)
    elseif b.t == "CodeBlock" then
      table.insert(out, b.text)
    end
  end
end

function Div(el)
  local style = el.attributes["custom-style"]
  if not style then return nil end
  if style == "Code Block" then
    -- Images embedded in Code Block paragraphs (diagram fragments) must
    -- not be swallowed by the text extraction: collect them separately
    -- and emit them as their own paragraph ahead of the code block.
    local images = {}
    pandoc.Div(el.content):walk{
      Image = function(im) table.insert(images, im) end }
    local lines = {}
    gather_code_lines(el.content, lines)
    local out = pandoc.Blocks{}
    if #images > 0 then
      out:insert(pandoc.Para(images))
    end
    local text = table.concat(lines, "\n")
    if text:match("%S") or #images == 0 then
      out:insert(pandoc.CodeBlock(text))
    end
    return out
  end
  -- Any other custom style: unwrap the Div and any top-level BlockQuote.
  local out = pandoc.Blocks{}
  for _, b in ipairs(el.content) do
    if b.t == "BlockQuote" then
      out:extend(b.content)
    else
      out:insert(b)
    end
  end
  return out
end

function Header(el)
  if #el.content == 0 or pandoc.utils.stringify(el.content):match("^%s*$") then
    return pandoc.Blocks{}          -- drop empty headings (docx debris)
  end
  el.attr = pandoc.Attr()           -- drop {#id .unnumbered} clutter
  return el
end

-- Merge adjacent CodeBlocks: the docx stores each transcript line as its
-- own Code Block paragraph; visually they are one block.
function Blocks(blocks)
  local out = pandoc.Blocks{}
  for _, b in ipairs(blocks) do
    local last = out[#out]
    if b.t == "CodeBlock" and last and last.t == "CodeBlock" then
      out[#out] = pandoc.CodeBlock(last.text .. "\n" .. b.text)
    else
      out:insert(b)
    end
  end
  return out
end
