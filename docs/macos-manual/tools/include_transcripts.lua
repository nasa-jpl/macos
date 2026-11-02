-- include_transcripts.lua -- a chapter line `<!-- transcript: NAME -->` becomes the generated file
-- src/transcripts/NAME.md (tools/regen_transcripts.py writes them: the example listings verbatim
-- and the CLI dialogs from tools/transcripts/*.jou; BRIEF_to_manual_regen.md, Dave 2026-10-19).
-- The marker is the only thing a chapter keeps.  A marker whose file is missing STOPS the build:
-- a listing cannot silently drop out of the manual.  Runs before number_sections.lua (the included
-- blocks carry no headings, so the numbering is unchanged).

local function include(name)
  local path = pandoc.path.join({ "src", "transcripts", name .. ".md" })
  local f = io.open(path, "r")
  if not f then
    error("include_transcripts: " .. path .. " is missing -- run `make transcripts`")
  end
  local text = f:read("a")
  f:close()
  local blocks = pandoc.read(text, "markdown").blocks
  return blocks:walk({ CodeBlock = wrap_code })
end

-- Listings are included verbatim, but the page is 104 columns of Courier New at the code size and
-- Word breaks a longer line wherever it falls.  Wrap a long line at the last blank before the limit
-- and indent the continuation under the line's own indent (Dave's hand edits of 2026-11-02 are the
-- model).  Only the rendering changes: the transcript files stay as the tools printed them.
local WRAP = 100
function wrap_code(el)
  local out = {}
  for line in (el.text .. "\n"):gmatch("(.-)\n") do
    local indent = line:match("^(%s*)")
    local cont = indent .. "    "
    local rest = line
    while #rest > WRAP do
      local cut = nil
      for i = WRAP, #indent + 20, -1 do
        if rest:sub(i, i) == " " then cut = i break end
      end
      if not cut then break end
      table.insert(out, rest:sub(1, cut - 1))
      rest = cont .. rest:sub(cut + 1)
    end
    table.insert(out, rest)
  end
  el.text = table.concat(out, "\n")
  return el
end

function RawBlock(el)
  if el.format == "html" then
    local name = el.text:match("^<!%-%-%s*transcript:%s*([%w_%-%.]+)%s*%-%->%s*$")
    if name then
      return include(name)
    end
  end
  return nil
end
