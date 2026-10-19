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
  return pandoc.read(text, "markdown").blocks
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
