-- pins.lua -- a number the prose quotes is written %pin:NAME% and becomes the value the chapter's
-- driver produced (src/transcripts/pins.tsv, written by tools/regen_mmacos_transcripts.py after the
-- pin held against tools/pins/<driver>.tsv).  An unknown NAME STOPS the build: a quoted number
-- cannot silently go stale or be typed.

local pins = nil

local function load()
  pins = {}
  local f = io.open(pandoc.path.join({ "src", "transcripts", "pins.tsv" }), "r")
  if not f then
    error("pins.lua: src/transcripts/pins.tsv is missing -- run `make transcripts`")
  end
  for line in f:lines() do
    local name, shown = line:match("^([%w_%.%-]+)\t(.*)$")
    if name then pins[name] = shown end
  end
  f:close()
end

function Str(el)
  if not el.text:find("%pin:", 1, true) then return nil end
  if not pins then load() end
  local out = el.text:gsub("%%pin:([%w_%.%-]+)%%", function(name)
    local v = pins[name]
    if not v then error("pins.lua: no pin named " .. name .. " -- run `make transcripts`") end
    return v
  end)
  return pandoc.Str(out)
end
