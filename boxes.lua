function Div(el)
  local c = el.classes[1]
  if (c == "retenir" or c == "attention") and FORMAT:match("latex") then
    return { pandoc.RawBlock('latex', '\\begin{' .. c .. '}'), el, pandoc.RawBlock('latex', '\\end{' .. c .. '}') }
  end
end
