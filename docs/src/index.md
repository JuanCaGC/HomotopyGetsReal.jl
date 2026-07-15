```@meta
CurrentModule = HomotopyGetsReal
```

# HomotopyGetsReal

```@eval
using Markdown
readme = read(joinpath(@__DIR__, "..", "..", "README.md"), String)
# Drop the leading H1 so Documenter's page title is not duplicated.
readme = replace(readme, r"^# HomotopyGetsReal\r?\n+" => ""; count = 1)
Markdown.parse(readme)
```

## Documentation pages

```@contents
Pages = ["api.md"]
Depth = 2
```
