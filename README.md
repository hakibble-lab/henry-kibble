# henry-kibble

R data analysis project.

## Layout

- `R/` — analysis scripts
- `data/` — input data (gitignored contents aside from `.gitkeep`)
- `output/` — generated results/plots

## Workflow

This project is developed jointly: Claude pushes changes to the
`claude/sleepy-pasteur-tcuj78` branch, and you pull them into RStudio.

```r
# one-time setup in RStudio: File > New Project > Version Control > Git
# paste the repo URL, or if already cloned:
usethis::pr_fetch()  # or plain git:
```

```sh
git fetch origin claude/sleepy-pasteur-tcuj78
git checkout claude/sleepy-pasteur-tcuj78
git pull
```

After that, `git pull` in RStudio's Git pane (or terminal) brings in the
latest changes.
