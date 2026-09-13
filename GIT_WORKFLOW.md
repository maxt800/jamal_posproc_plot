# Local version control

The repository is this folder, including CODEX_CONTEXT.md and jamal_dashboard.
The main branch starts with the validated v25.3 checkpoint (tag v25.3).
Source, tests, documentation and synthetic inputs are tracked. Generated reports,
caches, Python environments and release ZIPs are excluded. Fixture bytes are
preserved without Git newline conversion so their checksums remain reproducible.

From this folder in PowerShell:

```powershell
git status
git diff
git log --oneline --decorate
git diff v25.3 -- jamal_dashboard/distributions.js
```

To save future changes after testing:

```powershell
git add jamal_dashboard CODEX_CONTEXT.md
git commit -m "Describe the change"
```

This is a local repository. No remote is configured and nothing is published.

## Windows ownership check

Git was initialized by the Codex sandbox account. If your terminal reports
"dubious ownership", use this repository-specific, per-command option:

```powershell
git -c safe.directory=C:/Users/User/Documents/ChatGPT/CFD_Plotter_Results status
```

Use the same `-c` option before other Git subcommands when needed. It changes
neither folder permissions nor persistent Git configuration.
