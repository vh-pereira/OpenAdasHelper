from pathlib import Path

def find_repo_root(start=None):
    p = Path(start or __file__).resolve()
    for d in [p] + list(p.parents):
        if (d / '.git').exists() or (d / 'pyproject.toml').exists() or (d / 'README.md').exists():
            return d
    return Path.cwd()