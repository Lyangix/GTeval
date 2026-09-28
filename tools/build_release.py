#!/usr/bin/env python3
"""Export only the explicitly listed public source files; no directory globbing."""
import argparse
from pathlib import Path, PurePosixPath
from zipfile import ZIP_DEFLATED, ZipFile


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=root.parent / "github_release.zip")
    args = parser.parse_args()
    names = (root / "MANIFEST.txt").read_text(encoding="utf-8").splitlines()
    if not names or len(names) != len(set(names)):
        raise SystemExit("The public manifest must be nonempty and contain unique paths.")
    files = []
    for name in names:
        relative = PurePosixPath(name)
        if not name or relative.is_absolute() or ".." in relative.parts or "\\" in name:
            raise SystemExit("Invalid manifest path: " + name)
        path = root / name
        if path.is_symlink() or any(p.is_symlink() for p in path.parents if p != root.parent):
            raise SystemExit("Symlinks cannot be exported: " + name)
        if not path.is_file() or root not in path.resolve().parents:
            raise SystemExit("Missing or external manifest file: " + name)
        if name != ".gitignore" and path.suffix not in {".R", ".md", ".txt", ".py"}:
            raise SystemExit("Only public source and documentation may be exported: " + name)
        contents = path.read_bytes()
        if b"\x00" in contents or len(contents) > 1_000_000:
            raise SystemExit("Unexpected binary or oversized file: " + name)
        contents.decode("utf-8")
        files.append((name, contents))
    # Exclusive creation avoids silently overwriting an existing release archive.
    with ZipFile(args.output, "x", compression=ZIP_DEFLATED) as archive:
        for name, contents in files:
            archive.writestr("github_release/" + name, contents)
    print("Created", args.output, "with", len(files), "public source/documentation files.")


if __name__ == "__main__":
    main()
