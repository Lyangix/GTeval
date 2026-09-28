# Publishing this code-only repository

The intended Git repository root is this release folder, not the parent research
directory. The parent contains study inputs and generated analyses and should not
be uploaded as a whole. The release folder has no Git history and no dependency
on files outside its own directory.

## Make a portable archive

Python 3 is needed only for the optional archive exporter:

```sh
python3 tools/build_release.py
```

This creates `../github_release.zip` from exactly the files in `MANIFEST.txt`.
The exporter never recursively copies the research directory and never includes
`outputs/`. It rejects symlinks, files outside the release root, and binary inputs.
It refuses to overwrite an existing archive; use `--output /path/to/new-release.zip`
for another build. The manifest is the list of files to review before publishing.

Extract the archive into a separate working location. Run the example and tests
there to verify that no private file paths are needed:

```sh
cd github_release
Rscript --vanilla examples/run_simulation.R
Rscript --vanilla tests/run_tests.R
```

## Create the GitHub repository

Choose the repository name, author information, citation details, and a code
license before public distribution. No license or paper citation has been
invented for this draft. Add the chosen `LICENSE` and any citation file to both
the manifest and `.gitignore` allowlist if they should be exported and tracked.
If adding an extension not currently supported by the exporter, update its
allowed source-file extensions as part of that change.

Inside the extracted release folder:

```sh
git init -b main
git add .gitignore README.md MANIFEST.txt R docs examples tests tools
git diff --cached --stat
git diff --cached
git commit -m "Add PRS variance decomposition and transfer learning example"
```

After creating an empty GitHub repository, replace `OWNER` and `REPOSITORY` with
your chosen values:

```sh
git remote add origin https://github.com/OWNER/REPOSITORY.git
git push -u origin main
```

The supplied `.gitignore` ignores everything except the named public files.
Adding a new intended public file requires updating this allowlist and the
manifest. Private dosages, phenotypes, IDs, real PRS weights, fitted study models,
subject-level predictions, and generated study results do not belong in this
release. Ignoring files does not remove anything already committed, so begin
with the separate code-only folder and inspect the staged diff.

The archive exporter packages the listed contents; it cannot determine whether
someone later embeds real data inside a listed script or document. Keep examples
synthetic when extending this repository.
