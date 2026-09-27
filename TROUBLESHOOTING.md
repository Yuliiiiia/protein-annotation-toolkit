# Troubleshooting

Real issues hit while setting this project up on WSL (Windows Subsystem for
Linux), and how each was resolved. Kept here because these are common
enough on fresh WSL installs that they're worth documenting rather than
rediscovering.

## `unzip: command not found`

WSL's minimal base image doesn't ship `unzip`. Either install it:
```bash
sudo apt install unzip
```
or skip it entirely with Python (no install needed):
```bash
python3 -c "import zipfile; zipfile.ZipFile('protein-toolkit.zip').extractall()"
```

## `sudo` password not accepted / forgotten

If you don't know your WSL user's password, `sudo` is a dead end until you
reset it. From an **admin PowerShell** on the Windows side:
```powershell
wsl -u root
```
then inside that root shell:
```bash
passwd yourusername
exit
```

**Better fix, avoids sudo entirely:** if you don't strictly need system-wide
installs, use [Miniconda](https://docs.conda.io/) instead — it installs
everything under your home directory with no admin rights required at all.
This is the route this project's setup ended up taking:
```bash
cd ~
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
source $HOME/miniconda3/bin/activate
pip install biopython pytest requests
```

## `pip3: command not found` / `No module named pip`

A minimal `python3` install on WSL sometimes doesn't include `pip`. Try, in
order:
```bash
python3 -m pip --version        # check if it's there but not linked
python3 -m ensurepip --user     # try installing pip for your user only
```
If neither works, the Miniconda route above sidesteps the problem
completely — conda ships its own `pip`.

## `CondaToSNonInteractiveError: Terms of Service have not been accepted`

Newer conda versions require explicitly accepting channel ToS before
installing from `defaults`:
```bash
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
```

## `mkdssp: error while loading shared libraries: libboost_thread.so.1.73.0`

The `salilab` conda channel's DSSP build (3.0.0) links against an old Boost
version that isn't present. Switching to conda-forge's build (4.4.11) gets
past this specific error — but keep reading, because it runs into a worse,
unfixable one next.

## DSSP fails on `.cif` files: `Error while loading dictionary mmcif_pdbx ... Is a directory`, and it's not fixable

This one took a lot of digging, so the full story is worth recording.

Every conda-forge build of `mkdssp` (4.x) that was tried — including
setting the `LIBCIFPP_DATA_DIR` environment variable explicitly to the
directory containing the dictionary file (confirmed present and readable
on disk), and even passing the exact file path directly via
`--mmcif-dictionary` — still failed with the same
`basic_filebuf::underflow error reading the file: Is a directory` error.

Tracing this back to [conda-forge/dssp-feedstock#4](https://github.com/conda-forge/dssp-feedstock/pull/4)
confirmed it's not a local misconfiguration: it's a genuine, unresolved bug
in how conda's binary-relocation mechanism interacts with `libcifpp` (the
C++ library `mkdssp` uses to load its dictionary). Even the library's own
maintainer, debugging the exact same error in that thread, couldn't fix
it and closed the PR unresolved:

> "Seems like the hack conda uses to make packages relocatable is not
> working with libcifpp. Until that is fixed... this feedstock is not
> going to work."

Installing DSSP via `apt` (Debian/Ubuntu's system package, which bundles
the `libcifpp-data` package correctly) does work — confirmed by testing it
directly. But that requires `sudo`, which isn't always available (see the
sudo section above), and doesn't help on Windows/WSL setups without a
working system package manager for it.

**The actual fix: this project no longer uses `mkdssp` at all.**
`analyze.py`'s `run_dssp()` now uses [`pydssp`](https://github.com/ShintaroMinami/PyDSSP)
— a pure NumPy re-implementation of the DSSP hydrogen-bond algorithm,
installed via plain `pip install pydssp`. No external binary, no compiled
dictionary files, no relocation issues. It's also the same backend
`MDAnalysis` uses as a DSSP alternative, so it's not an obscure choice.
The one tradeoff: `pydssp` pulls in `torch` as a dependency (even though
only its NumPy backend is used here), so the install is a few hundred MB
larger than a typical lightweight package — worth knowing in advance if
disk space is tight.

## `git commit` fails: `Please tell me who you are` / `empty ident name`

First-time git use on a machine needs an identity set once:
```bash
git config --global user.email "you@example.com"
git config --global user.name "Your Name"
```

## `git push` fails: `Password authentication is not supported for Git operations`

GitHub no longer accepts account passwords for git over HTTPS. Generate a
**Personal Access Token** instead:

1. [github.com/settings/tokens](https://github.com/settings/tokens) →
   **Generate new token (classic)**
2. Check the **`repo`** scope only (nothing else is needed for push/pull)
3. Set an expiration (90 days is a reasonable default)
4. Copy the token immediately — it's shown once

At the `git push` username/password prompt, use your GitHub username and
paste the **token** as the password.

## `<placeholder>` text typed literally into a command

A few commands in setup instructions use `<your-username>` or similar as a
**placeholder** meaning "put your real value here" — the angle brackets
themselves are never meant to be typed. Bash interprets a leading `<` as
input redirection, which produces confusing "No such file or directory"
errors that have nothing to do with the actual command.
