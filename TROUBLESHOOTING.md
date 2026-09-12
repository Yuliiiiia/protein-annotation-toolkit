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
version that isn't present. Use the conda-forge build instead, which is
newer and more actively maintained:
```bash
conda remove dssp -y
conda install -c conda-forge dssp -y
```

## DSSP fails on `.cif` files: `Error while loading dictionary mmcif_pdbx ... Is a directory`

Even after installing conda-forge's DSSP (4.4.11), running it against an
mmCIF file can fail. Newer `mkdssp` builds parse mmCIF via `libcifpp`,
which expects local chemical-component-dictionary data that usually isn't
configured by default — DSSP itself works fine, it's specifically the
mmCIF dictionary lookup that's broken.

**This is now handled automatically by the toolkit** (`analyze.py`,
`run_dssp()`): if DSSP fails on the mmCIF input, the code converts the
already-parsed structure to a temporary legacy `.pdb` file (which doesn't
trigger the dictionary lookup) and retries once before giving up. No
action needed on your end — if DSSP is installed and working at all,
you should get real helix/sheet/coil percentages.

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
