# power.R needs only base R + stats; power.py needs numpy + scipy.
# For the joint-model fallback (README.md), use rWrapper so rPackages are
# visible to Rscript, e.g.
#   (pkgs.rWrapper.override { packages = with pkgs.rPackages; [ lme4 ]; })
{ pkgs ? import <nixpkgs> { } }:
pkgs.mkShell {
  packages = [
    pkgs.R
    (pkgs.python3.withPackages (ps: [ ps.numpy ps.scipy ]))
  ];
}
