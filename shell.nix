# power.R needs only base R + stats. For the joint-model fallback (README.md),
# use rWrapper so rPackages are visible to Rscript, e.g.
#   (pkgs.rWrapper.override { packages = with pkgs.rPackages; [ lme4 simr ]; })
{ pkgs ? import <nixpkgs> { } }:
pkgs.mkShell {
  packages = [ pkgs.R ];
}
