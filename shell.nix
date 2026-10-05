# numpy, scipy, matplotlib run the analyses; pandas and statsmodels run checks/.
{ pkgs ? import <nixpkgs> { } }:
pkgs.mkShell {
  packages = [ (pkgs.python3.withPackages (ps: [ ps.numpy ps.scipy ps.matplotlib ps.pandas ps.statsmodels ])) ];
}
