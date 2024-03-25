{ pkgs ? import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/b06025f1533a.tar.gz") {} }:

pkgs.mkShell {
  buildInputs = [
    pkgs.faust2jack
    pkgs.x42-plugins
  ];
}
