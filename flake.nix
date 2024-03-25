{
  description = "A flake that includes Faust";

  outputs = { self, nixpkgs }: {
    devShell.x86_64-linux = let
      pinnedPkgs = import (builtins.fetchTarball {
        url = "https://github.com/NixOS/nixpkgs/archive/b06025f1533a.tar.gz";
        sha256 = "sha256:1b8dim6xpcg3wyb0xa0w4h4m22npbzl2np822x4r7wiw5wnnzg5a";
      }) {
        system = "x86_64-linux";
        overlays = [ ];
        config = { };
      };
    in pinnedPkgs.mkShell {
      packages = [
        pinnedPkgs.faust2jack
        pinnedPkgs.x42-plugins
      ];
    };
  };
}
