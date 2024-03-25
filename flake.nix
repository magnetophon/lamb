{
  description = "A flake that includes Faust";

  outputs = { self, nixpkgs }: {
    # Define a development shell that includes faust and other necessary tools
    devShell.x86_64-linux = nixpkgs.legacyPackages.x86_64-linux.mkShell {
      buildInputs = [
        nixpkgs.legacyPackages.x86_64-linux.faust2jack
      ];
    };
  };
}
