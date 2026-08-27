{
  description = "AlterPCB Transmission Line Simulator";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/26.05";

  outputs = {
    self,
    nixpkgs,
  }: let
    # Generate outputs for each supported system. Nix commands select the
    # output matching the current host system.
    supportedSystems = [
      "aarch64-linux"
      "x86_64-linux"
    ];
    forAllSystems = nixpkgs.lib.genAttrs supportedSystems;
  in {
    # `nix build` uses this derivation
    # packages can be reused later by other derivations
    packages = forAllSystems (system: let
      pkgs = nixpkgs.legacyPackages.${system};
    in {
      default = pkgs.stdenv.mkDerivation {
        pname = "alterpcb-tlinesim";
        version = "0.0.0";

        src = self;

        # Tools used while building
        # `qt5.qmake` also installs a Nix hook that makes `qmake` the automatic configurePhase
        # `wrapQtAppsHook` wraps installed GUI executables during fixupPhase
        nativeBuildInputs = with pkgs; [
          pkg-config
          qt5.qmake
          qt5.wrapQtAppsHook
        ];

        # Libraries and headers needed to compile and run the program.
        buildInputs = with pkgs; [
          eigen
          qt5.qtbase
        ];

        # Default unpackPhase: copy the flake source tree (`src = self`) into a
        # temporary build directory.

        # preConfigure hook: enter the directory containing the qmake project
        # before qmake runs. The later buildPhase also runs in this directory.
        preConfigure = ''
          cd src
        '';

        # qmake configurePhase: generate the Makefile with PREFIX=$out.
        # This also compiles DATADIR as $out/share/alterpcb-tlinesim.

        # Implied buildPhase: run `make`, which compiles and links C++

        # Defined installPhase: copy the results into the Nix store
        installPhase = ''
          runHook preInstall

          install -Dm755 alterpcb-tlinesim "$out/bin/alterpcb-tlinesim"
          mkdir -p "$out/share/alterpcb-tlinesim"
          cp -r ../data/. "$out/share/alterpcb-tlinesim/"

          runHook postInstall
        '';

        # Implied fixupPhase: wrap the Qt executable so it can find Qt plugins

        meta = {
          description = "Field solver for PCB transmission lines";
          homepage = "https://github.com/alterpcb/alterpcb-tlinesim";
          license = pkgs.lib.licenses.gpl3Plus;
          mainProgram = "alterpcb-tlinesim";
          platforms = supportedSystems;
        };
      };
    });

    # An unqualified `nix run` selects this app definition. Its program path
    # references the package above, so Nix builds it if needed and runs it.
    apps = forAllSystems (system: {
      default = {
        type = "app";
        program = "${self.packages.${system}.default}/bin/alterpcb-tlinesim";
        meta = self.packages.${system}.default.meta;
      };
    });
  };
}
