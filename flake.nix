{
  description = "ddf-pipeline flake using uv2nix";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

    pyproject-nix = {
      url = "github:pyproject-nix/pyproject.nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    uv2nix = {
      url = "github:pyproject-nix/uv2nix";
      inputs.pyproject-nix.follows = "pyproject-nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    pyproject-build-systems = {
      url = "github:pyproject-nix/build-system-pkgs";
      inputs.pyproject-nix.follows = "pyproject-nix";
      inputs.uv2nix.follows = "uv2nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs =
    {
      self,
      nixpkgs,
      uv2nix,
      pyproject-nix,
      pyproject-build-systems,
      ...
    }:
    let
      inherit (nixpkgs) lib;

      # Load a uv workspace from a workspace root.
      # Uv2nix treats all uv projects as workspace projects.
      workspace = uv2nix.lib.workspace.loadWorkspace { workspaceRoot = ./.; };

      # Create package overlay from workspace.
      overlay = workspace.mkPyprojectOverlay {
        # Prefer prebuilt binary wheels as a package source.
        # Sdists are less likely to "just work" because of the metadata missing from uv.lock.
        # Binary wheels are more likely to, but may still require overrides for library dependencies.
        sourcePreference = "wheel"; # or sourcePreference = "sdist";
        # Optionally customise PEP 508 environment
        # environ = {
        #   platform_release = "5.10.65";
        # };
      };

      # Extend generated overlay with build fixups
      #
      # Uv2nix can only work with what it has, and uv.lock is missing essential metadata to perform some builds.
      # This is an additional overlay implementing build fixups.
      # See:
      # - https://pyproject-nix.github.io/uv2nix/FAQ.html
      pyprojectOverrides = _final: _prev: let
        withSetupTools = pkg: { "${pkg}" = _prev.${pkg}.overrideAttrs (old: {
          buildInputs = (old.buildInputs or [ ]) ++ _final.resolveBuildSystem { setuptools = [ ]; };
        });
                              };

      in {
        # Implement build fixups here.
        # Note that uv2nix is _not_ using Nixpkgs buildPythonPackage.
        # It's using https://pyproject-nix.github.io/pyproject.nix/build.html
        astlib = _prev.astlib.overrideAttrs(old: {
          NIX_CFLAGS_COMPILE = "-Wno-error=implicit-function-declaration";
          buildInputs = (old.buildInputs or [ ]) ++ _final.resolveBuildSystem { setuptools = [ ]; };
        });
        mpi4py = _prev.mpi4py .overrideAttrs (old: {
          buildInputs = (old.buildInputs or [ ]) ++ _final.resolveBuildSystem { setuptools = [ ]; }
                        ++ [ _prev.pkgs.openmpi ];
        });
        numba= _prev.numba.overrideAttrs (old: {
          buildInputs = (old.buildInputs or [ ]) ++ [ _prev.pkgs.onetbb ];
        });
        astro-pyxis = _prev.astro-pyxis.overrideAttrs (old: {
          buildInputs = (old.buildInputs or [ ]) ++ _final.resolveBuildSystem {
            setuptools = [ ];
            six = [ ];
          };
        });
        ddfacet = _prev.ddfacet.overrideAttrs (old: {
          buildInputs = (old.buildInputs or [ ]) ++ _final.resolveBuildSystem {
            # add all build-system.requires manually here (https://github.com/astral-sh/uv/issues/5190)
            scikit-build-core = [ ];
            cmake = [ ];
            numpy = [];
            pybind11 = [ ];
          };
        });
        #montblanc = _prev.montblanc.overrideAttrs (old: {
        #  buildInputs = (old.buildInputs or [ ]) ++ _final.resolveBuildSystem {
        #    setuptools = [ ];
        #    tensorflow = [ ];
        #  };
        #});
        #tensorflow-io-gcs-filesystem = _prev.tensorflow-io-gcs-filesystem.overrideAttrs (old: {
        #  buildInputs = (old.buildInputs or [ ]) ++ _final.resolveBuildSystem {
        #    setuptools = [ ];
        #  } ++ [
        #    _prev.pkgs.libtensorflow
        #  ];
        #});
        sharedarray = _prev.sharedarray.overrideAttrs (old: {
          buildInputs = (old.buildInputs or [ ]) ++ _final.resolveBuildSystem {
            setuptools = [ ];
            numpy = [ ];
          };
        });
        pyfftw = _prev.pyfftw.overrideAttrs (old: {
          NIX_LDFLAGS="-L${_prev.pkgs.fftw.out}/lib -L${_prev.pkgs.fftwFloat.out}/lib -L${_prev.pkgs.fftwLongDouble.out}/lib";
          NIX_CFLAGS_COMPILE="-I${_prev.pkgs.fftw.dev}/include -I${_prev.pkgs.fftwFloat.dev}/include -I${_prev.pkgs.fftwLongDouble.dev}/include";
          buildInputs = (old.buildInputs or [ ]) ++ _final.resolveBuildSystem {
            setuptools = [ ];
            numpy = [ ];
            scipy = [ ];
          }
          ++ [
            _prev.pkgs.fftw
            _prev.pkgs.fftwFloat
            _prev.pkgs.fftwLongDouble
            ];
        });
      }
        // (withSetupTools "astro-kittens")
        // (withSetupTools "kittens")
        // (withSetupTools "meqtrees-cattery")
        // (withSetupTools "owlcat")
        // (withSetupTools "donfig")
        // (withSetupTools "hypercube")
        // (withSetupTools "polygon3")
        // (withSetupTools "purr")
        // (withSetupTools "pymoresane")
        ;

      # This example is only using x86_64-linux
      pkgs = nixpkgs.legacyPackages.x86_64-linux;

      # Use Python 3.12 from nixpkgs
      python = pkgs.python310;

      # Construct package set
      pythonSet =
        # Use base package set from pyproject.nix builders
        (pkgs.callPackage pyproject-nix.build.packages {
          inherit python;
        }).overrideScope
          (
            lib.composeManyExtensions [
              pyproject-build-systems.overlays.default
              overlay
              pyprojectOverrides
            ]
          );

    in
    {
      # Package a virtual environment as our main application.
      #
      # Enable no optional dependencies for production build.
      packages.x86_64-linux.default = pythonSet.mkVirtualEnv "ddf-pipeline" workspace.deps.default;

      # Make hello runnable with `nix run`
      apps.x86_64-linux = {
        default = {
          type = "app";
          program = "${self.packages.x86_64-linux.default}/bin/hello";
        };
      };

      # This example provides two different modes of development:
      # - Impurely using uv to manage virtual environments
      # - Pure development using uv2nix to manage virtual environments
      devShells.x86_64-linux = {
        # It is of course perfectly OK to keep using an impure virtualenv workflow and only use uv2nix to build packages.
        # This devShell simply adds Python and undoes the dependency leakage done by Nixpkgs Python infrastructure.
        impure = pkgs.mkShell {
          packages = [
            python
            pkgs.uv
          ];
          env =
            {
              # Prevent uv from managing Python downloads
              UV_PYTHON_DOWNLOADS = "never";
              # Force uv to use nixpkgs Python interpreter
              UV_PYTHON = python.interpreter;
              ENVRC = "impure";
            }
            // lib.optionalAttrs pkgs.stdenv.isLinux {
              # Python libraries often load native shared objects using dlopen(3).
              # Setting LD_LIBRARY_PATH makes the dynamic library loader aware of libraries without using RPATH for lookup.
              LD_LIBRARY_PATH = lib.makeLibraryPath pkgs.pythonManylinuxPackages.manylinux1;
            };
          shellHook = ''
            unset PYTHONPATH
          '';
        };

        # This devShell uses uv2nix to construct a virtual environment purely from Nix, using the same dependency specification as the application.
        # The notable difference is that we also apply another overlay here enabling editable mode ( https://setuptools.pypa.io/en/latest/userguide/development_mode.html ).
        #
        # This means that any changes done to your local files do not require a rebuild.
        #
        # Note: Editable package support is still unstable and subject to change.
        uv2nix =
          let
            # Create an overlay enabling editable mode for all local dependencies.
            editableOverlay = workspace.mkEditablePyprojectOverlay {
              # Use environment variable
              root = "$REPO_ROOT";
              # Optional: Only enable editable for these packages
              members = [ "ddf-pipeline" ];
            };

            # Override previous set with our overrideable overlay.
            editablePythonSet = pythonSet.overrideScope (
              lib.composeManyExtensions [
                editableOverlay

                # Apply fixups for building an editable package of your workspace packages
                (final: prev: {
                  ddf-pipeline = prev.ddf-pipeline.overrideAttrs (old: {
                    # It's a good idea to filter the sources going into an editable build
                    # so the editable package doesn't have to be rebuilt on every change.
                    src = lib.fileset.toSource {
                      root = old.src;
                      fileset = lib.fileset.unions [
                        (old.src + "/pyproject.toml")
                        (old.src + "/README.md")
                        (old.src + "/LICENSE.md")
                        (old.src + "/examples")
                        (old.src + "/utils")
                        (old.src + "/scripts")
                        #(old.src + "/src/hello_world/__init__.py")
                      ];
                    };
                  });
                })
              ]
            );

            # Build virtual environment, with local packages being editable.
            #
            # Enable all optional dependencies for development.
            virtualenv = editablePythonSet.mkVirtualEnv "ddf-pipeline-dev" workspace.deps.default;

          in
          pkgs.mkShell {
            packages = [
              virtualenv
              pkgs.uv
            ];

            env = {
              # Don't create venv using uv
              UV_NO_SYNC = "1";

              # Force uv to use Python interpreter from venv
              UV_PYTHON = "${virtualenv}/bin/python";

              # Prevent uv from downloading managed Python's
              UV_PYTHON_DOWNLOADS = "never";

              ENVRC = "uv2nix";
              DDF_LOCAL_DEV = "1";
            };

            shellHook = ''
              # Undo dependency propagation by nixpkgs.
              unset PYTHONPATH

              # Get repository root using git. This is expanded at runtime by the editable `.pth` machinery.
              export REPO_ROOT=$(git rev-parse --show-toplevel)
            '';
          };
      };
    };
}
