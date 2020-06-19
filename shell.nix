let
  pkgs = import ./nixpkgs.nix;

  ml_libs =  pkgs.python37Packages ;
in
  pkgs.mkShell {
    name = "astro";
    buildInputs = with pkgs //  ml_libs; [
      python37
      numpy
      typeguard
      numba
      matplotlib
      pandas
      scipy
      scikitlearn
      h5py
      seaborn
      astropy
      astroquery
      beautifulsoup4
   ];
   shellHook = ''
   '';
}
