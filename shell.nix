let
  pkgs = import ./nixpkgs.nix;

  ml_libs =  pkgs.python37Packages ;
in
  pkgs.mkShell {
    name = "astro";
    buildInputs = with pkgs //  ml_libs; [
      joblib
      cython
      plotly
      jupyterlab
      python37
      virtualenv
      tqdm
      numpy
      typeguard
      numba
      matplotlib
      pandas
      scipy
      h5py
      seaborn
      pytest
      astropy
      beautifulsoup4
      R
      rPackages.rmarkdown
      rPackages.reticulate
   ];
   shellHook = ''
      export SOURCE_DATE_EPOCH=$(date +%s)   '';
}
