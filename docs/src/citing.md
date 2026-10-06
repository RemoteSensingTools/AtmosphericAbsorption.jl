# Citing

AtmosphericAbsorption.jl makes HITRAN data easy to use, so it is easy to forget that the
spectroscopic parameters behind every line-by-line cross-section come from the
[HITRAN database](https://hitran.org). If you publish results computed with HITRAN data
through this package, **please cite the HITRAN edition you used** (and HAPI where noted below).

## Which edition did I use?

- **Direct downloads** (`load_hitran`, `fetch_hitran`, `HitranPort(; edition)`,
  `load_hitran_nonvoigt`): hitran.org serves only the *current* edition, which is
  **HITRAN2024** at the time of writing. The `edition` keyword is a cache/provenance label;
  it does not select an older edition from the server.
- **Local `.par` files** (`HitranPort("file.par"; edition)`): cite the edition the file
  came from. Set the `edition` keyword to match so the provenance record is accurate.
- Each cached download has a `.meta.toml` sidecar recording the edition label, spectral
  window, source URL, SHA-256 checksum, and download date (see
  [Provenance and reproducibility](data_sources.md#provenance-and-reproducibility)). Use it to
  confirm which data a result was computed from.

## HITRAN references

**HITRAN2024**

Gordon IE, Rothman LS, Hargreaves RJ, Gomez FM, Bertin T, Hill C, Kochanov RV, Tan Y, Wcisło P, Makhnev VYu, Bernath PF, Birk M, Boudon V, Campargue A, Coustenis A, Drouin BJ, Gamache RR, Hodges JT, Jacquemart D, Mlawer EJ, Nikitin AV, Perevalov VI, Rotger M, Robert S, Tennyson J, Toon GC, Tran H, Tyuterev VG, Adkins EM, Barbe A, Bailey DM, Bielska K, Bizzocchi L, Blake TA, Bowesman CA, Cacciani P, Čermák P, Császár AG, Denis L, Egbert SC, Egorov O, Ermilov AYu, Fleisher AJ, Fleurbaey H, Foltynowicz A, Furtenbacher T, Germann M, Guest ER, Harrison JJ, Hartmann J-M, Hjältén A, Hu S-M, Huang X, Johnson TJ, Jóźwiak H, Kassi S, Khan MV, Kwabia-Tchana F, Lee TJ, Lisak D, Liu A-W, Lyulin OM, Malarich NA, Manceron L, Marinina AA, Massie ST, Mascio J, Medvedev ES, Meshkov VV, Mellau GCh, Melosso M, Mikhailenko SN, Mondelain D, Müller HSP, O'Donnell M, Owens A, Perrin A, Polyansky OL, Raston PL, Reed ZD, Rey M, Richard C, Rieker GB, Röske C, Sharpe SW, Starikova E, Stolarczyk N, Stolyarov AV, Sung K, Tamassia F, Terragni J, Ushakov VG, Vasilchenko S, Vispoel B, Vodopyanov KL, Wagner G, Wójtewicz S, Yurchenko SN, Zobov NF. The HITRAN2024 molecular spectroscopic database. *J Quant Spectrosc Radiat Transf* 2026;353:109807. [doi:10.1016/j.jqsrt.2026.109807](https://doi.org/10.1016/j.jqsrt.2026.109807)

**HITRAN2020**

Gordon IE, Rothman LS, Hargreaves RJ, Hashemi R, Karlovets EV, Skinner FM, Conway EK, Hill C, Kochanov RV, Tan Y, Wcisło P, Finenko AA, Nelson K, Bernath PF, Birk M, Boudon V, Campargue A, Chance KV, Coustenis A, Drouin BJ, Flaud J-M, Gamache RR, Hodges JT, Jacquemart D, Mlawer EJ, Nikitin AV, Perevalov VI, Rotger M, Tennyson J, Toon GC, Tran H, Tyuterev VG, Adkins EM, Baker A, Barbe A, Canè E, Császár AG, Dudaryonok A, Egorov O, Fleisher AJ, Fleurbaey H, Foltynowicz A, Furtenbacher T, Harrison JJ, Hartmann J-M, Horneman V-M, Huang X, Karman T, Karns J, Kassi S, Kleiner I, Kofman V, Kwabia-Tchana F, Lavrentieva NN, Lee TJ, Long DA, Lukashevskaya AA, Lyulin OM, Makhnev VYu, Matt W, Massie ST, Melosso M, Mikhailenko SN, Mondelain D, Müller HSP, Naumenko OV, Perrin A, Polyansky OL, Raddaoui E, Raston PL, Reed ZD, Rey M, Richard C, Tóbiás R, Sadiek I, Schwenke DW, Starikova E, Sung K, Tamassia F, Tashkun SA, Vander Auwera J, Vasilenko IA, Vigasin AA, Villanueva GL, Vispoel B, Wagner G, Yachmenev A, Yurchenko SN. The HITRAN2020 molecular spectroscopic database. *J Quant Spectrosc Radiat Transf* 2022;277:107949. [doi:10.1016/j.jqsrt.2021.107949](https://doi.org/10.1016/j.jqsrt.2021.107949)

**HITRAN2016**

Gordon IE, Rothman LS, Hill C, Kochanov RV, Tan Y, Bernath PF, Birk M, Boudon V, Campargue A, Chance KV, Drouin BJ, Flaud J-M, Gamache RR, Hodges JT, Jacquemart D, Perevalov VI, Perrin A, Shine KP, Smith M-AH, Tennyson J, Toon GC, Tran H, Tyuterev VG, Barbe A, Császár AG, Devi VM, Furtenbacher T, Harrison JJ, Hartmann J-M, Jolly A, Johnson TJ, Karman T, Kleiner I, Kyuberis AA, Loos J, Lyulin OM, Massie ST, Mikhailenko SN, Moazzen-Ahmadi N, Müller HSP, Naumenko OV, Nikitin AV, Polyansky OL, Rey M, Rotger M, Sharpe SW, Sung K, Starikova E, Tashkun SA, Vander Auwera J, Wagner G, Wilzewski J, Wcisło P, Yu S, Zak EJ. The HITRAN2016 molecular spectroscopic database. *J Quant Spectrosc Radiat Transf* 2017;203:3–69. [doi:10.1016/j.jqsrt.2017.06.038](https://doi.org/10.1016/j.jqsrt.2017.06.038)

## HAPI

The HITRAN Application Programming Interface (HAPI) is the reference implementation this
package is validated against. Its `pcqsdhc` routine is the basis of our speed-dependent and
Hartmann–Tran profiles, and non-Voigt parameters are fetched from the HITRANonline API. Please also cite HAPI
when you use those profiles or parameters:

Kochanov RV, Gordon IE, Rothman LS, Wcisło P, Hill C, Wilzewski JS. HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data. *J Quant Spectrosc Radiat Transf* 2016;177:15–30. [doi:10.1016/j.jqsrt.2016.03.005](https://doi.org/10.1016/j.jqsrt.2016.03.005)
