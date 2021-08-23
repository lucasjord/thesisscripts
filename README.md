Scripts used for data reduction in PhD Thesis. 
Please contact lucas.hyland@utas.edu.au or mreid@cfa.harvard.edu for help with relevent scripts or programmes.

=========================================
Chapter 3:

fit_geoblocks_tropos_v7a.f
	Geoblock fitting programme. Written by Prof Mark Reid. Ref Reid et al. 2009.

fit_parallax_multi_4d.f
	Programme for least squares parallax and proper motion fitting for multiple quasars. Credit to Prof Mark Reid.

=========================================
Chapter 4:

maser_amplitude_calibration.py:
	Script used to perform secondary amplitude calibration using the autocorrelated spectrum of the masers.

find_peaks_uv.py:
    Script to get peaks in scalar and baseline averaged xcorrelated spectrum of masers, then print the visibilities (real,imag) out.

lsq_fit_maser_uv.py:
	Script to lsq fit core--halo model to the ourput of 'find_peaks_uv.py' and print results to file.


=========================================
Chapter 6:

control_mv028.py
	Example mv* control script. Modified from BeSSeL version.

definitions_mv.py
	mv* slave script. Modified from original definitions_bessel_current.py written by Dr Andreas Brunthalar.

fit_geoblocks_w_baselines_v2b.f
	Modified 'fit_geoblocks_tropos_v7a.f' with additional parameters for telescope offsets in X,Y,Z.

quasar_amplitude_autocorrect.py:
	Script for quasar amplitude calibration over baselines and frequency.

fit_phase_plane_v6.py:
	Script for 2D phase plane fitting to mv* data, with approximate goodness of fit and basic phase wrap handling.

splatnimagr.py:
	Multiview master script. Collects and organised visibility data in AIPS, fring's data etc, then runs slave script 'fit_phase_plane_v6.py', formats output, reuploads it, images and jmfits images.












