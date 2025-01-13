DIALS 3.23.1 (2025-01-13)
=========================

No significant changes.


DIALS 3.23.0 (2025-01-08)
=========================

Features
--------

- ``dials.estimate_resolution``: Add resolution estimate based on limit of cc1/2 significance. (`#2580 <https://github.com/dials/dials/issues/2580>`_)
- Added initial support for a H5-based reflection table. This is turned off by default, but can be activated by setting the ``DIALS_USE_H5`` environment variable. (`#2677 <https://github.com/dials/dials/issues/2677>`_)
- ``dials.refine_error_model``: Allow grouped and individual error model refinements. (`#2685 <https://github.com/dials/dials/issues/2685>`_)
- ``dials.symmetry``: Enable oversampling in ``systematic_absences.method=fourier``, for screw-axis detection. (`#2701 <https://github.com/dials/dials/issues/2701>`_)
- ``dials.index`` and ``dials.ssx_index``: Add the CUDA-accelerated fast-feedback-indexer to DIALS as a lattice search algorithm. See https://github.com/paulscherrerinstitute/fast-feedback-indexer for more details. (`#2717 <https://github.com/dials/dials/issues/2717>`_)
- More DIALS tests are changed to use publicly-available data via ``dials-data`` (`#2758 <https://github.com/dials/dials/issues/2758>`_)
- ``dials.correlation_matrix``: Include reachability plot for visualisation of coordinate-clustering. (`#2778 <https://github.com/dials/dials/issues/2778>`_)
- Python 3.10 is now the minimum python supported. (`#2792 <https://github.com/dials/dials/issues/2792>`_)
- ``dials.image_viewer``: Add a line tool that gives start, end, and midpoints, plus a line profile graph. (`#2798 <https://github.com/dials/dials/issues/2798>`_)
- ``dials.search_beam_position``: Include three new projection methods. (`#2809 <https://github.com/dials/dials/issues/2809>`_)
- ``dials.import``: Add syntax for splitting a long scan into chunks e.g. for sequential data collections on small crystals. New syntax is either ``split=frames_per_block`` or ``split=start,end,frames_per_block``, if you wish to specify the range, where start and end are 1-indexed. (`#2816 <https://github.com/dials/dials/issues/2816>`_)


Bugfixes
--------

- Ensure logger info is always printed correctly in spot-finding. (`#2678 <https://github.com/dials/dials/issues/2678>`_)
- Fix ``beam.get_s0()`` errors when trying to use image viewer with polychromatic data. (`#2765 <https://github.com/dials/dials/issues/2765>`_)
- Fix symmetry not being applied consistently post indexing, when multiple crystals and multiple orientations are present. (`#2786 <https://github.com/dials/dials/issues/2786>`_)
- ``dials.generate_distortion_maps``: Correct definitions for creating elliptical distortion correction maps. (`#2787 <https://github.com/dials/dials/issues/2787>`_)
- ``dials.predict``: Now adds experiment identifiers to output reflections. (`#2788 <https://github.com/dials/dials/issues/2788>`_)
- ``dials.image_viewer``: Recover missing image scrollbar by increasing the panel height. (`#2795 <https://github.com/dials/dials/issues/2795>`_)
- ``dials.index``: fix inconsistent experiment IDs after indexing with ``index_assignment.method=local`` (`#2799 <https://github.com/dials/dials/issues/2799>`_)
- ``dials.index``: Fix minor bug in global index assignment. (`#2800 <https://github.com/dials/dials/issues/2800>`_)
- ``dials.predict``: Do not try to filter shadowed reflections, if no masker is provided. (`#2805 <https://github.com/dials/dials/issues/2805>`_)
- ``dials.cosym``: Use unbiased formula for standard error. (`#2810 <https://github.com/dials/dials/issues/2810>`_)
- ``dials.cosym``: Fix weighting of cosym objective function, use count weights by default. (`#2813 <https://github.com/dials/dials/issues/2813>`_)
- ``dials.image_viewer``: Do not allow integer controls to take values outside of their bounds. (`#2822 <https://github.com/dials/dials/issues/2822>`_)
- ``dials.refine_bravais_settings``: Select only the reflections used in refinement in the indexing step for the calculation, to reduce memory footprint with very high resolution data. (`#2824 <https://github.com/dials/dials/issues/2824>`_)
- ``dials.image_viewer``: Fix display bug with the trusted range mask, when viewing stacked images. (`#2830 <https://github.com/dials/dials/issues/2830>`_)
- ``dials.export``: Fixes to work with gemmi v0.7. (`#2834 <https://github.com/dials/dials/issues/2834>`_)
- ``dials.image_viewer``: Avoid error messages when quitting while one of the tools (ring, unit cell, etc.) is open. (`#2836 <https://github.com/dials/dials/issues/2836>`_)
- ``dials.integrate``: Allow command line option ``d_min=`` without requiring full qualification as ``prediction.d_min=``. (`#2837 <https://github.com/dials/dials/issues/2837>`_)
- More fixes for gemmi 0.7. (`#2839 <https://github.com/dials/dials/issues/2839>`_)


Misc
----

- `#2528 <https://github.com/dials/dials/issues/2528>`_, `#2579 <https://github.com/dials/dials/issues/2579>`_, `#2769 <https://github.com/dials/dials/issues/2769>`_, `#2776 <https://github.com/dials/dials/issues/2776>`_, `#2779 <https://github.com/dials/dials/issues/2779>`_, `#2780 <https://github.com/dials/dials/issues/2780>`_, `#2791 <https://github.com/dials/dials/issues/2791>`_, `#2794 <https://github.com/dials/dials/issues/2794>`_, `#2796 <https://github.com/dials/dials/issues/2796>`_, `#2801 <https://github.com/dials/dials/issues/2801>`_, `#2802 <https://github.com/dials/dials/issues/2802>`_, `#2808 <https://github.com/dials/dials/issues/2808>`_, `#2818 <https://github.com/dials/dials/issues/2818>`_, `#2821 <https://github.com/dials/dials/issues/2821>`_, `#2832 <https://github.com/dials/dials/issues/2832>`_, `#2833 <https://github.com/dials/dials/issues/2833>`_, `#2838 <https://github.com/dials/dials/issues/2838>`_


DIALS 3.22.0 (2024-10-15)
=========================

Features
--------

- Added scaling corrections for time-of-flight data. (`#2704 <https://github.com/dials/dials/issues/2704>`_)
- ``dials.generate_mask``: Allow combining masks. (`#2711 <https://github.com/dials/dials/issues/2711>`_)
- ``dials.image_viewer``: Added multi-panel polygon masking. (`#2735 <https://github.com/dials/dials/issues/2735>`_)
- ``dials.correlation_matrix``: Add dimension optimisation for intensity-based dataset clustering. (`#2743 <https://github.com/dials/dials/issues/2743>`_)
- ``dials.generate_distortion_maps``: Much faster map generation for ``mode=ellipse``. (`#2745 <https://github.com/dials/dials/issues/2745>`_)
- Switch ``bootstrap.py`` to use CMake by default. You can pass ``--libtbx`` to get the previous behaviour, if you need to develop cctbx. (`#2755 <https://github.com/dials/dials/issues/2755>`_)
- ``dials.correlation_matrix``: Add auto-selection and output of significant clusters and extended graphs in html output (`#2770 <https://github.com/dials/dials/issues/2770>`_)


Bugfixes
--------

- ``dials.index``: Ensure reflections associated with removed experiments are reset to unindexed. (`#2653 <https://github.com/dials/dials/issues/2653>`_)
- ``dials.find_spots``: Enable ``threshold.algorithm=radial_profile`` for polychromatic data. (`#2729 <https://github.com/dials/dials/issues/2729>`_)
- ``dials.image_viewer``: Fixed stacking of masks, and stacking of multiple experiments. (`#2730 <https://github.com/dials/dials/issues/2730>`_)
- ``dials.image_viewer``: Fix UI glitches in masking tools. (`#2731 <https://github.com/dials/dials/issues/2731>`_)
- ``dials.resolution_analysis``: Suppress output of potential scipy OptimizeWarning. (`#2737 <https://github.com/dials/dials/issues/2737>`_)
- ``dials.generate_distortion_maps``: Fix bug with ``mode=ellipse`` for detectors with oblong panels. (`#2740 <https://github.com/dials/dials/issues/2740>`_)
- ``dials.cosym``: Handle the case of very few reflections, leading to zero weights in weighted cosym analysis. (`#2741 <https://github.com/dials/dials/issues/2741>`_)
- ``dials.find_spots_server``: fix multiprocessing to use fork on macOS and on all platforms (and needed for Python 3.14) (`#2747 <https://github.com/dials/dials/issues/2747>`_)
- DIALS viewers: Tweak size of UI elements for GTK3. (`#2753 <https://github.com/dials/dials/issues/2753>`_)
- ``dials.integrate``: Ensure the polarization correction is only performed for X-ray data (not electron nor neutron). (`#2754 <https://github.com/dials/dials/issues/2754>`_)
- ``dials.image_viewer``: Fix broken "Save As" PNG functionality. (`#2759 <https://github.com/dials/dials/issues/2759>`_)
- API: Ensure correct results from the ``determine_reindex_operator_against_reference`` function when the test and reference miller sets are the same object. (`#2761 <https://github.com/dials/dials/issues/2761>`_)
- ``dials.index``: When indexing multiple imagesets, ensure a failure to index one does not result in failure for all. (`#2764 <https://github.com/dials/dials/issues/2764>`_)
- ``dials.image_viewer``: Fix errors when trying to view polychromatic data. (`#2765 <https://github.com/dials/dials/issues/2765>`_)
- ``dials.reciprocal_lattice_viewer``: Fix errors when trying to view polychromatic data. (`#2766 <https://github.com/dials/dials/issues/2766>`_)
- ``dials.geometry_viewer``: Fix errors when trying to view polychromatic experiments. (`#2773 <https://github.com/dials/dials/issues/2773>`_)


Improved Documentation
----------------------

- ``dials.find_spots``: Added second implementation of signal mask erosion, for documentation. (`#2665 <https://github.com/dials/dials/issues/2665>`_)


Deprecations and Removals
-------------------------

- ``bootstrap.py`` no longer accepts ``--prebuilt-cctbx``. This was only ever supported in CI. If you want to build a distribution with prebuilt cctbx, please use the ``--cmake`` flag and build instead. (`#2727 <https://github.com/dials/dials/issues/2727>`_)
- ``dials.find_spots_server``: No longer supported in Windows. (`#2747 <https://github.com/dials/dials/issues/2747>`_)


Misc
----

- `#2727 <https://github.com/dials/dials/issues/2727>`_, `#2736 <https://github.com/dials/dials/issues/2736>`_, `#2738 <https://github.com/dials/dials/issues/2738>`_, `#2742 <https://github.com/dials/dials/issues/2742>`_, `#2756 <https://github.com/dials/dials/issues/2756>`_, `#2767 <https://github.com/dials/dials/issues/2767>`_, `#2768 <https://github.com/dials/dials/issues/2768>`_


Dials 3.22.0 (2024-10-15)
=========================

Features
--------

- ``dials.find_spots``: add second implementation of signal mask erosion for documentation and clarify comments (`#2665 <https://github.com/dials/dials/issues/2665>`_)
- Added time-of-flight Lorentz, spherical absorption, and normalisation w.r.t empty and incident run scaling corrections. (`#2704 <https://github.com/dials/dials/issues/2704>`_)
- Allow combining masks using dials.generate_mask (`#2711 <https://github.com/dials/dials/issues/2711>`_)
- Implemented multi-panel polygon masking in dials.image_viewer (`#2735 <https://github.com/dials/dials/issues/2735>`_)
- ``dials.correlation_matrix``: Add dimension optimisation for intensity-based dataset clustering (`#2743 <https://github.com/dials/dials/issues/2743>`_)
- ``dials.generate_distortion_maps``: Much faster map generation for ``mode=ellipse``. (`#2745 <https://github.com/dials/dials/issues/2745>`_)
- Make CMake the default build mode of bootstrap. (`#2755 <https://github.com/dials/dials/issues/2755>`_)
- ``dials.correlation_matrix``: Add auto-selection and output of significant clusters and extended graphs in html output (`#2770 <https://github.com/dials/dials/issues/2770>`_)


Bugfixes
--------

- ``dials.index``: ensure reflections associated with removed experiments are reset to unindexed. (`#2653 <https://github.com/dials/dials/issues/2653>`_)
- Enable radial_profile spot finding algorithm for polychromatic data by removing call to beam s0. (`#2729 <https://github.com/dials/dials/issues/2729>`_)
- Fixed stacking of masks and stacking of multiple experiments (e.g. stills) in dials.image_viewer (#1512, #2724) (`#2730 <https://github.com/dials/dials/issues/2730>`_)
- ``dials.image_viewer``: Fix UI glitches in masking tools. (`#2731 <https://github.com/dials/dials/issues/2731>`_)
- ``dials.resolution_analysis``: Suppress output of potential scipy OptimizeWarning. (`#2737 <https://github.com/dials/dials/issues/2737>`_)
- ``dials.generate_distortion_maps``: fix bug with ``mode=ellipse`` for detectors with oblong panels. (`#2740 <https://github.com/dials/dials/issues/2740>`_)
- ``dials.cosym``: Handle the case of very few reflections leading to zero weights in weighted cosym analysis. (`#2741 <https://github.com/dials/dials/issues/2741>`_)
- ``dials.find_spots_server``: fix multiprocessing to use fork on macOS and on all platforms (and needed for Python 3.14) (`#2747 <https://github.com/dials/dials/issues/2747>`_)
- DIALS viewers: tweak size of UI elements. (`#2753 <https://github.com/dials/dials/issues/2753>`_)
- ``dials.integrate``: Ensure the polarization correction is only performed for X-ray data (not electron nor neutron). (`#2754 <https://github.com/dials/dials/issues/2754>`_)
- ``dials.image_viewer``: Fix broken "Save As" PNG functionality. (`#2759 <https://github.com/dials/dials/issues/2759>`_)
- Ensure correct results from the ``determine_reindex_operator_against_reference``
  function when the test and reference miller sets are the same object. (`#2761 <https://github.com/dials/dials/issues/2761>`_)
- ``dials.index``: when indexing multiple imagesets, ensure a failure to
  index one does not result in failure for all. (`#2764 <https://github.com/dials/dials/issues/2764>`_)
- Fix beam.get_s0() errors when trying to use image viewer with polychromatic data. (`#2765 <https://github.com/dials/dials/issues/2765>`_)
- Fix reciprocal lattice viewer not working with polychromatic data due to call to beam.get_s0(). (`#2766 <https://github.com/dials/dials/issues/2766>`_)
- Fix errors with loading polychromatic experiments into dials.geometry_viewer. (`#2773 <https://github.com/dials/dials/issues/2773>`_)


Deprecations and Removals
-------------------------

- ``bootstrap.py`` no longer accepts ``--prebuilt-cctbx``. This was only ever supported in CI. If you want to build a distribution with prebuilt cctbx, please use the ``--cmake`` flag and build instead. (`#2727 <https://github.com/dials/dials/issues/2727>`_)
- ``dials.find_spots_server``: No longer supported in Windows. (`#2747 <https://github.com/dials/dials/issues/2747>`_)


Misc
----

- `#2727 <https://github.com/dials/dials/issues/2727>`_, `#2736 <https://github.com/dials/dials/issues/2736>`_, `#2738 <https://github.com/dials/dials/issues/2738>`_, `#2742 <https://github.com/dials/dials/issues/2742>`_, `#2756 <https://github.com/dials/dials/issues/2756>`_, `#2767 <https://github.com/dials/dials/issues/2767>`_, `#2768 <https://github.com/dials/dials/issues/2768>`_


DIALS 3.21.1 (2024-08-23)
=========================

Bugfixes
--------

- CMake build can now be used for release branches. (`#2727 <https://github.com/dials/dials/issues/2727>`_)


DIALS 3.21.0 (2024-08-20)
=========================

Features
--------

- Add classes to support time-of-flight and Laue indexing and refinement. (`#2662 <https://github.com/dials/dials/issues/2662>`_)
- ``dials.symmetry``: Allow free selection of ``significance_level`` in the range [0,1]. (`#2696 <https://github.com/dials/dials/issues/2696>`_)
- ``dials.export``: Add support for exporting still data in mmcif format that GEMMI can read. (`#2709 <https://github.com/dials/dials/issues/2709>`_)


Bugfixes
--------

- ``dials.index``: Fix a ``pink_indexer`` error that caused failures for images with electron diffraction geometry. (`#2680 <https://github.com/dials/dials/issues/2680>`_)
- ``dials.correlation_matrix``: Correctly select datasets for output json after filtering, when used by multiplex. (`#2681 <https://github.com/dials/dials/issues/2681>`_)
- ``dials.index``: Avoid mm to px conversion when the ``refinement_protocol`` is set to do no refinement, as the required data are not available. (`#2687 <https://github.com/dials/dials/issues/2687>`_)
- Avoid deprecated ``matplotlib.cm.get_cmap`` calls (`#2688 <https://github.com/dials/dials/issues/2688>`_)
- Change Docker base image to rockylinux:8 as centos:7 is EOL (`#2689 <https://github.com/dials/dials/issues/2689>`_)
- ``dials.image_viewer``: Increase the maximum resolution for the ring tool radius. This was too small for new detectors. (`#2697 <https://github.com/dials/dials/issues/2697>`_)
- ``dials.index``: Fix potential crash in max_cell estimation when all spots are at the resolution of ice rings. (`#2699 <https://github.com/dials/dials/issues/2699>`_)
- ``dials.merge``: Use GEMMI to output merged MTZs, for consistency with ``dials.export`` (`#2700 <https://github.com/dials/dials/issues/2700>`_)
- ``dials.reciprocal_lattice_viewer``: Fix middle mouse drag to translate function. (`#2707 <https://github.com/dials/dials/issues/2707>`_)
- Performance improvement for selections from large reflection tables. For a table containing 165k experiment identifiers the speedup is 1000x (12 minutes per call). (`#2718 <https://github.com/dials/dials/issues/2718>`_)


Improved Documentation
----------------------

- Docstrings and type hints are added to the ``reindex_experiments`` and ``reindex_reflections`` functions to make it easier to use these outside the ``dials.reindex`` program. (`#2683 <https://github.com/dials/dials/issues/2683>`_)
- Improvements to the small molecule tutorial. (`#2692 <https://github.com/dials/dials/issues/2692>`_)


Deprecations and Removals
-------------------------

- `dials.import`: remove useless as_grid_scan option as misleading. To import as stills set scan.oscillation=0,0 (`#2454 <https://github.com/dials/dials/issues/2454>`_)


Misc
----

- `#2675 <https://github.com/dials/dials/issues/2675>`_, `#2702 <https://github.com/dials/dials/issues/2702>`_, `#2714 <https://github.com/dials/dials/issues/2714>`_, `#2715 <https://github.com/dials/dials/issues/2715>`_, `#2723 <https://github.com/dials/dials/issues/2723>`_


DIALS 3.20.0 (2024-06-19)
=========================

Features
--------

- New tool: ``dials.correlation_matrix``: A new command-line tool for correlation and cosine similarity clustering of multi-crystal datasets, independent of ``xia2.multiplex``. It provides HTML output, including clustering heatmaps, dendrograms and corresponding ``dials.cosym`` graphs. (`#2632 <https://github.com/dials/dials/issues/2632>`_)
- DIALS is now compatible with Python 3.12. (`#2651 <https://github.com/dials/dials/issues/2651>`_)
- ``dials.scale``: Add filtering options to default basic error model to allow error modelling of stills data. (`#2654 <https://github.com/dials/dials/issues/2654>`_)
- ``dials.cosym``: Add alternative weighting during cosym CC calculation with a new ``cc_weights=`` option. (`#2666 <https://github.com/dials/dials/issues/2666>`_)
- ``dials.correlation_matrix``: Add tables with cluster information to html output. (`#2671 <https://github.com/dials/dials/issues/2671>`_)
- New tool: ``dials.split_still_data`` for splitting dials-processed still data based on image number (e.g. dose series). (`#2672 <https://github.com/dials/dials/issues/2672>`_)


Bugfixes
--------

- ``dials.find_rotation_axis``: Correctly set the orientation of the rotation axis for a multi-axis goniometer. (`#2658 <https://github.com/dials/dials/issues/2658>`_)
- Fix ``dials.show`` beam checks for time of flight experiments. (`#2660 <https://github.com/dials/dials/issues/2660>`_)
- When masking, raise an error if ``d_min > d_max`` (where no spots would be found). (`#2664 <https://github.com/dials/dials/issues/2664>`_)
- ``dials.cosym``: Make function return structure correctly, in the recently added ``cc_weights=`` option. (`#2668 <https://github.com/dials/dials/issues/2668>`_)
- ``dials.find_rotation_axis``: Add reflection selection criteria, to avoid runs that use a very large amount of memory. (`#2670 <https://github.com/dials/dials/issues/2670>`_)
- ``dials.cosym``: For ``cc_weights=sigma``, ensure correct filtering based on min_pairs parameters. (`#2673 <https://github.com/dials/dials/issues/2673>`_)
- ``dials.cosym``: Fix to give more accurate cc calculation when running with a ``space_group=`` set. (`#2674 <https://github.com/dials/dials/issues/2674>`_)
- Fix ``flex.reflection_table.concat``, to not modify in-place. (`#2679 <https://github.com/dials/dials/issues/2679>`_)


Misc
----

- `#2613 <https://github.com/dials/dials/issues/2613>`_, `#2616 <https://github.com/dials/dials/issues/2616>`_, `#2641 <https://github.com/dials/dials/issues/2641>`_, `#2667 <https://github.com/dials/dials/issues/2667>`_


DIALS 3.19.1 (2024-05-23)
=========================

No significant changes.


DIALS 3.19.0 (2024-04-17)
=========================

Features
--------

- ``dials.predict``: Allow usage when image data are not available. (`#2553 <https://github.com/dials/dials/issues/2553>`_)
- Add ``TOFSpotFinder`` to tailor default params to time of flight experiments and add additional reflection table data. (`#2602 <https://github.com/dials/dials/issues/2602>`_)
- ``dials.ssx_index``: Allow use of sequences indexer, pink_indexer and low_res_spot_match indexing algorithms. (`#2614 <https://github.com/dials/dials/issues/2614>`_)
- ``dials.export``: Add option ``composition=`` to specify of asymmetric unit composition for SHELX ``.ins`` file output. (`#2623 <https://github.com/dials/dials/issues/2623>`_)


Bugfixes
--------

- Fix building on RHEL8 and other more recent distributions. (`#1465 <https://github.com/dials/dials/issues/1465>`_)
- ``dials.index``: Joint indexing is automatically set on for rotation data, off for still data. This can be overridden by explicit use of ``joint_indexing=``. (`#2605 <https://github.com/dials/dials/issues/2605>`_)
- ``dials.export`` ``format=shelx``: Increased precision of unit cell parameters and their estimated standard deviations written to ``.ins`` file. (`#2624 <https://github.com/dials/dials/issues/2624>`_)
- ``dials.ssx_index``: Don't combine detector models if individually refined. (`#2634 <https://github.com/dials/dials/issues/2634>`_)
- ``dials.scale``: Fix error in propagation of partiality variance, when scaling still data. (`#2642 <https://github.com/dials/dials/issues/2642>`_)
- ``dials.export``: Fix crash when exporting scaled still-shot data. (`#2646 <https://github.com/dials/dials/issues/2646>`_)


Deprecations and Removals
-------------------------

- API: ``array_family/flex_ext.py``: remove ``nthread``s parameter from ``extract_shoeboxes``, as it was never implemented. (`#2638 <https://github.com/dials/dials/issues/2638>`_)


Misc
----

- `#2617 <https://github.com/dials/dials/issues/2617>`_, `#2618 <https://github.com/dials/dials/issues/2618>`_, `#2619 <https://github.com/dials/dials/issues/2619>`_, `#2626 <https://github.com/dials/dials/issues/2626>`_, `#2633 <https://github.com/dials/dials/issues/2633>`_, `#2648 <https://github.com/dials/dials/issues/2648>`_


DIALS 3.18.1 (2024-03-26)
=========================

Bugfixes
--------

- wxPython was inadvertently reverted to 4.1 instead of 4.2 for python 3.9 and 3.10 releases. (`#2636 <https://github.com/dials/dials/issues/2636>`_)


DIALS 3.18.0 (2024-02-22)
=========================

Features
--------

- ``dials.slice_sequence``: Add ``exclude_images_multiple=`` option to split a scan at calibration images, as an alternative to excluding those images in ``dials.integrate``. (`#2565 <https://github.com/dials/dials/issues/2565>`_)
- ``dials.index``: Add ``indexing.method=pink_indexer`` for still images, based on `Gevorkov Y et al. pinkIndexer – a universal indexer for pink-beam X-ray and electron diffraction snapshots <https://doi.org/10.1107/S2053273319015559>`_. (`#2604 <https://github.com/dials/dials/issues/2604>`_)


Bugfixes
--------

- ``dials.index``: Correctly output imageset id when ``joint_indexing=False`` and ``max_lattices>1`` (`#2538 <https://github.com/dials/dials/issues/2538>`_)
- ``dials.index``: Ensure the unindexed_reflections table in the indexer classes is updated during refinement. (`#2545 <https://github.com/dials/dials/issues/2545>`_)
- When building a new developer installation without CMake, always defer to environmental flags when configuring compilers. This fixes building on MacOS Sonoma; for which the default libtbx configuration failed. (`#2546 <https://github.com/dials/dials/issues/2546>`_)
- ``dials.cosym``: Ensure correct reindexing of low symmetry data to the correct setting. (`#2548 <https://github.com/dials/dials/issues/2548>`_)
- ``dials.export``: Unmerged MTZ output now uses gemmi, fixing various longstanding issues. (`#2549 <https://github.com/dials/dials/issues/2549>`_)
- ``dials.symmetry/dials.cosym``: Avoid crashes when unable to find consistent symmetry during cell reduction routine. (`#2552 <https://github.com/dials/dials/issues/2552>`_)
- ``dials.sequence_to_stills``: Avoid creating an experiment for an image with no reflections on it. (`#2554 <https://github.com/dials/dials/issues/2554>`_)
- ``dials.reciprocal_lattice_viewer``: Fix toggling of ids on spotfinding output when viewing multiple files. (`#2566 <https://github.com/dials/dials/issues/2566>`_)
- ``dials.estimate_resolution``: Fix bug causing fitting failures for isigma, misigma curves. (`#2568 <https://github.com/dials/dials/issues/2568>`_)
- ``dials.reciprocal_lattice_viewer``: Ensure the correct panel for the beam centre is recorded when loading models. (`#2572 <https://github.com/dials/dials/issues/2572>`_)
- ``dials.reindex``: When reindexing against a reference, use the best cell when determining lattice symmetry. (`#2573 <https://github.com/dials/dials/issues/2573>`_)
- ``dials.cosym``: Filter out datasets that can't be mapped through a consistent best cell, rather than exiting. (`#2574 <https://github.com/dials/dials/issues/2574>`_)
- ``dials.export_bitmaps``: Fix for use of imageset_index parameter on image sets that don't start at image numbers of 1. (`#2591 <https://github.com/dials/dials/issues/2591>`_)
- Reorganize imports in ``export_mtz.py`` so that its gemmi-less objects can be used by other gemmi-less software. (`#2592 <https://github.com/dials/dials/issues/2592>`_)
- ``dials.image_viewer``: Fix for stills experiments viewing crash. (`#2594 <https://github.com/dials/dials/issues/2594>`_)
- ``dials.refine``: Clarify a misleading error message when there are no reflections to refine. (`#2607 <https://github.com/dials/dials/issues/2607>`_)
- Dependency update for https://github.com/cctbx/dxtbx/pull/697: adds ``natsort`` (`#2611 <https://github.com/dials/dials/issues/2611>`_)


Improved Documentation
----------------------

- ``dials.refine``: Update online documentation. (`#2551 <https://github.com/dials/dials/issues/2551>`_)


Misc
----

- `#2408 <https://github.com/dials/dials/issues/2408>`_, `#2547 <https://github.com/dials/dials/issues/2547>`_, `#2550 <https://github.com/dials/dials/issues/2550>`_, `#2556 <https://github.com/dials/dials/issues/2556>`_, `#2557 <https://github.com/dials/dials/issues/2557>`_, `#2563 <https://github.com/dials/dials/issues/2563>`_, `#2570 <https://github.com/dials/dials/issues/2570>`_, `#2576 <https://github.com/dials/dials/issues/2576>`_, `#2578 <https://github.com/dials/dials/issues/2578>`_, `#2586 <https://github.com/dials/dials/issues/2586>`_, `#2590 <https://github.com/dials/dials/issues/2590>`_, `#2599 <https://github.com/dials/dials/issues/2599>`_, `#2600 <https://github.com/dials/dials/issues/2600>`_, `#2603 <https://github.com/dials/dials/issues/2603>`_


Dials 3.18 (2024-02-22)
=======================

Features
--------

- ``dials.slice_sequence``: Add ``exclude_images_multiple`` option to split a scan at calibration images, as an alternative to excluding those images in ``dials.integrate``. See #1992 for details. (`#2565 <https://github.com/dials/dials/issues/2565>`_)
- ``dials.index``: Add ``indexing.method=pink_indexer`` for still images based on Gevorkov Y, Barty A, Brehm W, White TA, Tolstikova A, Wiedorn MO, et al. pinkIndexer – a universal indexer for pink-beam X-ray and electron diffraction snapshots. Acta Cryst A. 2020 Mar 1;76(2):121–31. (`#2604 <https://github.com/dials/dials/issues/2604>`_)


Bugfixes
--------

- ``dials.index``: Correctly output imageset id when joint_indexing=False and max_lattices>1 (`#2538 <https://github.com/dials/dials/issues/2538>`_)
- ``dials.index``: Ensure the unindexed_reflections table in the indexer classes is updated during refinement. (`#2545 <https://github.com/dials/dials/issues/2545>`_)
- Always use configured environment when bootstrapping CCTBX. This fixes building on MacOS Sonoma; for which the default libtbx configuration fails. (`#2546 <https://github.com/dials/dials/issues/2546>`_)
- ``dials.cosym``: Ensure correct reindexing of low symemtry data to the correct setting (`#2548 <https://github.com/dials/dials/issues/2548>`_)
- ``dials.export``: Unmerged MTZ output now uses gemmi, which fixes various issues (#1099, #1100, #2057, #2379, #2505) (`#2549 <https://github.com/dials/dials/issues/2549>`_)
- ``dials.refine``: Ensure online documentation is up to date. (`#2551 <https://github.com/dials/dials/issues/2551>`_)
- ``dials.symmetry/dials.cosym``: Avoid crashes when unable to find consistent symmetry during cell reduction routine. (`#2552 <https://github.com/dials/dials/issues/2552>`_)
- ``dials.sequence_to_stills``: Avoid creating an experiment for an image with no reflections on it. (`#2554 <https://github.com/dials/dials/issues/2554>`_)
- ``dials.reciprocal_lattice_viewer``: Fix toggling of ids on spotfinding output when viewing multiple files (`#2566 <https://github.com/dials/dials/issues/2566>`_)
- ``dials.estimate_resolution``: Fix bug causing fitting failures for isigma, misigma curves (`#2568 <https://github.com/dials/dials/issues/2568>`_)
- ``dials.reciprocal_lattice_viewer``: ensure the correct panel for the beam centre is recorded when loading models. (`#2572 <https://github.com/dials/dials/issues/2572>`_)
- ``dials.reindex``: When reindexing against a reference, use the best cell when determining lattice symmetry. (`#2573 <https://github.com/dials/dials/issues/2573>`_)
- ``dials.cosym``: Filter out datasets that can't be mapped through a consistent best cell, rather than exiting. (`#2574 <https://github.com/dials/dials/issues/2574>`_)
- ``dials.export_bitmaps``: Fix for use of imageset_index parameter on image sets that don't start at image numbers of 1. (`#2591 <https://github.com/dials/dials/issues/2591>`_)
- Reorganize imports in ``export_mtz.py`` so that its gemmi-less objects can be used by other gemmi-less software (`#2592 <https://github.com/dials/dials/issues/2592>`_)
- ``dials.image_viewer``: Fix for stills experiments viewing crash introduced in #2556 (`#2594 <https://github.com/dials/dials/issues/2594>`_)
- ``dials.refine``: clarify a misleading error message when there are no reflections to refine. (`#2607 <https://github.com/dials/dials/issues/2607>`_)


Misc
----

- `#2408 <https://github.com/dials/dials/issues/2408>`_, `#2547 <https://github.com/dials/dials/issues/2547>`_, `#2550 <https://github.com/dials/dials/issues/2550>`_, `#2556 <https://github.com/dials/dials/issues/2556>`_, `#2557 <https://github.com/dials/dials/issues/2557>`_, `#2563 <https://github.com/dials/dials/issues/2563>`_, `#2570 <https://github.com/dials/dials/issues/2570>`_, `#2576 <https://github.com/dials/dials/issues/2576>`_, `#2578 <https://github.com/dials/dials/issues/2578>`_, `#2586 <https://github.com/dials/dials/issues/2586>`_, `#2590 <https://github.com/dials/dials/issues/2590>`_, `#2599 <https://github.com/dials/dials/issues/2599>`_, `#2600 <https://github.com/dials/dials/issues/2600>`_


DIALS 3.17.0 (2023-11-03)
=========================

Features
--------

- Circular dependencies have been removed between DIALS and ``cctbx.xfel`` by using the new ``serialtbx``. (`#2404 <https://github.com/dials/dials/issues/2404>`_)
- ``dials.import``: Enable use of an image range selection when importing a still sequence. (`#2490 <https://github.com/dials/dials/issues/2490>`_)
- ``dials.estimate_resolution``: Improved CC½ fitting by using a weighted tanh fit. (`#2499 <https://github.com/dials/dials/issues/2499>`_)
- ``dials.find_spots`` and ``dials.integrate`` now accept ``exclude_images_multiple=N``, which will be expanded to the appropriate ``exclude_images=`` definition that excludes every Nth image from the data set. This is intended for the situation where a scan of diffraction images is regularly interrupted by calibration images. (`#2511 <https://github.com/dials/dials/issues/2511>`_)
- ``dials.ssx_integrate``: Add ``mosaicity_max_limit=`` setting, to control what mosaicity limit is considered unphysically large. (`#2518 <https://github.com/dials/dials/issues/2518>`_)
- ``dials.ssx_integrate``: Add ``max_cell_volume_change_fraction=`` parameter to catch highly divergent cell refinements (`#2521 <https://github.com/dials/dials/issues/2521>`_)


Bugfixes
--------

- ``dials.cosym``: Fix handling of refinement weights, for sparse data with a low number of common reflections. This improves the handling of serial still data. (`#2503 <https://github.com/dials/dials/issues/2503>`_)
- ``dials.reflection_viewer``: Fix a wsPython-related type error. (`#2506 <https://github.com/dials/dials/issues/2506>`_)
- ``dials.sequence_to_stills``: Now works with statically-refined experiments. (`#2516 <https://github.com/dials/dials/issues/2516>`_)
- ``dials.sequence_to_stills``: Generate identifiers for the split experiments. (`#2517 <https://github.com/dials/dials/issues/2517>`_)
- ``dials.sequence_to_stills``: No longer crashes when operating on a scan that does not start at image 1. (`#2520 <https://github.com/dials/dials/issues/2520>`_)
- ``dials.symmetry``: Print reindexed unit cell in log output. (`#2523 <https://github.com/dials/dials/issues/2523>`_)
- ``dials.index``: Correctly index multiple sweeps for multi-sweep indexed input. (`#2534 <https://github.com/dials/dials/issues/2534>`_)
- ``dials.ssx_integrate``: Fix bug in foreground mask calculation, affects datasets where wavelength is much bigger or smaller than 1 Å. (`#2535 <https://github.com/dials/dials/issues/2535>`_)
- ``dials.index``: Fix lattice counting when indexing with >1 sweep and >1 lattice (`#2539 <https://github.com/dials/dials/issues/2539>`_)
- ``dials.cosym/dials.scale``: Handle potential quasi_normalisation failures (for serial data). (`#2540 <https://github.com/dials/dials/issues/2540>`_)
- Bootstrap ``--cmake``: Fix configuration when CCP4 is already loaded. Previously, this could cause a failure as it tried to resolve the CCP4 copy of HDF5 in order to determine version. (`#2544 <https://github.com/dials/dials/issues/2544>`_)


Deprecations and Removals
-------------------------

- ``dials.ssx_integrate``: Remove duplicate named parameter ``output.shoeboxes=`` that has no effect (`#2537 <https://github.com/dials/dials/issues/2537>`_)


Misc
----

- `#2446 <https://github.com/dials/dials/issues/2446>`_, `#2465 <https://github.com/dials/dials/issues/2465>`_, `#2492 <https://github.com/dials/dials/issues/2492>`_, `#2493 <https://github.com/dials/dials/issues/2493>`_, `#2494 <https://github.com/dials/dials/issues/2494>`_, `#2501 <https://github.com/dials/dials/issues/2501>`_, `#2504 <https://github.com/dials/dials/issues/2504>`_, `#2508 <https://github.com/dials/dials/issues/2508>`_, `#2514 <https://github.com/dials/dials/issues/2514>`_, `#2524 <https://github.com/dials/dials/issues/2524>`_, `#2525 <https://github.com/dials/dials/issues/2525>`_


Dials 3.17 (2023-11-03)
=======================

Features
--------

- Remove circular dependencies between DIALS and ``cctbx.xfel`` by using the new ``serialtbx``. (`#2404 <https://github.com/dials/dials/issues/2404>`_)
- ``dials.import``: Enable use of an image range selection when importing a still sequence (`#2490 <https://github.com/dials/dials/issues/2490>`_)
- ``dials.estimate_resolution``: Improved cc1/2 fitting by using a weighted tanh fit (`#2499 <https://github.com/dials/dials/issues/2499>`_)
- ``dials.find_spots`` and ``dials.integrate`` now accept a parameter ``exclude_images_multiple=N``, which will be expanded to the appropriate ``exclude_images`` definition that excludes every Nth image from the data set. This is intended for the situation where a scan of diffraction images is regularly interrupted by calibration images. (`#2511 <https://github.com/dials/dials/issues/2511>`_)
- ``dials.ssx_integrate``: Allow setting of the mosaicity limit that is considered unphysically large, with mosaicity_max_limit option (`#2518 <https://github.com/dials/dials/issues/2518>`_)


Bugfixes
--------

- ``dials.cosym``: Fix a bug affecting refinement weights for sparse data with a low number of common reflections. This improves the handling of serial still data. (`#2503 <https://github.com/dials/dials/issues/2503>`_)
- ``dials.reflection_viewer``: Fix wxpython-related type error. (`#2506 <https://github.com/dials/dials/issues/2506>`_)
- The ``exclude_images_multiple`` parameter is adjusted to work correctly for scans that do not start with image number 1. (`#2514 <https://github.com/dials/dials/issues/2514>`_)
- ``dials.sequence_to_stills``: allow to work with statically-refined experiments. (`#2516 <https://github.com/dials/dials/issues/2516>`_)
- ``dials.sequence_to_stills``: Generate identifiers for the split experiments. (`#2517 <https://github.com/dials/dials/issues/2517>`_)
- ``dials.sequence_to_stills``: fix crash when operating on a scan that does not start at image 1. (`#2520 <https://github.com/dials/dials/issues/2520>`_)
- ``dials.ssx_integrate``: Add max_cell_volume_change_fraction parameter to catch highly divergent cell refinements (`#2521 <https://github.com/dials/dials/issues/2521>`_)
- ``dials.symmetry``: Print reindexed unit cell in log output (`#2523 <https://github.com/dials/dials/issues/2523>`_)
- ``dials.index``: Correctly index multiple sweeps for multi-sweep indexed input (`#2534 <https://github.com/dials/dials/issues/2534>`_)
- ``dials.ssx_integrate``: Fix bug in foreground mask calculation, affects datasets where wavelength is much bigger or smaller than 1. (`#2535 <https://github.com/dials/dials/issues/2535>`_)
- ``dials.ssx_integrate``: Remove duplicate named parameter output.shoeboxes that has no effect (`#2537 <https://github.com/dials/dials/issues/2537>`_)
- ``dials.cosym/dials.scale``: Handle potential quasi_normalisation failures (serial data) (`#2540 <https://github.com/dials/dials/issues/2540>`_)
- Bootstrap ``--cmake``: Fix configuration when CCP4 is already loaded. Previously, this could cause a failure as it tried to resolve the CCP4 copy of HDF5 in order to determine version. (`#2544 <https://github.com/dials/dials/issues/2544>`_)


Misc
----

- `#2446 <https://github.com/dials/dials/issues/2446>`_, `#2465 <https://github.com/dials/dials/issues/2465>`_, `#2492 <https://github.com/dials/dials/issues/2492>`_, `#2493 <https://github.com/dials/dials/issues/2493>`_, `#2494 <https://github.com/dials/dials/issues/2494>`_, `#2501 <https://github.com/dials/dials/issues/2501>`_, `#2504 <https://github.com/dials/dials/issues/2504>`_, `#2508 <https://github.com/dials/dials/issues/2508>`_, `#2524 <https://github.com/dials/dials/issues/2524>`_, `#2525 <https://github.com/dials/dials/issues/2525>`_


DIALS 3.16.1 (2023-09-05)
=========================

Bugfixes
--------

- ``dials.merge``: Fix potential for crash with ``r_free_flags.extend=True``, if there is no new flags to extend. (`#2491 <https://github.com/dials/dials/issues/2491>`_)


DIALS 3.16.0 (2023-08-14)
=========================

Features
--------

- ``dials.scale``: Absorption correction is now switched off if all scattering angles are small (cos(2ϑ) ~ 1). (`#2455 <https://github.com/dials/dials/issues/2455>`_)
- ``dials.export_bitmaps``: Set ``resolution_rings.d_spacings`` to manually specify the size of resolution rings. (`#2461 <https://github.com/dials/dials/issues/2461>`_)
- ``dials.index``: Experiment number is now shown when indexing multiple experiments. (`#2470 <https://github.com/dials/dials/issues/2470>`_)
- ``dials.export``: You can now specifying a wavelength tolerance when exporting multi-wavelength data to mtz. (`#2471 <https://github.com/dials/dials/issues/2471>`_)
- Support for building a CMake developer environment with ``bootstrap.py --cmake``. This should be considered experimental, for now. (`#2474 <https://github.com/dials/dials/issues/2474>`_)
- ``dials.merge``: Added stereographic projections to the html report, if merging more than one experiment. (`#2476 <https://github.com/dials/dials/issues/2476>`_)


Bugfixes
--------

- ``dials.combine_experiments``: No longer requires reflection data. (`#2323 <https://github.com/dials/dials/issues/2323>`_)
- ``dials.model_background``: Works again for data sets that do not start from image "1". (`#2405 <https://github.com/dials/dials/issues/2405>`_)
- Updated DIALS to imageset slicing properly using 0-based indices. Previously, the slice had to be accessed with the original image indices. (`#2411 <https://github.com/dials/dials/issues/2411>`_)
- ``dials.reindex``: Fix crash reindexing against a reference file, for more than one input file. (`#2439 <https://github.com/dials/dials/issues/2439>`_)
- ``dials.index``: ``indexing.refinement_protocol.mode=None`` now applies to both the scans and stills indexers. (`#2456 <https://github.com/dials/dials/issues/2456>`_)
- ``dials.index``: Now fails with a useful error message when there are too few spots for indexing. (`#2457 <https://github.com/dials/dials/issues/2457>`_)
- ``dials.ssx_integrate``: Address issue with ellipsoid angular models by replacing with new, correct, r-dependent models. (`#2463 <https://github.com/dials/dials/issues/2463>`_)
- ``dials.import``: Now ensures manual geometry updates are only applied once to each model. This ensures ``invert_rotation_axis=True`` will only invert the rotation axis once. (`#2469 <https://github.com/dials/dials/issues/2469>`_)
- ``dials.image_viewer``: Fix crash when the detector is rotated more than 90° from the beam direction. (`#2479 <https://github.com/dials/dials/issues/2479>`_)
- ``dials.reciprocal_lattice_viewer``: Handle cases where a beam centre cannot be calculated, such as a detector with a high 2θ swing angle. (`#2483 <https://github.com/dials/dials/issues/2483>`_)


Improved Documentation
----------------------

- Update data files documentation to refer to MessagePack rather than pickle format. (`#2444 <https://github.com/dials/dials/issues/2444>`_)


Deprecations and Removals
-------------------------

- ``dials.image_viewer``: The non-functional "Show mask" checkbox has been removed from the Mask tool. (`#2458 <https://github.com/dials/dials/issues/2458>`_)
- ``dials.image_viewer``: The basis vector display has been removed. (`#2484 <https://github.com/dials/dials/issues/2484>`_)


Misc
----

- `#2390 <https://github.com/dials/dials/issues/2390>`_, `#2433 <https://github.com/dials/dials/issues/2433>`_, `#2434 <https://github.com/dials/dials/issues/2434>`_, `#2435 <https://github.com/dials/dials/issues/2435>`_, `#2437 <https://github.com/dials/dials/issues/2437>`_, `#2438 <https://github.com/dials/dials/issues/2438>`_, `#2445 <https://github.com/dials/dials/issues/2445>`_, `#2462 <https://github.com/dials/dials/issues/2462>`_, `#2466 <https://github.com/dials/dials/issues/2466>`_, `#2472 <https://github.com/dials/dials/issues/2472>`_, `#2482 <https://github.com/dials/dials/issues/2482>`_, `#2489 <https://github.com/dials/dials/issues/2489>`_


DIALS 3.15.1 (2023-06-29)
=========================

Bugfixes
--------

- ``dials.export_bitmaps``: Fix the ``resolution_rings.fontsize=`` feature to work on Mac, and more reliably across platforms. (`#2441 <https://github.com/dials/dials/issues/2441>`_)


DIALS 3.15.0 (2023-06-13)
=========================

Features
--------

- ``dials.refine``: Independent groups of experiments are now refined separately by default. This may not then produce identical output, but results are believed to be of equivalent quality. The previous behaviour can be selected by using ``separate_independent_sets=False``. (`#2336 <https://github.com/dials/dials/issues/2336>`_)
- Python 3.11 is now supported as a bootstrap target. (`#2338 <https://github.com/dials/dials/issues/2338>`_)
- Support ``IOBS/SIGIOBS`` and ``FOBS/SIGFOBS`` columns when reading reference intensities from MTZ files. (`#2386 <https://github.com/dials/dials/issues/2386>`_)
- ``dials.export_bitmaps``: Optionally overlay ice rings. (`#2393 <https://github.com/dials/dials/issues/2393>`_)
- ``export_bitmaps`` API code has been separated out of CLI script into the new ``dials.util.export_bitmaps`` module. (`#2417 <https://github.com/dials/dials/issues/2417>`_)
- ``dials.reindex``: Allow reindexing against a reference file (mtz, pdb, cif) from which intensities can be extracted/generated (`#2420 <https://github.com/dials/dials/issues/2420>`_)
- ``dials.refine``: new option ``auto_reduction.ignore=True`` allows avoiding the parameter auto-reduction procedure for special cases. (`#2431 <https://github.com/dials/dials/issues/2431>`_)


Bugfixes
--------

- ``dials.export_bitmaps``: Correctly display the requested number of resolution rings. (`#2397 <https://github.com/dials/dials/issues/2397>`_)
- More directly support cctbx bootstrap environments without fixed-version dependencies. (`#2409 <https://github.com/dials/dials/issues/2409>`_)
- ``dials.scale``: Fix crash in html plotting for scaling and filtering when R-pim/R-merge values contain 'None' values (`#2410 <https://github.com/dials/dials/issues/2410>`_)
- ``dials.scale``, ``dials.cosym``: Use bulk solvent model when calculating intensities from a reference file. (`#2414 <https://github.com/dials/dials/issues/2414>`_)
- ``dials.integrate``: Fix multiprocessing on Windows with Python 3.11. (`#2426 <https://github.com/dials/dials/issues/2426>`_)
- ``dials.ssx_integrate``: Fix error in partiality calculation for ellipsoid integration algorithm. (`#2429 <https://github.com/dials/dials/issues/2429>`_)
- ``dials.import``: Don't print filename template for every image of an ``ImageSequence``. (`#2432 <https://github.com/dials/dials/issues/2432>`_)


Deprecations and Removals
-------------------------

- Moved ``dials::af::flex_table`` and associated tooling to dxtbx. Functions that are specific to ``reflection_table`` have been moved to ``reflection_table_suite.h``. (`#2113 <https://github.com/dials/dials/issues/2113>`_)
- ``dials.find_spots``: Remove ``force_2d=`` option, which did not work as expected. (`#2422 <https://github.com/dials/dials/issues/2422>`_)


Misc
----

- `#2402 <https://github.com/dials/dials/issues/2402>`_, `#2403 <https://github.com/dials/dials/issues/2403>`_, `#2419 <https://github.com/dials/dials/issues/2419>`_, `#2424 <https://github.com/dials/dials/issues/2424>`_, `#2425 <https://github.com/dials/dials/issues/2425>`_, `#2427 <https://github.com/dials/dials/issues/2427>`_, `#2428 <https://github.com/dials/dials/issues/2428>`_


DIALS 3.14.2 (2023-05-16)
=========================

No significant changes.


DIALS 3.14.1 (2023-04-20)
=========================

Bugfixes
--------

- ``dials.export_bitmaps``: Correctly scale resolution rings with image ``binning=`` set. (`#2392 <https://github.com/dials/dials/issues/2392>`_)
- Revert GTK-message-supression change that was causing ``reciprocal_lattice_viewer`` errors. (`#2395 <https://github.com/dials/dials/issues/2395>`_)


DIALS 3.14.0 (2023-04-12)
=========================

Features
--------

- ``dials.damage_analysis``: Add ``dose_group_size=`` option to allow generation of plots for subsets of the data based on dose. (`#1903 <https://github.com/dials/dials/issues/1903>`_)
- ``dials.import``/``dials.export``: The image template is now reported in the log file, for rotation data. (`#2278 <https://github.com/dials/dials/issues/2278>`_)
- Refinement tests now use publicly available data, though the ``dials-data`` package. (`#2343 <https://github.com/dials/dials/issues/2343>`_)
- ``dials.merge``: Now suggests a resolution limit, based on a fit to CC½. (`#2351 <https://github.com/dials/dials/issues/2351>`_)
- ``dials.merge``: Save merged half datasets to the MTZ file when ``additional_stats=True``. (`#2353 <https://github.com/dials/dials/issues/2353>`_)
- ``dials.export_bitmaps``: Option to render resolution rings. (`#2360 <https://github.com/dials/dials/issues/2360>`_)
- ``dials.rs_mapper``: Now works with multi-panel detectors. (`#2362 <https://github.com/dials/dials/issues/2362>`_)
- ``dials.ssx_integrate``: Add support for multi-panel detectors, for the ellipsoid model. (`#2368 <https://github.com/dials/dials/issues/2368>`_)
- ``dials.find_spots``: Add the ``exclude_images=`` option. Previously this was only on ``dials.import``. (`#2374 <https://github.com/dials/dials/issues/2374>`_)
- ``dials.merge``: Generate ``"FreeR_flag"`` column in output mtz, optionally copying from a reference mtz containing a ``FreeR_flag`` column. (`#2378 <https://github.com/dials/dials/issues/2378>`_)


Bugfixes
--------

- ``dials.refine``: Does not call ``sys.exit`` when running as a library function. (`#2325 <https://github.com/dials/dials/issues/2325>`_)
- ``dials.generate_mask``: Fix a failure when the direct beam intersected the exact centre of a pixel. (`#2329 <https://github.com/dials/dials/issues/2329>`_)
- Providing invalid PHIL parameters no longer causes a traceback. (`#2331 <https://github.com/dials/dials/issues/2331>`_)
- ``dials.index``: Handle possibility of Krivy-Gruber iteration limit exceeded in Niggli cell reduction. (`#2332 <https://github.com/dials/dials/issues/2332>`_)
- CPU allocation limits will now be correctly inherited from the Slurm cluster scheduler. (`#2335 <https://github.com/dials/dials/issues/2335>`_)
- ``dials.ssx_integrate``: Fix cases where the background and centroid algorithms choices were not being used with ellipsoid-model integration. (`#2337 <https://github.com/dials/dials/issues/2337>`_)
- ``dials.assign_experiment_identifiers`` now works with experiment lists without identifiers. (`#2342 <https://github.com/dials/dials/issues/2342>`_)
- Test data files are no longer patched to correct their trusted range values, and usages of the deprecated ``DataBlock`` have been removed. (`#2346 <https://github.com/dials/dials/issues/2346>`_)
- ``dials.ssx_integrate``: catch rare ``ValueError`` in ellipsoid profile model refinement (`#2349 <https://github.com/dials/dials/issues/2349>`_)
- ``dials.symmetry``: Handle cases when an experiment has no scan object. (`#2350 <https://github.com/dials/dials/issues/2350>`_)
- ``dials.ssx_index``: Ensure file paths are split correctly on Windows. (`#2358 <https://github.com/dials/dials/issues/2358>`_)
- Replace the deprecated ``imp`` library with ``importlib``. (`#2363 <https://github.com/dials/dials/issues/2363>`_)
- ``dials.rs_mapper``: Fix bug causing the image data to be accessed with out-of-bounds indices. (`#2364 <https://github.com/dials/dials/issues/2364>`_)
- ``dials.image_viewer``: Fix error when using newer versions of wxPython. (`#2365 <https://github.com/dials/dials/issues/2365>`_)
- DIALS viewers: Suppress terminal output saying "``Gtk-Message``". (`#2366 <https://github.com/dials/dials/issues/2366>`_)
- Prevent occasional failures testing constrained geometry refinement. (`#2367 <https://github.com/dials/dials/issues/2367>`_)
- Slicing a ``reflection_table`` does not lose the experiment identifiers map. (`#2371 <https://github.com/dials/dials/issues/2371>`_)
- ``dials.integrate``: Update ``exclude_images=`` syntax to match ``dials.import`` and ``dials.find_spots``. The older form is now deprecated. (`#2374 <https://github.com/dials/dials/issues/2374>`_)
- ``reflection_table.match_by_hkle``: Fix failure on Windows with pandas 2.0.0. (`#2383 <https://github.com/dials/dials/issues/2383>`_)


Improved Documentation
----------------------

- Clarified coordinate system conventions in https://dials.github.io/documentation/conventions.html (`#2373 <https://github.com/dials/dials/issues/2373>`_)


Misc
----

- `#2251 <https://github.com/dials/dials/issues/2251>`_, `#2324 <https://github.com/dials/dials/issues/2324>`_, `#2344 <https://github.com/dials/dials/issues/2344>`_, `#2345 <https://github.com/dials/dials/issues/2345>`_, `#2352 <https://github.com/dials/dials/issues/2352>`_, `#2356 <https://github.com/dials/dials/issues/2356>`_, `#2357 <https://github.com/dials/dials/issues/2357>`_, `#2361 <https://github.com/dials/dials/issues/2361>`_, `#2369 <https://github.com/dials/dials/issues/2369>`_, `#2372 <https://github.com/dials/dials/issues/2372>`_, `#2381 <https://github.com/dials/dials/issues/2381>`_, `#2385 <https://github.com/dials/dials/issues/2385>`_


DIALS 3.13.0 (2023-01-26)
=========================

Features
--------

- ``dev.dials.napari_rlv``: A reciprocal lattice viewer based on napari. This currently requires the ``napari`` module to be manually added into your DIALS installation. (`#2229 <https://github.com/dials/dials/issues/2229>`_)
- ``dials.stills_process``: Adds ``suppressed_logging=`` option, and minimal progress feedback. (`#2263 <https://github.com/dials/dials/issues/2263>`_)
- ``dials.refine``: Use an overall ``interval_width_degrees=`` parameter to set the default smoothness in scan-varying refinement for all models. This does not affect which models will be parameterised as scan-varying, which is controlled by their individual ``force_static=`` parameters. (`#2268 <https://github.com/dials/dials/issues/2268>`_)
- Use Python 3.10 by default when running bootstrap. (`#2272 <https://github.com/dials/dials/issues/2272>`_)
- ``dials.stills_process``: Added reflection subsampling. If ``reflection_subsampling.enable=True``, and an image fails to index, the reflections will be randomly subsampled, and indexing tried again. Reflections will be randomly subsampled in smaller amounts until a threshold is reached or the image successfully indexes. (`#2275 <https://github.com/dials/dials/issues/2275>`_)
- ``dials.merge``: Allow ``exclude_images=`` parameter, as in ``dials.scale``. (`#2280 <https://github.com/dials/dials/issues/2280>`_)
- ``dials.scale``: Add ability to use a precalculated analytical correction as part of scaling models. (`#2313 <https://github.com/dials/dials/issues/2313>`_)
- Add ``additional_stats`` option to generate R-split statistic for stills data in ``dials.merge`` and ``dials.scale``. (`#2314 <https://github.com/dials/dials/issues/2314>`_)


Bugfixes
--------

- ``dials.integrate``: Fail for negative profile-fitting intensity variance, not zero variance. (`#2271 <https://github.com/dials/dials/issues/2271>`_)
- ``dials.import``: Fix ``convert_stills_to_sequences=`` option for h5 data formats. (`#2273 <https://github.com/dials/dials/issues/2273>`_)
- Slightly better support when handling empty reflection tables. (`#2281 <https://github.com/dials/dials/issues/2281>`_)
- ``dials.scale``: Catch rare crash when making summary stats table. (`#2284 <https://github.com/dials/dials/issues/2284>`_)
- Fix wxPython 4.2.0 type error affecting RLV and geometry viewer. (`#2287 <https://github.com/dials/dials/issues/2287>`_)
- Allow ``reflection_table.remove_on_experiment_identifiers()`` for an empty ``flex.reflection_table``. (`#2298 <https://github.com/dials/dials/issues/2298>`_)
- ``dials.cluster_unit_cell``: Correctly report the number of lattices in each cluster, instead of the number of clusters. (`#2300 <https://github.com/dials/dials/issues/2300>`_)
- ``dials.image_viewer``: Fix error when using newer versions of wxPython. (`#2306 <https://github.com/dials/dials/issues/2306>`_)
- ``dials.combine_experiments``: Prevent default ``clustering.max_clusters=None`` from raising error when clustering. (`#2311 <https://github.com/dials/dials/issues/2311>`_)
- ``dials.scale``: Fix crash when scaling against reference, and some datasets were removed during scaling. (`#2312 <https://github.com/dials/dials/issues/2312>`_)
- ``dials.integrate``: Clearer error message when kapton ``absorption_correction.apply=True``, but ``algorithm=None``. (`#2317 <https://github.com/dials/dials/issues/2317>`_)


Improved Documentation
----------------------

- ``integration/corrections.h``: Improve comments pertaining to the polarization correction. (`#2274 <https://github.com/dials/dials/issues/2274>`_)


Misc
----

- `#2266 <https://github.com/dials/dials/issues/2266>`_, `#2267 <https://github.com/dials/dials/issues/2267>`_, `#2269 <https://github.com/dials/dials/issues/2269>`_, `#2282 <https://github.com/dials/dials/issues/2282>`_, `#2289 <https://github.com/dials/dials/issues/2289>`_, `#2290 <https://github.com/dials/dials/issues/2290>`_, `#2293 <https://github.com/dials/dials/issues/2293>`_, `#2301 <https://github.com/dials/dials/issues/2301>`_, `#2303 <https://github.com/dials/dials/issues/2303>`_, `#2307 <https://github.com/dials/dials/issues/2307>`_, `#2308 <https://github.com/dials/dials/issues/2308>`_, `#2309 <https://github.com/dials/dials/issues/2309>`_, `#2310 <https://github.com/dials/dials/issues/2310>`_, `#2319 <https://github.com/dials/dials/issues/2319>`_, `#2321 <https://github.com/dials/dials/issues/2321>`_


DIALS 3.12.1 (2022-12-05)
=========================

No significant changes.


DIALS 3.12.0 (2022-10-31)
=========================

Features
--------

- ``dials.stills_process``: Add ``known_orientations=`` feature, useful for specifying previous processing results to use when reprocessing data, or for multiple detectors in an experiment. (`#2110 <https://github.com/dials/dials/issues/2110>`_)
- ``dials.filter_reflections``: Added ``remove_by_index=`` option. (`#2201 <https://github.com/dials/dials/issues/2201>`_)
- ``dials.rs_mapper``: Perform the calculation in parallel to improve speed. (`#2238 <https://github.com/dials/dials/issues/2238>`_)
- ``dials.background``: Add ``n_checkpoints=`` parameter, to run the analysis on evenly-spaced images. (`#2240 <https://github.com/dials/dials/issues/2240>`_)
- `dials.image_viewer`: Change default brightness to 10. (`#2254 <https://github.com/dials/dials/issues/2254>`_)
- The ``dials.ssx_index`` and ``dials.ssx_integrate`` programs are now considered stable, so have dropped the ``dev.`` prefix. (`#2265 <https://github.com/dials/dials/issues/2265>`_)


Bugfixes
--------

- ``dials.stills_process``: Fix crashes for raster scans of stills. (`#2128 <https://github.com/dials/dials/issues/2128>`_)
- ``dials.image_viewer``: Fix intensity readout mouseover, and beamcenter half-pixel errors. (`#2194 <https://github.com/dials/dials/issues/2194>`_)
- **trusted_range** is now defined throughout DIALS consistently as the *inclusive* range between the minimum and maximum trusted values - valid pixels are those less than or equal to the maximum trusted value and greater than or equal to the minimum trusted value. (`#2195 <https://github.com/dials/dials/issues/2195>`_)
- ``dials.find_rotation_axis``: Recover a plot that was broken by previous changes. (`#2225 <https://github.com/dials/dials/issues/2225>`_)
- ``dials.scale``: Restore consistent dataset id output numbering when a dataset is excluded. (`#2246 <https://github.com/dials/dials/issues/2246>`_)
- ``dials.reciprocal_lattice_viewer``: Set the maximum number of experiments to show in the selector toggles (default 15). Datasets with a number larger than this will not be individually selectable. (`#2248 <https://github.com/dials/dials/issues/2248>`_)
- ``dials.ssx_integrate``: Gracefully handle assertion error for particular bad data. (`#2264 <https://github.com/dials/dials/issues/2264>`_)


Improved Documentation
----------------------

- Removed references to pickle and json from PHIL config help strings. (`#2208 <https://github.com/dials/dials/issues/2208>`_)


Deprecations and Removals
-------------------------

- The deprecated ``dials.algorithms.symmetry.cosym.target.Target.get_sym_ops()`` function has been removed. Please use the ``Target.sym_ops`` property instead. The ``nproc`` argument to the ``Target`` constructor has also been removed.
  The deprecated ``dials.util.mp.parallel_map`` function has been removed.
  The deprecated ``dials.command_line.dials_import.ImageImporter`` has been removed. Please use ``do_import`` instead.
  The command ``dials.rl_csv`` has been removed. Similar functionality is available with ``dials.export format=json``.
  The command ``dials.find_shared_models`` has been removed. Use ``dials.show show_shared_modules=True ...`` instead. (`#2212 <https://github.com/dials/dials/issues/2212>`_)


Misc
----

- `#2230 <https://github.com/dials/dials/issues/2230>`_


DIALS 3.11.2 (2022-09-27)
=========================

Bugfixes
--------

- ``dials.scale``: Fix bug in intensity combination scoring for multi-sweep datasets, affecting midpoint test values. (`#2199 <https://github.com/dials/dials/issues/2199>`_)


DIALS 3.11.1 (2022-09-02)
=========================

Bugfixes
--------

- Revert default installation to Python 3.9, to avoid WXPython incompatibilities. (`#2216 <https://github.com/dials/dials/issues/2216>`_)


DIALS 3.11.0 (2022-08-24)
=========================

Features
--------

- ``dials.scale``: Added the ``reference=`` option, with support for using a cif data file as a scaling reference. This replaces the old usages of ``target_mtz=`` and ``target_model=``, which are now deprecated. (`#2148 <https://github.com/dials/dials/issues/2148>`_)
- ``dials.algorithms.clustering.unit_cell`` API: Include the linkage matrix in the returned ``ClusteringResult`` object. (`#2152 <https://github.com/dials/dials/issues/2152>`_)
- ``dials.cosym``: Allow use of a reference dataset, to consistently index against. (`#2154 <https://github.com/dials/dials/issues/2154>`_)
- ``dials.find_rotation_axis``: Dramatically improve execution time by performing the search in parallel. (`#2160 <https://github.com/dials/dials/issues/2160>`_)
- ``dials.image_viewer``: Score tool now accepts keyboard entry. (`#2162 <https://github.com/dials/dials/issues/2162>`_)
- ``dials.merge``: Add merging statistics to output html report, and the PHIL option ``output.json=`` to write this this to json. (`#2178 <https://github.com/dials/dials/issues/2178>`_)
- ``dials.symmetry``: Add option to score systematic absences using fourier analysis. Thanks to Kevin Dalton for contributing this feature. (`#2184 <https://github.com/dials/dials/issues/2184>`_)
- Unit cell clustering reports: Add plots of unit cell angle distribution. (`#2197 <https://github.com/dials/dials/issues/2197>`_)
- Use Python 3.10 by default when running bootstrap. (`#2206 <https://github.com/dials/dials/issues/2206>`_)
- ``dials.merge``: include merging statistics in json output (`#2207 <https://github.com/dials/dials/issues/2207>`_)


Bugfixes
--------

- ``dials.find_rotation_axis``: Consistently report rotation axis angle at start and end of the search. (`#2145 <https://github.com/dials/dials/issues/2145>`_)
- ``dials.generate_mask``: Make the log filename consistent with the program name. (`#2147 <https://github.com/dials/dials/issues/2147>`_)
- ``dials.image_viewer``: Automatically complete an active polygon mask when saving mask parameters. (`#2150 <https://github.com/dials/dials/issues/2150>`_)
- ``dials.refine``: Correct an error mapping constraint definitions to models. (`#2155 <https://github.com/dials/dials/issues/2155>`_)
- ``dials.background``: Check for empty experiments. (`#2163 <https://github.com/dials/dials/issues/2163>`_)
- ``dials.estimate_gain``: Now works on multiple experiments from a single image sequence. (`#2164 <https://github.com/dials/dials/issues/2164>`_)
- ``dials.powder_calibrate``: More general determination of background mask removal. (`#2165 <https://github.com/dials/dials/issues/2165>`_)
- Fix DIALS version reporting for release builds. (`#2166 <https://github.com/dials/dials/issues/2166>`_)
- Logging messages from dxtbx are now included in program output. (`#2171 <https://github.com/dials/dials/issues/2171>`_)
- ``dev.dials.ssx_integrate``: Handle potential assertion errors from FastMCD. (`#2179 <https://github.com/dials/dials/issues/2179>`_)
- ``dials.symmetry``: Ensure data for systematic absences check is in the correct setting for non-conventional minimum cells. (`#2183 <https://github.com/dials/dials/issues/2183>`_)
- ``dials.scale``: Fix bug in intensity combination scoring for multi-sweep datasets, affecting midpoint test values. (`#2199 <https://github.com/dials/dials/issues/2199>`_)
- ``dev.dials.ssx_index``: Fix crash when generating plotting data where rmsd values are almost all equivalent. (`#2203 <https://github.com/dials/dials/issues/2203>`_)
- ``dev.dials.ssx_integrate``: Fix divergent refinement bug when n_macro_cycles > 1. (`#2204 <https://github.com/dials/dials/issues/2204>`_)


Deprecations and Removals
-------------------------

- ``dials.scale``: The ``target_mtz=`` and ``target_model=`` options are deprecated. Please use ``reference=`` instead for both uses. (`#2148 <https://github.com/dials/dials/issues/2148>`_)


Misc
----

- `#2115 <https://github.com/dials/dials/issues/2115>`_, `#2138 <https://github.com/dials/dials/issues/2138>`_, `#2143 <https://github.com/dials/dials/issues/2143>`_, `#2144 <https://github.com/dials/dials/issues/2144>`_, `#2169 <https://github.com/dials/dials/issues/2169>`_, `#2180 <https://github.com/dials/dials/issues/2180>`_, `#2181 <https://github.com/dials/dials/issues/2181>`_, `#2185 <https://github.com/dials/dials/issues/2185>`_, `#2188 <https://github.com/dials/dials/issues/2188>`_


DIALS 3.10.3 (2022-08-02)
=========================

Bugfixes
--------

- ``dials.image_viewer``: Fix error after loading images with the "Load" button. (`#2168 <https://github.com/dials/dials/issues/2168>`_)
- ``dials.merge``: Fix crash for P-1 datasets. (`#2175 <https://github.com/dials/dials/issues/2175>`_)
- ``dials.export format=mtz``: Handle shared experiment models when converting to cambridge frame (`#2182 <https://github.com/dials/dials/issues/2182>`_)


DIALS 3.10.2 (2022-07-13)
=========================

Bugfixes
--------

- Fix DIALS version reporting for release builds. (`#2166 <https://github.com/dials/dials/issues/2166>`_)


DIALS 3.10.1 (2022-07-12)
=========================

Bugfixes
--------

- Refinement internals: Terminate Levenberg-Marquardt refinement if the objective is not expected to decrease. (`#2135 <https://github.com/dials/dials/issues/2135>`_)
- ``dials.scale``: Avoid edge case crashes when all reflections filtered out during reflection selection for a dataset in multi-dataset scaling (`#2146 <https://github.com/dials/dials/issues/2146>`_)
- ``dials.cosym``: Fix bug whereby the change of basis op was applied twice to the crystal model (`#2151 <https://github.com/dials/dials/issues/2151>`_)
- ``dials.ssx_index``: Fix potential numpy crash when trying to generate histograms from one datapoint (`#2156 <https://github.com/dials/dials/issues/2156>`_)
- ``dials.image_viewer``: Allow viewing still collections where some images have no reflections. (`#2157 <https://github.com/dials/dials/issues/2157>`_)
- ``dev.dials.ssx_index``: Fix cluster reporting for only one image indexed, fix crash when all images filtered out (`#2159 <https://github.com/dials/dials/issues/2159>`_)


DIALS 3.10.0 (2022-06-09)
=========================

Features
--------

- The DIALS code package now uses ``src/`` layout. You will need to at least ``libtbx.refresh`` if on a development install. (`#2077 <https://github.com/dials/dials/issues/2077>`_)
- ``dials.cosym``: Enable skipping of unit cell clustering by setting ``unit_cell_clustering.threshold`` parameter to 0 or None. (`#2058 <https://github.com/dials/dials/issues/2058>`_)
- ``dials.export``: Add extra unmerged data categories for mmcif output, conforming to the latest mmcif_pdbx.dic. Make v5 the default (rather than v5_next). (`#2078 <https://github.com/dials/dials/issues/2078>`_)
- ``dials.export``: MTZ files are now exported with geometry in the Cambridge frame. (`#2054 <https://github.com/dials/dials/issues/2054>`_)
- ``dials.index``: Performance improvements for serial indexing (``dials.stills-process``, ``dev.dials.ssx_index``). (`#2035 <https://github.com/dials/dials/issues/2035>`_)
- ``dials.merge``: New implementation of the French & Wilson (1978) algorithm for correction of negative intensities when estimating amplitudes. This implementation makes use of the standardized median as an M-estimator for the average intensity of resolution bins, which makes the procedure robust against the presence of very negative intensities. (`#2100 <https://github.com/dials/dials/issues/2100>`_)
- ``dials.powder_calibrate``: Add detector distance calibration. (`#2075 <https://github.com/dials/dials/issues/2075>`_)
- ``dials.refine``: New ``separate_images`` option performs outlier rejection on each image independently. (`#2036 <https://github.com/dials/dials/issues/2036>`_)
- ``dials.refine``: Parallelise outlier rejection to reduce overall run times. (`#1427 <https://github.com/dials/dials/issues/1427>`_)
- ``dials.refine``: Use sparse storage for scan-varying runs, reducing memory requirements and run times. (`#2022 <https://github.com/dials/dials/issues/2022>`_)
- ``dials.scale``: Allow use of a pdb model to calculate target intensities for scaling (phil option ``target_model``). (`#2053 <https://github.com/dials/dials/issues/2053>`_)
- ``dials.stills_process``: Validate command-line arguments to prevent confusion when there is a typo. (`#2106 <https://github.com/dials/dials/issues/2106>`_)
- ``dev.dials.ssx_index``, ``dev.dials.ssx_integrate``: Add option ``output.nuggets=``, which can be used to specify a directory to which in-process results are stored. (`#2114 <https://github.com/dials/dials/issues/2114>`_)
- Add a CMake build of DIALS. (`#2096 <https://github.com/dials/dials/issues/2096>`_)


Bugfixes
--------

- ``dials.cluster_unit_cell``: Modify test so that it runs on Windows. (`#2027 <https://github.com/dials/dials/issues/2027>`_)
- ``dials.cluster_unit_cells``: Correctly handle cases with only one input crystal. (`#2120 <https://github.com/dials/dials/issues/2120>`_)
- ``dials.combine_experiments``: Exit with a helpful error message, if experiments have the same identifiers. (`#2069 <https://github.com/dials/dials/issues/2069>`_)
- ``dials.export``: Fix crash for exporting ssx data. (`#2126 <https://github.com/dials/dials/issues/2126>`_)
- ``dials.find_bad_pixels``: Output a mask file, as expected in the phil scope. Remove unimplemented png output option. (`#2122 <https://github.com/dials/dials/issues/2122>`_)
- ``dials.image_viewer``: Fix downstream SEGV involving wxPython. (`#2134 <https://github.com/dials/dials/issues/2134>`_)
- ``dials.image_viewer``: Fixed user setting for ``show_beam_center=`` being overridden. (`#2103 <https://github.com/dials/dials/issues/2103>`_)
- ``dials.image_viewer``: the ``n_iqr`` value for ``radial_profile`` thresholding is now correctly handled. (`#2116 <https://github.com/dials/dials/issues/2116>`_)
- ``dials.import``: When trying to import a missing file, say which file was missing. Print a warning if trying to import with a wildcard and no files were found. (`#1863 <https://github.com/dials/dials/issues/1863>`_)
- ``dials.refine_bravais_settings``: Ensure that the reported reindexing operators correctly map the input symmetry to the given Bravais settings, regardless of whether the input symmetry was a primitive or non-primitive setting. (`#2105 <https://github.com/dials/dials/issues/2105>`_)
- ``dials.reindex``: Fail with a helpful error message when attempting to reindex to a left-handed cell. (`#1779 <https://github.com/dials/dials/issues/1779>`_)
- ``dials.scale``: If scaling against a target, do targeted outlier rejection. (`#2052 <https://github.com/dials/dials/issues/2052>`_)
- ``dials.scale``: Respect user supplied high resolution limit in summary table. (`#2118 <https://github.com/dials/dials/issues/2118>`_)
- ``dials.show``: Avoid crash when an experiment does not have an imageset. (`#2056 <https://github.com/dials/dials/issues/2056>`_)
- ``dev.dials.ssx_index``: Fix reporting of results for h5 files, skip indexing of an image if fewer than ``min_spots=`` strong spots. (default 10). (`#2055 <https://github.com/dials/dials/issues/2055>`_)
- ``dev.dials.ssx_index``: Handle case where the input ``strong.refl`` file has no spots for some images. (`#2039 <https://github.com/dials/dials/issues/2039>`_)
- ``dev.dials.ssx_integrate``: Correctly handle input data containing multiple imagesets. (`#2124 <https://github.com/dials/dials/issues/2124>`_)
- Handle reflection simulation case where test reflections could fail to generate. (`#2094 <https://github.com/dials/dials/issues/2094>`_)


Improved Documentation
----------------------

- ``dials.model_background``: Add help message and test program. (`#1109 <https://github.com/dials/dials/issues/1109>`_)
- Remove reference to ``nproc`` in tutorials where it is not needed. (`#2030 <https://github.com/dials/dials/issues/2030>`_)
- Update DPF3 part 2 tutorial, to match current output. (`#2030 <https://github.com/dials/dials/issues/2030>`_)
- Add ``dials.two_theta_refine`` to documentation (`#2061 <https://github.com/dials/dials/issues/2061>`_)


Deprecations and Removals
-------------------------

- ``dials.find_shared_models`` has been retired. The command will now redirect users to ``dials.show``, in combination with the ``show_shared_models=True`` option. This command stub will be removed in a future version. (`#1070 <https://github.com/dials/dials/issues/1070>`_)
- ``dials.integrate``: The unused ``background.algorithm=median`` has been removed. (`#2066 <https://github.com/dials/dials/issues/2066>`_)
- The API interface ``dials.command_line.dials_import.ImageImporter`` is now deprecated. Please use ``...dials_import.do_import`` instead. (`#2080 <https://github.com/dials/dials/issues/2080>`_)


Misc
----

- `#1973 <https://github.com/dials/dials/issues/1973>`_, `#2037 <https://github.com/dials/dials/issues/2037>`_, `#2038 <https://github.com/dials/dials/issues/2038>`_, `#2041 <https://github.com/dials/dials/issues/2041>`_, `#2043 <https://github.com/dials/dials/issues/2043>`_, `#2047 <https://github.com/dials/dials/issues/2047>`_, `#2051 <https://github.com/dials/dials/issues/2051>`_, `#2062 <https://github.com/dials/dials/issues/2062>`_, `#2065 <https://github.com/dials/dials/issues/2065>`_, `#2070 <https://github.com/dials/dials/issues/2070>`_, `#2071 <https://github.com/dials/dials/issues/2071>`_, `#2073 <https://github.com/dials/dials/issues/2073>`_, `#2074 <https://github.com/dials/dials/issues/2074>`_, `#2079 <https://github.com/dials/dials/issues/2079>`_, `#2081 <https://github.com/dials/dials/issues/2081>`_, `#2082 <https://github.com/dials/dials/issues/2082>`_, `#2083 <https://github.com/dials/dials/issues/2083>`_, `#2086 <https://github.com/dials/dials/issues/2086>`_, `#2087 <https://github.com/dials/dials/issues/2087>`_, `#2088 <https://github.com/dials/dials/issues/2088>`_, `#2089 <https://github.com/dials/dials/issues/2089>`_, `#2091 <https://github.com/dials/dials/issues/2091>`_, `#2092 <https://github.com/dials/dials/issues/2092>`_, `#2093 <https://github.com/dials/dials/issues/2093>`_, `#2095 <https://github.com/dials/dials/issues/2095>`_, `#2099 <https://github.com/dials/dials/issues/2099>`_, `#2101 <https://github.com/dials/dials/issues/2101>`_, `#2102 <https://github.com/dials/dials/issues/2102>`_, `#2104 <https://github.com/dials/dials/issues/2104>`_, `#2112 <https://github.com/dials/dials/issues/2112>`_, `#2119 <https://github.com/dials/dials/issues/2119>`_, `#2131 <https://github.com/dials/dials/issues/2131>`_, `#2133 <https://github.com/dials/dials/issues/2133>`_


DIALS 3.8.6 (2022-06-07)
========================

- Resolve xia2 installation issues for downstream packaging.


DIALS 3.8.5 (2022-06-01)
========================

Features
--------

- ``dials.stills_process``: validate command-line arguments to prevent confusion when there is a typo (`#2106 <https://github.com/dials/dials/issues/2106>`_)


Bugfixes
--------

- ``dials.show``: Fix display of unknown vector columns. (`#2048 <https://github.com/dials/dials/issues/2048>`_)
- ``dials.image_viewer``: Fixed user setting for ``show_beam_center=`` being overridden. (`#2103 <https://github.com/dials/dials/issues/2103>`_)


DIALS 3.9.2 (2022-05-09)
========================

Bugfixes
--------

- ``dials.show``: Fix display of unknown vector columns. (`#2048 <https://github.com/dials/dials/issues/2048>`_)


DIALS 3.8.4 (2022-04-01)
========================

Bugfixes
--------

- ``dials.scale``: Fix crash when a dataset is filtered out during the scaling process (issue #2045). (`#2045 <https://github.com/dials/dials/issues/2045>`_)


DIALS 3.9.1 (2022-03-31)
========================

Bugfixes
--------

- ``dials.scale``: Fix crash when a dataset is filtered out during the scaling process (issue #2045). (`#2045 <https://github.com/dials/dials/issues/2045>`_)


DIALS 3.9.0 (2022-03-14)
========================

Features
--------

- Bootstrap now allows creating a Python 3.10 environment. This should be considered experimental at this stage, and may fail because not all our dependencies have Python 3.10 support yet. (`#1866 <https://github.com/dials/dials/issues/1866>`_)
- ``dials.export``: Add SHELX ``.hkl`` file output. (`#1925 <https://github.com/dials/dials/issues/1925>`_)
- ``dials.background`` now writes to a log file. (`#1948 <https://github.com/dials/dials/issues/1948>`_)
- ``dials.cluster_unit_cell``: Add option ``output.clusters=True/False`` to generate output files for each cluster generated by splitting the dendrogram at the given ``threshold``. (`#1950 <https://github.com/dials/dials/issues/1950>`_)
- Add ``reflection_table.match_by_hkle`` method, to match reflections between tables that have the same miller index and entering flags. (`#1951 <https://github.com/dials/dials/issues/1951>`_)
- Add experimental ``dev.dials.ssx_integrate`` script for profile modelling and integration of SSX data, including `ellipsoid` profile modelling. (`#1974 <https://github.com/dials/dials/issues/1974>`_)
- ``dials.reindex``: Allow reindexing using multi-crystal reference data files. (`#1977 <https://github.com/dials/dials/issues/1977>`_)
- add flex.reflection_table.concat method, to concatenate a list of reflection tables, including handling their experiment identifiers and ids (`#1994 <https://github.com/dials/dials/issues/1994>`_)
- For data reduction programs, allow exclude_images option to take a single multi-sweep command, e.g. ``exclude_images=0:100:120,1:150:180`` (`#1996 <https://github.com/dials/dials/issues/1996>`_)
- Enable spot-finding threshold algorithms to use information about experimental models. (`#2001 <https://github.com/dials/dials/issues/2001>`_)
- ``dials.find_spots``: Added new ``spotfinder.threshold=radial_profile``
  threshold algorithm. This calculates an average background in 2θ shells,
  and identifies peak pixels at a user-controllable level above the
  background. This simple method is particularly appropriate for cases
  with strong rotationally-symmetric background, such as electron
  diffraction images. An optional blurring function helps to suppress
  noise peaks and to join split spots. (`#2009 <https://github.com/dials/dials/issues/2009>`_)
- ``dials.export``: Add `PETS 2 <http://pets.fzu.cz/>`_ exporting with ``format=pets``. This is used for processing electron diffraction data. (`#2014 <https://github.com/dials/dials/issues/2014>`_)
- New tool - ``dials.powder_calibrate`` to help calibrate the geometry of an electron powder pattern. (`#2016 <https://github.com/dials/dials/issues/2016>`_)
- Colours in plots:  Some of our plots and figures were still using the Matplotlib colour map Jet.  Matplotlib moved some time ago to using `the perceptually uniform colour map Viridis <https://matplotlib.org/stable/users/prev_whats_new/dflt_style_changes.html#colormap>`_, which was designed to be relatively colourblind-friendly, as its default.  In this version of DIALS, we too have moved to use Viridis for our Plotly plots.  With thanks to `Ammaar Saeed (ammsa23) <https://github.com/ammsa23>`_ for this change. (`#2026 <https://github.com/dials/dials/issues/2026>`_)


Bugfixes
--------

- ``dials.refine`` could in some rare cases introducing an unphysical beam polarization vector. (`#1939 <https://github.com/dials/dials/issues/1939>`_)
- ``dials.scale``: Fixes to properly handle partiality of ssx data (`#1965 <https://github.com/dials/dials/issues/1965>`_)
- This round includes modifications to enable Kapton absorption correction
  for higher angles of rotation of the Kapton tape. (`#1968 <https://github.com/dials/dials/issues/1968>`_)
- ``ThreadPool`` and ``Socket`` resources were not properly closed in ``dials.find_spots_client``. (`#1976 <https://github.com/dials/dials/issues/1976>`_)
- Fix regression in installer python 2 compatibility. (`#1990 <https://github.com/dials/dials/issues/1990>`_)
- Avoid using the ``uuid`` module, to avoid mpi errors on Python 3.8. (`#2000 <https://github.com/dials/dials/issues/2000>`_)
- ``dials.refine``: Avoid spike in memory usage while saving reflections. (`#2024 <https://github.com/dials/dials/issues/2024>`_)


Improved Documentation
----------------------

- Improved the "small molecule" tutorial by including the symmetry determination, scaling and export to e.g. SHELX format output. (`#1900 <https://github.com/dials/dials/issues/1900>`_)
- Add a new multi-crystal tutorial, discussing the analysis of Br-lysozyme microcrystal data with ``dials.cosym`` and ``xia2.multiplex``. (`#1960 <https://github.com/dials/dials/issues/1960>`_)
- Update the betalactamase tutorial to use automatic scan-varying refinement. (`#1971 <https://github.com/dials/dials/issues/1971>`_)
- Remove outdated developer tutorial. (`#2003 <https://github.com/dials/dials/issues/2003>`_)
- Documentation font has changed to sans-serif. (`#2010 <https://github.com/dials/dials/issues/2010>`_)


Deprecations and Removals
-------------------------

- Starting with this release DIALS requires a minimum Python version of 3.8. Bootstrap no longer allows the creation of Python 3.7 environments. (`#1866 <https://github.com/dials/dials/issues/1866>`_)


Misc
----

- `#1576 <https://github.com/dials/dials/issues/1576>`_, `#1930 <https://github.com/dials/dials/issues/1930>`_, `#1953 <https://github.com/dials/dials/issues/1953>`_, `#1966 <https://github.com/dials/dials/issues/1966>`_, `#1967 <https://github.com/dials/dials/issues/1967>`_, `#1972 <https://github.com/dials/dials/issues/1972>`_, `#1975 <https://github.com/dials/dials/issues/1975>`_, `#1978 <https://github.com/dials/dials/issues/1978>`_, `#1986 <https://github.com/dials/dials/issues/1986>`_, `#1989 <https://github.com/dials/dials/issues/1989>`_, `#1995 <https://github.com/dials/dials/issues/1995>`_, `#1999 <https://github.com/dials/dials/issues/1999>`_, `#2013 <https://github.com/dials/dials/issues/2013>`_, `#2015 <https://github.com/dials/dials/issues/2015>`_, `#2017 <https://github.com/dials/dials/issues/2017>`_, `#2018 <https://github.com/dials/dials/issues/2018>`_, `#2019 <https://github.com/dials/dials/issues/2019>`_, `#2020 <https://github.com/dials/dials/issues/2020>`_, `#2021 <https://github.com/dials/dials/issues/2021>`_, `#2023 <https://github.com/dials/dials/issues/2023>`_


DIALS 3.8.3 (2022-02-22)
========================

Bugfixes
--------

- ``dials.cosym``: Fix crash for edge case of a full dataset being excluded by the resolution filter (`#1993 <https://github.com/dials/dials/issues/1993>`_)
- ``dials.two_theta_refine``: Fix crash if running on scaled data with ``exclude_datasets=``. (`#2006 <https://github.com/dials/dials/issues/2006>`_)
- Fix downloads failing on MacOS with "426 Upgrade Required" (`#2012 <https://github.com/dials/dials/issues/2012>`_)


DIALS 3.8.2 (2022-02-07)
========================

No significant changes.


DIALS 3.8.1 (2022-01-25)
========================

Features
--------

- It is now possible to explicitly bootstrap all supported python versions. (`#1988 <https://github.com/dials/dials/issues/1988>`_)


Bugfixes
--------

- ``dials.refine``: Fix crash from recording of parameter correlations in ``history.json``. (`#1923 <https://github.com/dials/dials/issues/1923>`_)
- Correct reporting of phi angle in scan-varying model plots. (`#1929 <https://github.com/dials/dials/issues/1929>`_)
- ``dials.sequence_to_stills``: correct error in the crystal model for scans that do not start from image 1. (`#1933 <https://github.com/dials/dials/issues/1933>`_)
- `dials.show`: only show format class if meaningful (i.e. not Format or FormatMultiImage) (`#1981 <https://github.com/dials/dials/issues/1981>`_)


Improved Documentation
----------------------

- Update some Cosym and resolution-related PHIL descriptions. (`#1969 <https://github.com/dials/dials/issues/1969>`_)


DIALS 3.8.0 (2022-01-11)
========================

Features
--------

- ``dials.indexed_as_integrated``: manipulate an indexed reflection file to look as if it were summation integrated. This simply takes the spot intensities that have been indexed, assigns a resolution and sets the summation integrated flag. Allows symmetry analysis using intensities and scaling on indexed data for very rapid feedback data processing. (`#1912 <https://github.com/dials/dials/issues/1912>`_)
- All command line programs now allow passing -h argument multiple times to increase verbosity. (`#1920 <https://github.com/dials/dials/issues/1920>`_)
- Add ``Dockerfile`` to automatically build and push images on new releases. (`#1936 <https://github.com/dials/dials/issues/1936>`_)
- Add experimental ``dev.dials.ssx_index`` script to index a block of ssx data. (`#1955 <https://github.com/dials/dials/issues/1955>`_)
- Add alias ``dials.rlv`` for ``dials.reciprocal_lattice_viewer``, and ``dials.rbs`` for ``dials.refine_bravais_settings``. (`#1959 <https://github.com/dials/dials/issues/1959>`_)


Bugfixes
--------

- ``dials.import``: Support arbitrary P1 known unit cells. (`#1880 <https://github.com/dials/dials/issues/1880>`_)
- Utility fixes: Better handle cases of missing partiality data in reflection files. Treat missing resolution values as ``0.0`` instead of being empty. (`#1911 <https://github.com/dials/dials/issues/1911>`_)
- ``dials.scale``: Add missing "expids_and_image_ranges" information to the json output. This is required for some of the plots in ``dials.report`` output. (`#1913 <https://github.com/dials/dials/issues/1913>`_)
- Copy docker entrypoint script with exec permissions. (`#1940 <https://github.com/dials/dials/issues/1940>`_)
- Don't line-wrap command-line help messages (`#1954 <https://github.com/dials/dials/issues/1954>`_)
- `dials.check_indexing_symmetry`: correctly handle d_max parameter if left at default value when d_min set. (`#1957 <https://github.com/dials/dials/issues/1957>`_)
- Fixed bug that prevented the ability to plot absorption end of max and min due to Kapton (`#1962 <https://github.com/dials/dials/issues/1962>`_)


Improved Documentation
----------------------

- Updated MyD88 tutorial to make use of the new ``dials.find_rotation_axis`` command. (`#1885 <https://github.com/dials/dials/issues/1885>`_)
- Document the dxtbx convention for representing the goniostat rotation operator :math:`\mathbf{R}` on `the conventions page. <https://dials.github.io/documentation/conventions.html#the-dxtbx-goniometer-model>`_ of the online documentation. (`#1917 <https://github.com/dials/dials/issues/1917>`_)
- Update tutorial for DUI 2021.11.1. (`#1938 <https://github.com/dials/dials/issues/1938>`_)


Deprecations and Removals
-------------------------

- The `cosym nproc=` and ``dials.util.parallel_map`` warnings have been made more visible. (`#1909 <https://github.com/dials/dials/issues/1909>`_)


Misc
----

- `#1907 <https://github.com/dials/dials/issues/1907>`_, `#1908 <https://github.com/dials/dials/issues/1908>`_, `#1910 <https://github.com/dials/dials/issues/1910>`_, `#1928 <https://github.com/dials/dials/issues/1928>`_, `#1947 <https://github.com/dials/dials/issues/1947>`_


DIALS 3.7.2 (2021-12-02)
========================

Features
--------

- ``dials.integrate``: When determining available memory, take into account ``MemoryProvisioned`` from HTCondor machine ad if the ``_CONDOR_JOB_AD`` environment variable is set.
  ``nproc=auto``: Take into account ``CpusProvisioned`` from HTCondor machine ad. (`#1943 <https://github.com/dials/dials/issues/1943>`_)


Bugfixes
--------

- Read ``_CONDOR_JOB_AD`` not ``_CONDOR_MACHINE_AD`` (`#1945 <https://github.com/dials/dials/issues/1945>`_)


DIALS 3.7.1 (2021-11-17)
========================

Bugfixes
--------

- ``dials.export``: No longer allow (erroneous) MTZ export for multiple experiments with multiple space groups. (`#1915 <https://github.com/dials/dials/issues/1915>`_)
- ``dials.export``: No longer fails for XDS_ASCII and SADABS export with ``intensity=auto``. (`#1926 <https://github.com/dials/dials/issues/1926>`_)
- ``dials.report``: Fix broken json output option. Include more graphs in json output. (`#1932 <https://github.com/dials/dials/issues/1932>`_)


DIALS 3.7.0 (2021-11-01)
========================

Features
--------

- Bootstrap support for MacOS M1 platforms. (`#1841 <https://github.com/dials/dials/issues/1841>`_)
- New ``dials.find_rotation_axis`` program optimises the rotation axis from strong spot positions prior to indexing. (`#1884 <https://github.com/dials/dials/issues/1884>`_)
- ``dials.import``: Allow importing templates with no template characters. (`#1840 <https://github.com/dials/dials/issues/1840>`_)
- ``dials.stills_process``: Performance improvements in Kapton absorption correction and in rare cases of highly mosaic crystals. (`#1846 <https://github.com/dials/dials/issues/1846>`_)
- ``dials.image_viewer`` Coordinates are now given in fast, slow order. (`#1849 <https://github.com/dials/dials/issues/1849>`_)
- ``dials.image_viewer``: Crystal basis vectors are now shown in the same colour as their predictions. (`#1855 <https://github.com/dials/dials/issues/1855>`_)
- ``dials.image_viewer``: Add the option to display the rotation axis (`#1856 <https://github.com/dials/dials/issues/1856>`_)
- ``dials.image_viewer``: Draw resolution rings for curved detectors. (`#1899 <https://github.com/dials/dials/issues/1899>`_)
- ``dials.import``: Unhandled files are now by default ignored. This means that e.g. ``*.log`` files alongside images will no longer prevent a successful import. Set ``ignore_unhandled=False`` to restore the previous behaviour. (`#1881 <https://github.com/dials/dials/issues/1881>`_)
- ``dials.scale``: Allow fixing of a particular correction with e.g. ``physical.correction.fix=absorption``. (`#1883 <https://github.com/dials/dials/issues/1883>`_)
- Installer now accepts a ``--raw-prefix`` option to use the target destination directly, instead of in a ``dials-X.Y`` subdirectory. (`#1896 <https://github.com/dials/dials/issues/1896>`_)


Bugfixes
--------

- ``dials.compute_delta_cchalf``: Fix crash when only using passing dataset/group. (`#1892 <https://github.com/dials/dials/issues/1892>`_)
- ``dials.find_bad_pixels``: Pixel coordinates are now reported in row-major order, and mask value is now set to 16, which corresponds internally to "noisy pixel". (`#1876 <https://github.com/dials/dials/issues/1876>`_)
- ``dials.find_rotation_axis``: removed unused parameter ``optimise={True|False}``. (`#1898 <https://github.com/dials/dials/issues/1898>`_)
- ``dials.report``: Don't show otherwise empty sections. (`#1875 <https://github.com/dials/dials/issues/1875>`_)


Improved Documentation
----------------------

- Improvements to 3DED tutorials. (`#1850 <https://github.com/dials/dials/issues/1850>`_)
- SARS-CoV-2 main protease tutorial: process in C2 rather than I2 setting for consistency with published structures. (`#1854 <https://github.com/dials/dials/issues/1854>`_)
- Removed outdated lysozyme nanocrystals tutorial. (`#1877 <https://github.com/dials/dials/issues/1877>`_)
- Add an associated projects page to the website. (`#1893 <https://github.com/dials/dials/issues/1893>`_)


Deprecations and Removals
-------------------------

- Bootstrap no longer allows creating Python 3.6 environments. (`#1852 <https://github.com/dials/dials/issues/1852>`_)
- ``dials.util.mp``: deprecate ``parallel_map()`` function and remove previously deprecated ``preserve_exception_message=`` parameters. (`#1860 <https://github.com/dials/dials/issues/1860>`_)


Misc
----

- `#1851 <https://github.com/dials/dials/issues/1851>`_, `#1853 <https://github.com/dials/dials/issues/1853>`_, `#1862 <https://github.com/dials/dials/issues/1862>`_, `#1865 <https://github.com/dials/dials/issues/1865>`_, `#1867 <https://github.com/dials/dials/issues/1867>`_, `#1869 <https://github.com/dials/dials/issues/1869>`_, `#1882 <https://github.com/dials/dials/issues/1882>`_, `#1887 <https://github.com/dials/dials/issues/1887>`_, `#1888 <https://github.com/dials/dials/issues/1888>`_, `#1889 <https://github.com/dials/dials/issues/1889>`_, `#1891 <https://github.com/dials/dials/issues/1891>`_, `#1894 <https://github.com/dials/dials/issues/1894>`_, `#1902 <https://github.com/dials/dials/issues/1902>`_


DIALS 3.6.2 (2021-09-21)
========================

Bugfixes
--------

- ``dials.reciprocal_lattice_viewer``: In cases with multiple lattices, "Crystal Frame" now aligns all crystal frames, rather than just the first. Unindexed reflections are no longer shown in this mode. (`#1868 <https://github.com/dials/dials/issues/1868>`_)


DIALS 3.6.1 (2021-09-06)
========================

No significant changes.


DIALS 3.6.0 (2021-08-16)
========================

This is the last release to support Python 3.6. Future releases will require a
minimum of Python 3.7.

Features
--------

- DIALS bootstrap now creates a Python 3.9 environment by default (`#1735 <https://github.com/dials/dials/issues/1735>`_)
- New program: ``dials.reference_profile_viewer`` for viewing reference profiles dumped by ``dials.integrate`` when using the ``debug.reference.output=True`` option. (`#1759 <https://github.com/dials/dials/issues/1759>`_)
- ``dials.combine_experiments``: Unindexed reflections are now included in the combined output (`#1760 <https://github.com/dials/dials/issues/1760>`_)
- ``dials.image_viewer``: Image overlays are now accumulated over stacks of images (`#1750 <https://github.com/dials/dials/issues/1750>`_)
- ``dials.image_viewer``: Allow control of the basis vector scale from the settings window (`#1780 <https://github.com/dials/dials/issues/1780>`_)
- ``dials.image_viewer``: Better colour choice for text overlays. Labels will now be light grey on Black, or Dark grey on White. The previous settings were sometimes hard to read on narrow-contrast images. (`#1781 <https://github.com/dials/dials/issues/1781>`_)
- ``dials.merge``: Include DANO/SIGDANO columns in output .mtz when ``anomalous=True`` and ``truncate=True`` (`#1809 <https://github.com/dials/dials/issues/1809>`_)
- ``dials.reciprocal_lattice_viewer``: Show resolution on the "nearest point" label (`#1770 <https://github.com/dials/dials/issues/1770>`_)
- ``dials.reciprocal_lattice_viewer`` now shows the path to the reflections in the title bar (`#1771 <https://github.com/dials/dials/issues/1771>`_)
- ``dials.reciprocal_lattice_viewer``: The default marker size now scaled automatically based on the data density (`#1773 <https://github.com/dials/dials/issues/1773>`_)
- ``dials.scale``: Always enable absorption correction if the ``absorption_level=`` parameter is set. Previously it was only enabled for sweeps >= 60° or if ``absorption_correction=True``. (`#1793 <https://github.com/dials/dials/issues/1793>`_)
- ``dials.scale``: Allow a shared absorption correction between sweeps if using the physical model, with the option ``share.absorption=True``. Extra absorption correction plots have also been added; and multiple sweeps are now aligned to the same reference frame. (`#1811 <https://github.com/dials/dials/issues/1811>`_)
- API: ``...scaling_library.scaled_data_as_miller_array`` now sets wavelength in the returned ``miller.array`` (`#1808 <https://github.com/dials/dials/issues/1808>`_)
- ``reflection_table.match()`` now returns ``flex.size_t`` index arrays, instead of ``flex.int``. (`#1784 <https://github.com/dials/dials/issues/1784>`_)
- New bootstrap option: ``--conda`` to install with miniconda instead of micromamba. (`#1730 <https://github.com/dials/dials/issues/1730>`_)


Bugfixes
--------

- ``dials.combine_experiments``: Correctly preserve mapping to images. This affects ``dials.image_viewer`` and ``dial.reciprocal_lattice_viewer``. (`#1093 <https://github.com/dials/dials/issues/1093>`_)
- ``dials.compute_delta_cchalf``: Unwarranted precision in the output has been reduced (`#1751 <https://github.com/dials/dials/issues/1751>`_)
- ``dials.find_spots``: Fix counting of imagesets in histogram output (`#1827 <https://github.com/dials/dials/issues/1827>`_)
- ``dials.image_viewer``: Add buttons to clear unit cell and generic ring display (`#1777 <https://github.com/dials/dials/issues/1777>`_)
- ``dials.image_viewer``: Fix various minor behavioural bugs in the spot-finding and image type controls. (`#1796 <https://github.com/dials/dials/issues/1796>`_)
- ``dials.import``: Fail gracefully when `#` is missing from template. (`#1840 <https://github.com/dials/dials/issues/1840>`_)
- ``dials.integrate``: change default filename of debug reference profile to ``reference_profiles.pickle``. (`#1747 <https://github.com/dials/dials/issues/1747>`_)
- ``dials.integrate``: Change default configuration so that unintegrated reflections are not retained. This helps reduce memory usage of downstream tools. Set ``output_unintegrated_reflections=True`` to restore the previous behaviour. (`#1753 <https://github.com/dials/dials/issues/1753>`_)
- ``dials.integrate``: ensure imageset_ids are always output. Affected use of image viewer, reciprocal lattice viewer on multi-sweep data. (`#1762 <https://github.com/dials/dials/issues/1762>`_)
- ``dials.reciprocal_lattice_viewer``: When starting with ``black_background=False``, ensure the rotation axis and beam vector are displayed in black. (`#1540 <https://github.com/dials/dials/issues/1540>`_)
- ``dials.reciprocal_lattice_viewer``: More robust beam centre control that works for multiple panel detectors (`#1842 <https://github.com/dials/dials/issues/1842>`_)
- ``dials.refine_bravais_settings``: correctly report mI Bravais settings (`#1825 <https://github.com/dials/dials/issues/1825>`_)
- ``dials.split_experiments``: Update the imageset_id column in the output reflection files. (`#1792 <https://github.com/dials/dials/issues/1792>`_)
- Don't fail ``bootstrap.py`` if a submodule is missing a reference (`#1834 <https://github.com/dials/dials/issues/1834>`_)

- Correctly handle reflection ``imageset_id`` column in ``dials.scale``, ``dials.cosym``, and ``dials.symmetry``. (`#1763 <https://github.com/dials/dials/issues/1763>`_)

Improved Documentation
----------------------

- ``dials.anvil_correction``: Made a small improvement to the developer documentation. (`#1788 <https://github.com/dials/dials/issues/1788>`_)
- Fix help string for ``best_monoclinic_beta=`` parameter (for ``dials.cosym``, ``dials.refine_bravais_settings`` and ``dials.symmetry``) (`#1833 <https://github.com/dials/dials/issues/1833>`_)
- Added a new tutorial on 3DED/MicroED data processing. (`#1837 <https://github.com/dials/dials/issues/1837>`_)
- Add a "Getting started" page to the documentation on the website. (`#1844 <https://github.com/dials/dials/issues/1844>`_)
- Add a tutorial on processing small molecule 3DED data. (`#1847 <https://github.com/dials/dials/issues/1847>`_)


Deprecations and Removals
-------------------------

- The previously deprecated ``dials.resolutionizer`` command has been removed. Please use ``dials.estimate_resolution`` instead. (`#1330 <https://github.com/dials/dials/issues/1330>`_)
- The previously deprecated ``dials.refine`` parameter ``trim_scan_edges`` has been removed. Please use ``scan_margin=...`` instead. (`#1374 <https://github.com/dials/dials/issues/1374>`_)
- The previously deprecated ``Spotfinder()()`` interface has been removed. Please use ``Spotfinder().find_spots()`` instead. (`#1484 <https://github.com/dials/dials/issues/1484>`_)
- The previously deprecated ``dials.util.masking.MaskGenerator`` has been removed. Please use ``dials.util.masking.generate_mask`` instead. (`#1569 <https://github.com/dials/dials/issues/1569>`_)
- The bootstrap option ``--mamba`` has become the default and will be removed in the future. (`#1730 <https://github.com/dials/dials/issues/1730>`_)
- ``dials.anvil_correction``:  Drop compatibility support for SciPy < 1.4 (`#1787 <https://github.com/dials/dials/issues/1787>`_)


Misc
----

- `#1746 <https://github.com/dials/dials/issues/1746>`_, `#1733 <https://github.com/dials/dials/issues/1733>`_, `#1752 <https://github.com/dials/dials/issues/1752>`_, `#1755 <https://github.com/dials/dials/issues/1755>`_, `#1756 <https://github.com/dials/dials/issues/1756>`_, `#1764 <https://github.com/dials/dials/issues/1764>`_, `#1767 <https://github.com/dials/dials/issues/1767>`_, `#1772 <https://github.com/dials/dials/issues/1772>`_, `#1783 <https://github.com/dials/dials/issues/1783>`_, `#1789 <https://github.com/dials/dials/issues/1789>`_, `#1791 <https://github.com/dials/dials/issues/1791>`_, `#1794 <https://github.com/dials/dials/issues/1794>`_, `#1795 <https://github.com/dials/dials/issues/1795>`_, `#1799 <https://github.com/dials/dials/issues/1799>`_, `#1802 <https://github.com/dials/dials/issues/1802>`_, `#1804 <https://github.com/dials/dials/issues/1804>`_, `#1806 <https://github.com/dials/dials/issues/1806>`_, `#1807 <https://github.com/dials/dials/issues/1807>`_, `#1812 <https://github.com/dials/dials/issues/1812>`_, `#1816 <https://github.com/dials/dials/issues/1816>`_, `#1817 <https://github.com/dials/dials/issues/1817>`_, `#1823 <https://github.com/dials/dials/issues/1823>`_, `#1830 <https://github.com/dials/dials/issues/1830>`_, `#1835 <https://github.com/dials/dials/issues/1835>`_, `#1836 <https://github.com/dials/dials/issues/1836>`_, `#1839 <https://github.com/dials/dials/issues/1839>`_


DIALS 3.5.4 (2021-07-27)
========================

Bugfixes
--------

- ``dials.stills_process``: Fix case where imagesets and experiment filenames could potentially disagree (`#1814 <https://github.com/dials/dials/issues/1814>`_)
- ``dials.scale``: Fix incorrect output files, for targeted scaling with more than one target dataset. (`#1815 <https://github.com/dials/dials/issues/1815>`_)
- ``dials.image_viewer``: Fix opening datasets with ``load_models=False`` (`#1818 <https://github.com/dials/dials/issues/1818>`_)


DIALS 3.5.3 (2021-07-12)
========================

Bugfixes
--------

- ``dials.image_viewer``: Fix the ``basis_vector_scale=`` parameter. (`#1769 <https://github.com/dials/dials/issues/1769>`_)


DIALS 3.5.2 (2021-06-28)
========================

Bugfixes
--------

- ``dials.image_viewer``: Fix display of spotfinding intermediates (threshold, dispersion, etc) when viewing multiple still experiments (`#1734 <https://github.com/dials/dials/issues/1734>`_)
- ``dials.image_viewer``: Stacking images no longer gives incorrect results for multi-sweep data beyond the first sweep (`#1758 <https://github.com/dials/dials/issues/1758>`_)


DIALS 3.5.1 (2021-06-14)
========================

No significant changes.


DIALS 3.5.0 (2021-05-27)
========================

Features
--------

- ``dials.integrate``: Avoid crash when data is too large to process in memory, by splitting into subsets (`#1392 <https://github.com/dials/dials/issues/1392>`_)
- New bootstrap options: ``--mamba`` to install with `micromamba <https://github.com/mamba-org/mamba#micromamba>`_, and ``--clean`` to remove installation caches immediately after completion. (`#1676 <https://github.com/dials/dials/issues/1676>`_)
- ``dials.find_spots_server``: Faster filtering of reflections by resolution (`#1680 <https://github.com/dials/dials/issues/1680>`_)
- ``dials.scale``: Add option ``error_model.grouping=`` to control refinement of either individual or grouped error models during scaling (`#1684 <https://github.com/dials/dials/issues/1684>`_)
- ``dials.scale``: Added ``physical.absorption_level=[low|medium|high]`` option for automatic setting of suitable absorption correction parameters. (`#1688 <https://github.com/dials/dials/issues/1688>`_)
- ``dials.cosym``: Significantly faster calculation of Rij matrix of pairwise correlation coefficients (`#1693 <https://github.com/dials/dials/issues/1693>`_)
- ``dials.sort_reflections`` and ``dials.merge_reflection_lists`` are now available without a ``dev.`` prefix. (`#1703 <https://github.com/dials/dials/issues/1703>`_)
- New command: ``dials.find_bad_pixels`` to identify pixels which are identified as signal in >= 50% of images (`#1710 <https://github.com/dials/dials/issues/1710>`_)
- ``dials.image_viewer``: Add selector to choose between a new default "image" and traditional "lab" coordinate frames. "image" frame attempts to align the fast/slow axes of the detector panels to screen x and y coordinates, so overall detector rotations will mostly be invisible. "lab" frame is the previous projection, where rotated detectors will appear rotated.
  ``dials.export_bitmaps``: Gained this same ``projection=`` option. (`#1716 <https://github.com/dials/dials/issues/1716>`_)
- ``dials.find_spots`` and ``dials.integrate``: `nproc=` now works with N > 1 on Windows. (`#1724 <https://github.com/dials/dials/issues/1724>`_)


Bugfixes
--------

- Fix rare crash in symmetry calculations when no resolution limit could be calculated (`#1641 <https://github.com/dials/dials/issues/1641>`_)
- ``dials.report``: Add units of pixels / images to centroid difference histograms (`#1677 <https://github.com/dials/dials/issues/1677>`_)
- ``dials.refine``: Scan-varying refinement failed when ``trim_scan_to_observations=False`` was used. (`#1686 <https://github.com/dials/dials/issues/1686>`_)
- ``dials.spot_counts_per_image``: Show an explicit error if given data that isn't spotfinding output (i.e. unindexed reflections/experiments). (`#1690 <https://github.com/dials/dials/issues/1690>`_)
- ``dials.integrate``: Improved background model variance calculation for integrating detectors. (`#1692 <https://github.com/dials/dials/issues/1692>`_)
- ``dials.stills_process``: improve processing performance by preventing re-reading of image data (`#1705 <https://github.com/dials/dials/issues/1705>`_)
- ``dials.background``: Correctly identify signal pixels for integrating detectors, and respect pre-calculated masks. (`#1726 <https://github.com/dials/dials/issues/1726>`_)
- ``dials.integrate``: Fixed bug in memory-use calculation for multi-sweep integration runs (`#1728 <https://github.com/dials/dials/issues/1728>`_)


Improved Documentation
----------------------

- Remove remaining 'master' references in the documentation. (`#1632 <https://github.com/dials/dials/issues/1632>`_)


Deprecations and Removals
-------------------------

- The previously deprecated ``dials.util.masking.MaskGenerator`` now prints a user warning. Please use ``dials.util.masking.generate_mask`` instead. (`#1643 <https://github.com/dials/dials/issues/1643>`_)
- ``dials.cosym``: Remove clustering code as this is no longer a necessary part of determination of symmetry or reindexing operations, and serves no useful purporse. (`#1647 <https://github.com/dials/dials/issues/1647>`_)
- ``dials.cosym``: ``nproc=`` parameter is deprecated. The algorithm is much faster on single cores. (`#1693 <https://github.com/dials/dials/issues/1693>`_)
- The pytest option ``--runslow`` was retired. The tests that it triggered will now always run. (`#1695 <https://github.com/dials/dials/issues/1695>`_)
- ``dev.dials.csv`` has been deprecated. Similar functionality is available with ``dials.export format=json``. (`#1708 <https://github.com/dials/dials/issues/1708>`_)
- ``dials.util.mp``: The ``preserve_exception_message`` argument has been deprecated. (`#1722 <https://github.com/dials/dials/issues/1722>`_)


Misc
----

- `#1631 <https://github.com/dials/dials/issues/1631>`_, `#1633 <https://github.com/dials/dials/issues/1633>`_, `#1648 <https://github.com/dials/dials/issues/1648>`_, `#1649 <https://github.com/dials/dials/issues/1649>`_, `#1652 <https://github.com/dials/dials/issues/1652>`_, `#1661 <https://github.com/dials/dials/issues/1661>`_, `#1672 <https://github.com/dials/dials/issues/1672>`_, `#1673 <https://github.com/dials/dials/issues/1673>`_, `#1674 <https://github.com/dials/dials/issues/1674>`_, `#1675 <https://github.com/dials/dials/issues/1675>`_, `#1676 <https://github.com/dials/dials/issues/1676>`_, `#1678 <https://github.com/dials/dials/issues/1678>`_, `#1679 <https://github.com/dials/dials/issues/1679>`_, `#1687 <https://github.com/dials/dials/issues/1687>`_, `#1696 <https://github.com/dials/dials/issues/1696>`_, `#1697 <https://github.com/dials/dials/issues/1697>`_, `#1698 <https://github.com/dials/dials/issues/1698>`_, `#1701 <https://github.com/dials/dials/issues/1701>`_, `#1706 <https://github.com/dials/dials/issues/1706>`_, `#1707 <https://github.com/dials/dials/issues/1707>`_, `#1711 <https://github.com/dials/dials/issues/1711>`_, `#1713 <https://github.com/dials/dials/issues/1713>`_, `#1717 <https://github.com/dials/dials/issues/1717>`_, `#1718 <https://github.com/dials/dials/issues/1718>`_, `#1720 <https://github.com/dials/dials/issues/1720>`_


DIALS 3.4.3 (2021-04-20)
========================

Bugfixes
--------

- ``dials.scale``: Fix crash when full-matrix minimisation is unsuccessful due to indeterminate normal equations. (`#1653 <https://github.com/dials/dials/issues/1653>`_)
- ``dials.scale``: Fix crash when no reflections remain after initial filtering. (`#1654 <https://github.com/dials/dials/issues/1654>`_)
- ``dials.export``: Fix error observed with ``format=mmcif`` for narrow sweeps with low symmetry (`#1656 <https://github.com/dials/dials/issues/1656>`_)
- Fix image numbering inconsistency in ascii histogram of per-image spot counts (`#1660 <https://github.com/dials/dials/issues/1660>`_)
- ``dials.find_spots_server``: Significant performance improvement for HDF5 grid scans. (`#1665 <https://github.com/dials/dials/issues/1665>`_)


DIALS 3.4.2 (2021-04-12)
========================

Bugfixes
--------

- Log messages from spot finding and integration no longer ignore logging level when using ``nproc > 1``. This mainly affects usage of dials from outside contexts. (`#1645 <https://github.com/dials/dials/issues/1645>`_)


DIALS 3.4.1 (2021-04-01)
========================

Features
--------

- ``dials.cosym``: Significantly faster via improved computation of functional, gradients and curvatures (`#1639 <https://github.com/dials/dials/issues/1639>`_)
- ``dials.integrate``: Added parameter ``valid_foreground_threshold=``, to require a minimum fraction of valid pixels before profile fitting is attempted (`#1640 <https://github.com/dials/dials/issues/1640>`_)


Bugfixes
--------

- ``dials.cosym``: Cache cases where Rij is undefined, rather than recalculating each time. This can have significant performance benefits when handling large numbers of sparse data sets. (`#1634 <https://github.com/dials/dials/issues/1634>`_)
- ``dials.cosym``: Fix factor of 2 error when calculating target weights (`#1635 <https://github.com/dials/dials/issues/1635>`_)
- ``dials.cosym``: Fix broken ``engine=scipy`` option (`#1636 <https://github.com/dials/dials/issues/1636>`_)
- ``dials.integrate``: Reject reflections with a high number of invalid pixels, which were being integrated since 3.4.0. This restores better merging statistics, and prevents many reflections being incorrect profiled as zero-intensity. (`#1640 <https://github.com/dials/dials/issues/1640>`_)


DIALS 3.4.0 (2021-03-15)
========================

Features
--------

- ``dials.integrate``: Profile-fitting improvements; Profile fitting will now be attempted on
  reflections with masked pixels, and the number of reflections qualifying for profile-fitting on
  multi-panel detectors has dramatically increased. (`#1297 <https://github.com/dials/dials/issues/1297>`_)
- ``dials.import``: When using ``reference_models=``, individual components of the model can be excluded with ``use_beam_reference=``, ``use_gonio_reference=`` and ``use_detector_reference=``. (`#1371 <https://github.com/dials/dials/issues/1371>`_)
- ``flex.reflection_table.match`` can now match reflections with configurable
  distance and scaling between any 3-vector column in the reflection tables. The
  default is still ``"xyzobs.px.value"``. (`#1398 <https://github.com/dials/dials/issues/1398>`_)
- ``dials.background``: Add option ``output.plot=`` to save an image to
  disk, instead of displaying interactively. Image files can now also be
  used directly. (`#1537 <https://github.com/dials/dials/issues/1537>`_)
- ``dials.import``: The default ``tolerance.scan.oscillation=`` is increased to
  3% of the image width, in order to accommodate electron diffraction datasets
  with poor rotation stages. (`#1543 <https://github.com/dials/dials/issues/1543>`_)
- ``dials.background``: Add support for multiple imagesets (`#1554 <https://github.com/dials/dials/issues/1554>`_)
- dials.estimate_resolution: reject Wilson outliers to minimise effect of spurious observations from e.g. ice rings on the resulting resolution estimates (`#1580 <https://github.com/dials/dials/issues/1580>`_)
- ``dials.cosym``: Use numpy in place of flex for large parts of cosym analysis (`#1581 <https://github.com/dials/dials/issues/1581>`_)
- ``dials.cosym``: Add option to use scipy `L-BFGS-B <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html>` minimization engine (``minimization.engine=scipy``) (`#1581 <https://github.com/dials/dials/issues/1581>`_)
- New masking parameter ``disable_parallax_correction=False``. Set to ``True`` to speed up generation of resolution masks by disabling parallax correction (this is only likely to have significant effect when spotfinding is spread across many independent processes). (`#1590 <https://github.com/dials/dials/issues/1590>`_)
- ``dials.image_viewer``: New parameter ``basis_vector_scale=`` to adjust the length of the basis vector overlay (`#1598 <https://github.com/dials/dials/issues/1598>`_)
- ``dials.merge``: add option to set wavelength_tolerance for MAD datasets (`#1609 <https://github.com/dials/dials/issues/1609>`_)
- ``dials.reciprocal_lattice_viewer``: Added an option to label the reciprocal lattice point nearest the centre (`#1614 <https://github.com/dials/dials/issues/1614>`_)
- ``dials.scale``: An additional outlier rejection based on normalised intensities has been added (`#1627 <https://github.com/dials/dials/issues/1627>`_)


Bugfixes
--------

- ``dials.image_viewer``: Fix various display issues relating to viewing still images (`#1463 <https://github.com/dials/dials/issues/1463>`_)
- ``dials.background``: Fix crash when writing output plot with bad display configuration (`#1550 <https://github.com/dials/dials/issues/1550>`_)
- ``dials.scale``: Fix issue of error model not always being carried through after
  the profile/summation intensity combination step. (`#1566 <https://github.com/dials/dials/issues/1566>`_)
- Fail bootstrap step if the git checkout fails in a non-interactive or non-posix environment (`#1572 <https://github.com/dials/dials/issues/1572>`_)
- Fixes working towards direct support of Windows builds:

  * Fix build errors by ensuring conda environment is correctly set up. (`#1575 <https://github.com/dials/dials/issues/1575>`_)
  * Fix importing using paths with wildcards (`#1583 <https://github.com/dials/dials/issues/1583>`_)
  * Fix ``dials.*`` commands crashing when unicode output is directed to a file (`#1602 <https://github.com/dials/dials/issues/1602>`_)
  * Fix some type-related test failures (`#1608 <https://github.com/dials/dials/issues/1608>`_)


Improved Documentation
----------------------

- Describe how to fix gltbx build failures for development installations on non-RHEL distributions (`#1561 <https://github.com/dials/dials/issues/1561>`_)
- Replace references to ``.pickle`` with reflections / ``.refl`` in docstrings (`#1619 <https://github.com/dials/dials/issues/1619>`_)
- Add documentation for ``dials.filter_reflections`` to the website. (`#1625 <https://github.com/dials/dials/issues/1625>`_)


Deprecations and Removals
-------------------------

- Remove previously deprecated ``use_trusted_range=`` parameter from masking configuration (`#1156 <https://github.com/dials/dials/issues/1156>`_)
- The main development branch of dials was renamed from 'master' to 'main'. (`#1546 <https://github.com/dials/dials/issues/1546>`_)
- ``dials.background``: The ``plot=`` parameter to interactively display the background plot has
  been removed. Use ``output.plot=`` to save to file instead. (`#1554 <https://github.com/dials/dials/issues/1554>`_)
- Remove ``*.o`` files from the DIALS installer package (`#1564 <https://github.com/dials/dials/issues/1564>`_)
- ``dials.util.masking.MaskGenerator`` is deprecated in favour of ``dials.util.masking.generate_mask`` (`#1569 <https://github.com/dials/dials/issues/1569>`_)


Misc
----

- `#1530 <https://github.com/dials/dials/issues/1530>`_, `#1531 <https://github.com/dials/dials/issues/1531>`_, `#1532 <https://github.com/dials/dials/issues/1532>`_, `#1534 <https://github.com/dials/dials/issues/1534>`_, `#1535 <https://github.com/dials/dials/issues/1535>`_, `#1536 <https://github.com/dials/dials/issues/1536>`_, `#1542 <https://github.com/dials/dials/issues/1542>`_, `#1567 <https://github.com/dials/dials/issues/1567>`_, `#1570 <https://github.com/dials/dials/issues/1570>`_, `#1571 <https://github.com/dials/dials/issues/1571>`_, `#1588 <https://github.com/dials/dials/issues/1588>`_, `#1593 <https://github.com/dials/dials/issues/1593>`_, `#1597 <https://github.com/dials/dials/issues/1597>`_, `#1599 <https://github.com/dials/dials/issues/1599>`_, `#1600 <https://github.com/dials/dials/issues/1600>`_, `#1601 <https://github.com/dials/dials/issues/1601>`_, `#1603 <https://github.com/dials/dials/issues/1603>`_, `#1604 <https://github.com/dials/dials/issues/1604>`_, `#1613 <https://github.com/dials/dials/issues/1613>`_, `#1620 <https://github.com/dials/dials/issues/1620>`_, `#1621 <https://github.com/dials/dials/issues/1621>`_, `#1624 <https://github.com/dials/dials/issues/1624>`_, `#1626 <https://github.com/dials/dials/issues/1626>`_, `#1630 <https://github.com/dials/dials/issues/1630>`_


DIALS 3.3.4 (2021-03-05)
========================

Bugfixes
--------

- ``dials.import``: Selecting a subset of images with ``image_range=`` now works on stills (`#1592 <https://github.com/dials/dials/issues/1592>`_)
- `dials.search_beam_centre`: Dramatically improved execution time for large data sets (`#1612 <https://github.com/dials/dials/issues/1612>`_)
- ``dials.reindex``: Write ``.refl`` file output in the default
  "MessagePack" format for better compatibility with downstream programs (`#1616 <https://github.com/dials/dials/issues/1616>`_)
- ``dials.scale``: Fix rare memory crash from infinite loop, that could
  occur with very bad quality datasets (`#1622 <https://github.com/dials/dials/issues/1622>`_)


Improved Documentation
----------------------

- ``dials.refine``: More informative error message when reflections have weights of zero (`#1584 <https://github.com/dials/dials/issues/1584>`_)


DIALS 3.3.3 (2021-02-15)
========================

No changes to core DIALS in 3.3.3.


DIALS 3.3.2 (2021-02-01)
========================

Bugfixes
--------

- Remove unnecessary call to ``imageset.get_raw_data()`` while generating
  masks. This was causing performance issues when spotfinding. (`#1449 <https://github.com/dials/dials/issues/1449>`_)
- ``dials.export``: Allow data with either summation or profile fitted
  intensities to be exported. Previously, both were (erroneously)
  required to be present. (`#1556 <https://github.com/dials/dials/issues/1556>`_)
- ``dials.scale``: Fix crash if only summation intensities present and ``intensity_choice=combine`` (`#1557 <https://github.com/dials/dials/issues/1557>`_)
- Fix unicode logging errors on Windows (`#1565 <https://github.com/dials/dials/issues/1565>`_)


DIALS 3.3.1 (2021-01-18)
========================

Features
--------

- ``dials.index``: More verbose debug logs when rejecting crystal models that are inconsistent with input symmetry (`#1538 <https://github.com/dials/dials/issues/1538>`_)


Bugfixes
--------

- ``dials.stills_process``: Fix spotfinding error "Failed to remap experiment IDs" (`#1180 <https://github.com/dials/dials/issues/1180>`_)
- Improved spotfinding performance for HDF5 when using a single processor. (`#1539 <https://github.com/dials/dials/issues/1539>`_)


DIALS 3.3.0 (2021-01-04)
========================

Features
--------

- DIALS is now using `GEMMI <https://gemmi.readthedocs.io/>`_. (`#1266 <https://github.com/dials/dials/issues/1266>`_)
- Upgrade ``h5py`` requirement to 3.1+ for SWMR-related functionality. (`#1495 <https://github.com/dials/dials/issues/1495>`_)
- Added support for small integer types to DIALS flex arrays. (`#1488 <https://github.com/dials/dials/issues/1488>`_)
- ``dials.estimate_resolution``: Only use cc_half in default resolution analysis. (`#1492 <https://github.com/dials/dials/issues/1492>`_)
- ``dials.export``: Allow on-the-fly bzip2 or gzip compression for mmCIF
  output, because unmerged mmCIF reflection files are large. (`#1480 <https://github.com/dials/dials/issues/1480>`_)
- ``dials.find_spots`` and ``dials.integrate`` both now have ``nproc=Auto`` by
  default, which uses the number of allowed/available cores detected. (`#1441 <https://github.com/dials/dials/issues/1441>`_)
- ``dials.merge``: Report ``<dF/s(dF)>``, if ``anomalous=True``. An html report
  is also generated to plot this statistic. (`#1483 <https://github.com/dials/dials/issues/1483>`_)
- ``dials.scale``: Apply a more realistic initial error model, or load the
  existing error model, if rescaling. (`#1526 <https://github.com/dials/dials/issues/1526>`_)
- ``dials.stills_process``: allow using different saturation cutoffs for
  indexing and integration. Useful for using saturated reflections for indexing
  while still rejecting them during integration. (`#1473 <https://github.com/dials/dials/issues/1473>`_)


Bugfixes
--------

- Internal: Logging metadata is now preserved when running spotfinding and
  integration across multiple processes. (`#1484 <https://github.com/dials/dials/issues/1484>`_)
- Fix NXmx behaviour with h5py 3.1. (`#1523 <https://github.com/dials/dials/issues/1523>`_)
- ``dials.cosym``: Choose the cluster containing the most identity reindexing
  ops by default. Under some circumstances, particularly in the case of
  approximate pseudosymmetry, the previous behaviour could result in reindexing
  operators being chosen that weren't genuine indexing ambiguities, instead
  distorting the input unit cells. (`#1514 <https://github.com/dials/dials/issues/1514>`_)
- ``dials.estimate_resolution``: Handle very low multiplicity datasets without
  crashing, and better error handling. (`#1494 <https://github.com/dials/dials/issues/1494>`_)
- ``dials.export``,``dials.two_theta_refine``: Updates to mmcif output to
  conform to latest pdb dictionaries (v5). (`#1528 <https://github.com/dials/dials/issues/1528>`_)
- ``dials.find_spots``: fix crash when ``nproc=Auto``. (`#1019 <https://github.com/dials/dials/issues/1019>`_)
- ``dials.image_viewer``: Fix crash on newer wxPython versions. (`#1476 <https://github.com/dials/dials/issues/1476>`_)
- ``dials.index``: Fix configuration error when there is more than one lattice
  search indexing method. (`#1515 <https://github.com/dials/dials/issues/1515>`_)
- ``dials.merge``: Fix incorrect output of SigF, N+, N- in ``merged.mtz``. (`#1522 <https://github.com/dials/dials/issues/1522>`_)
- ``dials.reciprocal_lattice_viewer``: Fix error opening with wxPython 4.1+. (`#1511 <https://github.com/dials/dials/issues/1511>`_)
- ``dials.scale``: fix issues for some uses of multi-crystal rescaling if ``full_matrix=False``. (`#1479 <https://github.com/dials/dials/issues/1479>`_)


Improved Documentation
----------------------

- Update information on how to care for an existing development environment,
  and remove outdated information. (`#1472 <https://github.com/dials/dials/issues/1472>`_)
- Each of the available indexing strategies in ``dials.index`` now has some
  help text explaining how it works. You can view this help by calling
  ``dials.index -c -a1 -e1`` and looking for ``method`` under ``indexing``. (`#1519 <https://github.com/dials/dials/issues/1519>`_)
- Include ``__init__`` methods in autodoc generated library documentation. (`#1520 <https://github.com/dials/dials/issues/1520>`_)
- ``dials.estimate_resolution``: Improved documentation. (`#1493 <https://github.com/dials/dials/issues/1493>`_)


Deprecations and Removals
-------------------------

- ``dials.algorithms.spot_finding.finder.SpotFinder``: Use of ``__call__`` to
  run spotfinding has been deprecated in favor of ``SpotFinder.find_spots(experiments)``. (`#1484 <https://github.com/dials/dials/issues/1484>`_)


Misc
----

- `#1469 <https://github.com/dials/dials/issues/1469>`_, `#1481 <https://github.com/dials/dials/issues/1481>`_,
  `#1484 <https://github.com/dials/dials/issues/1484>`_, `#1487 <https://github.com/dials/dials/issues/1487>`_,
  `#1491 <https://github.com/dials/dials/issues/1491>`_, `#1496 <https://github.com/dials/dials/issues/1496>`_,
  `#1497 <https://github.com/dials/dials/issues/1497>`_, `#1498 <https://github.com/dials/dials/issues/1498>`_,
  `#1499 <https://github.com/dials/dials/issues/1499>`_, `#1500 <https://github.com/dials/dials/issues/1500>`_,
  `#1501 <https://github.com/dials/dials/issues/1501>`_, `#1514 <https://github.com/dials/dials/issues/1514>`_.


DIALS 3.2.3 (2020-12-07)
========================

Bugfixes
--------

- ``dials.slice_sequence``: Fix crash using ``block_size=`` option (`#1502 <https://github.com/dials/dials/issues/1502>`_)
- ``dials.scale``: Fix broken ``exclude_images=`` option (`#1509 <https://github.com/dials/dials/issues/1509>`_)


DIALS 3.2.2 (2020-11-23)
========================

Bugfixes
--------

- Fix case where ``dials.stills_process`` could swallow error messages
- ``dials.cosym``: Fix non-determinism. Repeat runs will now give identical results. (`#1490 <https://github.com/dials/dials/issues/1490>`_)
- Developers: Fix precommit installation failure on MacOS (`#1489 <https://github.com/dials/dials/issues/1490>`_)


DIALS 3.2.1 (2020-11-09)
========================

3.2 Branch releases will now use a fixed conda environment. This release
is the first to use the same versions of all dependencies as 3.2.0.

Bugfixes
--------

- ``dials.symmetry``, ``dials.cosym`` and ``dials.two_theta_refine``: Lowered
  default partiality_threshold from ``0.99`` to to ``0.4``. The previous
  default could occasionally result in too many reflections being rejected for
  particularly narrow wedges. (`#1470 <https://github.com/dials/dials/issues/1470>`_)
- ``dials.stills_process`` Improve performance when using MPI by avoiding
  unnecessary log file writing (`#1471 <https://github.com/dials/dials/issues/1471>`_)
- ``dials.scale``: Fix scaling statistics output of r_anom data. (`#1478 <https://github.com/dials/dials/issues/1478>`_)


DIALS 3.2.0 (2020-10-27)
========================

Features
--------

- DIALS development environments are now running Python 3.8 by default.  (`#1373 <https://github.com/dials/dials/issues/1373>`_)
- Add a scaled flag to the reflection table. Indicates which reflections are
  good after the scaling process.  (`#1377 <https://github.com/dials/dials/issues/1377>`_)
- Python warnings are now highlighted on the console log and written to log files  (`#1401 <https://github.com/dials/dials/issues/1401>`_)
- Exit error messages from commands will now be colourized  (`#1420 <https://github.com/dials/dials/issues/1420>`_)
- Change the way ``dials.integrate`` splits data into blocks, to reduce
  unnecessary data reads, increasing performance up to 35% in some cases  (`#1396 <https://github.com/dials/dials/issues/1396>`_)
- Add ``dials.util.mp.available_cores`` function  (`#1430 <https://github.com/dials/dials/issues/1430>`_)
- ``dials.refine``: Trimming scans to observations for scan-varying refinement can
  now be turned off, using the parameter ``trim_scan_to_observations=False``  (`#1374 <https://github.com/dials/dials/issues/1374>`_)
- ``dials.refine``: Change default to ``separate_panels=False``. This speeds up
  outlier rejection for multi-panel detectors. For metrology refinement this
  should be set to ``True``  (`#1424 <https://github.com/dials/dials/issues/1424>`_)
- ``dials.merge``: Add best_unit_cell option. If the best_unit_cell option is set
  in ``dials.scale``, this will now propagate to the merged mtz output file.  (`#1444 <https://github.com/dials/dials/issues/1444>`_)
- DIALS bootstrap now allow creating a Python 3.9 environment  (`#1452 <https://github.com/dials/dials/issues/1452>`_)
- DIALS now uses pytype for limited static type checking. We hope that this
  will, over time, improve code quality.  (`#1364 <https://github.com/dials/dials/issues/1364>`_)
- ``dials.stills_process``: Added ``process_percent=`` to restrict processing
  to a subset of data, sync reference geometry instead of overwriting it and
  handle composite spotfinding modes.  (`#1409 <https://github.com/dials/dials/issues/1409>`_)


Bugfixes
--------

- ``dials.stills_process``: Prevent memory usage getting too high by clearing the
  imageset cache during processing.  (`#1412 <https://github.com/dials/dials/issues/1412>`_)
- ``dials.find_spots_server``: Return HTTP 500 instead of 200 when running fails  (`#1443 <https://github.com/dials/dials/issues/1443>`_)
- ``dials.find_spots_server``: Fix multiprocessing-related crash on macOS with Python3.8  (`#1447 <https://github.com/dials/dials/issues/1447>`_)
- ``dials.integrate``: Fix failures when building with GCC 9  (`#1456 <https://github.com/dials/dials/issues/1456>`_)
- ``dials.image_viewer``: Fix deprecation warnings  (`#1462 <https://github.com/dials/dials/issues/1462>`_)
- ``dials.index``: When using local index assignment, take into account phi in
  nearest neighbour analysis. This can significantly improve indexing rates in
  some cases with scans > 360°  (`#1459 <https://github.com/dials/dials/issues/1459>`_)
- ``dials.reindex``: Show an error instead of crashing for bad reindex operations.  (`#1282 <https://github.com/dials/dials/issues/1282>`_)

Deprecations and Removals
-------------------------

- dials.refine: the parameter ``trim_scan_edges`` is renamed ``scan_margin``
  and the former name is deprecated  (`#1374 <https://github.com/dials/dials/issues/1374>`_)
- The developer command ``dev.dials.show_test_failure_reasons`` was removed.  (`#1436 <https://github.com/dials/dials/issues/1436>`_)
- Remove clipper sources from new development installations  (`#1437 <https://github.com/dials/dials/issues/1437>`_)


Misc
----

- `#1175 <https://github.com/dials/dials/issues/1175>`_, `#1337 <https://github.com/dials/dials/issues/1337>`_,
  `#1354 <https://github.com/dials/dials/issues/1354>`_, `#1379 <https://github.com/dials/dials/issues/1379>`_,
  `#1381 <https://github.com/dials/dials/issues/1381>`_, `#1400 <https://github.com/dials/dials/issues/1400>`_,
  `#1406 <https://github.com/dials/dials/issues/1406>`_, `#1416 <https://github.com/dials/dials/issues/1416>`_,
  `#1423 <https://github.com/dials/dials/issues/1423>`_, `#1426 <https://github.com/dials/dials/issues/1426>`_,
  `#1432 <https://github.com/dials/dials/issues/1432>`_, `#1433 <https://github.com/dials/dials/issues/1433>`_,
  `#1435 <https://github.com/dials/dials/issues/1435>`_, `#1446 <https://github.com/dials/dials/issues/1446>`_,
  `#1454 <https://github.com/dials/dials/issues/1454>`_, `#1466 <https://github.com/dials/dials/issues/1466>`_,
  `#1468 <https://github.com/dials/dials/issues/1468>`_


DIALS 3.1.4 (2020-10-12)
========================

No changes to core DIALS in 3.1.4.


DIALS 3.1.3 (2020-09-28)
========================

Bugfixes
--------

- ``dials.integrate``: fix integrator=3d_threaded crash if njobs > 1 (`#1410 <https://github.com/dials/dials/issues/1410>`_)
- ``dials.integrate``: Check for and show error message if shoebox data is missing (`#1421 <https://github.com/dials/dials/issues/1421>`_)
- ``dials.refine``: Avoid crash for experiments with zero reflections if the
  `auto_reduction.action=remove` option was active (`#1417 <https://github.com/dials/dials/issues/1417>`_)

Improved Documentation
----------------------

- ``dials.merge``: improve help message by adding usage examples (`#1413 <https://github.com/dials/dials/issues/1413>`_)
- ``dials.refine``: More helpful error message when too few reflections (`#1431 <https://github.com/dials/dials/issues/1431>`_)


DIALS 3.1.2 (2020-09-14)
========================

Features
--------

- ``dials.stills_process``: Add parameter ``max_images=`` to limit the number
  of processed images

Bugfixes
--------

- ``dials.integrate``: fix crash when run with integrator=3d_threaded (`#1404 <https://github.com/dials/dials/issues/1404>`_)
- ``dials.integrate``: Minor performance improvements (`#1399 <https://github.com/dials/dials/issues/1399>`_)
- ``dials.stills_process``: MPI performance improvements for large datasets
- ``dials.stills_process``: Fix error when using split logs


DIALS 3.1.1 (2020-09-01)
========================

Bugfixes
--------

- ``dials.scale``: Prevent discarding of resolution limits in rare cases, which
  could cause incorrect symmetry determination, and worse results. (`#1378 <https://github.com/dials/dials/issues/1378>`_)
- ``dials.cosym``: filter out experiments with inconsistent unit cells (`#1380 <https://github.com/dials/dials/issues/1380>`_)
- Internally slicing experiments now works if image range doesn't start at 1 (`#1383 <https://github.com/dials/dials/issues/1383>`_)
- Restore missing I/sigma(I) resolution estimate log output (`#1384 <https://github.com/dials/dials/issues/1384>`_)
- ``dials.image_viewer``: "Save As" button now works, for single panels
- Fix developer ``libtbx.precommit`` installation error (`#1375 <https://github.com/dials/dials/issues/1375>`_)


DIALS 3.1 (2020-08-17)
======================

Features
--------

- Supports Python 3.7 and 3.8. Python 3.6 remains the default. (`#1236 <https://github.com/dials/dials/issues/1236>`_)
- Switch DIALS environment to use conda compilers. For development environments,
  a new ``dials`` script, located above the build directory, replaces the
  existing 'setpaths'-family of scripts. This means that all commands within
  the conda environment will now be available. (`#1235 <https://github.com/dials/dials/issues/1235>`_)
- New command: ``dials.missing_reflections`` to identify connected regions of
  missing reflections in the asymmetric unit. (`#1285 <https://github.com/dials/dials/issues/1285>`_)
- Improvements to image stacking in ``dials.image_viewer``:
  - add pull-down selector for stacking mode
  - add modes for mean and max
  - add command-line selection for stacking mode
  - rename ``sum_images`` command-line option to ``stack_images`` (`#1302 <https://github.com/dials/dials/issues/1302>`_)
- Reduce volume of output in ``dials.integrate``; histograms and other less
  important information only shows in debug output. Pass the ``-vv`` option
  to restore the previous behaviour (`#1319 <https://github.com/dials/dials/issues/1319>`_)
- ``dials.integrate``: Experimental feature: Specifying
  ``output_unintegrated_reflections=False`` discards unintegrated data from
  output reflection file, for smaller output and faster post-processing (`#1343 <https://github.com/dials/dials/issues/1343>`_)
- Rename ``dials.resolutionizer`` command to ``dials.estimate_resolution``,
  and includes a html report. Writing png plot output is now turned off by
  default (passing ``plot=True`` will restore this behaviour). (`#1330 <https://github.com/dials/dials/issues/1330>`_)
- ``dials.scale`` now separates anomalous pairs during error model analysis (`#1332 <https://github.com/dials/dials/issues/1332>`_)
- ``dials.background``: Add parameter ``corrected=`` to optionally use
  pedestal-and-gain corrected data (`#1348 <https://github.com/dials/dials/issues/1348>`_)
- ``dials.combine_experiments``: Add option ``output.max_reflections_per_experiment=``,
  to reject experiments with too many reflections (`#1369 <https://github.com/dials/dials/issues/1369>`_)


Bugfixes
--------

- ``dials.image_viewer``: Shoeboxes are now shown rotated with rotated detector panels. (`#1189 <https://github.com/dials/dials/issues/1189>`_)
- ``dials.index``: In multi-lattice indexing, ensure that reflections where
  refinement fails are flagged as unindexed. (`#1350 <https://github.com/dials/dials/issues/1350>`_)
- ``dials.scale``: Reflections excluded from scaling are no longer permanently
  excluded from any subsequent ``dials.scale`` jobs. (`#1275 <https://github.com/dials/dials/issues/1275>`_)
- ``dials.scale``: When using ``intensity_choice=combine`` (the default), don't
  exclude reflections that only have one of summed or profiled intensities
  available, but not both. (`#1300 <https://github.com/dials/dials/issues/1300>`_)
- ``dials.split_experiments``: Don't generate extra leading zeros in the output
  filename when not required e.g. ``output_09.expt`` -> ``output_9.expt`` (`#1316 <https://github.com/dials/dials/issues/1316>`_)
- ``dials.plot_reflections``: Fix invisible white spots on white background. (`#1346 <https://github.com/dials/dials/issues/1346>`_)


Deprecations and Removals
-------------------------

- ``dials.find_spots``: Deprecate ``spotfinder.filter.use_trusted_range=`` (`#1156 <https://github.com/dials/dials/issues/1156>`_)
- ``setpaths.sh`` and related scripts in newly created DIALS development
  environments are made obsolete and will no longer work. (`#1235 <https://github.com/dials/dials/issues/1235>`_)
- ``dials.show``: Remove ``show_image_statistics=`` parameter. Use
  ``image_statistics.show_raw=`` for equivalent output (`#1306 <https://github.com/dials/dials/issues/1306>`_)
- Log files will omit timings unless the relevant dials program was run with ``-v`` (`#1313 <https://github.com/dials/dials/issues/1313>`_)

Misc
----

- `#1184 <https://github.com/dials/dials/issues/1184>`_, `#1216 <https://github.com/dials/dials/issues/1216>`_, `#1288 <https://github.com/dials/dials/issues/1288>`_, `#1312 <https://github.com/dials/dials/issues/1312>`_, `#1320 <https://github.com/dials/dials/issues/1320>`_, `#1322 <https://github.com/dials/dials/issues/1322>`_, `#1325 <https://github.com/dials/dials/issues/1325>`_, `#1328 <https://github.com/dials/dials/issues/1328>`_, `#1352 <https://github.com/dials/dials/issues/1352>`_, `#1365 <https://github.com/dials/dials/issues/1365>`_, `#1366 <https://github.com/dials/dials/issues/1366>`_, `#1370 <https://github.com/dials/dials/issues/1370>`_


DIALS 3.0.4 (2020-07-20)
========================

- ``dials.scale``: Allow usage of ``mode=image_group`` with ``filtering.method=deltacchalf`` when
  only providing a single data set (`#1334 <https://github.com/dials/dials/issues/1334>`_)
- ``dials.import``: When using a template and specifying an image_range, missing images outside of
  the range will not cause a failure (`#1333 <https://github.com/dials/dials/issues/1333>`_)
- ``dials.stills_process``: Show better error message in specific spotfinding failure case (`#1180 <https://github.com/dials/dials/issues/1180>`_)


DIALS 3.0.3 (2020-07-06)
========================

Features
--------

- Developer tool: On posix systems, sending SIGUSR2 to DIALS commands will now print a stack trace (`#1277 <https://github.com/dials/dials/issues/1277>`_)

Bugfixes
--------
- HTML reports: Plot bin centres instead bin minimum for d_min line plots vs. resolution (`#1323 <https://github.com/dials/dials/issues/1323>`_)
- ``dials.export``: Fix inconsistency in mtz export when given a non-reference (e.g. I2 or primitive) setting (`#1279 <https://github.com/dials/dials/issues/1279>`_)
- ``dials.refine_bravais_settings``: Fix crash with large (>2gb) reflection tables and reduce memory use (`#1274 <https://github.com/dials/dials/issues/1274>`_)
- ``dials.scale``: Fix bug in outlier rejection code causing misidentification of outliers (with outlier_rejection=standard).
- ``dials.scale``: Fix outlier rejection formula to avoid overconfidence in spuriously low values


DIALS 3.0.2 (2020-06-23)
========================

Bugfixes
--------

- Fix crash in scaling error model handling (`#1243 <https://github.com/dials/dials/issues/1243>`_)


DIALS 3.0.1 (2020-06-11)
========================

Features
--------

- dials.reciprocal_lattice_viewer: Add an option to show lattice(s) in the crystal rather than laboratory frame. (`#1259 <https://github.com/dials/dials/issues/1259>`_)
- Support for mtz project_name in export and scaling

Bugfixes
--------

- dials.reciprocal_lattice_viewer: fix multiple experiment view for integrated data (`#1284 <https://github.com/dials/dials/issues/1284>`_)


DIALS 3.0 (2020-05-22)
======================

Features
--------

- Show more useful output when crashing in C++ code (`#659 <https://github.com/dials/dials/issues/659>`_)
- dials.image_viewer: for the unit cell tool, rename parameters for consistency and add a new show_hkl option to filter displayed powder rings to select only those of interest. (`#1192 <https://github.com/dials/dials/issues/1192>`_)
- In dials.integrate: changed the background box size multiplier to be a parameter (sigma_b_multiplier) - setting to small values significantly reduces memory requirements. (`#1195 <https://github.com/dials/dials/issues/1195>`_)
- dials.image_viewer: add an overlaying showing pixels marked as strong by the spot-finding operations. That is, the pixels picked out by the "threshold" image. (`#1200 <https://github.com/dials/dials/issues/1200>`_)
- dials.scale report file was renamed from scaling.html to dials.scale.html
  dials.symmetry report file was renamed from dials-symmetry.html to dials.symmetry.html (`#1202 <https://github.com/dials/dials/issues/1202>`_)
- dials.report output file was renamed from dials-report.html to dials.report.html (`#1206 <https://github.com/dials/dials/issues/1206>`_)
- dials.image_viewer: faster navigation between different image types. (`#1213 <https://github.com/dials/dials/issues/1213>`_)
- Crystal model now has a new recalculated_unit_cell attribute. This allows it to store
  a post-refined unit cell (e.g. from dials.two_theta_refine) in addition to that from
  traditional geometry refinement (which was used for prediction). Downstream programs
  such as dials.scale and dials.export will now use the recalculated unit cell
  where appropriate. (`#1214 <https://github.com/dials/dials/issues/1214>`_)
- New best_monoclinic_beta parameter for dials.refine_bravais_settings and dials.symmetry.
  Setting this to False will ensure that C2 is selected in preference to I2, where I2
  would lead to a less oblique cell (i.e. smaller beta angle). (`#1226 <https://github.com/dials/dials/issues/1226>`_)
- New scaling model, model=dose_decay, implementing a shared exponential decay component for multicrystal experiments (`#1183 <https://github.com/dials/dials/issues/1183>`_)


Bugfixes
--------

- Avoid empty "Unable to handle" messages on failed dials.import (`#600 <https://github.com/dials/dials/issues/600>`_)
- Functions from dials.export now raise exceptions on errors rather than exit. This improves their use elsewhere (such as in dials.scale). (`#1205 <https://github.com/dials/dials/issues/1205>`_)
- Ensure dials.index chooses the C2 setting with the smallest beta angle (`#1217 <https://github.com/dials/dials/issues/1217>`_)
- Fix propagation of best_unit_cell and application of resolution cutoffs in dials.scale and export_mtz.
  Add a new mtz.best_unit_cell parameter to dials.export (`#1248 <https://github.com/dials/dials/issues/1248>`_)
- Make some of the DIALS tools furthest downstream (``dials.scale``, ``dials.symmetry``, ``dials.merge`` and ``dials.report``) more robust in the case of very few reflections. (`#1263 <https://github.com/dials/dials/issues/1263>`_)


Misc
----

- `#1221 <https://github.com/dials/dials/issues/1221>`_


DIALS 2.2 (2020-03-15)
======================

Features
--------

- dials.image_viewer: Add a choice between displaying the raw or the corrected image. (`#634 <https://github.com/dials/dials/issues/634>`_)
- Automatically choose between the `simple` and `glm` background determination
  algorithms, depending on whether the detector appears to be integrating or
  counting. (`#706 <https://github.com/dials/dials/issues/706>`_)
- Allow adjustment of font size for overlay text, such as Miller indices and
  resolution ring values. (`#1074 <https://github.com/dials/dials/issues/1074>`_)
- Keep goniometer and scan objects in indexing of still data, if provided in input (`#1076 <https://github.com/dials/dials/issues/1076>`_)
- Experimental: ``dials.image_viewer`` can be remotely controlled via a
  ZeroMQ endpoint with the ``zmq_endpoint`` PHIL parameter. Initially,
  the viewer can be commanded to load new images. This requires the
  (optional) ``pyzmq``package. (`#1085 <https://github.com/dials/dials/issues/1085>`_)
- Programs now generate a unique identifier for each experiment created, and reflection tables are linked via the experiment_identifiers map (`#1086 <https://github.com/dials/dials/issues/1086>`_)
- Introduce `dials.anvil_correction` to correct the absorption of the incident and diffracted X-ray beam by the diamond anvils in a pressure cell.
  Call `dials.anvil_correction` on the output of `dials.integrate` and then proceed to use post-integration tools as normal, just as though the sample had been measured in air. (`#1090 <https://github.com/dials/dials/issues/1090>`_)
- Map of detector efficiency for photon counting detectors as a function of
  detector position added to report, based on the qe value applied at the end
  of integration. (`#1108 <https://github.com/dials/dials/issues/1108>`_)
- Significantly reduce the amount of memory required to write .refl output files (`#1115 <https://github.com/dials/dials/issues/1115>`_)
- Add maximum_trusted_value=N option to spot finding to temporarily allow override of trusted range, e.g. to find overloaded spots in spot finding. (`#1157 <https://github.com/dials/dials/issues/1157>`_)
- array_family.flex interface has changed: background and centroid algorithms are
  set via public properties. Instead of flex.strategy use functools.partial with
  the same signature. as_miller_array() raises KeyError instead of Sorry.
  .extract_shoeboxes() lost its verbosity parameter, use log levels instead. (`#1158 <https://github.com/dials/dials/issues/1158>`_)
- dials.stills_process now supports imagesets of length > 1 (e.g. grid scans) (`#1174 <https://github.com/dials/dials/issues/1174>`_)


Bugfixes
--------

- Fixed prediction on images numbered zero, so integrating works correctly. (`#1128 <https://github.com/dials/dials/issues/1128>`_)
- Fix an issue (`#1097 <https://github.com/dials/dials/issues/1097>`_) whereby aggregating small numbers of reflections into resolution bins could sometimes result in empty bins and consequent errors. (`#1130 <https://github.com/dials/dials/issues/1130>`_)
- Ensure that restraints are ignored for parameterisations that are anyway fixed (`#1142 <https://github.com/dials/dials/issues/1142>`_)
- Fix dials.search_beam_centre to ensure that the correct detector models are
  output when multiple detector models are present in the input.
  Fix dials.search_beam_centre n_macro_cycles option (previously it was starting
  from the original geometry every macro cycle). (`#1145 <https://github.com/dials/dials/issues/1145>`_)
- dials.find_spots_server no longer slows down 3x when using resolution filters (`#1170 <https://github.com/dials/dials/issues/1170>`_)


Misc
----

- `#932 <https://github.com/dials/dials/issues/932>`_, `#1034 <https://github.com/dials/dials/issues/1034>`_, `#1050 <https://github.com/dials/dials/issues/1050>`_, `#1077 <https://github.com/dials/dials/issues/1077>`_


DIALS 2.1 (2019-12-12)
======================

Features
--------

- We now fully support Python 3 environments.
- MessagePack is now the default reflection table file format. Temporarily, the
  environment variable ``DIALS_USE_PICKLE`` can be used to revert to the previous
  pickle-based format, however this will be removed in a future version. (`#986 <https://github.com/dials/dials/issues/986>`_)
- new option for dials.show 'show_shared_models=True' displays which beam, crystal, and detector models are used across experiments (`#996 <https://github.com/dials/dials/issues/996>`_)
- Import still image sequence as N experiments dereferencing into one image set
  rather than one experiment. (`#1014 <https://github.com/dials/dials/issues/1014>`_)
- Add `reflection_table.get` method for defaulted column access (`#1031 <https://github.com/dials/dials/issues/1031>`_)


Bugfixes
--------

- Don't use -2 to indicate masked pixels, except for DECTRIS detectors where this
  is to be expected. (`#536 <https://github.com/dials/dials/issues/536>`_)
- No longer show pixels that are above the trusted range upper bound as
  "saturated" on the "variance" image. (`#846 <https://github.com/dials/dials/issues/846>`_)
- Correctly account for scan-varying crystals while providing a scan range to
  dials.integrate (`#962 <https://github.com/dials/dials/issues/962>`_)
- Ensure that generated masks do not include pixels that are overloaded on a few
  images, but only pixels that are always outside the trusted range. (`#978 <https://github.com/dials/dials/issues/978>`_)
- Rewritten parameter auto-reduction code for dials.refine provides finer-grained
  fixing of individual parameters rather than whole parameterisations and
  correctly takes constrained parameters into account (`#990 <https://github.com/dials/dials/issues/990>`_)
- Fix output of predictions in dials.refine.
  A recently-introduced bug meant that the updated predictions weren't
  being copied to the output reflections file. (`#991 <https://github.com/dials/dials/issues/991>`_)
- Allow scan-varying refinement where either the crystal cell or
  orientation is fixed. (`#999 <https://github.com/dials/dials/issues/999>`_)
- Respect batch= option to dials.symmetry - can reduce time taken for finding
  the symmetry for large data sets. (`#1000 <https://github.com/dials/dials/issues/1000>`_)
- Scan-varying refinement no longer fails when the scan is wider than the
  observed reflections (e.g. when the crystal has died). Instead, the scan
  is first trimmed to match the range of the diffraction. (`#1025 <https://github.com/dials/dials/issues/1025>`_)
- If convert_sequences_to_stills then delete the goniometer and scan. (`#1035 <https://github.com/dials/dials/issues/1035>`_)
- Correctly account for scan-varying crystals in dials.slice_sequence (`#1040 <https://github.com/dials/dials/issues/1040>`_)
- Eliminate systematic absences before applying change of basis op to minimum
  cell in dials.symmetry. (`#1064 <https://github.com/dials/dials/issues/1064>`_)


Improved Documentation
----------------------

- Add "Extending DIALS" page to developer documentation (`#893 <https://github.com/dials/dials/issues/893>`_)


Deprecations and Removals
-------------------------

- The command dials.analyse_output was removed.
  Its replacement, dials.report, will give you more useful output. (`#1009 <https://github.com/dials/dials/issues/1009>`_)


Misc
----

- `#983 <https://github.com/dials/dials/issues/983>`_, `#1004 <https://github.com/dials/dials/issues/1004>`_


DIALS 2.0 (2019-10-23)
======================

Features
--------

- Support exporting multi-dataset and still experiments to XDS_ASCII (`#637 <https://github.com/dials/dials/issues/637>`_)
- Replace default spotfinder with improved dispersion algorithm (`#758 <https://github.com/dials/dials/issues/758>`_)
- ``dials.report`` now displays oscillation data with units and more significant figures (`#896 <https://github.com/dials/dials/issues/896>`_)
- A new program, ``dials.sequence_to_stills`` to create split a sequence into a
  separate still Experiment for every scan point in the sequence, splitting
  reflections as necessary. (`#917 <https://github.com/dials/dials/issues/917>`_)
- Moved ``dials.export format=best`` to ``dials.export_best`` as that one needed
  access to the format object, the rest do not, and having ``dials.export`` work
  in the general case seems like a better idea... (`#921 <https://github.com/dials/dials/issues/921>`_)
- Unified logging output for dials programs - logs are no longer split into .log
  and .debug.log. Use -v to get debug output. (`#923 <https://github.com/dials/dials/issues/923>`_)
- New command ``dials.resolutionizer`` (replaces ``xia2.resolutionizer``). Add support for ``expt``/``refl``
  in ``dials.resolutionizer``. (`#933 <https://github.com/dials/dials/issues/933>`_)
- Changed the selection of reflections used for determination of the reflection
  profile parameters in integration. Now uses reflections which were previously
  used in refinement rather than all reflections, resulting in a speed
  improvement for large data sets and a negligible difference in the quality
  of the integrated results. (`#942 <https://github.com/dials/dials/issues/942>`_)
- ``dials.image_viewer`` now allows the choice between
  ``dispersion_extended`` (new default) and ``dispersion`` (old default)
  thresholding algorithms for investigating the effect of different
  spot-finding parameters. (`#948 <https://github.com/dials/dials/issues/948>`_)
- ``dials.rs_mapper`` now respects masked regions of images (including
  the trusted range mask). (`#955 <https://github.com/dials/dials/issues/955>`_)


Bugfixes
--------

- Fix and reinstate normalisation option in ``dials.option`` (`#919 <https://github.com/dials/dials/issues/919>`_)


Misc
----

- `#795 <https://github.com/dials/dials/issues/795>`_, `#862 <https://github.com/dials/dials/issues/862>`_, `#895 <https://github.com/dials/dials/issues/895>`_, `#915 <https://github.com/dials/dials/issues/915>`_, `#924 <https://github.com/dials/dials/issues/924>`_
