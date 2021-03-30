Introduction
============

``profitroll`` is a high-level Object-Oriented Python package which simulates the atmospheric dynamics near the tropopause as well as water-vapor satellite images. It is a perfectly suited sandbox for meteorological algorithm implementation as its architecture is not problem dependent. Indeed, the architecture handles any 2D time-problem. I/O are managed through the use of NetCDF files for the sake of generality. The user is then free to implement his own scientific algorithm on top of the architecture.

The current version of ``profitroll`` is focused on simulating the evolution of a stratospheric intrusion into the tropopause. In particular, few assumptions are made:
- Neglect effects from the lower boundary (Earth surface)
- Quasigeostrophic theory
- f-plane (uniform Coriolis coefficient)
- Boussinesq approximation

Specifically, optimised built-in methods are provided to simulate satellite observations in the water-vapor channel.

Motivation
**********

This programme is the result of an academic project in collaboration with the "Centre National de Recherche en Météorologie". It has been identified as an excellent tool for conducting preliminary / 1st order studies of the tropopause's dynamics at very low computational cost. To provide an order of magnitude, running single simulation takes less than five minutes on a mid-range laptop, thanks to various code optimisations. The programme is also meant to be extended by whoever wants to add his own scientific methods. In that respect, ``profitroll`` was made completely open source.

Limitations
***********

- ``profitroll`` is currently restricted to quasigeostrophic dynamics and one would have to further implement new scientific methods to account for more various typical phenomena in meteorology. Yet the programme architecture is flexible enough to build third-party classes and methods on top of the existing version.

- The use of Python 3 is highly recommended. The compatibility with Python 2 is not garanteed.

- While interactive plot features different colormaps (user's choice), videos are currently restricted to a single colormap. An easy fix consists in diretly modifying the source code of the corresponding python file.