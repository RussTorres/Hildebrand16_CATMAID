# Hildebrand16
Simplified CATMAID analysis examples for [Hildebrand et al. public dataset](https://neurodata.io/data/hildebrand16/).


Originally developed as part of the catmaid_utils package for [Brett Graham](https://github.com/braingram), [Wei-Chung Lee](https://github.com/wclee), and [David G.C. Hildebrand](https://github.com/davidhildebrand)
For use with [CATMAID](https://github.com/catmaid/CATMAID) software by [Albert Cardona](https://github.com/acardona), [Stephan Saalfeld](https://github.com/axtimwalde), [Andrew Champion](https://github.com/aschampion), and [Tom Kazimiers](https://github.com/tomka).
Developed in collaboration with Brett Graham, [Logan Thomas](https://github.com/Lathomas42), [Rui Zheng](https://github.com/rui14), [Brendan Shanny](https://github.com/brenshanny), [Alex Coda](https://github.com/alexcoda), and David G.C. Hildebrand.

## Pathlength and morphology metrics
representing skeleton morphologies as directed acyclic graphs allows simple analysis of neuronal morphology using the [networkx](https://github.com/networkx/networkx) library.
Important analysis functionality has been included in the [morphology module](catmaid_analysis/algorithms/morphology.py), while additional graph-based functionality can be built from graph representations in the [neuron class.](catmaid_analysis/neuron.py)

## Kalman-filter based neurite smoothing
Residuals due to tracing and alignment errors can add up significantly in near-isotropic EM reconstruction datasets.  Quality of reconstructions can be improved using a kalman smoother on individual neurites (implemented using the [pykalman](https://github.com/pykalman/pykalman) library.)
See our [example Jupyter notebook.](notebooks/smoothed_pathlength.ipynb)

## Copy the Hildebrand16 annotations to your own CATMAID database
It is possible to use [a Python script](sql/import_dataset.py) to add Hildebrand16 reconstructions to another [configured](http://catmaid.readthedocs.io/en/stable/installation.html) CATMAID project.  Although primary keys are not retained, other properties (names, annotations, creation and edition data, etc.) are, and could be used to link datasets.

## Further Expansion
CATMAID has functionality beyond the simple skeleton GET requests used in this example allowing for intricate graph-based connectivity and morphology analysis.  Please look at the [CATMAID swagger doc]('http://hildebrand16.neurodata.io/catmaid/apis') for usage examples.
