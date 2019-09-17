mismatch_crispri
================

Author: John S. Hawkins [really@gmail.com]

Python scripts to facilitate the design of single-mismatch CRISPRi guide
libraries and analyze the results of pooled growth experiments using those
libraries.

Introduction
------------

The basic anticipated workflow is:

[optional] train predictive model

->

design guides

->

count sequenced samples

->

compute gammas

->

convert to relative fitness

->

analyze relative fitness

Training a Predictive Model
---------------------------

For most users, this step is not necessary, as we have included our trained
model in this code, under gfpdata/.  If you wish to re-train the model
yourself to accommodate a different data source, an update to the machine
learning libraries, etc., you may do so by calling:

::

    ./train_linear_model.py

Designing Guides
----------------

The first step for most users will be to design a guide library.  First you
will need to locate/generate the following prerequisites:

* TARGETFILE: generally the output of
  https://www.github.com/traeki/sgrna_design, this is a tab-separated-value
  file containing annotated targets for a genome.

* LOCIFILE: a file containing one locus_tag per line, matching the locus_tag
  column of the TARGETFILE.  The code will only design guides for the listed
  subset of locus_tag entries in TARGETFILE.

then call

::

    ./choose_guides -h

to see the call signature.  --families, --n, and --divide_evenly control the
size and diversity of the library.  --outfile specifies the destination for
the designed guides.

When run successfully, this code will output a list of guides for each
targeted locus_tag containing a range of predicted knockdowns, in
tab-separated-value format.


Counting Sequenced Samples
--------------------------

Once you have constructed your library and sequenced one or more samples, you
will need to generate a file for each sample containing two tab separated
columns, "variant", and "count", which indicate the raw frequency with which
each guide was seen in that sample.  For convenience we provide
count_guides.py as is. It can be altered to suit your needs, but is somewhat
dependent on the particulars of the method by which your sequencing library is
prepared.  At minimum, you will need to set the --reverse flag depending on
the orientation of your reads relative to the sequences provided via
--guide_set.  You may also need to alter the internal base sequence the code
is matching, depending on your DNA construct.  Call

::

    ./count_guides.py -h

to see rough call signature, with the aforementioned caveats.

::

    ./count_guides.py

with no arguments will run the script on the sample files in testdata/ as a
demonstration.


Computing Gamma Values
----------------------

Once guide frequency counts have been generated (by count_guides.py or
otherwise), you will need to create a configuration file for each replicate
group for which you wish to compute gamma values.  See testdata/config.tsv for
an example.  This file should specify, for each replicate, the replicate
identifier, pre-growth countfile, and post-growth countfile.  The config should
be placed in the same directory as the count files it names, and should be
called "config.tsv"

To see how to run the script once these configuration files and directories
have been constructed, use compute_gammas.py.  Run

::

    ./compute_gammas.py -h

for call signature information, and call

::

    ./compute_gammas.py

with no arguments to compute gammas for the example data in testdata/.

This script generates a tsv file with the annotated gamma scores, and also a
file (using the same name, with .mean included prior to the final suffix) with
the replicate-averaged gamma values for each guide variant.


Converting to Relative Fitness
------------------------------

For most purposes, it is easier to reason about relative fitness (growth rate
as a fraction of wildtype growth rate) than about gamma.  Relative fitness is
just (1 + γ), so the computation is simple, but for convenience we have a
script that converts a tsv file from gamma to relative fitness.

::

    ./gamma_to_relfit.py -h

for call signature information, and

::

    ./gamma_to_relfit.py --gammafile <infile> --relfitfile <outfile>

to convert a file.  kvf_by_gene.py, below, assumes this has been done, and we
recommend you do this as a matter of course before any other downstream
analysis.  (Indeed, we may eventually change the code base to use this metric
by default.)


Analyze Fitness
---------------

Analysis will depend heavily on application.  We provide kvf_by_gene.py for
simple visualization of prediction vs outcome, broken down by locus_tag.

As usual,

::

    ./kfv_by_gene.py -h

gives usage information, and

::

    ./kfv_by_gene.py

with no arguments applies the script to the sample data in testdata/.
