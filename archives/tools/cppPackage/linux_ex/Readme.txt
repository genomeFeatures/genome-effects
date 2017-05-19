Nonnegative matrix factorization for motif identification program package.
Release 0.1

Joel Graber lab, The Jackson Laboratory, Bar Harbor, ME
contact: joel.graber@jax.org

-----------------------------------------------------------------------------

Copyright (c) 2008 The Jackson Laboratory

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE. 

-----------------------------------------------------------------------------

This package includes executables only that run a series of programs 
designed to analyze sets of sequences aligned on a common feature or 
functional site, such as a polaydenylation or splice site.

Source code will be released at a later date, but as of this release, we are
including only linux binary executables.

The following packages/libraries must also be available on any machine on which this 
set of programs will be run:

1) Weblogo:  available at http://weblogo.berkeley.edu/  
The location of the weblogo scripts must be hard-coded into the main controlling
script "nmf_pipeline.prl"

2) Gnuplot, version 4.0 or later.  The executable is expected to be in the path such
that a call to "gnuplot" will successfully start the program.  This must be installed 
with PNG as an output option.  The "Arial" font types are also expected to be available,
typically in the ghostscript fonts directory.  (needed for the line plot labels)

3) The Gnu Scientific Library (GSL), specifically the dynamic libraries libgsl.so and 
libgslcblas.so.  These can be obtained from www.gnu.org/software/gsl/

If you have problems making the programs run, please contact Joel Graber at joel.graber@jax.org.

-----------------------------------------------------------------------------

Running the programs. 

While the programs can be run individually, this documentation assumes that everything 
will be controlled from the master perl script nmf_pipeline.prl, and provides detailed 
descriptions for the input parameters for this program only.  All other program inputs and
outputs are set automatically when run from nmf_pipeline.prl.

Output files will be named with a user provided prefix, and will also include automatically
generated portions that include k-mer and window size.  If run to completion, the last 
program in the analysis will generate an html formatted summary of the result.  It is highly 
recommended that this be used.  The raw output is pretty ugly, especially if you're not
used to it.

Mandatory arguments:

-i filename      : designates the input file of sequences.  This is assumed to be fasta-
                   formatted, with all sequences of identical length and already aligned
                   on their common feature.

-a integer (1-11): starting point in the analysis

-o string        : output file prefix  (all generated files will start with this)

-L integer       : length of each sequence in the input fasta file


Optional arguments:
-p integer (1-11): stopping point in the analysis  (default=11, generation of the html summary)
-d integer       : if set to anything other than 0, it will cause print out of commands, but will
                   not execute any of them.
-k integer       : k-mer size to count (default=4)
-w integer       : window size in the positional word count (default=5)
-r integer       : number of patterns to generate (default=8)
-z integer       : zero point for referencing positions in the sequence (default=0)
-L integer       : maximum allowed length of input lines
-t float         : tolerance threshold for convergence of nmf; smaller values are stricter convergence
                   (default = 1.0e-6)

-smN integer     : maximum number of k-mers to use in the nmf (default=all)

-nmfI integer    : maximum number of iterations per NMF restart (default=1000);
-nmfS integer    : number of NMF restarts

-xLab string     : label for the zero point in the sequence alignments (default="3'-processing site")

-bmS integer     : number of solutions for each motif in conversion of word list to motif (default=300)
-bmMinW integer  : smallest motif size to consider (default=4)
-bmMaxW integer  : largest motif size to consider (default=8)
-bmP integer     : number of iterations without improved model for termination of motif search (default=250)
-bmI integer     : maximum number of iterations to attempt in each solution (default=2000)
-bmx integer     : minimum overlap between k-mers and the existing motif model (controls the extent
                   that the k-mer can extend into the background, default=1)
-bmff float      : buildmotif renormalization factor,  Power to which the weighted word list elements 
                   are raised prior to motif build (empirically found to be best at about 1.4 for large 
                   data sets (>1000) and up to 2.0 for smaller sets (<50))  (default=1.0) 

-mmMinIC float   : trimming threshold more motifs.  We trim the final motif eliminating positions at
                   either beginning or end with information content less than this threshold.  (default=0.2)
-mmMinW integer  : minimum allowed size of the final motif after trimming

-nEx integer     : number of examples of each motif to generate for the file that is used for the sequence
                   logo generation

Examination of the nmf_pipeline.prl script will reveal additional variables, but at this point, we advise against
changing them.

Good luck!

-----------------------------------------------------------------------------

Programs run in the analysis (-a and -p options control where to start and stop in these)
(1) fa2Stats: collections mono- and di-nucleotide frequencies from the input file
(2) WindowCount: performs the PWC analysis
(3) nmfSmoothAndPseudoCounts: smooths the PWC output rows and adds pseudocounts
(4) nnmf: performs the NMF decomposition
(5) nmfSortMatrix: reorders the output of nnmf by their peak position density
(6) nmfWplot.prl: generates a script to plot the positioning matrices (piped into gnuplot)
(7) nmfBuildMotifs: performs the MCMC motif generation from the nmf W matrix
(8) nmfMotifsToModels: converts the format ofthe motifs to that expected by pwmToExamples
(9) pwmToExamples: generates example files sampled from each motif for feeding into weblogo
(10) weblogo: generates sequence logo representations of all motifs (weblogo must be separately installed)
(11) nmfWebPage: formats all output into an html formatted summary


