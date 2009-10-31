DAWG VERSION 1.2-RELEASE

Copyright (c) (2004-2009) Reed A. Cartwright - All rights reserved.

DESCRIPTION

Dawg is an application that will simulate nucleotide evolution with
gaps.

ABSTRACT

DNA Assembly with Gaps (Dawg) is an application designed to simulate the
evolution of recombinant DNA sequences in continuous time based on the
robust general time reversible model with gamma and invariant rate
heterogeneity and a novel length-dependent model of gap formation. The
application accepts phylogenies in Newick format and can return the
sequence of any node, allowing for the exact evolutionary history to be
recorded at the discretion of users. Dawg records the gap history of
every lineage to produce the true alignment in the output. Many options
are available to allow users to customize their simulations and results.

Many tools and procedures exist for reconstructing alignments and
phylogenies and estimating evolutionary parameters from extant data.
True phylogenies and alignments are known in very rare instances. In the
absence of known data with true phylogenies, we are left with using
simulations to test the accuracy of such procedures. Proper simulation
of sequence evolution should involve both nucleotide substitution and
indel formation. However, existing tools for simulating sequence
evolution either do not include indels, like Seq-gen or evolver, or
include a rather inexact model of indel formation, like Rose. I
developed Dawg to fill in these gaps.

CONTACT

racartwr@ncsu.edu or reed@scit.us

Reed A. Cartwright, PhD
Postdoctoral Research Associate
Department of Genetics
Bioinformatics Research Center
North Carolina State University
Campus Box 7566
Raleigh, NC 27695-7566

Most work was done while I was a PhD student:

Department of Genetics
University of Georgia
Athens, GA

REFERENCE

Cartwright, R.A. (2005) DNA Assembly With Gaps (Dawg): Simulating Sequence
Evolution. Bioinformatics 21 (Suppl. 3): iii31-iii38

LICENSE

See COPYING for license information.

DOWNLOAD

Dawg can be downloaded from <http://scit.us/projects/dawg/>.

INSTALLATION

See Dawg's website for binary packages for Windows, Mac OSX, and other
systems.  Alternatively, you can compile Dawg from the source.  Dawg
requires CMake 2.6 (http://www.cmake.org/) to build it from sources.  Many
Unix-like operating systems can install CMake through their package
systems.  Extract the Dawg source code and issue the following commands in
the extracted directory:

    cmake .
    make
    make install

The '-G' option to cmake is used to specify different build systems, e.g. Unix
Makefiles versus KDevelop3 project.  The '-D' option to cmake can be used to
set different cmake variables from the command line:

    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr .
    make
    make install

This will build an optimized version of Dawg and install it to '/usr/bin'.
To specify your own build flags you need to set the environment variables
CFLAGS and LDFLAGS as necessary.  Then specify

    cmake -DCMAKE_BUILD_TYPE= .
    
See CMake's manual for additional information.

If you want to build the source code on Windows you will need to install Flex
and Bison from the Gnuwin32 project <http://gnuwin32.sourceforge.net/>, and make
sure that they are in your path.  You can then run the CMake GUI interface.  If
you would prefer to run the command line version, then open up a command console
through the Visual Studio tools shortcut (or similar shortcut).  This will add
the required compiler programs to your command console environment.  After
changing to the source code directory issue the following commands:

    cmake -G "NMake Makefiles" .
    nmake

If successful, you should find dawg.exe in the "src" directory.

If you are trying to compile Dawg on a UNIX machine that does not have CMake
installed, and you can't install it from a package, then you may need to install
it locally.  After downloading and extracting CMake in your home directory,
change to its directory and issue the following commands.

  ./configure --prefix=$HOME
  make
  make install
  
If "make" fails, try using "gmake" instead.

EXAMPLES

example0.dawg - minimal
example1.dawg - typical usage
example2.dawg - simple indel formation
example3.dawg - robust indel formation
example4.dawg - recombination 

COMMAND LINE USAGE

dawg -[scubvhqew?] [-o outputfile] file1 [file2...]  
  -s: process files serially [default]
  -c: process files combined together
  -u: unbuffered output
  -b: buffered output [default]
  -q: disable error and warning reports (quiet)
  -e: enable error reports [default]
  -w: enable warning reports [default]
  -v: display version information
  -h: display help information
  -?: same as -h
  -o outputfile: override ouput filename in the configuration file

  Dawg will read stdin if filename is "-".

FILE FORMAT

The file format takes a series of statements in the form of "[name]" or 
"name = value", where "name" is alphanumeric and value can be a string,
number, Boolean, tree, or vector of values.  The former specifies a heading,
which can simplify variable assignment.  A single variable is equivalent to
a vector of a single entry.  

When using headings, the following statements are equivalent:

  Out.Block.Head = "A Comment"
  Out.Block.Tail = "B Comment"

  [Out.Block]
  Head = "A Comment"
  Tail = "B Comment"

  Out.Block.Head = "A Comment"
  [Out.Block]
  Tail = "B Comment"

  [Out.Block]
  Head = "A Comment"
  []
  Out.Block.Tail = "B Comment"

  [Out]
  Block.Head = "A Comment"
  [Out.Block]
  Tail = "B Comment"

  [Out]
  Block.Head = "A Comment"
  [.Block]
  Tail = "B Comment"

Values can be specified via the following syntaxes.

string:  "[char-sequence]"
         '[char-sequence]'
         """[multi-line char-sequence]""" (removes initial and final newlines)
         '''[multi-line char-sequence]''' (preserves initial and final newlines)
number:  [sign]digits[.digits][(e|E)[sign]digits]
boolean: true|false
tree:    Newick Format
vector:  { value, value, ...}

OPTIONS

  Name            Type            Description
--------------------------------------------------------------------------
  Tree             VT  phylogeny
  TreeScale        N   coefficient to scale branch lengths by
  Sequence         VS  root sequences 
  Length           VN  length of generated root sequences
  Rates            VVN rate of evolution of each root nucleotide
  Model            S   model of evolution: GTR|JC|K2P|K3P|HKY|F81|F84|TN
  Freqs            VN  nucleotide (ACGT) frequencies 
  Params           VN  parameters for the model of evolution
  Width            N   block width for indels and recombination
  Scale            VN  block position scales
  Gamma            VN  coefficients of variance for rate heterogenity
  Alpha            VN  shape parameters
  Iota             VN  proportions of invariant sites
  GapModel         VS  models of indel formation: NB|PL|US
  Lambda           VN  rates of indel formation
  GapParams        VVN parameter for the indel model
  Reps             N   number of data sets to output
  File             S   output file 
  Format           S   output format: Fasta|Nexus|Phylip|Clustal
  GapSingleChar    B   output gaps as a single character
  GapPlus          B   distinguish insertions from deletions in alignment
  KeepFlank        N   undeletable flanking regions N nucs from sequence
  KeepEmpty        B   preserve empty columns in final alignment
  LowerCase        B   output sequences in lowercase
  Translate        B   translate outputted sequences to amino acids
  Seed             VN  pseudo-random-number-generator seed (integers)
  Out.Block.Head   S   string to insert at the start of the output
  Out.Block.Tail   S   string to insert at the end of the output
  Out.Block.Before S   string to insert before a sequence set in the output
  Out.Block.After  S   string to insert after a sequence set in the output
  Out.Subst        B   do variable substitution in Out.Block.*

DEFAULTS

  TreeScale = 1.0
  Length = 100
  Model = "JC"
  Freqs = {0.25,0.25,0.25,0.25}
  Params = {1.0,1.0,1.0,1.0,1.0,1.0}
  Width = 1
  Scale = 1.0
  Gamma = 0.0
  Iota =  0.0
  GapModel = "US"
  GapParams = 1.0
  Reps = 1
  Format = "Fasta"
  GapSingleChar = false
  GapPlus = false
  LowerCase = false
  Translate = false
  Out.Subst = true

VARIABLE SUBSTITUTION

If Out.Subst is true (the default), then Dawg will preform variable substitution
in any Out.Block that it outputs.  Currently three variables are supported.
  %r is replaced by the current dataset number
  %R is replaced by the total dataset number
  %% is replaced by a percent sign.

OUTPUT FILE

Dawg can automatically detect the format of the output file based on its extension.
Supported extensions and their formats are:

  Clustal: aln, poo, txt, out, Clustal
  Fasta:   fas, Fasta
  Nexus:   nex, Nexus
  Phylip:  phy, Phylip

Dawg also supports the filename format of "ext:file" to output to "file" with
the format specified by extension "ext".  That way one can use "nex:-" to output
to stdout in Nexus format.

NOTES

The meaning of the "Params" vector is different for each substitution model.
  GTR: Substitution rates A-C, A-G, A-T, C-G, C-T, G-T
  JC:  Ignored
  K2P: Transition rate, Transversion rate
  K3P: Alpha (Transitions), Beta (A-T & G-C), Gamma (A-C & G-T)
  HKY: Transition rate, Transversion rate
  F81: Ignored
  F84: Kappa
  TN:  Alpha1 (A-G), Alpha2 (C-T), Beta (Transversions)

Parameter "Freqs" is ignored by the models "JC", "K2P", and "K3P".

If "Lambda" is a single value, then it specifies the rate of indel formation,
e.g. "Lambda = 0.1" is the same as "Lambda = {0.05, 0.05}".  The first
parameter is the insertion rate and the second parameter is the deletion rate.

The first parameter of "GapModel" specifies the distribution model of 
insertion sizes. The second parameter specifies the distribution model of 
deletion sizes.  If only one parameter is given it is the model for both 
insertions and deletions.

The first parameter of "GapParams" is a vector specifying the parameters for the
gap model of insertions.  Likewise the second parameter is a vector specifying 
the parameters for the gap model of deletions.  If "GapParams" is not a vector 
of vectors, then it specifies the vector of parameters for both insertions and 
deletions.

The meaning of the GapParams vector is different for each gap model.
  US: The distribution of gap sizes.
  NB: The number of failures (r), the probability of success (q).
  PL: The rate parameter (a), the maximum gap size.

To create a recombinant tree, you may need to specifically describe and label
the inner nodes at which the recombination events occur.  See example4.dawg.

Gamma takes precedence over Alpha.

Sequence takes precedence over Length.

If Out.Block.* is the name of a file, the code is read from that file.

The following vector parameters have a size of "Width": "Scale", "Alpha", 
"Gamma", and "Iota".  If their size is less than width then the first value in 
the vector will be used to fill in the rest of the values, e.g. "Scale = 1.0" 
is the same as "Scale = {1.0,1.0,1.0}" when "Width = 3".
