DAWG VERSION 2-CURRENT

Copyright (c) (2004-2010) Reed A.  Cartwright - All rights reserved.

DESCRIPTION

Dawg is an application that will simulate sequence evolution with gaps.

ABSTRACT

DNA Assembly with Gaps (Dawg) is an application designed to simulate the
evolution of recombinant DNA sequences in continuous time based on the robust
general time reversible model with gamma and invariant rate heterogeneity and a
novel length-dependent model of gap formation.  The application accepts
phylogenies in Newick format and can return the sequence of any node, allowing
for the exact evolutionary history to be recorded at the discretion of users.
Dawg records the gap history of every lineage to produce the true alignment in
the output.  Many options are available to allow users to customize their
simulations and results.

Many tools and procedures exist for reconstructing alignments and phylogenies
and estimating evolutionary parameters from extant data.  True phylogenies and
alignments are known in very rare instances.  In the absence of known data with
true phylogenies, we are left with using simulations to test the accuracy of
such procedures.  Proper simulation of sequence evolution should involve both
nucleotide substitution and indel formation.  However, existing tools for
simulating sequence evolution either do not include indels, like Seq-gen or
evolver, or include a rather inexact model of indel formation, like Rose.  I
developed Dawg to fill in these gaps.

CONTACT

racartwright@uh.edu or reed@scit.us

REFERENCE

Cartwright, R.A. (2005) DNA Assembly With Gaps (Dawg): Simulating Sequence
Evolution.  Bioinformatics 21 (Suppl. 3): iii31-iii38

LICENSE

See copying.txt for license information.

DOWNLOAD

Dawg can be downloaded from <http://scit.us/projects/dawg/>.

PREREQUISITES

If installing from source you, need to ensure that the development libraries of
Boost <http://www.boost.org/> and GSL <http://www.gnu.org/software/gsl/> are
installed and findable on your machine.  Dawg can compile with older versions of
Boost, but requires a recent version of Boost's Spirit library.  If the version
of Boost on your machine doesn't comes with a recent version of Spirit, just
download a recent version of the boost library and copy the header files for
spirit into [dawg source]/src/include/boost/spirit.  That way they will be found
by the compiler before older versions of spirit.

If you have to install Boost locally on a Unix machine, the following works for
me:

  cd boost_1_44_0
  ./bootstrap.sh --prefix=$HOME
  ./bjam install
  cd ../dawg-current
  cmake -DBOOST_ROOT=$HOME . make

Change as needed.

INSTALLATION

See Dawg's website for binary packages for Windows, Mac OSX, and other systems.
Alternatively, you can compile Dawg from the source.  Dawg requires CMake 2.6
(http://www.cmake.org/) to build it from sources.  Many Unix-like operating
systems can install CMake through their package systems.  Extract the Dawg
source code and issue the following commands in the extracted directory:

  cmake .
  make
  make install

The '-G' option to cmake is used to specify different build systems, e.g.  Unix
Makefiles versus KDevelop3 project.  The '-D' option to cmake can be used to set
different cmake variables from the command line:

  cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr .
  make
  make install

This will build an optimized version of Dawg and install it to '/usr/bin'.  To
specify your own build flags you need to set the environment variables CFLAGS
and LDFLAGS as necessary.  Then specify

  cmake -DCMAKE_BUILD_TYPE= .

See CMake's manual for additional information.

If you would prefer to run the command line version, then open up a command
console through the Visual Studio tools shortcut (or similar shortcut).  This
will add the required compiler programs to your command console environment.
After changing to the source code directory issue the following commands:

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

COMMAND LINE USAGE

dawg [options] trick.dawg

Use "dawg --help" for help information.  Dawg will read stdin if filename is
"-".

Use "dawg --help-trick" for help regarding input files.

DESCRIPTION OF SIMULATION

Dawg splits the simulation of a single replicate into jobs based on the tree and
sequence section.  Different areas of the tree and different areas of the
sequence can have different evolutionary models.  By varying the phylogeny in
different parts of the sequence you can produce recombinant sequences.

At the first level, the root sequence is split into a series of "segments",
which evolve independently, except for deletions which may span many segments.
Next the phylogeny of each segment is split into multiple sections allowing
different evolutionary models for different parts of the tree.  This all can be
easily controlled from the input file format.

INPUT FILE FORMAT

Dawg is controlled by a series of input files, referred to as tricks.  Tricks
contain a series of Sections which define different models for different
sequence or tree regions:

  [[SectionA]]

  Parameter.A = valueA
  Parameter.B = valueB, valueC

  [[SectionB]]

  Parameter.C = valueZ
  Parameter.B = valueD, valueDD

Double square brackets define a new section.  At the start of a trick file there
is an implied [[_initial_]] section header, allowing you to skip specifying one
if it is unneeded.

By default sections inherit the values of the section above it.  Section
_initial_ inherits from an implied _default_ section.  In the example above,
SectionA inherits the values from _initial_, and SectionB from SectionA.  This
makes it easy to specify a new section without having to write out everything.
For SectionB, Parameter.A has valueA, just like in SectionA.

You can also change the default behavior of inheriting from the previous
section:

  [[SectionC = SectionA]]

Here SectionC will inherit the values from SectionA, not SectionB.  If you use a
blank header, [[]], the name of the section will be generated automatically.

Another shortcut is to use parameter headers:

  [[SectionA]]
  
  Foo.Bar = A
  Foo.Par = B
  
  [[SectionB]]
  
  [Foo]
  Bar = A
  Par = B

Both sections produce identical results.  Parameters in Dawg are named such that
the trick files can be simplified using headers.  Some parameter headers are
also special:

 []           # clears the current header
 [.part]      # add a new part to the current header
 [..part]     # replace last part of the header
 [....part]   # replace last two parts
 
If 'part' is blank it simply deletes the last part of the current header.

Each parameter line in a trick file contains a parameter id, equals, and a list
of strings, separated by commas.  Ids can contain one or more numbers, letters,
dashes, underscores, and periods: e.g. [A-Za-z0-9._-]+.  There are four
different ways to specify strings.

Bare Strings contain a series of non-space characters except, excluding ,#"[]=()

Tree Strings lack spaces and start with '(' and end with ';'.

Quoted Strings are a series of printable characters between double quotation
marks.  Newlines are not acceptable

Triple Quoted Strings can contain any character between two sets of three double
quotes.

An example of all four string types:

  AList = Bare_String, (Tree:0.1,String:0.1);,
  "Quoted String",
  """Double
           Quoted
                 String"""
 
Comments start with '#' and go to the end of the line.

OUTPUT FILE

Dawg can automatically detect the format of the output file based on its
extension.  Supported formats and their extensions are:

  Clustal: aln
  Fasta:   fasta, fas, fsa
  Nexus:   nexus, nex
  Phylip:  phylip, phy
  Poo:     poo

Dawg also supports the filename format of "ext:file" to output to "file" with
the format specified by extension "ext".  That way one can use "nex:-" to output
to stdout in Nexus format.  Partial matches of extensions are allowed.

If the --split option is on, the each replicate will be saved to its own file,
based on the filename in the output option.

NOTES

The meaning of the "Params" vector is different for each substitution model.

  GTR: Substitution rates A-C, A-G, A-T, C-G, C-T, G-T
  JC: Ignored
  K2P: Transition rate, Transversion rate
  K3P: Alpha (Transitions), Beta (A-T & G-C), Gamma (A-C & G-T)
  HKY: Transition rate, Transversion rate
  F81: Ignored
  F84: Kappa
  TN: Alpha1 (A-G), Alpha2 (C-T), Beta (Transversions)

Parameter "Freqs" is ignored by the models "JC", "K2P", and "K3P".

If "Lambda" is a single value, then it specifies the rate of indel formation,
e.g. "Lambda = 0.1" is the same as "Lambda = {0.05, 0.05}".  The first parameter
is the insertion rate and the second parameter is the deletion rate.

The first parameter of "GapModel" specifies the distribution model of insertion
sizes.  The second parameter specifies the distribution model of deletion sizes.
If only one parameter is given it is the model for both insertions and
deletions.

The first parameter of "GapParams" is a vector specifying the parameters for the
gap model of insertions.  Likewise the second parameter is a vector specifying
the parameters for the gap model of deletions.  If "GapParams" is not a vector
of vectors, then it specifies the vector of parameters for both insertions and
deletions.

The meaning of the GapParams vector is different for each gap model.  US: The
distribution of gap sizes.  NB: The number of failures (r), the probability of
success (q).  PL: The rate parameter (a), the maximum gap size.

To create a recombinant tree, you may need to specifically describe and label
the inner nodes at which the recombination events occur.

Gamma takes precedence over Alpha.

Sequence takes precedence over Length.

If Out.Block.* is the name of a file, the code is read from that file.

The following vector parameters have a size of "Width": "Scale", "Alpha",
"Gamma", and "Iota".  If their size is less than width then the first value in
the vector will be used to fill in the rest of the values, e.g. "Scale = 1.0" is
the same as "Scale = {1.0,1.0,1.0}" when "Width = 3".

