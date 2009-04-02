"COMMAND LINE USAGE\n" \
"  dawg -[scubvh?] file1 [file2...]\n" \
"  -s: process files serially [default]\n" \
"  -c: process files combined together\n" \
"  -u: unbuffered output\n" \
"  -b: buffered output [default]\n" \
"  -v: display version information\n" \
"  -h: display help information\n" \
"  -?: same as -h\n" \
"\n" \
"  Dawg will read stdin if filename is \"-\".\n" \
"\n" \
"FILE FORMAT\n" \
"  The file format takes a series of statements in the form of \"name = value,\"\n" \
"  where \"name\" is alphanumeric and value can be a string, number, boolean,\n" \
"  tree, or vector of values.  A single variable is equivalent to a vector of\n" \
"  a single entry.\n" \
"\n" \
"  string:  \"[char-sequence]\"\n" \
"           <<EOF [several lines] EOF\n" \
"  number:  [sign]digits[.digits][(e|E)[sign]digits]\n" \
"  boolean: true|false\n" \
"  tree:    Newick Format\n" \
"  vector:  { value, value, ...}\n" \
"\n" \
"OPTIONS\n" \
"  Name          Type            Description\n" \
"--------------------------------------------------------------------------\n" \
"  Tree           VT  phylogeny\n" \
"  TreeScale      N   coefficient to scale branch lengths by\n" \
"  Sequence       VS  root sequences\n" \
"  Length         VN  length of generated root sequences\n" \
"  Rates          VVN rate of evolution of each root nucleotide\n" \
"  Model          S   model of evolution: GTR|JC|K2P|K3P|HKY|F81|F84|TN\n" \
"  Freqs          VN  nucleotide (ACGT) frequencies\n" \
"  Params         VN  parameters for the model of evolution\n" \
"  Width          N   block width for indels and recombination\n" \
"  Scale          VN  block position scales\n" \
"  Gamma          VN  coefficients of variance for rate heterogenity\n" \
"  Alpha          VN  shape parameters\n" \
"  Iota           VN  proportions of invariant sites\n" \
"  GapModel       VS  models of indel formation: NB|PL|US\n" \
"  Lambda         VN  rates of indel formation\n" \
"  GapParams      VVN parameter for the indel model\n" \
"  Reps           N   number of data sets to output\n" \
"  File           S   output file\n" \
"  Format         S   output format: Fasta|Nexus|Phylip|Clustal\n" \
"  GapSingleChar  B   output gaps as a single character\n" \
"  GapPlus        B   distinguish insertions from deletions in alignment\n" \
"  EmptyCol       B   preserve empty columns in final alignment\n" \
"  LowerCase      B   output sequences in lowercase\n" \
"  Translate      B   translate outputed sequences to amino acids\n" \
"  NexusCode      S   text or file to include between datasets in Nexus format\n" \
"  Seed           VN  PRNG seed (integers)\n" \

