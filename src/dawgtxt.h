"COMMAND LINE USAGE\n" \
"  dawg -[scvh?] file1 [file2...]\n" \
"  -s: process files serially\n" \
"  -c: process files combined together\n" \
"  -u: unbuffered output\n" \
"  -v: display version information\n" \
"  -h: display help information\n" \
"  -?: same as -h\n" \
"\n" \
"  Dawg will read stdin if filename is \"-\".\n" \
"\n" \
"FILE FORMAT\n" \
"  The file format takes a series of statements in the form of\n" \
"  \"name = value,\" where \"name\" is alphanumeric and value can\n" \
"  be a string, number, boolean, tree, or vector of values.\n" \
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
"  Sequence       VS  root sequences \n" \
"  Length         VN  length of generated root sequences\n" \
"  Rates          VVN rate of evolution of each root nucleotide\n" \
"  Model          S   model of evolution: GTR|JC|K2P|K3P|HKY|F81|F84|TN\n" \
"  Freqs          VN  nucleotide (ACGT) frequencies \n" \
"  Params         VN  parameters for the model of evolution\n" \
"  Width          N   block width for indels and recombination\n" \
"  Scale          VN  block position scales\n" \
"  Gamma          VN  coefficients of variance for rate heterogenity\n" \
"  Alpha          VN  shape parameters\n" \
"  Iota           VN  proportions of invariant sites\n" \
"  GapModel       VS  models of indel formation: NB|PL|US\n" \
"  Lambda         VN  rates of indel formation\n" \
"  GapParams      VN  parameter for the indel model\n" \
"  Reps           N   number of data sets to output\n" \
"  File           S   output file \n" \
"  Format         S   output format: Fasta|Nexus|Phylip\n" \
"  GapSingleChar  B   output gaps as a single character\n" \
"  GapPlus        B   distinguish insertions from deletions in alignment\n" \
"  LowerCase      B   output sequences in lowercase\n" \
"  NexusCode      S   text or file to include between datasets in Nexus format\n" \
"  Seed           VN  PRNG seed (integers)\n" \

